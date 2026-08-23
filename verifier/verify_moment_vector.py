#!/usr/bin/env python3
"""Recompute the frozen order-10 moment vector in exact rational arithmetic.

This does not trust the cached ``mom_term_exact.pkl`` values.  It rebuilds every
coordinate from the regenerated order-9 flag matrices, the regenerated 9->10
deletion lift, and the public manifest-Gram atoms.  Equality with the frozen
vector is checked as an equality of exact GMP rationals.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import pickle
import sys
import time
from fractions import Fraction
from math import comb, prod
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

try:
    from gmpy2 import mpq
except ModuleNotFoundError:  # Exact, portable, and slower.
    from exact_arithmetic import mpq


_WORK = {}


def compute_atom_chunk(task):
    """Compute a disjoint group of Gram atoms in a forked worker."""
    chunk_number, atom_indices = task
    n_order10 = _WORK["n_order10"]
    subtotal = [mpq(0)] * n_order10
    for atom_index in atom_indices:
        label = _WORK["labels"][atom_index]
        vector = _WORK["vectors"][atom_index]
        weight = _WORK["weights"][atom_index]
        integer_matrices, denominators, safe_bound = _WORK["label_info"][label]
        integer_vector = np.asarray(
            [int(round(float(value) * 10**6)) for value in vector], dtype=np.int64
        )
        if safe_bound >= 2**63:
            raise OverflowError(f"int64 Gram bound failed for label {label}: {safe_bound}")
        numerators = np.einsum(
            "i,tij,j->t",
            integer_vector,
            integer_matrices,
            integer_vector,
            dtype=np.int64,
            optimize=True,
        )
        if weight:
            for row, numerator in enumerate(numerators):
                if not numerator:
                    continue
                scaled_value = weight * mpq(
                    int(numerator), 10**12 * int(denominators[row])
                )
                for col, coefficient in _WORK["deletion_by_row"][row]:
                    subtotal[col] += coefficient * scaled_value
    return chunk_number, len(atom_indices), subtotal


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--public-anc", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    flagsdp = args.flagsdp.resolve()
    public_anc = args.public_anc.resolve()
    sys.path.insert(0, str(flagsdp))

    started = time.time()
    with (flagsdp / "cache_n9.pkl").open("rb") as handle:
        cache = pickle.load(handle)
    lift = np.load(flagsdp / "c5lift_cache.npz", allow_pickle=True)
    with (public_anc / "moment_gram_w.pkl").open("rb") as handle:
        gram = pickle.load(handle)
    with (public_anc / "mom_term_exact.pkl").open("rb") as handle:
        expected_raw = pickle.load(handle)

    states = cache["states"]
    n_order9 = len(states)
    n_order10 = int(lift["nJ"])
    if n_order9 != 1897 or n_order10 != 12172:
        raise RuntimeError(
            f"catalogue dimensions differ: order9={n_order9}, order10={n_order10}"
        )

    deletion = csr_matrix(
        (lift["Dval"], (lift["Drow"], lift["Dcol"])),
        shape=(n_order9, n_order10),
    ).tocoo()
    deletion_by_row = [[] for _ in range(n_order9)]
    for row, col, value in zip(deletion.row, deletion.col, deletion.data):
        deletion_by_row[int(row)].append(
            (int(col), mpq(int(round(float(value) * 10)), 10))
        )

    label_info = {}
    for label, _size, sigma, _flags, free_size, _flat, integer_matrices in cache["moments"]:
        root_size = sigma[0]
        denominators = np.asarray(
            [
                int(
                    prod(order - i for i in range(root_size))
                    * comb(order - root_size, free_size) ** 2
                )
                if order - root_size >= free_size
                else 1
                for order, _adjacency in states
            ],
            dtype=np.int64,
        )
        label_info[label] = (
            np.asarray(integer_matrices, dtype=np.int64),
            denominators,
            None,
        )

    support = list(gram["support"])
    support_weights = [gram["w"][index] for index in support]
    labels = gram["atoms_lab"]
    vectors = gram["atoms_vv"]
    if not (len(support) == len(labels) == len(vectors)):
        raise RuntimeError("manifest-Gram support arrays have inconsistent lengths")

    for label, (matrix_stack, denominators, _unused) in list(label_info.items()):
        label_vectors = [
            vector
            for atom_label, vector in zip(labels, vectors)
            if atom_label == label
        ]
        max_vector_entry = max(
            abs(int(round(float(value) * 10**6)))
            for vector in label_vectors
            for value in vector
        )
        max_matrix_mass = max(int(np.abs(matrix).sum()) for matrix in matrix_stack)
        safe_bound = max_vector_entry * max_vector_entry * max_matrix_mass
        label_info[label] = (matrix_stack, denominators, safe_bound)

    rational_weights = []
    for weight in support_weights:
        fraction = Fraction(float(weight)).limit_denominator(10**6)
        rational_weights.append(mpq(fraction.numerator, fraction.denominator))
    if any(weight < 0 for weight in rational_weights):
        raise RuntimeError("a rationalized manifest-Gram weight is negative")

    _WORK.update(
        n_order10=n_order10,
        labels=labels,
        vectors=vectors,
        weights=rational_weights,
        label_info=label_info,
        deletion_by_row=deletion_by_row,
    )
    jobs = max(1, min(args.jobs, len(support)))
    chunks = [list(range(offset, len(support), jobs)) for offset in range(jobs)]
    tasks = [(number, indices) for number, indices in enumerate(chunks) if indices]
    actual = [mpq(0)] * n_order10
    if jobs == 1:
        results = map(compute_atom_chunk, tasks)
        pool = None
    else:
        context = mp.get_context("fork")
        pool = context.Pool(processes=jobs)
        results = pool.imap_unordered(compute_atom_chunk, tasks)
    completed = 0
    try:
        for chunk_number, atom_count, subtotal in results:
            completed += atom_count
            for index, value in enumerate(subtotal):
                actual[index] += value
            print(
                f"chunk {chunk_number + 1}/{len(tasks)} complete; "
                f"atoms {completed}/{len(support)}; elapsed={time.time() - started:.1f}s",
                flush=True,
            )
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    expected = [mpq(numerator, denominator) for numerator, denominator in expected_raw]
    mismatches = [index for index, pair in enumerate(zip(actual, expected)) if pair[0] != pair[1]]
    print(f"order-9 states: {n_order9}")
    print(f"order-10 states: {n_order10}")
    print(f"manifest-Gram atoms: {len(support)}")
    print(f"nonnegative rational weights: {all(weight >= 0 for weight in rational_weights)}")
    print(f"exact coordinate matches: {n_order10 - len(mismatches)}/{n_order10}")
    if mismatches:
        first = mismatches[0]
        print(f"first mismatch at state {first}")
        print(f"actual   = {actual[first]}")
        print(f"expected = {expected[first]}")
        return 1
    print("PASS: public exact moment vector was independently reconstructed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
