#!/usr/bin/env python3
"""Turn a feasible floating fixed-Gram dual into an exact rational certificate.

All nonnegative multipliers except rho use one decimal denominator.  The two
envelope-leg equalities/inequalities are repaired over the integers.  Finally
rho is defined as the exact maximum of the 12,172 rational coordinate loads,
so dual feasibility is a theorem of the emitted integer data rather than a
floating-point tolerance statement.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import pickle
import time
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

try:
    import gmpy2
except ModuleNotFoundError:  # Exact, portable, and slower.
    import exact_arithmetic as gmpy2


K7_DENOMINATOR = 10 * 25 * 181440
COMMON_FUNCTIONAL_DENOMINATOR = K7_DENOMINATOR
DIGIT_BASE = 100_000


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--public-anc", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--dual", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--denominator", type=int, default=10**12)
    parser.add_argument("--target-n", type=int, default=41)
    parser.add_argument("--threads", type=int, default=8)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_generator(args):
    spec = importlib.util.spec_from_file_location("fixed_gram_generator", args.generator)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    namespace = argparse.Namespace(
        flagsdp=args.flagsdp,
        public_anc=args.public_anc,
        state=args.state,
        max_iterations=0,
        threads=args.threads,
        k7_restarts=1,
        k7_keep=1,
        horn_keep=1,
        pool_cap=700,
        slack_drop=1e-6,
        resume=True,
        crossover=False,
        solver="ipm",
    )
    return module.CertificateGenerator(namespace)


def build_exact_functional_matrix(generator):
    """Rows are descriptor functionals scaled to denominator 45,360,000."""
    row_indices = []
    row_values = []
    indptr = [0]
    started = time.time()
    for row, descriptor in enumerate(generator.descriptors):
        numerator, denominator = generator.functional(descriptor)
        denominator = int(denominator)
        if COMMON_FUNCTIONAL_DENOMINATOR % denominator:
            raise RuntimeError(f"functional denominator {denominator} does not divide common one")
        factor = COMMON_FUNCTIONAL_DENOMINATOR // denominator
        nonzero = np.flatnonzero(numerator)
        values = np.asarray(numerator[nonzero], dtype=np.int64)
        if factor != 1:
            if values.size and int(np.max(np.abs(values))) > np.iinfo(np.int64).max // factor:
                raise OverflowError("scaled functional entry does not fit int64")
            values = values * factor
        row_indices.append(nonzero.astype(np.int32))
        row_values.append(values)
        indptr.append(indptr[-1] + len(nonzero))
        if (row + 1) % 500 == 0:
            print(
                f"exact rows {row + 1}/{len(generator.descriptors)}; nnz={indptr[-1]}",
                flush=True,
            )
    indices = np.concatenate(row_indices) if row_indices else np.empty(0, dtype=np.int32)
    values = np.concatenate(row_values) if row_values else np.empty(0, dtype=np.int64)
    matrix = csr_matrix(
        (values, indices, np.asarray(indptr, dtype=np.int64)),
        shape=(len(generator.descriptors), generator.n10),
        dtype=np.int64,
    )
    matrix.check_format(full_check=True)
    print(
        f"exact functional matrix: shape={matrix.shape} nnz={matrix.nnz} "
        f"build={time.time() - started:.1f}s",
        flush=True,
    )
    return matrix


def exact_transpose_product(matrix, multipliers):
    """Compute A^T y exactly via base-B digits, with proven int64 headroom."""
    multipliers = np.asarray(multipliers, dtype=np.int64)
    if np.any(multipliers < 0):
        raise RuntimeError("negative descriptor multiplier")
    absolute_column_sums = np.asarray(abs(matrix).sum(axis=0)).ravel().astype(np.int64)
    max_column_sum = int(absolute_column_sums.max(initial=0))
    digit_bound = max_column_sum * (DIGIT_BASE - 1)
    if digit_bound >= np.iinfo(np.int64).max:
        raise OverflowError(
            f"digit product could overflow int64: bound={digit_bound}"
        )

    remaining = multipliers.copy()
    pieces = []
    while np.any(remaining):
        digit = remaining % DIGIT_BASE
        piece = np.asarray(matrix.T @ digit, dtype=np.int64).ravel()
        pieces.append(piece)
        remaining //= DIGIT_BASE
    if not pieces:
        pieces = [np.zeros(matrix.shape[1], dtype=np.int64)]
    print(
        f"exact sparse product: digits={len(pieces)} base={DIGIT_BASE} "
        f"int64_bound={digit_bound}",
        flush=True,
    )
    return pieces


def combine_digits_at(pieces, index):
    value = 0
    scale = 1
    for piece in pieces:
        value += int(piece[index]) * scale
        scale *= DIGIT_BASE
    return value


def main():
    args = parse_args()
    denominator = int(args.denominator)
    target_n = int(args.target_n)
    if denominator <= 0 or denominator > 10**15:
        raise ValueError("denominator must lie in [1, 10^15]")
    if target_n <= 0:
        raise ValueError("target-n must be positive")
    generator = load_generator(args)
    with args.dual.open("rb") as handle:
        dual = pickle.load(handle)
    if dual.get("model_status") not in ("Optimal", "Unknown"):
        raise RuntimeError(f"unsupported source dual status {dual.get('model_status')}")
    source_violation = max(
        float(dual.get("lower_violation", float("inf"))),
        float(dual.get("upper_violation", float("inf"))),
        float(dual.get("multiplier_violation", float("inf"))),
        float(dual.get("equality_error", float("inf"))),
    )
    if source_violation > 1e-7:
        raise RuntimeError(f"source dual is not feasible enough: {source_violation}")
    if len(dual["descriptors"]) != len(generator.descriptors):
        raise RuntimeError("source dual and state descriptor counts differ")
    for left, right in zip(dual["descriptors"], generator.descriptors):
        if generator.descriptor_key(left) != generator.descriptor_key(right):
            raise RuntimeError("source dual and state descriptor order differ")

    values = np.asarray(dual["values"], dtype=float)
    n_descriptors = len(generator.descriptors)
    if len(values) != 6 + n_descriptors:
        raise RuntimeError("unexpected source multiplier vector length")

    # Static upper rows: density high, density low, fixed Gram, K7 leg, K8 leg.
    # Retain the density multipliers: later pools can make a band endpoint
    # active even though they were negligible in the n<=48 certificate.
    density_high_num = max(0, int(np.rint(values[0] * denominator)))
    density_low_num = max(0, int(np.rint(values[1] * denominator)))
    moment_num = max(0, int(np.rint(values[2] * denominator)))
    k7_leg_num = int(np.rint(values[3] * denominator))
    k7_leg_num = min(max(k7_leg_num, 0), denominator)
    k8_leg_num = denominator - k7_leg_num
    descriptor_nums = np.rint(values[5 : 5 + n_descriptors] * denominator).astype(np.int64)
    descriptor_nums = np.maximum(descriptor_nums, 0)

    # Repair the two envelope legs exactly, root by root.  Adding to a row has
    # zero RHS cost; any coordinate cost is paid exactly by rho below.
    repair_total = 0
    repair_count = 0
    for kind, target, nroots in (
        ("k7", k7_leg_num, generator.n7),
        ("k8", k8_leg_num, generator.n8),
    ):
        by_root = [[] for _ in range(nroots)]
        for index, descriptor in enumerate(generator.descriptors):
            if descriptor["kind"] == kind:
                by_root[int(descriptor["root"])].append(index)
        for root, indices in enumerate(by_root):
            if not indices:
                raise RuntimeError(f"no {kind} descriptor for root {root}")
            current = sum(int(descriptor_nums[index]) for index in indices)
            if current < target:
                # Repair the already-largest multiplier to avoid creating a
                # new tiny support entry.
                chosen = max(indices, key=lambda index: values[5 + index])
                increment = target - current
                descriptor_nums[chosen] += increment
                repair_total += increment
                repair_count += 1

    functional_matrix = build_exact_functional_matrix(generator)
    pieces = exact_transpose_product(functional_matrix, descriptor_nums)
    with (args.public_anc / "mom_term_exact.pkl").open("rb") as handle:
        moment_raw = pickle.load(handle)
    if len(moment_raw) != generator.n10:
        raise RuntimeError("exact moment vector has the wrong dimension")

    D = gmpy2.mpz(denominator)
    L = gmpy2.mpz(COMMON_FUNCTIONAL_DENOMINATOR)
    rho = None
    argmax = None
    loads = []
    started = time.time()
    for state, (moment_numerator, moment_denominator) in enumerate(moment_raw):
        functional_numerator = combine_digits_at(pieces, state)
        load = (
            gmpy2.mpq(moment_num * int(moment_numerator), denominator * int(moment_denominator))
            + gmpy2.mpq(functional_numerator, D * L)
            - gmpy2.mpq(
                (density_high_num - density_low_num)
                * int(generator.density_numerators[state]),
                D * 45,
            )
        )
        loads.append(load)
        if rho is None or load > rho:
            rho = load
            argmax = state
    if rho is None:
        raise RuntimeError("empty coordinate set")
    minimum_residual = min(rho - load for load in loads)
    if minimum_residual < 0:
        raise AssertionError("rho maximum construction failed")

    objective = (
        rho
        + gmpy2.mpq(density_high_num * 3197, D * 10000)
        - gmpy2.mpq(density_low_num * 2486, D * 10000)
        - gmpy2.mpq(2 * k8_leg_num, 25 * denominator)
    )
    threshold = gmpy2.mpq(2, 25 * target_n * target_n)
    margin = threshold - objective

    # Independent integer checks of all non-q coordinates.
    if k7_leg_num + k8_leg_num != denominator:
        raise AssertionError("eta equality repair failed")
    for kind, target, nroots in (
        ("k7", k7_leg_num, generator.n7),
        ("k8", k8_leg_num, generator.n8),
    ):
        sums = [0] * nroots
        for multiplier, descriptor in zip(descriptor_nums, generator.descriptors):
            if descriptor["kind"] == kind:
                sums[int(descriptor["root"])] += int(multiplier)
        if min(sums) < target:
            raise AssertionError(f"{kind} envelope repair failed")

    payload = {
        "format": "erdos23-fixed-gram-exact-dual-v1",
        "multiplier_denominator": denominator,
        "static_numerators": [
            density_high_num,
            density_low_num,
            moment_num,
            k7_leg_num,
            k8_leg_num,
        ],
        "descriptor_numerators": [int(value) for value in descriptor_nums],
        "rho": (int(gmpy2.numer(rho)), int(gmpy2.denom(rho))),
        "objective": (int(gmpy2.numer(objective)), int(gmpy2.denom(objective))),
        "global_closure": bool(objective <= 0),
        "target_n": target_n,
        "threshold_target": (int(gmpy2.numer(threshold)), int(gmpy2.denom(threshold))),
        "margin_target": (int(gmpy2.numer(margin)), int(gmpy2.denom(margin))),
        "argmax_state": int(argmax),
        "minimum_coordinate_residual": (
            int(gmpy2.numer(minimum_residual)),
            int(gmpy2.denom(minimum_residual)),
        ),
        "descriptors": generator.descriptors,
        "source": {
            "dual_sha256": sha256(args.dual),
            "state_sha256": sha256(args.state),
            "moment_exact_sha256": sha256(args.public_anc / "mom_term_exact.pkl"),
            "floating_objective": float(dual["objective"]),
            "floating_status": dual["model_status"],
        },
        "repair": {
            "roots_repaired": repair_count,
            "total_numerator_added": repair_total,
        },
    }
    with args.output.open("wb") as handle:
        pickle.dump(payload, handle, protocol=4)
    if args.json_output is not None:
        args.json_output.write_text(json.dumps(payload, separators=(",", ":")) + "\n")
    print(
        f"exact rho argmax={argmax}; coordinate_check={time.time() - started:.1f}s; "
        f"repairs={repair_count} total={repair_total}/{denominator}",
        flush=True,
    )
    print(f"exact objective={objective} ~= {float(objective):+.15e}", flush=True)
    print(f"global closure (objective <= 0): {objective <= 0}", flush=True)
    print(f"n={target_n} threshold={threshold} ~= {float(threshold):+.15e}", flush=True)
    print(f"exact margin={margin} ~= {float(margin):+.15e}", flush=True)
    print(f"saved {args.output}", flush=True)
    if args.json_output is not None:
        print(f"saved {args.json_output}", flush=True)
    if margin <= 0:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
