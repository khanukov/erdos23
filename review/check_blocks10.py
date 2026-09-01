#!/usr/bin/env python3
"""Independent check of blocks10.c on a sample of states.

For each sampled state G and each type size s in the requested levels, this
script recomputes the moment blocks M_sigma(G) in plain Python -- its own
canonical forms, its own flag indexing -- and compares them with the C helper's
blocks for q = indicator of G through indexing-free invariants: the sorted
eigenvalues, the total entry sum, and the sorted row sums.  Any mismatch is a
bug in one of the two implementations.
"""

from __future__ import annotations

import argparse
import itertools
import pickle
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from blocks10 import Blocks10  # noqa: E402

NV = 10


def induced(adj, order):
    """adjacency bitmasks of the subgraph on `order`, relabelled 0..k-1."""
    k = len(order)
    out = [0] * k
    for i in range(k):
        for j in range(k):
            if i != j and (adj[order[i]] >> order[j]) & 1:
                out[i] |= 1 << j
    return tuple(out)


def canonical_type(sub):
    """min over all relabellings of the adjacency tuple, and the relabellings attaining it."""
    k = len(sub)
    best, attaining = None, []
    for perm in itertools.permutations(range(k)):
        moved = [0] * k
        for i in range(k):
            for j in range(k):
                if (sub[i] >> j) & 1:
                    moved[perm[i]] |= 1 << perm[j]
        moved = tuple(moved)
        if best is None or moved < best:
            best, attaining = moved, [perm]
        elif moved == best:
            attaining.append(perm)
    return best, attaining


def canonical_flag(sub, s):
    """canonical form of a flag: roots 0..s-1 fixed, free vertices permuted."""
    k = len(sub)
    f = k - s
    best = None
    for perm in itertools.permutations(range(f)):
        full = list(range(s)) + [s + p for p in perm]
        moved = [0] * k
        for i in range(k):
            for j in range(k):
                if (sub[i] >> j) & 1:
                    moved[full[i]] |= 1 << full[j]
        moved = tuple(moved)
        if best is None or moved < best:
            best = moved
    return best


def python_blocks(adj, s):
    """dict: canonical type -> dict (F1, F2) -> count, over all ordered root tuples."""
    f = (NV - s) // 2
    blocks = {}
    for subset in itertools.combinations(range(NV), s):
        sub = induced(adj, subset)
        canon, attaining = canonical_type(sub)
        rest = [v for v in range(NV) if v not in subset]
        block = blocks.setdefault(canon, {})
        for perm in attaining:            # every ordering of the subset that induces exactly `canon`
            roots = [None] * s
            for i in range(s):
                roots[perm[i]] = subset[i]
            roots = tuple(roots)
            assert induced(adj, roots) == canon
            for A in itertools.combinations(rest, f):
                B = tuple(v for v in rest if v not in A)
                F1 = canonical_flag(induced(adj, roots + A), s)
                F2 = canonical_flag(induced(adj, roots + B), s)
                block[(F1, F2)] = block.get((F1, F2), 0) + 1
    return blocks


def dense(block):
    flags = sorted({F for pair in block for F in pair})
    index = {F: i for i, F in enumerate(flags)}
    M = np.zeros((len(flags), len(flags)))
    for (F1, F2), c in block.items():
        M[index[F1], index[F2]] += c
    return M


def invariants(M):
    """indexing-free signature: (entry sum, nonzero row sums, nonzero eigenvalues), rounded."""
    Ms = (M + M.T) / 2
    eig = np.linalg.eigvalsh(Ms)
    rows = M.sum(axis=1)
    return (round(float(M.sum()), 6),
            tuple(np.round(np.sort(rows[np.abs(rows) > 1e-9]), 6).tolist()),
            tuple(np.round(np.sort(eig[np.abs(eig) > 1e-7]), 6).tolist()))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuilt", type=Path, required=True)
    parser.add_argument("--work", type=Path, default=Path(".work"))
    parser.add_argument("--levels", default="0246")
    parser.add_argument("--samples", type=int, default=12)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()
    started = time.time()
    with args.rebuilt.open("rb") as handle:
        states10 = pickle.load(handle)["states10"]
    blocks = Blocks10(states10, args.work, args.levels, verbose=False)
    rng = np.random.default_rng(args.seed)
    sample = sorted(set(rng.choice(len(states10), args.samples, replace=False).tolist()) | {0, len(states10) - 1})
    checked = 0
    for st in sample:
        adj = states10[st]
        q = np.zeros(len(states10))
        q[st] = 1.0
        c_blocks = {key: M for key, M in blocks.matrices(q)}
        for s in [int(c) for c in args.levels]:
            py = python_blocks(adj, s)
            norm = blocks.normaliser[s]
            # C blocks of this level, nonzero ones, keyed by type index; Python keyed by canonical adjacency
            c_level = [(key, M * norm) for key, M in c_blocks.items() if key[0] == s and M.any()]
            py_level = [dense(block) for block in py.values()]
            if len(c_level) != len(py_level):
                raise SystemExit(f"state {st} s={s}: {len(c_level)} nonzero C blocks vs {len(py_level)} Python blocks")
            # C blocks carry every flag of the type, Python only those that occur; zero rows and
            # zero eigenvalues are dropped from the signatures, so the multisets must coincide
            c_inv = sorted(invariants(M) for _k, M in c_level)
            p_inv = sorted(invariants(M) for M in py_level)
            for c_sig, p_sig in zip(c_inv, p_inv):
                if c_sig != p_sig:
                    raise SystemExit(f"MISMATCH at state {st} s={s}: {c_sig[0]} vs {p_sig[0]}, "
                                     f"{len(c_sig[1])} vs {len(p_sig[1])} rows")
                checked += 1
        print(f"state {st:5d}: all levels agree  [{time.time() - started:.0f}s]", flush=True)
    print(f"\n{checked} nonzero blocks compared over {len(sample)} states, signatures rounded to 1e-6")
    print("VERDICT: the C helper's order-10 blocks agree with the independent Python enumeration")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
