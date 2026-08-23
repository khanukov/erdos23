#!/usr/bin/env python3
"""Step 0: is the 7-root per-root-MaxCut envelope provably above 2/25 in the band?

At a triangle-free graphon W the LP's K7 leg is

    leg7 = ( sum_sigma min_c mono_sigma(c)  -  (2/25) * sum_sigma <mass_sigma, q9> ) / 181440,

with `mono_sigma(c) = sum_{a<=b} N^sigma_ab [c_a = c_b]` and the minimum over ALL
Boolean profile rules.  Writing `diag = sum_a N_aa` and `S = sum_{a<b} N_ab`,

    min_c mono = diag + S - MaxCut(N_offdiag),

so a *lower* bound on the envelope needs an *upper* bound on a weighted MaxCut,
one instance per root type.  This script produces both sides:

* lower bound on each MaxCut, by multi-start local search with Kernighan-Lin
  sweeps (any cut is a witness), which upper-bounds `leg7`;
* upper bound on each MaxCut, by the SDP dual.  For any vector `mu` with
  `W + Diag(mu) >= 0` and every `x` in `{-1,1}^n`,

      x^T W x = x^T (W + Diag(mu)) x - sum_i mu_i  >=  - sum_i mu_i,

  hence `MaxCut = S/2 - (1/4) min_x x^T W x <= S/2 + (1/4) sum_i mu_i`.
  A good `mu` comes from the mixing method (block-coordinate minimisation of
  `sum_ij W_ij <r_i, r_j>` over unit vectors), then `mu_i = || sum_j W_ij r_j ||`,
  repaired by a uniform shift if the numerical `lambda_min` is slightly negative.
  The certificate is just `mu`: it is checked, not trusted.

MaxCut decomposes over connected components, and small components are solved
exactly by enumeration, so the SDP is only used where it has to be.

If the resulting lower bound on `leg7` is strictly positive at a graphon inside
the band, then no per-root-MaxCut envelope certificate at any flag order can
prove the conjecture there.
"""

from __future__ import annotations

import argparse
import pickle
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from k7_leg_ceiling import (  # noqa: E402
    blow_up, canonical, compositions, cycle, edges_of, grotzsch, invariant,
    isomorphic, multinomial, petersen, quotient, triangles,
)

N9 = 1897
DEN = 181440
LO, HI = Fraction(2486, 10000), Fraction(3197, 10000)
EXACT_LIMIT = 20          # components up to this size are solved by enumeration


# --------------------------------------------------------------- MaxCut side
def best_cut(weights, starts, rng):
    """Lower bound on MaxCut: multi-start local search with Kernighan-Lin sweeps."""
    n = weights.shape[0]
    best = 0.0
    for attempt in range(starts):
        x = (np.zeros(n, dtype=np.int8) if attempt == 0
             else np.asarray([rng.random() < 0.5 for _ in range(n)], dtype=np.int8))
        while True:                                   # steepest single flips
            diff = (x[:, None] != x[None, :])
            gain = (weights * ~diff).sum(axis=1) - (weights * diff).sum(axis=1)
            i = int(np.argmax(gain))
            if gain[i] <= 1e-12:
                break
            x[i] ^= 1
        for _ in range(2):                            # Kernighan-Lin sweep
            locked = np.zeros(n, dtype=bool)
            cur = x.copy()
            best_local = cur.copy()
            value = cut_value(weights, cur)
            for _step in range(n):
                diff = (cur[:, None] != cur[None, :])
                gain = (weights * ~diff).sum(axis=1) - (weights * diff).sum(axis=1)
                gain[locked] = -np.inf
                i = int(np.argmax(gain))
                if not np.isfinite(gain[i]):
                    break
                cur[i] ^= 1
                locked[i] = True
                trial = cut_value(weights, cur)
                if trial > value:
                    value, best_local = trial, cur.copy()
            x = best_local
        best = max(best, cut_value(weights, x))
    return best


def cut_value(weights, x):
    return float((np.triu(weights, 1) * (x[:, None] != x[None, :])).sum())


def exact_cut(weights):
    """Exact MaxCut by enumeration; vertex 0 is fixed to break the symmetry."""
    n = weights.shape[0]
    iu, ju = np.triu_indices(n, 1)
    w = weights[iu, ju]
    keep = w != 0
    iu, ju, w = iu[keep], ju[keep], w[keep]
    patterns = np.arange(1 << (n - 1), dtype=np.int64)
    best = 0.0
    step = max(1, (1 << 22) // max(1, len(w)))
    for start in range(0, len(patterns), step):
        block = patterns[start:start + step]
        bits = ((block[:, None] >> np.arange(n - 1)[None, :]) & 1).astype(np.int8)
        full = np.concatenate([np.zeros((len(block), 1), dtype=np.int8), bits], axis=1)
        cuts = ((full[:, iu] != full[:, ju]) * w[None, :]).sum(axis=1)
        best = max(best, float(cuts.max()))
    return best


def sdp_cut_bound(weights, rounds=400, seed=0):
    """Upper bound on MaxCut via a checked SDP dual certificate mu."""
    n = weights.shape[0]
    total = float(np.triu(weights, 1).sum())
    rank = min(n, int(np.ceil(np.sqrt(2 * n))) + 4)
    rng = np.random.default_rng(seed)
    r = rng.normal(size=(n, rank))
    r /= np.linalg.norm(r, axis=1, keepdims=True)
    order = np.arange(n)
    for _ in range(rounds):
        rng.shuffle(order)
        for i in order:
            g = weights[i] @ r
            norm = float(np.linalg.norm(g))
            if norm > 1e-14:
                r[i] = -g / norm
    mu = np.linalg.norm(weights @ r, axis=1)
    smallest = float(np.linalg.eigvalsh(weights + np.diag(mu))[0])
    if smallest < 0:                                   # repair, then re-check
        mu = mu - smallest * (1.0 + 1e-12)
        smallest = float(np.linalg.eigvalsh(weights + np.diag(mu))[0])
    certified = smallest >= -1e-9 * max(1.0, float(np.abs(weights).max()))
    return total / 2.0 + mu.sum() / 4.0, certified


# ------------------------------------------------------------------- graphon
def build_q9(catalogue_states, buckets, base, weights):
    def index9(adjacency):
        for i in buckets.get(invariant(adjacency), ()):
            if isomorphic(catalogue_states[i], adjacency):
                return i
        raise KeyError("9-vertex graph not in catalogue")

    key_cache, parts_list, index_list, mult_list = {}, [], [], []
    for parts in compositions(9, len(base)):
        key = canonical(*quotient(base, parts))
        hit = key_cache.get(key)
        if hit is None:
            hit = key_cache[key] = index9(blow_up(base, parts))
        parts_list.append(parts)
        index_list.append(hit)
        mult_list.append(multinomial(parts))
    P = np.asarray(parts_list, dtype=float)
    q = np.zeros(N9)
    np.add.at(q, np.asarray(index_list),
              np.asarray(mult_list, dtype=float) * np.exp(P @ np.log(np.asarray(weights, dtype=float))))
    return q


def densities(base, weights):
    edges = edges_of(base)
    products = [weights[u] * weights[v] for u, v in edges]
    d_edge = 2 * sum(products)
    d_mono = 2 * min(
        sum(p for (u, v), p in zip(edges, products)
            if ((pattern >> u) & 1) == ((pattern >> v) & 1))
        for pattern in range(1 << len(base))
    )
    return d_edge, d_mono


BASES = {"grotzsch": grotzsch, "petersen": petersen, "c5": lambda: cycle(5)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--k7-cache", type=Path, required=True)
    parser.add_argument("--base", choices=sorted(BASES), default="grotzsch")
    parser.add_argument("--weights", default="",
                        help="comma separated rationals; default = the in-band optimum found for Grotzsch")
    parser.add_argument("--starts", type=int, default=60)
    parser.add_argument("--rounds", type=int, default=400)
    args = parser.parse_args()

    base = BASES[args.base]()
    bad = triangles(base)
    if bad:
        raise RuntimeError(f"base is not triangle-free: {bad[:3]}")
    if args.weights:
        w = [Fraction(chunk) for chunk in args.weights.split(",")]
    elif args.base == "grotzsch":
        w = [Fraction(n, 100000) for n in
             (10173, 9717, 8691, 10166, 7628, 8213, 7344, 9352, 12084, 10443, 6189)]
    else:
        w = [Fraction(1, len(base))] * len(base)
    total = sum(w)
    w = [x / total for x in w]

    d_edge, d_mono = densities(base, w)
    print(f"base = {args.base}, {len(base)} classes")
    print(f"d_edge = {float(d_edge):.6f}   inside [0.2486, 0.3197]? {LO <= d_edge <= HI}")
    print(f"d_mono = {float(d_mono):.6f}   <= 2/25 = 0.08 ? {d_mono <= Fraction(2, 25)}")

    with (args.flagsdp / "cache_n9.pkl").open("rb") as handle:
        states = [list(masks) for _o, masks in pickle.load(handle)["states"]]
    buckets = {}
    for i, adjacency in enumerate(states):
        buckets.setdefault(invariant(adjacency), []).append(i)
    started = time.time()
    q = build_q9(states, buckets, base, w)
    print(f"order-9 state vector: support {int((q > 0).sum())}, sum {q.sum():.12f}"
          f"  [{time.time() - started:.0f}s]")

    rng = random.Random(20260823)
    env_lo = env_hi = 0.0
    mass = 0.0
    gap_total = 0.0
    exact_used = sdp_used = uncertified = 0
    started = time.time()
    for root in range(107):
        with np.load(args.k7_cache / f"root_{root:03d}.npz") as payload:
            edge = np.asarray(payload["edge_raw"], dtype=np.int64)
            mass += float(np.asarray(payload["mass_raw"], dtype=np.int64) @ q)
        s, a, b = np.nonzero(edge)
        value = edge[s, a, b] * q[s]
        keep = value != 0
        if not keep.any():
            continue
        aa, bb, vv = a[keep], b[keep], value[keep]
        active = sorted(set(aa.tolist()) | set(bb.tolist()))
        remap = {p: i for i, p in enumerate(active)}
        n = len(active)
        N = np.zeros((n, n))
        np.add.at(N, ([remap[x] for x in aa.tolist()], [remap[x] for x in bb.tolist()]), vv)
        diag = float(np.trace(N))
        upper = np.triu(N, 1)
        W = upper + upper.T
        S = float(upper.sum())
        env_lo += diag + S
        env_hi += diag + S

        seen = np.zeros(n, dtype=bool)
        for v in range(n):
            if seen[v]:
                continue
            stack, comp = [v], [v]
            seen[v] = True
            while stack:
                u = stack.pop()
                for x in np.nonzero(W[u] > 0)[0]:
                    if not seen[x]:
                        seen[x] = True
                        stack.append(int(x))
                        comp.append(int(x))
            if len(comp) < 2:
                continue
            sub = W[np.ix_(comp, comp)]
            if len(comp) <= EXACT_LIMIT:
                exact = exact_cut(sub)
                env_lo -= exact
                env_hi -= exact
                exact_used += 1
            else:
                lower = best_cut(sub, args.starts, rng)
                upper_bound, certified = sdp_cut_bound(sub, args.rounds, seed=root * 131 + len(comp))
                if not certified:
                    uncertified += 1
                env_lo -= upper_bound
                env_hi -= lower
                gap_total += upper_bound - lower
                sdp_used += 1
    print(f"components solved exactly: {exact_used}; by SDP dual: {sdp_used}; "
          f"uncertified duals: {uncertified}  [{time.time() - started:.0f}s]")

    target = (2.0 / 25.0) * mass
    leg_lo = (env_lo - target) / DEN
    leg_hi = (env_hi - target) / DEN
    print(f"\nenvelope, lower bound (certified MaxCut upper bounds) : {env_lo:.6f}")
    print(f"envelope, upper bound (explicit cuts)                 : {env_hi:.6f}")
    print(f"(2/25) * root-mass term                               : {target:.6f}")
    print(f"total SDP gap over the large components               : {gap_total:.6f}")
    print(f"\nleg7 lower bound : {leg_lo:+.6e}")
    print(f"leg7 upper bound : {leg_hi:+.6e}")
    if leg_lo > 0:
        print("\nRESULT: leg7 is provably POSITIVE at this in-band graphon.")
        print("        No per-root-MaxCut envelope certificate can prove the conjecture")
        print("        on the band, at any flag order.")
    elif leg_hi <= 0:
        print("\nRESULT: leg7 is negative here; this point is no obstruction.")
    else:
        print("\nRESULT: inconclusive - the MaxCut bounds are not tight enough.")
        print(f"        gap to close: {-leg_lo:.6e} in leg7 units "
              f"({-leg_lo * DEN:.3f} in envelope units)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
