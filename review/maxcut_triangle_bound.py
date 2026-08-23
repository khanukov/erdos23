#!/usr/bin/env python3
"""Upper bounds on weighted MaxCut from the SDP strengthened by triangle inequalities.

`envelope_maxcut_bound.py` stalls because the plain Goemans-Williamson dual is
4-6% loose on these instances and the Step-0 decision needs about 1%.  The
standard fix is to add the triangle inequalities of the cut polytope: for
`X = x x^T` with `x` in `{-1,1}^n` and any triple `i<j<k`,

    X_ij + X_ik + X_jk >= -1,   -X_ij - X_ik + X_jk >= -1,   and the two others.

Dualising them with multipliers `lam >= 0` gives, for `W' = W - sum_t lam_t A_t`,

    min_x x^T W x  >=  min_{X psd, diag 1} <W', X>  -  sum_t lam_t
                   >=  - sum_i mu_i - sum_t lam_t     whenever W' + Diag(mu) >= 0,

so `(lam, mu)` is a checkable certificate for any `lam >= 0`, and
`MaxCut <= S/2 + (sum mu + sum lam) / 4`.  The inner SDP is solved by the mixing
method, `mu` is read off its stationarity condition and repaired by a uniform
shift, and `lam` is driven by projected subgradient ascent on the Lagrangian,
with violated triples separated from the current `X = R R^T` each round.

Every bound returned is verified by an explicit eigenvalue check of
`W' + Diag(mu)`; an unverified bound is reported as such and never used.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np


def mixing(Wp, rounds, rng, r0=None):
    """Block-coordinate minimisation of sum_ij W'_ij <r_i, r_j> over unit rows."""
    n = Wp.shape[0]
    rank = min(n, int(np.ceil(np.sqrt(2 * n))) + 4)
    if r0 is not None and r0.shape == (n, rank):
        r = r0.copy()
    else:
        r = rng.normal(size=(n, rank))
        r /= np.linalg.norm(r, axis=1, keepdims=True)
    order = np.arange(n)
    for _ in range(rounds):
        rng.shuffle(order)
        for i in order:
            g = Wp[i] @ r
            norm = float(np.linalg.norm(g))
            if norm > 1e-14:
                r[i] = -g / norm
    return r


def inner_certificate(Wp, r):
    """mu with W' + Diag(mu) psd, from stationarity, repaired by a uniform shift."""
    mu = np.linalg.norm(Wp @ r, axis=1)
    smallest = float(np.linalg.eigvalsh(Wp + np.diag(mu))[0])
    if smallest < 0:
        mu = mu - smallest * (1.0 + 1e-12)
        smallest = float(np.linalg.eigvalsh(Wp + np.diag(mu))[0])
    tol = 1e-9 * max(1.0, float(np.abs(Wp).max()))
    return mu, smallest >= -tol


SIGNS = ((1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1))


def separate(X, limit, rng):
    """Most violated triangle inequalities at X, sampled when n is large."""
    n = X.shape[0]
    if n <= 60:
        idx = np.triu_indices(n, 1)
        triples = [(i, j, k) for i in range(n) for j in range(i + 1, n) for k in range(j + 1, n)]
    else:
        pool = 40000
        a = rng.integers(0, n, size=pool)
        b = rng.integers(0, n, size=pool)
        c = rng.integers(0, n, size=pool)
        ok = (a != b) & (b != c) & (a != c)
        triples = [tuple(sorted(t)) for t in zip(a[ok].tolist(), b[ok].tolist(), c[ok].tolist())]
        triples = list(dict.fromkeys(triples))
    found = []
    for (i, j, k) in triples:
        xij, xik, xjk = X[i, j], X[i, k], X[j, k]
        for s in SIGNS:
            value = s[0] * xij + s[1] * xik + s[2] * xjk
            if value < -1.0 - 1e-6:
                found.append((value, i, j, k, s))
                break
    found.sort()
    return found[:limit]


def triangle_bound(W, rounds=60, outer=60, cuts_per_round=None, seed=0,
                   incumbent=None, verbose=False):
    """Certified upper bound on MaxCut(W); returns (bound, certified, info).

    `incumbent` is a known cut value (a lower bound on MaxCut).  It fixes the
    Polyak target for the subgradient ascent: the Lagrangian value can never
    exceed `2S - 4 * MaxCut`, so `2S - 4 * incumbent` is a safe over-estimate of
    the optimum and makes the step size self-scaling.
    """
    n = W.shape[0]
    total = float(np.triu(W, 1).sum())
    rng = np.random.default_rng(seed)
    if incumbent is None:
        incumbent = 0.0
    target = 2.0 * total - 4.0 * incumbent          # >= min_x x^T W x
    active = {}                                     # (i,j,k,signs) -> multiplier
    best = np.inf
    best_certified = False
    r = None
    for it in range(outer):
        Wp = W.copy()
        for (i, j, k, s), lam in active.items():
            if lam <= 0:
                continue
            Wp[i, j] -= 0.5 * s[0] * lam; Wp[j, i] -= 0.5 * s[0] * lam
            Wp[i, k] -= 0.5 * s[1] * lam; Wp[k, i] -= 0.5 * s[1] * lam
            Wp[j, k] -= 0.5 * s[2] * lam; Wp[k, j] -= 0.5 * s[2] * lam
        r = mixing(Wp, rounds, rng, r)
        mu, ok = inner_certificate(Wp, r)
        lam_sum = float(sum(v for v in active.values() if v > 0))
        value = -float(mu.sum()) - lam_sum           # lower bound on min_x x^T W x
        bound = total / 2.0 - value / 4.0            # upper bound on MaxCut
        if ok and bound < best:
            best, best_certified = bound, True
        X = r @ r.T
        np.fill_diagonal(X, 1.0)

        limit = cuts_per_round or max(3 * n, 120)
        for (_v, i, j, k, s) in separate(X, limit, rng):
            active.setdefault((i, j, k, s), 0.0)
        if not active:
            break
        keys = list(active)
        grad = {}
        for key in keys:
            i, j, k, s = key
            slack = s[0] * X[i, j] + s[1] * X[i, k] + s[2] * X[j, k] + 1.0
            grad[key] = -slack                       # > 0 exactly when violated
        norm2 = sum(g * g for g in grad.values())
        if norm2 <= 1e-18:
            break
        step = 1.6 * max(0.0, target - value) / norm2
        for key in keys:
            active[key] = max(0.0, active[key] + step * grad[key])
            if active[key] <= 1e-14:
                del active[key]
        if verbose and (it + 1) % 15 == 0:
            print(f"    iter {it + 1:3d}: bound {bound:.6f} (best {best:.6f}), "
                  f"cuts {len(active)}", flush=True)
    return best, best_certified, {"active_cuts": len(active)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--instances", type=Path, required=True,
                        help="directory written by export_maxcut_instances.py")
    parser.add_argument("--outer", type=int, default=45)
    parser.add_argument("--rounds", type=int, default=60)
    parser.add_argument("--only", type=str, default="",
                        help="comma separated file names to process")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    manifest = json.loads((args.instances / "manifest.json").read_text())
    wanted = set(args.only.split(",")) if args.only else None
    results = {}
    improved_total = plain_total = 0.0
    started = time.time()
    for entry in manifest["instances"]:
        if wanted and entry["file"] not in wanted:
            continue
        lines = (args.instances / entry["file"]).read_text().split("\n")
        n, m = (int(x) for x in lines[0].split())
        W = np.zeros((n, n))
        for line in lines[1:]:
            if not line.strip():
                continue
            i, j, w = line.split()
            W[int(i) - 1, int(j) - 1] = W[int(j) - 1, int(i) - 1] = float(w)
        W /= entry["weight_scale"]
        lower = entry["local_search_lower_bound"]
        plain = entry["sdp_dual_upper_bound"]
        bound, ok, info = triangle_bound(W, args.rounds, args.outer, seed=n * 17,
                                         incumbent=lower, verbose=(n >= 60))
        bound = min(bound, plain)
        results[entry["file"]] = {"vertices": n, "lower": lower, "plain": plain,
                                  "triangle": bound, "certified": bool(ok)}
        improved_total += bound - lower
        plain_total += plain - lower
        print(f"{entry['file']:>20} n={n:4d}  lower {lower:12.6f}  plain {plain:12.6f}"
              f"  +triangle {bound:12.6f}  gap {bound - lower:9.6f}"
              f"  ({100 * (1 - (bound - lower) / max(1e-12, plain - lower)):5.1f}% closed)"
              f"  cuts {info['active_cuts']}", flush=True)
    settled = manifest["settled_part"]
    target = manifest["root_mass_term"]
    den = manifest["denominator"]
    leg_lo = (settled - sum(r["triangle"] for r in results.values())
              - sum(e["sdp_dual_upper_bound"] for e in manifest["instances"]
                    if e["file"] not in results) - target) / den
    print(f"\ntotal gap: plain {plain_total:.4f} -> with triangles {improved_total:.4f}"
          f"   [{time.time() - started:.0f}s]")
    print(f"budget to prove leg7 > 0: {manifest['leg7_if_local_search_is_exact'] * den:.4f}")
    print(f"leg7 lower bound now: {leg_lo:+.6e}")
    print("VERDICT: leg7 provably positive." if leg_lo > 0 else
          "VERDICT: still short; remaining gap "
          f"{-leg_lo * den:.4f} envelope units.")
    if args.out:
        args.out.write_text(json.dumps(results, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
