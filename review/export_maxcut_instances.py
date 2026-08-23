#!/usr/bin/env python3
"""Export the open MaxCut instances behind Step 0 in rudy format.

`envelope_maxcut_bound.py` reduces "is the 7-root envelope provably above 2/25
inside the band?" to a finite list of weighted MaxCut problems, one connected
component per root type.  Components up to 20 vertices are settled here by
enumeration; the rest are only bracketed, because a plain Goemans-Williamson
style SDP dual is 4-6% loose and the decision needs about 1%.

This script writes the unsettled components to disk in the `rudy` format that
BiqMac, BiqCrunch and most branch-and-bound MaxCut codes read:

    n m
    i j w          (1-indexed, one line per edge)

Weights are scaled to integers by a per-instance factor recorded in
`manifest.json`, together with the local-search lower bound and the SDP upper
bound, so a solver's answer can be dropped straight back into the accounting:

    envelope = sum over roots of ( diag + S - MaxCut )
    leg7     = ( envelope - (2/25) * root_mass ) / 181440

`manifest.json` also records `diag`, `S`, `root_mass` and the budget, so the
verdict follows from the exact MaxCut values alone.
"""

from __future__ import annotations

import argparse
import json
import pickle
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from envelope_maxcut_bound import (  # noqa: E402
    BASES, best_cut, build_q9, densities, exact_cut, sdp_cut_bound,
)
from k7_leg_ceiling import invariant, isomorphic, triangles  # noqa: E402

DEN = 181440
SCALE = 10 ** 12


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--k7-cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--base", choices=sorted(BASES), default="grotzsch")
    parser.add_argument("--min-size", type=int, default=21,
                        help="export components at least this large")
    args = parser.parse_args()

    base = BASES[args.base]()
    if triangles(base):
        raise RuntimeError("base is not triangle-free")
    w = [Fraction(n, 100000) for n in
         (10173, 9717, 8691, 10166, 7628, 8213, 7344, 9352, 12084, 10443, 6189)] \
        if args.base == "grotzsch" else [Fraction(1, len(base))] * len(base)
    total = sum(w)
    w = [x / total for x in w]
    d_edge, d_mono = densities(base, w)

    with (args.flagsdp / "cache_n9.pkl").open("rb") as handle:
        states = [list(masks) for _o, masks in pickle.load(handle)["states"]]
    buckets = {}
    for i, adjacency in enumerate(states):
        buckets.setdefault(invariant(adjacency), []).append(i)
    q = build_q9(states, buckets, base, w)

    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20260823)
    entries = []
    settled = 0.0          # diag + S - MaxCut, over everything already exact
    mass = 0.0
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
        upper = np.triu(N, 1)
        W = upper + upper.T
        settled += float(np.trace(N)) + float(upper.sum())

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
            if len(comp) < args.min_size:
                settled -= exact_cut(sub)
                continue
            lower = best_cut(sub, 60, rng)
            bound, certified = sdp_cut_bound(sub, 600, seed=root * 131 + len(comp))
            name = f"root{root:03d}_n{len(comp):03d}.rudy"
            iu, ju = np.triu_indices(len(comp), 1)
            weights = sub[iu, ju]
            mask = weights > 0
            iu, ju, weights = iu[mask], ju[mask], weights[mask]
            scaled = np.rint(weights * SCALE).astype(np.int64)
            with (args.out / name).open("w") as handle:
                handle.write(f"{len(comp)} {len(scaled)}\n")
                for i, j, x in zip(iu, ju, scaled):
                    handle.write(f"{i + 1} {j + 1} {int(x)}\n")
            entries.append({
                "file": name, "root": root, "vertices": len(comp), "edges": int(len(scaled)),
                "weight_scale": SCALE,
                "local_search_lower_bound": lower,
                "sdp_dual_upper_bound": bound,
                "sdp_certified": bool(certified),
            })
    target = (2.0 / 25.0) * mass
    manifest = {
        "base": args.base,
        "weights": [f"{x.numerator}/{x.denominator}" for x in w],
        "d_edge": float(d_edge), "d_mono": float(d_mono),
        "in_band": bool(Fraction(2486, 10000) <= d_edge <= Fraction(3197, 10000)),
        "settled_part": settled,
        "root_mass_term": target,
        "denominator": DEN,
        "note": ("leg7 = (settled_part - sum of the exported MaxCut values - root_mass_term)"
                 " / denominator; leg7 > 0 kills the 7-root envelope route"),
        "leg7_if_local_search_is_exact":
            (settled - sum(e["local_search_lower_bound"] for e in entries) - target) / DEN,
        "leg7_lower_bound_from_sdp":
            (settled - sum(e["sdp_dual_upper_bound"] for e in entries) - target) / DEN,
        "instances": sorted(entries, key=lambda e: -e["vertices"]),
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"exported {len(entries)} instances to {args.out}  [{time.time() - started:.0f}s]")
    print(f"largest: {max(e['vertices'] for e in entries)} vertices")
    print(f"leg7 if the local search is exact : "
          f"{manifest['leg7_if_local_search_is_exact']:+.6e}")
    print(f"leg7 lower bound from the SDP dual: "
          f"{manifest['leg7_lower_bound_from_sdp']:+.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
