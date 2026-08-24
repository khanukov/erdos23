#!/usr/bin/env python3
"""Rebuild the 8-rooted edge decomposition from scratch, trusting nothing.

`u8_decomp.pkl` is broken (INDEPENDENT_REVIEW.md section 9) and Aut-averaging
patches it only partway, so the decisive question — is a *correct* `U8` below
`2/25` inside the band? — cannot be answered from the shipped tables.  This
script builds the decomposition itself:

1. the order-10 catalogue, by extending each of the 1,897 order-9 triangle-free
   graphs with one vertex whose neighbourhood is an independent set, deduplicated
   by isomorphism (must come to exactly 12,172);
2. its identification with the certificate's state numbering, through the
   vertex-deletion profiles of `c5lift_cache.npz`, which separate all 12,172
   states (must be a bijection);
3. the 410 order-8 root types, by the same extension from the 107 order-7 types;
4. for every order-10 state and every edge between two vertices, the root
   `R` = the other eight, its type `tau`, an explicit isomorphism `R -> tau`,
   and the profiles of the two endpoints under it.

Summing over *all* isomorphisms `R -> tau` rather than picking one is what the
shipped data gets wrong; that sum is `|Aut(tau)|` copies of the Aut-orbit
average, so the decomposition is built once with one isomorphism and then
averaged over `Aut(tau)` computed from the type itself.

The acid test is the balanced `C5` blow-up: there `d_mono = U7 = 2/25` and
`d_mono <= U8 <= U7`, so a correct `U8` must be exactly `2/25`.  Everything the
script prints is worthless if that check fails.
"""

from __future__ import annotations

import argparse
import itertools
import pickle
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
from envelope_maxcut_bound import best_cut, exact_cut  # noqa: E402
from k7_leg_ceiling import (  # noqa: E402
    compositions, cycle, edges_of, grotzsch, multinomial, petersen, triangles,
)

N9, N10 = 1897, 12172


def degrees(adjacency):
    return [bin(mask).count("1") for mask in adjacency]


def invariant(adjacency):
    deg = degrees(adjacency)
    return (
        len(adjacency), sum(deg) // 2, tuple(sorted(deg)),
        tuple(sorted(
            tuple(sorted(deg[w] for w in range(len(adjacency)) if (mask >> w) & 1))
            for mask in adjacency)),
    )


def isomorphism(a, b):
    """An explicit isomorphism a -> b as a list, or None."""
    n = len(a)
    da, db = degrees(a), degrees(b)
    if sorted(da) != sorted(db):
        return None
    order = sorted(range(n), key=lambda v: -da[v])
    image = [-1] * n

    def extend(depth, used):
        if depth == n:
            return True
        u = order[depth]
        for w in range(n):
            if (used >> w) & 1 or db[w] != da[u]:
                continue
            if all(((a[u] >> order[j]) & 1) == ((b[w] >> image[order[j]]) & 1)
                   for j in range(depth)):
                image[u] = w
                if extend(depth + 1, used | (1 << w)):
                    return True
                image[u] = -1
        return False

    return list(image) if extend(0, 0) else None


def automorphisms(adjacency):
    n = len(adjacency)
    deg = degrees(adjacency)
    found, image = [], [-1] * n

    def extend(depth, used):
        if depth == n:
            found.append(tuple(image))
            return
        for w in range(n):
            if (used >> w) & 1 or deg[w] != deg[depth]:
                continue
            if all(((adjacency[depth] >> j) & 1) == ((adjacency[w] >> image[j]) & 1)
                   for j in range(depth)):
                image[depth] = w
                extend(depth + 1, used | (1 << w))
                image[depth] = -1

    extend(0, 0)
    return found


class Catalogue:
    """Triangle-free graphs on n vertices, with isomorphism lookup."""

    def __init__(self, graphs):
        self.graphs = graphs
        self.buckets = {}
        for i, adjacency in enumerate(graphs):
            self.buckets.setdefault(invariant(adjacency), []).append(i)

    def find(self, adjacency):
        for i in self.buckets.get(invariant(adjacency), ()):
            mapping = isomorphism(adjacency, self.graphs[i])
            if mapping is not None:
                return i, mapping
        return None, None

    def add(self, adjacency):
        i, _ = self.find(adjacency)
        if i is None:
            i = len(self.graphs)
            self.graphs.append(adjacency)
            self.buckets.setdefault(invariant(adjacency), []).append(i)
        return i


def extend_all(graphs):
    """Every triangle-free graph on n+1 vertices obtained by adding one vertex."""
    catalogue = Catalogue([])
    for adjacency in graphs:
        n = len(adjacency)
        for mask in range(1 << n):
            ok = True
            for i in range(n):
                if not (mask >> i) & 1:
                    continue
                if adjacency[i] & mask:
                    ok = False
                    break
            if not ok:
                continue
            catalogue.add([adjacency[i] | (((mask >> i) & 1) << n) for i in range(n)] + [mask])
    return catalogue.graphs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--starts", type=int, default=25)
    parser.add_argument("--cache", type=Path, default=None,
                        help="where to store/reuse the rebuilt decomposition")
    args = parser.parse_args()

    started = time.time()
    with (args.flagsdp / "cache_n9.pkl").open("rb") as handle:
        states9 = [list(masks) for _o, masks in pickle.load(handle)["states"]]
    cat9 = Catalogue(states9)
    print(f"order-9 catalogue: {len(states9)} graphs")

    lift = np.load(args.flagsdp / "c5lift_cache.npz", allow_pickle=True)
    counts = np.rint(np.asarray(lift["Dval"]) * 10).astype(np.int64)
    columns = csr_matrix((counts, (lift["Drow"], lift["Dcol"])),
                         shape=(N9, N10), dtype=np.int64).tocsc()
    profile_to_state = {}
    for j in range(N10):
        s, e = columns.indptr[j], columns.indptr[j + 1]
        profile_to_state[tuple(sorted(zip(columns.indices[s:e].tolist(),
                                          columns.data[s:e].tolist())))] = j
    if len(profile_to_state) != N10:
        raise RuntimeError("deletion profiles do not separate the order-10 states")

    if args.cache and args.cache.exists():
        with args.cache.open("rb") as handle:
            payload = pickle.load(handle)
        states10, roots8, pairs = payload["states10"], payload["roots8"], payload["pairs"]
        print(f"reloaded rebuilt decomposition from {args.cache}")
    else:
        graphs10 = extend_all(states9)
        print(f"order-10 catalogue rebuilt: {len(graphs10)} graphs "
              f"(expected {N10})  [{time.time() - started:.0f}s]")
        if len(graphs10) != N10:
            raise RuntimeError("order-10 catalogue size mismatch")

        states10 = [None] * N10
        for adjacency in graphs10:
            profile = {}
            for v in range(10):
                keep = [u for u in range(10) if u != v]
                position = {u: i for i, u in enumerate(keep)}
                reduced = [0] * 9
                for u in keep:
                    for w in keep:
                        if (adjacency[u] >> w) & 1:
                            reduced[position[u]] |= 1 << position[w]
                i, _ = cat9.find(reduced)
                profile[i] = profile.get(i, 0) + 1
            state = profile_to_state[tuple(sorted(profile.items()))]
            if states10[state] is not None:
                raise RuntimeError("two rebuilt graphs map to the same state")
            states10[state] = adjacency
        if any(a is None for a in states10):
            raise RuntimeError("some order-10 states were never produced")
        print(f"identified all {N10} states with the certificate numbering "
              f"[{time.time() - started:.0f}s]")

        graphs7 = extend_all(extend_all(extend_all(extend_all(extend_all(extend_all([[0]]))))))
        if len(graphs7) != 107:
            raise RuntimeError(f"order-7 catalogue is {len(graphs7)}, expected 107")
        roots8 = extend_all(graphs7)
        if len(roots8) != 410:
            raise RuntimeError(f"order-8 catalogue is {len(roots8)}, expected 410")
        cat8 = Catalogue(roots8)
        print(f"root types rebuilt: {len(graphs7)} on 7 vertices, {len(roots8)} on 8"
              f"  [{time.time() - started:.0f}s]")

        pairs = [dict() for _ in range(410)]          # root -> {(a,b): {state: count}}
        for state, adjacency in enumerate(states10):
            for u in range(10):
                for v in range(u + 1, 10):
                    if not (adjacency[u] >> v) & 1:
                        continue
                    rest = [x for x in range(10) if x != u and x != v]
                    position = {x: i for i, x in enumerate(rest)}
                    induced = [0] * 8
                    for x in rest:
                        for y in rest:
                            if (adjacency[x] >> y) & 1:
                                induced[position[x]] |= 1 << position[y]
                    tau, mapping = cat8.find(induced)
                    prof_u = tuple(sorted(mapping[position[x]] for x in rest
                                          if (adjacency[u] >> x) & 1))
                    prof_v = tuple(sorted(mapping[position[x]] for x in rest
                                          if (adjacency[v] >> x) & 1))
                    for key in ((prof_u, prof_v), (prof_v, prof_u)):
                        bucket = pairs[tau].setdefault(key, {})
                        bucket[state] = bucket.get(state, 0) + 1
            if (state + 1) % 2000 == 0:
                print(f"  decomposed {state + 1}/{N10} states "
                      f"[{time.time() - started:.0f}s]", flush=True)
        if args.cache:
            with args.cache.open("wb") as handle:
                pickle.dump({"states10": states10, "roots8": roots8, "pairs": pairs}, handle)
            print(f"cached to {args.cache}")

    total_pairs = sum(sum(sum(b.values()) for b in root.values()) for root in pairs)
    print(f"rebuilt decomposition: {total_pairs} ordered (edge, root) incidences")

    # ---- Aut-averaged matrices, per root -------------------------------
    aut = [automorphisms(root) for root in roots8]
    print(f"|Aut| over the 410 types: min {min(len(a) for a in aut)}, "
          f"max {max(len(a) for a in aut)}")

    def u8_at(q, starts, seed=3):
        rng = random.Random(seed)
        total = 0.0
        for tau in range(410):
            entries = pairs[tau]
            if not entries:
                continue
            profiles = sorted({p for key in entries for p in key},
                              key=lambda p: (len(p), p))
            closure = list(profiles)
            index = {p: i for i, p in enumerate(closure)}
            for g in aut[tau]:
                for p in profiles:
                    moved = tuple(sorted(g[x] for x in p))
                    if moved not in index:
                        index[moved] = len(closure)
                        closure.append(moved)
            width = len(closure)
            raw = np.zeros((width, width))
            for (a, b), bucket in entries.items():
                value = sum(count * q[state] for state, count in bucket.items())
                if value:
                    raw[index[a], index[b]] += value
            perms = np.asarray(
                [[index[tuple(sorted(g[x] for x in p))] for p in closure] for g in aut[tau]],
                dtype=np.int64)
            rows, cols = np.nonzero(raw)
            if len(rows):
                flat = (perms[:, rows] * width + perms[:, cols]).ravel()
                spread = np.tile(raw[rows, cols], len(perms))
                matrix = np.bincount(flat, weights=spread,
                                     minlength=width * width).reshape(width, width)
            else:
                matrix = np.zeros((width, width))
            matrix /= len(aut[tau]) * 90.0
            matrix = (matrix + matrix.T) / 2.0
            diag = float(np.trace(matrix))
            upper = np.triu(matrix, 1)
            S = float(upper.sum())
            W = upper + upper.T
            active = np.nonzero(W.sum(axis=1) + np.diag(matrix))[0]
            if len(active) > 1:
                sub = W[np.ix_(active, active)]
                cut = exact_cut(sub) if len(active) <= 18 else best_cut(sub, starts, rng)
            else:
                cut = 0.0
            total += diag + 2.0 * (S - cut)
        return total

    # ---- q10 of a blow-up graphon --------------------------------------
    def q10_of(base, weights):
        h = len(base)
        index9 = {}
        for parts in compositions(9, h):
            blocks, idx = [], 0
            for size in parts:
                blocks.append(list(range(idx, idx + size)))
                idx += size
            adjacency = [0] * idx
            for x in range(h):
                for y in range(x + 1, h):
                    if (base[x] >> y) & 1:
                        for p in blocks[x]:
                            for r in blocks[y]:
                                adjacency[p] |= 1 << r
                                adjacency[r] |= 1 << p
            i, _ = cat9.find(adjacency)
            index9[parts] = i
        q = np.zeros(N10)
        mass = 0.0
        for parts in compositions(10, h):
            profile = {}
            for i in range(h):
                if parts[i]:
                    reduced = list(parts)
                    reduced[i] -= 1
                    j = index9[tuple(reduced)]
                    profile[j] = profile.get(j, 0) + parts[i]
            state = profile_to_state[tuple(sorted(profile.items()))]
            weight = multinomial(parts) * float(np.prod(np.power(weights, parts)))
            q[state] += weight
            mass += weight
        return q / mass

    def densities(base, weights):
        edges = edges_of(base)
        products = [weights[u] * weights[v] for u, v in edges]
        d_edge = 2 * sum(products)
        d_mono = 2 * min(
            sum(p for (u, v), p in zip(edges, products)
                if ((pattern >> u) & 1) == ((pattern >> v) & 1))
            for pattern in range(1 << len(base)))
        return d_edge, d_mono

    print("\n--- acid test: the balanced C5 blow-up, where U8 must be exactly 2/25 ---")
    c5 = cycle(5)
    q = q10_of(c5, np.full(5, 0.2))
    d_edge, d_mono = densities(c5, np.full(5, 0.2))
    value = u8_at(q, args.starts)
    print(f"d_edge = {d_edge:.6f}   d_mono = {d_mono:.6f}")
    print(f"rebuilt U8 = {value:.9f}   target 2/25 = 0.080000000   "
          f"error {value - 0.08:+.3e}")
    if abs(value - 0.08) > 5e-4:
        print("\nREBUILD REJECTED: the acid test fails, nothing below is trustworthy.")
        return 1
    print("rebuild validated.")

    print("\n--- U8 inside the band ---")
    bases = [("Grotzsch", grotzsch(),
              np.asarray([10173, 9717, 8691, 10166, 7628, 8213, 7344, 9352,
                          12084, 10443, 6189], dtype=float)),
             ("Petersen", petersen(), np.full(10, 1.0))]
    print(f"{'graphon':>12} {'d_edge':>9} {'d_mono':>9} {'U8':>11} {'U8-2/25':>11}"
          f"  {'valid':>6} {'useful':>7}")
    worst = -1.0
    for name, base, weights in bases:
        if triangles(base):
            raise RuntimeError(f"{name} is not triangle-free")
        weights = weights / weights.sum()
        d_edge, d_mono = densities(base, weights)
        q = q10_of(base, weights)
        value = u8_at(q, args.starts)
        worst = max(worst, value)
        print(f"{name:>12} {d_edge:>9.6f} {d_mono:>9.6f} {value:>11.6f} "
              f"{value - 0.08:>+11.6f}  {str(value >= d_mono - 1e-9):>6} "
              f"{str(value < 0.08):>7}")
    print(f"\nlargest U8 seen in the band: {worst:.6f}   2/25 = 0.080000")
    print("The 8-root envelope clears 2/25 at these points."
          if worst < 0.08 else
          "The 8-root envelope does NOT clear 2/25 here; the whole envelope family fails.")
    print(f"\ntotal {time.time() - started:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
