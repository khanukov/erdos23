#!/usr/bin/env python3
"""How negative can the LP's K7 envelope leg get inside the band?

The K7 leg is a constraint of the LP:  eta <= leg7(q)  at every feasible point,
with

    leg7(q) = (1/10) * ( sum_sigma min_c <same_{sigma,c}, q> - (2/25) sum_sigma <mass_sigma, q> ),

the minimum running over the K7 rules the certificate actually carries.  Once
the Horn and K8 rows are repaired, every in-band triangle-free graphon is a
feasible point, so a certificate in the published 7-root framework (all dual
weight on the K7 leg, no K8 leg) obeys

    delta  >=  max{ leg7(q(W)) : W triangle-free, d_edge(W) in [LO, HI] }.

This script lower-bounds that maximum by hill-climbing the class weights of
blow-up graphons.  Only the order-9 state vector is needed, because contracting
the 9->10 deletion lift against q10 is the same as contracting q9 with a factor
of ten.

Blow-up state vectors are identified against the order-9 catalogue through the
twin-reduced weighted quotient, which is a complete isomorphism invariant for
blow-ups, so each distinct 9-vertex graph costs one isomorphism test globally.
"""

from __future__ import annotations

import argparse
import itertools
import json
import pickle
import random
import time
from math import factorial
from pathlib import Path

import numpy as np

LO, HI = 0.2486, 0.3197
DEN = 181440
N9 = 1897


# ------------------------------------------------------------ graph helpers
def degrees(adjacency):
    return [bin(mask).count("1") for mask in adjacency]


def invariant(adjacency):
    deg = degrees(adjacency)
    return (
        len(adjacency),
        sum(deg) // 2,
        tuple(sorted(deg)),
        tuple(sorted(
            tuple(sorted(deg[w] for w in range(len(adjacency)) if (mask >> w) & 1))
            for mask in adjacency
        )),
    )


def isomorphic(a, b):
    n = len(a)
    da, db = degrees(a), degrees(b)
    if sorted(da) != sorted(db):
        return False
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

    return extend(0, 0)


def triangles(adjacency):
    n = len(adjacency)
    return [(i, j, k) for i, j, k in itertools.combinations(range(n), 3)
            if (adjacency[i] >> j) & 1 and (adjacency[j] >> k) & 1 and (adjacency[i] >> k) & 1]


def edges_of(adjacency):
    n = len(adjacency)
    return [(u, v) for u in range(n) for v in range(u + 1, n) if (adjacency[u] >> v) & 1]


def quotient(base, parts):
    """Twin-reduced weighted quotient of base[parts]: a complete blow-up invariant."""
    live = [i for i in range(len(parts)) if parts[i]]
    neighbourhood = {i: frozenset(j for j in live if (base[i] >> j) & 1) for i in live}
    groups = {}
    for i in live:
        groups.setdefault(neighbourhood[i], []).append(i)
    reps = sorted(groups.values(), key=lambda g: g[0])
    weight = [sum(parts[i] for i in g) for g in reps]
    size = len(reps)
    adjacency = [0] * size
    for x in range(size):
        for y in range(size):
            if x != y and (base[reps[x][0]] >> reps[y][0]) & 1:
                adjacency[x] |= 1 << y
    return adjacency, weight


def canonical(adjacency, weight):
    size = len(adjacency)
    deg = degrees(adjacency)
    colour = [
        (weight[v], deg[v],
         tuple(sorted((weight[w], deg[w]) for w in range(size) if (adjacency[v] >> w) & 1)))
        for v in range(size)
    ]
    order = sorted(range(size), key=lambda v: colour[v])
    cells, start = [], 0
    for i in range(1, size + 1):
        if i == size or colour[order[i]] != colour[order[start]]:
            cells.append(order[start:i])
            start = i
    best = None
    for choice in itertools.product(*[itertools.permutations(c) for c in cells]):
        perm = [v for cell in choice for v in cell]
        key = (
            tuple(weight[v] for v in perm),
            tuple(1 if (adjacency[perm[x]] >> perm[y]) & 1 else 0
                  for x in range(size) for y in range(x + 1, size)),
        )
        if best is None or key < best:
            best = key
    return best


def blow_up(base, parts):
    blocks, index = [], 0
    for size in parts:
        blocks.append(list(range(index, index + size)))
        index += size
    adjacency = [0] * index
    for u in range(len(parts)):
        for v in range(u + 1, len(parts)):
            if (base[u] >> v) & 1:
                for x in blocks[u]:
                    for y in blocks[v]:
                        adjacency[x] |= 1 << y
                        adjacency[y] |= 1 << x
    return adjacency


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def multinomial(parts):
    value = factorial(sum(parts))
    for part in parts:
        value //= factorial(part)
    return value


def cycle(n):
    adjacency = [0] * n
    for u in range(n):
        v = (u + 1) % n
        adjacency[u] |= 1 << v
        adjacency[v] |= 1 << u
    return adjacency


def petersen():
    adjacency = [0] * 10
    for u, v in [(0,1),(1,2),(2,3),(3,4),(4,0),(5,7),(7,9),(9,6),(6,8),(8,5),
                 (0,5),(1,6),(2,7),(3,8),(4,9)]:
        adjacency[u] |= 1 << v
        adjacency[v] |= 1 << u
    return adjacency


def grotzsch():
    adjacency = [0] * 11
    pairs = [(i, (i + 1) % 5) for i in range(5)]
    for i in range(5):
        pairs += [(5 + i, (i - 1) % 5), (5 + i, (i + 1) % 5), (10, 5 + i)]
    for u, v in pairs:
        adjacency[u] |= 1 << v
        adjacency[v] |= 1 << u
    return adjacency


# --------------------------------------------------------------------- main
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--k7-cache", type=Path, required=True)
    parser.add_argument("--sample", type=int, default=45,
                        help="random 9-vertex triangle-free bases to add to the named ones")
    parser.add_argument("--restarts", type=int, default=4)
    parser.add_argument("--steps", type=int, default=150)
    args = parser.parse_args()

    with (args.flagsdp / "cache_n9.pkl").open("rb") as handle:
        states = [list(masks) for _o, masks in pickle.load(handle)["states"]]
    buckets = {}
    for i, adjacency in enumerate(states):
        buckets.setdefault(invariant(adjacency), []).append(i)

    def index9(adjacency):
        for i in buckets.get(invariant(adjacency), ()):
            if isomorphic(states[i], adjacency):
                return i
        raise KeyError("9-vertex graph not in catalogue")

    started = time.time()
    certificate = json.loads(args.certificate.read_text())
    pool, mass = {}, np.zeros((107, N9))
    for root in range(107):
        with np.load(args.k7_cache / f"root_{root:03d}.npz") as payload:
            edge = np.asarray(payload["edge_raw"], dtype=np.int64)
            mass[root] = np.asarray(payload["mass_raw"], dtype=np.int64)
        rows = []
        for descriptor in certificate["descriptors"]:
            if descriptor["kind"] != "k7" or int(descriptor["root"]) != root:
                continue
            rule = np.asarray(descriptor["rule"], dtype=np.uint8)
            same = np.triu(np.equal.outer(rule, rule)).astype(np.int64)
            rows.append(np.tensordot(edge, same, axes=([1, 2], [0, 1])))
        pool[root] = np.asarray(rows, dtype=np.float64)
    mass_total = mass.sum(axis=0)
    print(f"K7 pool prepared in {time.time() - started:.1f}s; "
          f"rules per root: {min(v.shape[0] for v in pool.values())}"
          f"..{max(v.shape[0] for v in pool.values())}")

    def leg7(q9):
        envelope = sum(float(pool[root].dot(q9).min()) for root in range(107))
        return (10.0 / DEN) * (envelope - (2.0 / 25.0) * float(mass_total.dot(q9))) / 10.0

    key_cache = {}

    def table(base):
        parts_list, index_list, weight_list = [], [], []
        for parts in compositions(9, len(base)):
            key = canonical(*quotient(base, parts))
            hit = key_cache.get(key)
            if hit is None:
                hit = key_cache[key] = index9(blow_up(base, parts))
            parts_list.append(parts)
            index_list.append(hit)
            weight_list.append(multinomial(parts))
        return (np.asarray(parts_list, dtype=np.float64),
                np.asarray(index_list, dtype=np.int64),
                np.asarray(weight_list, dtype=np.float64))

    def q9_of(tab, w):
        parts, index, mult = tab
        values = mult * np.exp(parts @ np.log(np.maximum(np.asarray(w, dtype=float), 1e-300)))
        q = np.zeros(N9)
        np.add.at(q, index, values)
        return q

    def densities(base, w):
        edges = edges_of(base)
        products = [w[u] * w[v] for u, v in edges]
        d_edge = 2 * sum(products)
        d_mono = 2 * min(
            sum(p for (u, v), p in zip(edges, products)
                if ((pattern >> u) & 1) == ((pattern >> v) & 1))
            for pattern in range(1 << len(base))
        )
        return d_edge, d_mono

    def search(base, seed):
        size = len(base)
        tab = table(base)
        rng = random.Random(seed)
        best = (-9.9, None, None)
        for attempt in range(args.restarts):
            if attempt == 0:
                w = np.full(size, 1.0 / size)
            else:
                w = np.asarray([rng.random() + 1e-3 for _ in range(size)])
                w /= w.sum()
            d_edge, _ = densities(base, w)
            value = leg7(q9_of(tab, w)) if LO <= d_edge <= HI else -9.9
            scale = 0.3
            for step in range(args.steps):
                i, j = rng.randrange(size), rng.randrange(size)
                if i == j:
                    continue
                trial = w.copy()
                shift = scale * rng.random() * trial[i]
                trial[i] -= shift
                trial[j] += shift
                d_edge, _ = densities(base, trial)
                if not LO <= d_edge <= HI:
                    continue
                candidate = leg7(q9_of(tab, trial))
                if candidate > value:
                    w, value = trial, candidate
                if step % 30 == 29:
                    scale *= 0.6
            if value > best[0]:
                best = (value, w.copy(), densities(base, w))
        return best

    named = [("C5", cycle(5)), ("C7", cycle(7)), ("C9", cycle(9)),
             ("Petersen", petersen()), ("Grotzsch", grotzsch())]
    for name, adjacency in named:
        bad = triangles(adjacency)
        if bad:
            raise RuntimeError(f"{name} is not triangle-free: {bad[:3]}")
    rng = random.Random(7)
    sample = rng.sample(range(N9), args.sample)
    candidates = named + [(f"tf9#{i}", states[i]) for i in sample if edges_of(states[i])]

    results = []
    started = time.time()
    for n, (name, base) in enumerate(candidates):
        value, w, dens = search(base, hash(name) & 0xFFFF)
        if value > -9:
            results.append((value, name, dens[0], dens[1]))
        if (n + 1) % 5 == 0:
            print(f"  {n + 1}/{len(candidates)} bases, {time.time() - started:.0f}s", flush=True)
    results.sort(reverse=True)

    print(f"\n{'leg7':>14} {'base':>10} {'d_edge':>8} {'d_mono':>8}")
    for value, name, d_edge, d_mono in results[:12]:
        print(f"{value:>14.6e} {name:>10} {d_edge:>8.4f} {d_mono:>8.4f}")
    top = results[0]
    print(f"\nmax leg7 found inside the band : {top[0]:+.6e}   (at {top[1]})")
    print(f"certificate objective delta*   : {-9.878886951679021e-04:+.6e}")
    print(f"=> a repaired certificate in the 7-root framework cannot go below "
          f"{top[0]:+.3e},")
    print(f"   so the claimed delta* is out of reach there by a factor of "
          f"{9.878886951679021e-04 / abs(top[0]):.0f}."
          if top[0] > -9.878886951679021e-04 else
          "   the claimed delta* is not excluded by this scan.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
