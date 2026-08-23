#!/usr/bin/env python3
"""Exact order-10 state vector of a blow-up graphon.

The certificate's LP lives on the 12,172 triangle-free order-10 states.  To
audit whether a row is a *valid* inequality one needs the state vector `q` of a
concrete triangle-free graphon, computed exactly rather than sampled.

For a base graph `H` on `k` vertices and rational class weights `w`, the
associated step graphon `W` samples ten independent points; the resulting
order-10 induced subgraph is the blow-up `H[a]` where `a` is the multiset of
class multiplicities.  Hence

    q[j] = sum over compositions a of 10 into k parts with H[a] ~ state j
           of  multinomial(10; a) * prod_i w_i^{a_i},

which is exact rational arithmetic over at most C(19, k-1) compositions.

Identifying `H[a]` with its index in the order-10 catalogue does not require
that catalogue: the 9->10 deletion lift `c5lift_cache.npz` stores, for every
order-10 state, the multiset of its ten vertex-deleted order-9 subgraphs, and
those profiles are pairwise distinct (checked here).  Matching the profile of
`H[a]` against the lift therefore pins the index, using only the order-9
catalogue from `cache_n9.pkl`.

Dependencies are exactly those of the repository's own replay (numpy, scipy);
the small isomorphism test is included rather than pulled from networkx.
"""

from __future__ import annotations

import argparse
import pickle
from fractions import Fraction
from math import factorial
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

N_ORDER9 = 1897
N_ORDER10 = 12172


# --------------------------------------------------------------------------
# tiny graph utilities: a graph is a list of adjacency bitmasks
# --------------------------------------------------------------------------
def degrees(adjacency):
    return [bin(mask).count("1") for mask in adjacency]


def invariant(adjacency):
    """Cheap isomorphism-invariant used only to bucket candidates."""
    deg = degrees(adjacency)
    edges = sum(deg) // 2
    neighbour_degrees = sorted(
        tuple(sorted(deg[w] for w in range(len(adjacency)) if (mask >> w) & 1))
        for mask in adjacency
    )
    return (len(adjacency), edges, tuple(sorted(deg)), tuple(neighbour_degrees))


def isomorphic(a, b):
    """Backtracking isomorphism test; exact, and fast enough at this size."""
    n = len(a)
    if n != len(b):
        return False
    deg_a, deg_b = degrees(a), degrees(b)
    if sorted(deg_a) != sorted(deg_b):
        return False
    order = sorted(range(n), key=lambda v: -deg_a[v])
    image = [-1] * n
    used = 0

    def extend(depth):
        nonlocal used
        if depth == n:
            return True
        u = order[depth]
        for w in range(n):
            if (used >> w) & 1 or deg_b[w] != deg_a[u]:
                continue
            if all(
                ((a[u] >> order[j]) & 1) == ((b[w] >> image[order[j]]) & 1)
                for j in range(depth)
            ):
                image[u] = w
                used |= 1 << w
                if extend(depth + 1):
                    return True
                used &= ~(1 << w)
                image[u] = -1
        return False

    return extend(0)


def delete_vertex(adjacency, v):
    n = len(adjacency)
    keep = [u for u in range(n) if u != v]
    position = {u: i for i, u in enumerate(keep)}
    out = [0] * (n - 1)
    for u in keep:
        for w in keep:
            if (adjacency[u] >> w) & 1:
                out[position[u]] |= 1 << position[w]
    return out


def blow_up(base_edges, multiplicities):
    """Blow up a base graph (given by its edge list) by the given class sizes."""
    blocks, index = [], 0
    for size in multiplicities:
        blocks.append(list(range(index, index + size)))
        index += size
    adjacency = [0] * index
    for u, v in base_edges:
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


# --------------------------------------------------------------------------
class Catalogue:
    """Order-9 catalogue plus the 9->10 deletion-profile index."""

    def __init__(self, flagsdp: Path):
        with (flagsdp / "cache_n9.pkl").open("rb") as handle:
            cache = pickle.load(handle)
        self.states9 = [list(masks) for _order, masks in cache["states"]]
        if len(self.states9) != N_ORDER9:
            raise RuntimeError("unexpected order-9 catalogue size")
        self.dedge9 = np.asarray(cache["dedge"], dtype=float)
        self.buckets9 = {}
        for index, adjacency in enumerate(self.states9):
            self.buckets9.setdefault(invariant(adjacency), []).append(index)

        lift = np.load(flagsdp / "c5lift_cache.npz", allow_pickle=True)
        if int(lift["nJ"]) != N_ORDER10:
            raise RuntimeError("unexpected order-10 catalogue size")
        counts = np.rint(np.asarray(lift["Dval"]) * 10).astype(np.int64)
        if not np.array_equal(counts.astype(float) / 10.0, np.asarray(lift["Dval"])):
            raise RuntimeError("deletion lift is not integral over denominator ten")
        self.deletion = csr_matrix(
            (counts, (lift["Drow"], lift["Dcol"])),
            shape=(N_ORDER9, N_ORDER10),
            dtype=np.int64,
        )
        columns = self.deletion.tocsc()
        self.by_profile = {}
        for j in range(N_ORDER10):
            start, end = columns.indptr[j], columns.indptr[j + 1]
            key = tuple(
                sorted(zip(columns.indices[start:end].tolist(), columns.data[start:end].tolist()))
            )
            if key in self.by_profile:
                raise RuntimeError("deletion profiles do not separate order-10 states")
            self.by_profile[key] = j

    def index9(self, adjacency):
        for candidate in self.buckets9.get(invariant(adjacency), ()):
            if isomorphic(self.states9[candidate], adjacency):
                return candidate
        raise KeyError("order-9 graph is absent from the catalogue")

    def index10(self, adjacency):
        profile = {}
        for v in range(len(adjacency)):
            i = self.index9(delete_vertex(adjacency, v))
            profile[i] = profile.get(i, 0) + 1
        return self.by_profile[tuple(sorted(profile.items()))]

    def order10_edge_counts(self):
        """Exact edge count of every order-10 state, via the deletion lift."""
        edge9 = np.rint(self.dedge9 * 36).astype(np.int64)
        total = np.asarray(self.deletion.T @ edge9).ravel()
        if np.any(total % 8):
            raise RuntimeError("order-10 deletion edge sum is not divisible by eight")
        return total // 8


def blow_up_state_vector(catalogue: Catalogue, base_edges, n_classes, weights):
    """Exact order-10 state vector of the weighted blow-up graphon of a base graph."""
    if len(weights) != n_classes or sum(weights) != 1:
        raise ValueError("weights must be a rational probability vector over the classes")
    q = [Fraction(0)] * N_ORDER10
    memo = {}
    for parts in compositions(10, n_classes):
        adjacency = blow_up(base_edges, parts)
        key = invariant(adjacency)
        seen = memo.setdefault(key, [])
        # The invariant only buckets; membership is settled by an exact
        # isomorphism test, so distinct graphs never share a cached index.
        for representative, cached in seen:
            if isomorphic(representative, adjacency):
                state = cached
                break
        else:
            state = catalogue.index10(adjacency)
            seen.append((adjacency, state))
        weight = Fraction(multinomial(parts))
        for part, w in zip(parts, weights):
            weight *= w ** part
        q[state] += weight
    if sum(q) != 1:
        raise RuntimeError("state vector does not sum to one")
    return q


def blow_up_densities(base_edges, n_classes, weights):
    """Exact d_edge and d_mono of a weighted blow-up graphon.

    d_mono is the minimum over *measurable* two-colourings, but the mono-edge
    mass is multilinear in the fraction of each class placed on one side, so the
    minimum is attained at a class-constant colouring (this is exactly the
    argument behind the blow-up identity bip(G[t]) = t^2 bip(G)).  Enumerating
    the 2^k class-constant colourings is therefore exact.
    """
    d_edge = 2 * sum(weights[u] * weights[v] for u, v in base_edges)
    best = None
    for pattern in range(1 << n_classes):
        mono = sum(
            weights[u] * weights[v]
            for u, v in base_edges
            if ((pattern >> u) & 1) == ((pattern >> v) & 1)
        )
        if best is None or mono < best:
            best = mono
    return d_edge, 2 * best


CYCLE5 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
PETERSEN = [
    (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
    (5, 7), (7, 9), (9, 6), (6, 8), (8, 5),
    (0, 5), (1, 6), (2, 7), (3, 8), (4, 9),
]
BASES = {"c5": (CYCLE5, 5), "petersen": (PETERSEN, 10)}


def parse_weights(text, n_classes):
    parts = [Fraction(chunk) for chunk in text.split(",")] if text else [Fraction(1)] * n_classes
    if len(parts) != n_classes:
        raise ValueError(f"expected {n_classes} weights")
    total = sum(parts)
    return [part / total for part in parts]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--base", choices=sorted(BASES), default="c5")
    parser.add_argument("--weights", default="", help="comma separated rationals, e.g. 6,1,1,1,1")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    base_edges, n_classes = BASES[args.base]
    weights = parse_weights(args.weights, n_classes)
    catalogue = Catalogue(args.flagsdp.resolve())
    q = blow_up_state_vector(catalogue, base_edges, n_classes, weights)

    edges10 = catalogue.order10_edge_counts()
    support = [j for j, value in enumerate(q) if value]
    sampled_d_edge = sum(Fraction(int(edges10[j]), 45) * q[j] for j in support)
    d_edge, d_mono = blow_up_densities(base_edges, n_classes, weights)
    if sampled_d_edge != d_edge:
        raise RuntimeError("state vector and closed-form edge density disagree")
    print(f"base={args.base} weights={[str(w) for w in weights]}")
    print(f"support={len(support)} states, sum={sum(q)}")
    print(f"d_edge = 2|E|/N^2 = {d_edge} = {float(d_edge)}  (matches the state vector)")
    print(f"d_mono = 2 bip/N^2 = {d_mono} = {float(d_mono)}")
    with args.out.open("wb") as handle:
        pickle.dump(
            {
                "base": args.base,
                "weights": [(w.numerator, w.denominator) for w in weights],
                "support": support,
                "values": [(q[j].numerator, q[j].denominator) for j in support],
                "d_edge": (d_edge.numerator, d_edge.denominator),
                "d_mono": (d_mono.numerator, d_mono.denominator),
            },
            handle,
        )
    print(f"saved {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
