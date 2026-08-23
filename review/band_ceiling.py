#!/usr/bin/env python3
"""How negative can a *valid* band certificate's objective possibly be?

Whatever rows a certificate uses, "the LP is a relaxation" means that at every
genuine triangle-free graphon W the point (q(W), eta = d_mono(W) - 2/25, ...)
is feasible.  The dual objective bounds eta from above at every feasible point,
so for every W in the band

    delta  >=  d_mono(W) - 2/25,

and therefore

    delta  >=  sup{ d_mono(W) : W triangle-free, d_edge(W) in [LO, HI] } - 2/25.

That is a ceiling no amount of extra rows, higher flag order or better Gram
blocks can get past.  This script lower-bounds the supremum by searching a wide
family of triangle-free graphons:

* weighted blow-ups of every triangle-free graph on 9 vertices (all 1,897 of
  them, read straight out of `cache_n9.pkl`), plus C5, C7, C9, C11, the
  Petersen graph and the Grotzsch graph;
* dilutions `theta * W`, which stay triangle-free (t(K3, theta W) = theta^3
  t(K3, W) = 0) and scale both densities, so a graphon above the band can be
  pulled back into it along its own ray.

For a blow-up, both densities have closed forms: `d_edge = 2 sum_{uv in E} w_u w_v`
and `d_mono = 2 min over the 2^h class-constant colourings`, the minimum over
*measurable* colourings being attained at a class-constant one by multilinearity.
No flag data is needed, so the scan is cheap.

Every base graph is asserted triangle-free before use.
"""

from __future__ import annotations

import argparse
import itertools
import pickle
import random
from fractions import Fraction
from pathlib import Path

import numpy as np

LO = Fraction(2486, 10000)
HI = Fraction(3197, 10000)
LO_F, HI_F = float(LO), float(HI)


def triangles(adjacency):
    n = len(adjacency)
    return [
        (i, j, k)
        for i, j, k in itertools.combinations(range(n), 3)
        if (adjacency[i] >> j) & 1 and (adjacency[j] >> k) & 1 and (adjacency[i] >> k) & 1
    ]


def edges_of(adjacency):
    n = len(adjacency)
    return [(u, v) for u in range(n) for v in range(u + 1, n) if (adjacency[u] >> v) & 1]


def cycle(n):
    adjacency = [0] * n
    for u in range(n):
        v = (u + 1) % n
        adjacency[u] |= 1 << v
        adjacency[v] |= 1 << u
    return adjacency


def petersen():
    adjacency = [0] * 10
    for u, v in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
                 (5, 7), (7, 9), (9, 6), (6, 8), (8, 5),
                 (0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]:
        adjacency[u] |= 1 << v
        adjacency[v] |= 1 << u
    return adjacency


def grotzsch():
    """Mycielskian of C5: u_i joins the *neighbours* of v_i, w joins every u_i."""
    adjacency = [0] * 11
    pairs = [(i, (i + 1) % 5) for i in range(5)]
    for i in range(5):
        pairs += [(5 + i, (i - 1) % 5), (5 + i, (i + 1) % 5), (10, 5 + i)]
    for u, v in pairs:
        adjacency[u] |= 1 << v
        adjacency[v] |= 1 << u
    return adjacency


class Blowup:
    """Densities of the weighted blow-up graphon of one base graph."""

    def __init__(self, name, adjacency):
        bad = triangles(adjacency)
        if bad:
            raise ValueError(f"{name} is not triangle-free: {bad[:3]}")
        self.name = name
        self.size = len(adjacency)
        self.edges = edges_of(adjacency)
        bits = (np.arange(1 << self.size)[:, None] >> np.arange(self.size)[None, :]) & 1
        self.same = np.asarray(
            [bits[:, u] == bits[:, v] for u, v in self.edges]
        ).T.astype(np.float64)

    def densities(self, weights):
        products = np.asarray([weights[u] * weights[v] for u, v in self.edges])
        return 2.0 * products.sum(), 2.0 * float((self.same @ products).min())

    def exact_densities(self, weights):
        d_edge = 2 * sum(weights[u] * weights[v] for u, v in self.edges)
        best = None
        for pattern in range(1 << self.size):
            mono = sum(
                weights[u] * weights[v]
                for u, v in self.edges
                if ((pattern >> u) & 1) == ((pattern >> v) & 1)
            )
            if best is None or mono < best:
                best = mono
        return d_edge, 2 * best

    def in_band(self, weights):
        """Best in-band d_mono on this ray; dilution pulls a dense point back in."""
        d_edge, d_mono = self.densities(weights)
        if d_edge < LO_F - 1e-12:
            return -1.0, d_edge, 1.0
        theta = min(1.0, HI_F / d_edge)
        return d_mono * theta, d_edge * theta, theta

    def optimise(self, seed, restarts=6, steps=260):
        rng = random.Random(seed)
        best = (-1.0, None)
        for attempt in range(restarts):
            if attempt == 0:
                weights = np.full(self.size, 1.0 / self.size)
            else:
                weights = np.asarray([rng.random() + 1e-3 for _ in range(self.size)])
                weights /= weights.sum()
            value = self.in_band(weights)[0]
            scale = 0.35
            for step in range(steps):
                i, j = rng.randrange(self.size), rng.randrange(self.size)
                if i == j:
                    continue
                trial = weights.copy()
                shift = scale * rng.random() * trial[i]
                trial[i] -= shift
                trial[j] += shift
                candidate = self.in_band(trial)[0]
                if candidate > value:
                    weights, value = trial, candidate
                if step % 60 == 59:
                    scale *= 0.55
            if value > best[0]:
                best = (value, weights.copy())
        return best


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--top", type=int, default=12)
    args = parser.parse_args()

    with (args.flagsdp / "cache_n9.pkl").open("rb") as handle:
        states = pickle.load(handle)["states"]

    bases = []
    for name, adjacency in [("C5", cycle(5)), ("C7", cycle(7)), ("C9", cycle(9)),
                            ("C11", cycle(11)), ("Petersen", petersen()),
                            ("Grotzsch", grotzsch())]:
        bases.append(Blowup(name, adjacency))
    print(f"named bases verified triangle-free: {[b.name for b in bases]}")
    for i, (_order, masks) in enumerate(states):
        blow = Blowup(f"tf9#{i}", list(masks))
        if blow.edges:
            bases.append(blow)
    print(f"bases scanned: {len(bases)} (all triangle-free graphs on 9 vertices + named)\n")

    results = []
    for base in bases:
        value, weights = base.optimise(seed=hash(base.name) & 0xFFFF)
        if value <= 0:
            continue
        d_edge, d_mono = base.densities(weights)
        theta = min(1.0, HI_F / d_edge)
        results.append((value, base, weights, d_edge * theta, d_mono / d_edge, theta))
    results.sort(key=lambda row: -row[0])

    print(f"{'d_mono':>9} {'base':>10} {'h':>3} {'m':>3} {'d_edge':>8} "
          f"{'theta':>7} {'bip/m':>7}")
    for value, base, _w, d_edge, ratio, theta in results[: args.top]:
        print(f"{value:>9.6f} {base.name:>10} {base.size:>3} {len(base.edges):>3} "
              f"{d_edge:>8.4f} {theta:>7.4f} {ratio:>7.4f}")

    best_value, best_base, best_w, _de, _r, best_theta = results[0]
    print(f"\nbest in-band d_mono found : {best_value:.9f}")
    print(f"best bip/m ratio seen     : {max(row[4] for row in results):.9f}"
          f"   (C5 and Petersen both give exactly 1/5)")
    print(f"2/25                      : {0.08:.9f}")
    print(f"\n=> every valid certificate must satisfy delta >= {best_value - 0.08:+.6e}")
    print(f"   the certificate's claimed delta* is        {-9.878886951679021e-04:+.6e}")
    print(f"   claimed value excluded by this ceiling?    "
          f"{best_value - 0.08 > -9.878886951679021e-04}")
    print(f"\n   The whole slack of the conjecture inside the band is therefore at most")
    print(f"   2/25 - {best_value:.6f} = {0.08 - best_value:.6f}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
