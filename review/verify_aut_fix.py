#!/usr/bin/env python3
"""Root cause of the K8/Horn defect, and a check that averaging over Aut repairs it.

`u8_decomp*.pkl` records exactly 90 contributions per order-10 state — that is
45 free pairs times 2 orderings, i.e. **one** labelling of each 8-element root
subset.  The K7 cache, by contrast, is normalised by 9P7 = 181,440: it sums over
*all* ordered injections of the 7 root vertices.  That is the difference between
the two families, and it is why only the K8 side misbehaves.

With a single labelling per subset, the labelling that a canonical form assigns
depends on how the subset is presented inside the state, hence on the two free
vertices.  Writing `lambda(R,u,v) = alpha . lambda_0` with `alpha in Aut(tau)`,
the recorded matrix is

    M = E[ e_{alpha a_0} e_{alpha b_0}^T ],       alpha = alpha(R,u,v),

instead of the true rooted pair density `M_0 = E_R[1(tau) p p^T]`.  Averaging
the recorded pairs over `Aut(tau)` removes the arbitrary `alpha` outright:

    (1/|Aut|) sum_beta  e_{beta alpha a_0} e_{beta alpha b_0}^T
        = (1/|Aut|) sum_gamma e_{gamma a_0} e_{gamma b_0}^T,

so the symmetrised matrix is `(1/|Aut|) sum_gamma gamma M_0 gamma^T`, a
nonnegative mixture of PSD matrices.  The same substitution repairs the K8
per-root MaxCut rows: after symmetrisation the row is an average over `gamma` of
`<h_{tau, c . gamma}, q>`, and each of those *is* an average over labelled roots
of the mono-edge density of one fixed two-colouring, hence at least `d_mono`
times the root mass.  Summing over root types then gives `U8 >= d_mono`, which
is exactly the property that fails today.

This script checks the mechanism numerically: it recovers each root graph,
computes `Aut(tau)`, symmetrises, and re-tests positive semidefiniteness.

Caveat: the root graphs are not stored in `u8_decomp_all.pkl`, so they are
recovered from the recorded profiles (a pair is an edge iff no recorded profile
contains both).  That recovery is exact for every root whose non-edges all
co-occur in some recorded profile; root 409, the densest type, comes back with
spurious edges and therefore a possibly too small Aut group.  It is not one of
the failing roots, and the conclusion does not depend on it.  Regenerating the
decomposition with the K7 convention (sum over all 10P8 ordered injections) is
the real fix; this script only demonstrates that the diagnosis is right.
"""

from __future__ import annotations

import argparse
import itertools
import pickle
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

N_ORDER10 = 12172
N_ROOTS = 410
ROOT_SIZE = 8


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--graphon", type=Path, required=True)
    return parser.parse_args()


def sorted_profiles(decomposition):
    return {
        root: tuple(
            sorted(
                {tuple(p) for p in decomposition["Rprofiles"][root]},
                key=lambda p: (len(p), p),
            )
        )
        for root in range(N_ROOTS)
    }


def recover_root(profiles):
    """A pair is a non-edge as soon as some recorded profile contains both."""
    together = set()
    for profile in profiles:
        for i, j in itertools.combinations(profile, 2):
            together.add((i, j))
    adjacency = [0] * ROOT_SIZE
    for i, j in itertools.combinations(range(ROOT_SIZE), 2):
        if (i, j) not in together:
            adjacency[i] |= 1 << j
            adjacency[j] |= 1 << i
    return adjacency


def automorphisms(adjacency):
    degree = [bin(mask).count("1") for mask in adjacency]
    found, image = [], [-1] * ROOT_SIZE

    def extend(depth, used):
        if depth == ROOT_SIZE:
            found.append(tuple(image))
            return
        for w in range(ROOT_SIZE):
            if (used >> w) & 1 or degree[w] != degree[depth]:
                continue
            if all(
                ((adjacency[depth] >> j) & 1) == ((adjacency[w] >> image[j]) & 1)
                for j in range(depth)
            ):
                image[depth] = w
                extend(depth + 1, used | (1 << w))
                image[depth] = -1

    extend(0, 0)
    return found


def main() -> int:
    args = parse_args()
    with (args.flagsdp / "u8_decomp_all.pkl").open("rb") as handle:
        decomposition = pickle.load(handle)
    profiles = sorted_profiles(decomposition)

    per_state = {len(entry) for entry in decomposition["decomp"]}
    print(f"contributions per order-10 state: {sorted(per_state)}"
          f"   (45 free pairs x 2 orderings = one labelling per 8-subset)")
    print("the K7 cache instead normalises by 9P7 = 181440, i.e. all ordered injections\n")

    roots = [recover_root(profiles[r]) for r in range(N_ROOTS)]
    groups = [automorphisms(adjacency) for adjacency in roots]
    print(f"|Aut(tau)|: min {min(len(g) for g in groups)}, "
          f"max {max(len(g) for g in groups)}, "
          f"nontrivial for {sum(1 for g in groups if len(g) > 1)} of {N_ROOTS} root types")

    with args.graphon.open("rb") as handle:
        payload = pickle.load(handle)
    q = np.zeros(N_ORDER10)
    for state, (numerator, denominator) in zip(payload["support"], payload["values"]):
        q[state] = numerator / denominator

    index = [{p: i for i, p in enumerate(profiles[r])} for r in range(N_ROOTS)]
    rows = [[] for _ in range(N_ROOTS)]
    columns = [[] for _ in range(N_ROOTS)]
    for state, contributions in enumerate(decomposition["decomp"]):
        for root, profile_a, profile_b in contributions:
            root = int(root)
            width = len(profiles[root])
            rows[root].append(
                index[root][tuple(profile_a)] * width + index[root][tuple(profile_b)]
            )
            columns[root].append(state)

    before, after = [], []
    for root in range(N_ROOTS):
        width = len(profiles[root])
        matrix = coo_matrix(
            (np.ones(len(rows[root])), (rows[root], columns[root])),
            shape=(width * width, N_ORDER10),
        ).tocsr()
        recorded = (matrix @ q).reshape(width, width) / 90.0
        if recorded.sum() <= 0:
            continue
        if np.linalg.eigvalsh((recorded + recorded.T) / 2)[0] < -1e-12 * max(
            1.0, float(np.abs(recorded).max())
        ):
            before.append(root)

        # extend the index to the Aut-closure, then average over Aut(tau)
        closure = list(profiles[root])
        position = {p: i for i, p in enumerate(closure)}
        for permutation in groups[root]:
            for profile in profiles[root]:
                moved = tuple(sorted(permutation[x] for x in profile))
                if moved not in position:
                    position[moved] = len(closure)
                    closure.append(moved)
        size = len(closure)
        padded = np.zeros((size, size))
        padded[:width, :width] = recorded
        averaged = np.zeros((size, size))
        for permutation in groups[root]:
            target = np.asarray(
                [position[tuple(sorted(permutation[x] for x in p))] for p in closure]
            )
            averaged[np.ix_(target, target)] += padded
        averaged /= len(groups[root])
        if np.linalg.eigvalsh((averaged + averaged.T) / 2)[0] < -1e-12 * max(
            1.0, float(np.abs(averaged).max())
        ):
            after.append(root)

    name = f"{payload['base']} blow-up, weights " + ",".join(
        str(Fraction(n, d)) for n, d in payload["weights"]
    )
    print(f"\nrooted pair-density matrices at: {name}")
    print(f"  non-PSD before Aut-averaging : {len(before)}  {before}")
    print(f"  non-PSD after  Aut-averaging : {len(after)}  {after}")
    print(f"  |Aut| of the failing types   : "
          f"{ {r: len(groups[r]) for r in before} }")
    if before and not after:
        print("\n  Every failing type has a nontrivial automorphism group, and averaging")
        print("  over it restores positive semidefiniteness: the defect is the arbitrary")
        print("  choice of one labelling per 8-subset, exactly as diagnosed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
