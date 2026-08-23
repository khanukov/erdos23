#!/usr/bin/env python3
"""Audit the certificate's row families against concrete triangle-free graphons.

`verify_exact_fixed_gram_dual.py` checks that the stored multipliers are an
exactly feasible dual of one specific finite LP.  It cannot check the step that
turns a negative dual objective into a statement about bip(G): that step needs
every row of the LP to be a *valid* inequality, and needs the LP optimum to
dominate d_mono(W) - 2/25.  Both are mathematical claims about the row
semantics, and both are testable by evaluating the rows the verifier itself
reconstructs at the exact order-10 state vector of a genuine triangle-free
graphon.

This script reuses the verifier's own row reconstruction verbatim and reports:

  0. how much dual weight the two edge-density band rows carry;
  1. the sign of the fixed-Gram moment row;
  2. the signs of all rooted-Horn rows;
  3. positive semidefiniteness of the rooted pair-density matrices that the
     rooted-Horn rows are built from (a necessary condition for any genuine
     rooted pair density, since it is a nonnegative mixture of p p^T);
  4. the K7 and K8 per-root MaxCut envelopes against the true d_mono;
  5. the primal LP value eta at that graphon against the certificate objective.

Everything is exact rational arithmetic.
"""

from __future__ import annotations

import argparse
import itertools
import json
import pickle
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

N_ORDER9 = 1897
N_ORDER10 = 12172
N_K7_ROOTS = 107
N_K8_ROOTS = 410
K7_DENOMINATOR = 10 * 25 * 181440


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--public-anc", type=Path, required=True)
    parser.add_argument("--k7-cache", type=Path, required=True)
    parser.add_argument(
        "--graphon",
        type=Path,
        action="append",
        required=True,
        help="state vector produced by graphon_state_vector.py; repeatable",
    )
    return parser.parse_args()


def profile_classes(k, adjacency):
    out = []
    for size in range(k + 1):
        for subset in itertools.combinations(range(k), size):
            if all(not ((adjacency[a] >> b) & 1) for a in subset for b in subset if a < b):
                out.append(frozenset(subset))
    return out


def sorted_profiles(decomposition):
    return {
        root: tuple(
            sorted(
                {tuple(profile) for profile in decomposition["Rprofiles"][root]},
                key=lambda profile: (len(profile), profile),
            )
        )
        for root in range(N_K8_ROOTS)
    }


def pair_matrices(decomposition, profiles, indices):
    rows = [[] for _ in range(N_K8_ROOTS)]
    columns = [[] for _ in range(N_K8_ROOTS)]
    for state, contributions in enumerate(decomposition):
        for root, profile_a, profile_b in contributions:
            root = int(root)
            width = len(profiles[root])
            rows[root].append(
                indices[root][tuple(profile_a)] * width + indices[root][tuple(profile_b)]
            )
            columns[root].append(state)
    matrices = []
    for root in range(N_K8_ROOTS):
        width = len(profiles[root])
        matrix = coo_matrix(
            (np.ones(len(rows[root]), dtype=np.int16), (rows[root], columns[root])),
            shape=(width * width, N_ORDER10),
            dtype=np.int16,
        ).tocsr()
        matrix.sum_duplicates()
        matrices.append(matrix)
    return matrices


class Model:
    """Exactly the objects `verify_exact_fixed_gram_dual.py` reconstructs."""

    def __init__(self, flagsdp: Path, public_anc: Path, k7_cache: Path):
        with (flagsdp / "cache_n9.pkl").open("rb") as handle:
            cache = pickle.load(handle)
        lift = np.load(flagsdp / "c5lift_cache.npz", allow_pickle=True)
        counts = np.rint(np.asarray(lift["Dval"]) * 10).astype(np.int64)
        self.deletion_t = csr_matrix(
            (counts, (lift["Drow"], lift["Dcol"])),
            shape=(N_ORDER9, N_ORDER10),
            dtype=np.int64,
        ).T.tocsr()

        self.k7 = []
        for root in range(N_K7_ROOTS):
            with np.load(k7_cache / f"root_{root:03d}.npz", allow_pickle=False) as payload:
                adjacency = tuple(int(x) for x in payload["adjacency"])
                edge_raw = np.asarray(payload["edge_raw"], dtype=np.uint32)
                mass_raw = np.asarray(payload["mass_raw"], dtype=np.uint32)
            self.k7.append((adjacency, edge_raw, mass_raw, profile_classes(7, adjacency)))

        with (flagsdp / "u8_decomp.pkl").open("rb") as handle:
            edge_decomposition = pickle.load(handle)
        with (flagsdp / "u8_decomp_all.pkl").open("rb") as handle:
            all_decomposition = pickle.load(handle)
        self.edge_profiles = sorted_profiles(edge_decomposition)
        self.all_profiles = sorted_profiles(all_decomposition)
        edge_indices = [
            {p: i for i, p in enumerate(self.edge_profiles[r])} for r in range(N_K8_ROOTS)
        ]
        self.all_indices = [
            {p: i for i, p in enumerate(self.all_profiles[r])} for r in range(N_K8_ROOTS)
        ]
        self.edge_pairs = pair_matrices(
            edge_decomposition["decomp"], self.edge_profiles, edge_indices
        )
        self.all_pairs = pair_matrices(
            all_decomposition["decomp"], self.all_profiles, self.all_indices
        )
        with (public_anc / "mom_term_exact.pkl").open("rb") as handle:
            self.moment = pickle.load(handle)

    def k7_numerator(self, root, rule):
        _adjacency, edge_raw, mass_raw, _classes = self.k7[root]
        rule = np.asarray(rule, dtype=np.uint8)
        same = np.triu(np.equal.outer(rule, rule)).astype(np.uint32)
        same_raw = np.tensordot(edge_raw, same, axes=([1, 2], [0, 1])).astype(np.int64)
        order9 = 25 * same_raw - 2 * mass_raw.astype(np.int64)
        return np.asarray(self.deletion_t @ order9).ravel().astype(np.int64)

    def k7_mono_numerator(self, root, rule):
        """The same row without the -2/25 root-mass term (mono-edge part only)."""
        _adjacency, edge_raw, _mass_raw, _classes = self.k7[root]
        rule = np.asarray(rule, dtype=np.uint8)
        same = np.triu(np.equal.outer(rule, rule)).astype(np.uint32)
        same_raw = np.tensordot(edge_raw, same, axes=([1, 2], [0, 1])).astype(np.int64)
        return np.asarray(self.deletion_t @ same_raw).ravel().astype(np.int64)

    def k7_mass_numerator(self, root):
        _adjacency, _edge_raw, mass_raw, _classes = self.k7[root]
        return np.asarray(self.deletion_t @ mass_raw.astype(np.int64)).ravel().astype(np.int64)

    def k8_numerator(self, root, sides):
        sides = np.asarray(sides, dtype=np.uint8)
        coefficients = np.equal.outer(sides, sides).ravel().astype(np.int64)
        return np.asarray(self.edge_pairs[root].T @ coefficients).ravel().astype(np.int64)

    def horn_numerator(self, root, cycle):
        cycle = tuple(tuple(profile) for profile in cycle)
        indices = [self.all_indices[root][profile] for profile in cycle]
        width = len(self.all_profiles[root])
        coefficients = np.zeros((width, width), dtype=np.int64)
        coefficients[np.ix_(indices, indices)] = 1
        for i in range(5):
            coefficients[indices[i], indices[(i + 1) % 5]] -= 4
        return np.asarray(self.all_pairs[root].T @ coefficients.ravel()).ravel().astype(np.int64)


class Point:
    """A graphon, presented by its exact order-10 state vector."""

    def __init__(self, path: Path):
        with path.open("rb") as handle:
            payload = pickle.load(handle)
        self.name = f"{payload['base']} blow-up, weights " + ",".join(
            str(Fraction(n, d)) for n, d in payload["weights"]
        )
        self.support = np.asarray(payload["support"], dtype=np.int64)
        self.values = [Fraction(n, d) for n, d in payload["values"]]
        self.d_edge = Fraction(*payload["d_edge"])
        self.d_mono = Fraction(*payload["d_mono"])

    def dot(self, numerator):
        total = Fraction(0)
        for value, weight in zip(numerator[self.support].tolist(), self.values):
            if value:
                total += value * weight
        return total


def report(model: Model, certificate, point: Point):
    descriptors = certificate["descriptors"]
    denominator = int(certificate["multiplier_denominator"])
    static = [int(v) for v in certificate["static_numerators"]]
    multipliers = certificate["descriptor_numerators"]
    objective = Fraction(*certificate["objective"])
    band_low = Fraction(2486, 10000)
    band_high = Fraction(3197, 10000)

    print("=" * 78)
    print(f"graphon: {point.name}")
    print(f"  d_edge = {point.d_edge} = {float(point.d_edge)}"
          f"   inside the certificate band [{band_low}, {band_high}]? "
          f"{band_low <= point.d_edge <= band_high}")
    print(f"  d_mono = {point.d_mono} = {float(point.d_mono)}   (2/25 = 0.08)")

    moment_value = sum(
        Fraction(int(self_n), int(self_d)) * weight
        for (self_n, self_d), weight in zip(
            (model.moment[s] for s in point.support.tolist()), point.values
        )
    )
    print("\n[1] fixed-Gram moment row  <m,q> >= 0 :",
          moment_value >= 0, f"({float(moment_value):+.9e})")

    u7_mono, u7_mass, u8, horn_negative, horn_min = {}, {}, {}, [], None
    for row, descriptor in enumerate(descriptors):
        kind = descriptor["kind"]
        root = int(descriptor["root"])
        if kind == "k7":
            value = point.dot(model.k7_mono_numerator(root, descriptor["rule"])) / 181440
            u7_mono[root] = value if root not in u7_mono else min(u7_mono[root], value)
            if root not in u7_mass:
                u7_mass[root] = point.dot(model.k7_mass_numerator(root)) / 181440
        elif kind == "k8":
            value = point.dot(model.k8_numerator(root, descriptor["sides"])) / 90
            u8[root] = value if root not in u8 else min(u8[root], value)
        else:
            value = point.dot(model.horn_numerator(root, descriptor["cycle"])) / 90
            horn_min = value if horn_min is None else min(horn_min, value)
            if value < 0:
                horn_negative.append((row, root, value, int(multipliers[row])))

    horn_total = sum(1 for d in descriptors if d["kind"] == "horn")
    horn_weight = sum(int(multipliers[i]) for i, d in enumerate(descriptors) if d["kind"] == "horn")
    print(f"\n[2] rooted-Horn rows  <horn,q> >= 0 : "
          f"{len(horn_negative)} of {horn_total} are STRICTLY NEGATIVE")
    print(f"    most negative value: {float(horn_min):+.9e}")
    if horn_negative:
        horn_negative.sort(key=lambda item: item[2])
        bad_weight = sum(item[3] for item in horn_negative)
        print(f"    dual weight carried by the invalid rows: "
              f"{bad_weight}/{horn_weight} = {100.0 * bad_weight / horn_weight:.1f}% of the Horn weight")
        print(f"    {'row':>6} {'root':>5} {'value':>16} {'multiplier x 10^12':>22}")
        for row, root, value, weight in horn_negative[:8]:
            print(f"    {row:>6} {root:>5} {float(value):>16.9f} {weight:>22}")
        print(f"    roots involved: {sorted({r for _, r, _, _ in horn_negative})}")

    envelope7 = sum(u7_mono.values())
    mass7 = sum(u7_mass.values())
    envelope8 = sum(u8.values())
    print(f"\n[3] K7 per-root MaxCut envelope, sum over the 107 roots of the pool minimum")
    print(f"    mono part / root-mass part = {envelope7 / mass7} = {float(envelope7 / mass7)}"
          f"   vs d_mono = {float(point.d_mono)}")
    print(f"\n[4] K8 per-root MaxCut envelope U8, sum over the 410 roots of the pool minimum")
    print(f"    U8 = {envelope8} = {float(envelope8)}"
          f"   vs d_mono = {float(point.d_mono)}   ->  U8 >= d_mono ? "
          f"{envelope8 >= point.d_mono}")

    leg7 = (envelope7 - Fraction(2, 25) * mass7) / 10
    leg8 = Fraction(-2, 25) + envelope8
    eta = min(leg7, leg8)
    target = point.d_mono - Fraction(2, 25)
    print(f"\n[5] LP value at this graphon; every leg must dominate d_mono - 2/25"
          f" = {float(target):+.9e}")
    print(f"    K7 leg (dual weight {Fraction(static[3], denominator)}) = "
          f"{float(leg7):+.9e}   valid here? {leg7 >= target}")
    print(f"    K8 leg (dual weight {Fraction(static[4], denominator)}) = "
          f"{float(leg8):+.9e}   valid here? {leg8 >= target}")
    print(f"    eta = min(legs)              = {float(eta):+.9e}   "
          f"dominates d_mono - 2/25? {eta >= target}")
    print(f"    certificate objective delta* = {float(objective):+.9e}")


def psd_report(model: Model, point: Point):
    """A genuine rooted pair density is a nonnegative mixture of p p^T, hence PSD."""
    print("\n" + "=" * 78)
    print(f"[6] rooted pair-density matrices from u8_decomp_all.pkl at: {point.name}")
    dense_q = np.zeros(N_ORDER10)
    for state, value in zip(point.support.tolist(), point.values):
        dense_q[state] = float(value)
    failures = []
    for root in range(N_K8_ROOTS):
        width = len(model.all_profiles[root])
        matrix = (model.all_pairs[root] @ dense_q).reshape(width, width) / 90.0
        mass = matrix.sum()
        if mass <= 0:
            continue
        scale = max(1.0, float(np.abs(matrix).max()))
        if np.abs(matrix - matrix.T).max() > 1e-13 * scale:
            failures.append((root, float("nan"), mass, "asymmetric"))
            continue
        smallest = float(np.linalg.eigvalsh((matrix + matrix.T) / 2)[0])
        if smallest < -1e-12 * scale:
            failures.append((root, smallest, mass, "indefinite"))
    print(f"    root types with a NEGATIVE eigenvalue: {len(failures)} of {N_K8_ROOTS}")
    failures.sort(key=lambda item: item[1])
    print(f"    {'root':>5} {'min eigenvalue':>18} {'root mass':>14}")
    for root, smallest, mass, _why in failures[:8]:
        print(f"    {root:>5} {smallest:>18.9e} {mass:>14.9f}")
    if failures:
        exact_psd_witness(model, point, failures[0][0])


def exact_psd_witness(model: Model, point: Point, root: int):
    """Re-run the worst failure in exact rational arithmetic."""
    width = len(model.all_profiles[root])
    block = model.all_pairs[root][:, point.support].toarray().astype(object)
    exact = [[Fraction(0)] * width for _ in range(width)]
    for a in range(width):
        for b in range(width):
            total = Fraction(0)
            for entry, weight in zip(block[a * width + b].tolist(), point.values):
                if entry:
                    total += int(entry) * weight
            exact[a][b] = total / 90
    approx = np.array([[float(x) for x in row] for row in exact])
    vector = np.linalg.eigh((approx + approx.T) / 2)[1][:, 0]
    rational = [Fraction(int(round(v * 10 ** 6)), 10 ** 6) for v in vector]
    form = sum(
        rational[a] * exact[a][b] * rational[b] for a in range(width) for b in range(width)
    )
    norm = sum(v * v for v in rational)
    print(f"\n    exact rational witness for root {root} (width {width}):")
    print(f"      x^T M x = {form}")
    print(f"              = {float(form):+.9e}   with x^T x = {float(norm):.6f}")
    print(f"      M is provably not positive semidefinite: {form < 0}")


def main() -> int:
    args = parse_args()
    if args.certificate.suffix.lower() == ".json":
        certificate = json.loads(args.certificate.read_text())
    else:
        with args.certificate.open("rb") as handle:
            certificate = pickle.load(handle)
    if certificate.get("format") != "erdos23-fixed-gram-exact-dual-v1":
        raise RuntimeError("unsupported certificate format")

    denominator = int(certificate["multiplier_denominator"])
    static = [int(v) for v in certificate["static_numerators"]]
    print("[0] dual weight on the two edge-density band rows")
    print(f"    y(d_edge <= 0.3197) = {static[0]}/{denominator}")
    print(f"    y(d_edge >= 0.2486) = {static[1]}/{denominator}")
    if static[0] == 0 and static[1] == 0:
        print("    Both are ZERO: the certificate never uses the middle band, so its")
        print("    bound applies verbatim to every triangle-free graphon, including the")
        print("    C5 blow-up at d_edge = 0.4 where d_mono = 2/25 exactly.")

    model = Model(args.flagsdp.resolve(), args.public_anc.resolve(), args.k7_cache.resolve())
    points = [Point(path) for path in args.graphon]
    for point in points:
        report(model, certificate, point)
    psd_report(model, points[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
