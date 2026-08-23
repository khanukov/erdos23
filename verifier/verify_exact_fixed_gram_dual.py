#!/usr/bin/env python3
"""Independent exact replay of an Erdős 23 fixed-Gram dual certificate.

This verifier does not import the LP generator or trust a floating matrix.  It
rebuilds every K7, K8 and Horn functional from exact combinatorial caches,
recomputes all 12,172 dual residuals with rational arithmetic, and checks both
the stored finite-order target and the global nonpositive-objective closure.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import pickle
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

try:
    import gmpy2
except ModuleNotFoundError:  # Exact, portable, and slower.
    import exact_arithmetic as gmpy2


K7_DENOMINATOR = 10 * 25 * 181440
COMMON_DENOMINATOR = K7_DENOMINATOR
DIGIT_BASE = 100_000


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--public-anc", type=Path, required=True)
    parser.add_argument("--k7-cache", type=Path, required=True)
    parser.add_argument("--target-n", type=int)
    return parser.parse_args()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_pair_matrices(decomposition, profiles_by_root, profile_indices, nstates):
    nroots = len(profiles_by_root)
    rows = [[] for _ in range(nroots)]
    columns = [[] for _ in range(nroots)]
    for state, contributions in enumerate(decomposition):
        for root, profile_a, profile_b in contributions:
            root = int(root)
            width = len(profiles_by_root[root])
            index = (
                profile_indices[root][tuple(profile_a)] * width
                + profile_indices[root][tuple(profile_b)]
            )
            rows[root].append(index)
            columns[root].append(state)
    matrices = []
    for root in range(nroots):
        width = len(profiles_by_root[root])
        values = np.ones(len(rows[root]), dtype=np.int16)
        matrix = coo_matrix(
            (values, (rows[root], columns[root])),
            shape=(width * width, nstates),
            dtype=np.int16,
        ).tocsr()
        matrix.sum_duplicates()
        matrices.append(matrix)
    return matrices


def descriptor_key(descriptor):
    kind = descriptor["kind"]
    if kind == "k7":
        return (kind, int(descriptor["root"]), tuple(descriptor["rule"]))
    if kind == "k8":
        return (kind, int(descriptor["root"]), tuple(descriptor["sides"]))
    if kind == "horn":
        return (kind, int(descriptor["root"]), tuple(tuple(x) for x in descriptor["cycle"]))
    raise ValueError(f"unknown descriptor kind {kind}")


def profile_classes(k, adjacency):
    """Enumerate independent root subsets without importing the LP generator."""
    classes = []
    for size in range(k + 1):
        for subset in itertools.combinations(range(k), size):
            if all(
                not ((adjacency[a] >> b) & 1)
                for a in subset
                for b in subset
                if a < b
            ):
                classes.append(frozenset(subset))
    return classes


def main():
    args = parse_args()
    flagsdp = args.flagsdp.resolve()
    public_anc = args.public_anc.resolve()
    k7_cache = args.k7_cache.resolve()
    if args.certificate.suffix.lower() == ".json":
        certificate = json.loads(args.certificate.read_text())
    else:
        with args.certificate.open("rb") as handle:
            certificate = pickle.load(handle)
    if certificate.get("format") != "erdos23-fixed-gram-exact-dual-v1":
        raise RuntimeError("unsupported certificate format")
    target_n = int(args.target_n or certificate.get("target_n", 41))
    if target_n <= 0:
        raise RuntimeError("target n must be positive")
    if "target_n" in certificate and int(certificate["target_n"]) != target_n:
        raise RuntimeError("requested target n differs from certificate metadata")
    descriptors = certificate["descriptors"]
    if len({descriptor_key(item) for item in descriptors}) != len(descriptors):
        raise RuntimeError("duplicate descriptor in certificate")
    D = int(certificate["multiplier_denominator"])
    static = [int(value) for value in certificate["static_numerators"]]
    multipliers = np.asarray(certificate["descriptor_numerators"], dtype=np.int64)
    if D <= 0 or len(static) != 5 or len(multipliers) != len(descriptors):
        raise RuntimeError("malformed multiplier block")
    if min(static) < 0 or np.any(multipliers < 0):
        raise RuntimeError("negative upper-row multiplier")
    if static[3] + static[4] != D:
        raise RuntimeError("eta-coordinate equality is not exact")

    with (flagsdp / "cache_n9.pkl").open("rb") as handle:
        cache = pickle.load(handle)
    lift = np.load(flagsdp / "c5lift_cache.npz", allow_pickle=True)
    n9 = len(cache["states"])
    n10 = int(lift["nJ"])
    if (n9, n10) != (1897, 12172):
        raise RuntimeError(f"unexpected state dimensions {(n9, n10)}")
    deletion_counts = np.rint(np.asarray(lift["Dval"]) * 10).astype(np.int64)
    if not np.array_equal(deletion_counts.astype(float) / 10.0, np.asarray(lift["Dval"])):
        raise RuntimeError("deletion lift is not exact over denominator 10")
    deletion_count_t = csr_matrix(
        (deletion_counts, (lift["Drow"], lift["Dcol"])),
        shape=(n9, n10),
        dtype=np.int64,
    ).T.tocsr()
    # Recover edge counts without comparing two binary-float evaluation paths.
    # Each edge of an order-10 graph survives exactly eight vertex deletions.
    edge9 = np.rint(np.asarray(cache["dedge"], dtype=float) * 36).astype(np.int64)
    if not np.allclose(edge9.astype(float) / 36.0, cache["dedge"], atol=2e-15, rtol=0):
        raise RuntimeError("order-9 edge-count recovery failed")
    deletion_edge_sum = np.asarray(deletion_count_t @ edge9).ravel().astype(np.int64)
    if np.any(deletion_edge_sum % 8):
        raise RuntimeError("order-10 deletion edge sum is not divisible by eight")
    density_numerators = deletion_edge_sum // 8

    manifest = json.loads((k7_cache / "manifest.json").read_text())
    if (
        manifest.get("format") != "erdos23-k7-exact-count-cache-v1"
        or int(manifest.get("denominator", -1)) != 181440
        or int(manifest.get("n_states", -1)) != n9
        or int(manifest.get("n_roots", -1)) != 107
    ):
        raise RuntimeError("invalid K7 exact-count manifest")
    k7 = []
    for root in range(107):
        with np.load(k7_cache / f"root_{root:03d}.npz", allow_pickle=False) as payload:
            if int(payload["root"]) != root:
                raise RuntimeError(f"K7 cache root mismatch at {root}")
            adjacency = tuple(int(x) for x in payload["adjacency"])
            edge_raw = np.asarray(payload["edge_raw"], dtype=np.uint32)
            mass_raw = np.asarray(payload["mass_raw"], dtype=np.uint32)
        classes = profile_classes(7, adjacency)
        if edge_raw.shape != (n9, len(classes), len(classes)) or mass_raw.shape != (n9,):
            raise RuntimeError(f"K7 cache shape mismatch at {root}")
        k7.append((adjacency, edge_raw, mass_raw, classes))

    with (flagsdp / "u8_decomp.pkl").open("rb") as handle:
        edge_decomp = pickle.load(handle)
    with (flagsdp / "u8_decomp_all.pkl").open("rb") as handle:
        all_decomp = pickle.load(handle)
    if int(edge_decomp["nR"]) != 410 or int(all_decomp["nR"]) != 410:
        raise RuntimeError("unexpected number of U8 roots")
    edge_profiles = {
        root: tuple(
            sorted(
                {tuple(profile) for profile in edge_decomp["Rprofiles"][root]},
                key=lambda profile: (len(profile), profile),
            )
        )
        for root in range(410)
    }
    all_profiles = {
        root: tuple(
            sorted(
                {tuple(profile) for profile in all_decomp["Rprofiles"][root]},
                key=lambda profile: (len(profile), profile),
            )
        )
        for root in range(410)
    }
    edge_indices = [
        {profile: index for index, profile in enumerate(edge_profiles[root])}
        for root in range(410)
    ]
    all_indices = [
        {profile: index for index, profile in enumerate(all_profiles[root])}
        for root in range(410)
    ]
    edge_pairs = exact_pair_matrices(
        edge_decomp["decomp"], edge_profiles, edge_indices, n10
    )
    all_pairs = exact_pair_matrices(
        all_decomp["decomp"], all_profiles, all_indices, n10
    )

    root_sums7 = [0] * 107
    root_sums8 = [0] * 410
    rows_indices = []
    rows_values = []
    indptr = [0]
    started = time.time()
    for row, (multiplier, descriptor) in enumerate(zip(multipliers, descriptors)):
        kind = descriptor["kind"]
        root = int(descriptor["root"])
        if kind == "k7":
            if not 0 <= root < 107:
                raise RuntimeError("K7 root out of range")
            _adjacency, edge_raw, mass_raw, classes = k7[root]
            rule = np.asarray(descriptor["rule"], dtype=np.uint8)
            if rule.shape != (len(classes),) or np.any(rule > 1):
                raise RuntimeError(f"invalid K7 rule at descriptor {row}")
            same = np.triu(np.equal.outer(rule, rule)).astype(np.uint32)
            same_raw = np.tensordot(edge_raw, same, axes=([1, 2], [0, 1])).astype(np.int64)
            order9 = 25 * same_raw - 2 * mass_raw.astype(np.int64)
            numerator = np.asarray(deletion_count_t @ order9).ravel().astype(np.int64)
            factor = 1
            root_sums7[root] += int(multiplier)
        elif kind == "k8":
            if not 0 <= root < 410:
                raise RuntimeError("K8 root out of range")
            sides = np.asarray(descriptor["sides"], dtype=np.uint8)
            if sides.shape != (len(edge_profiles[root]),) or np.any(sides > 1):
                raise RuntimeError(f"invalid K8 sides at descriptor {row}")
            coefficients = np.equal.outer(sides, sides).ravel().astype(np.int64)
            numerator = np.asarray(edge_pairs[root].T @ coefficients).ravel().astype(np.int64)
            factor = COMMON_DENOMINATOR // 90
            root_sums8[root] += int(multiplier)
        elif kind == "horn":
            if not 0 <= root < 410:
                raise RuntimeError("Horn root out of range")
            cycle = tuple(tuple(profile) for profile in descriptor["cycle"])
            if len(cycle) != 5 or len(set(cycle)) != 5:
                raise RuntimeError(f"invalid Horn cycle at descriptor {row}")
            try:
                indices = [all_indices[root][profile] for profile in cycle]
            except KeyError as error:
                raise RuntimeError(f"invalid Horn profile at descriptor {row}") from error
            width = len(all_profiles[root])
            coefficients = np.zeros((width, width), dtype=np.int64)
            coefficients[np.ix_(indices, indices)] = 1
            for index in range(5):
                coefficients[indices[index], indices[(index + 1) % 5]] -= 4
            numerator = np.asarray(all_pairs[root].T @ coefficients.ravel()).ravel().astype(np.int64)
            factor = COMMON_DENOMINATOR // 90
        else:
            raise RuntimeError(f"unknown descriptor kind {kind}")
        nonzero = np.flatnonzero(numerator)
        scaled = numerator[nonzero]
        if factor != 1:
            scaled = scaled * factor
        rows_indices.append(nonzero.astype(np.int32))
        rows_values.append(np.asarray(scaled, dtype=np.int64))
        indptr.append(indptr[-1] + len(nonzero))
        if (row + 1) % 500 == 0:
            print(f"replayed rows {row + 1}/{len(descriptors)}; nnz={indptr[-1]}", flush=True)

    if min(root_sums7) < static[3]:
        raise RuntimeError("a K7 root multiplier sum is below its leg")
    if min(root_sums8) < static[4]:
        raise RuntimeError("a K8 root multiplier sum is below its leg")
    exact_matrix = csr_matrix(
        (
            np.concatenate(rows_values),
            np.concatenate(rows_indices),
            np.asarray(indptr, dtype=np.int64),
        ),
        shape=(len(descriptors), n10),
        dtype=np.int64,
    )
    exact_matrix.check_format(full_check=True)
    print(
        f"replayed exact matrix: shape={exact_matrix.shape} nnz={exact_matrix.nnz} "
        f"build={time.time() - started:.1f}s",
        flush=True,
    )

    absolute_column_sums = np.asarray(abs(exact_matrix).sum(axis=0)).ravel().astype(np.int64)
    bound = int(absolute_column_sums.max(initial=0)) * (DIGIT_BASE - 1)
    if bound >= np.iinfo(np.int64).max:
        raise OverflowError(f"digit multiplication bound exceeds int64: {bound}")
    remaining = multipliers.copy()
    pieces = []
    while np.any(remaining):
        digit = remaining % DIGIT_BASE
        pieces.append(np.asarray(exact_matrix.T @ digit, dtype=np.int64).ravel())
        remaining //= DIGIT_BASE
    if not pieces:
        pieces = [np.zeros(n10, dtype=np.int64)]

    def combined(index):
        total = 0
        scale = 1
        for piece in pieces:
            total += int(piece[index]) * scale
            scale *= DIGIT_BASE
        return total

    # Guard the digit decomposition itself against an implementation or
    # overflow mistake by replaying representative columns as Python-integer
    # sparse dot products (including the tight coordinate 46).
    audit_states = sorted({0, 1, 17, 46, 107, 1897, 4096, n10 - 2, n10 - 1})
    for state in audit_states:
        column = exact_matrix.getcol(state).tocoo()
        direct = sum(
            int(multipliers[int(row)]) * int(value)
            for row, value in zip(column.row, column.data)
        )
        if direct != combined(state):
            raise RuntimeError(f"digit-product audit failed at state {state}")

    moment_path = public_anc / "mom_term_exact.pkl"
    if certificate.get("source", {}).get("moment_exact_sha256") != file_sha256(moment_path):
        raise RuntimeError("exact moment vector hash mismatch")
    with moment_path.open("rb") as handle:
        moment = pickle.load(handle)
    if len(moment) != n10:
        raise RuntimeError("exact moment vector dimension mismatch")

    rho = gmpy2.mpq(*certificate["rho"])
    loads = []
    for state, (moment_numerator, moment_denominator) in enumerate(moment):
        load = (
            gmpy2.mpq(static[2] * int(moment_numerator), D * int(moment_denominator))
            + gmpy2.mpq(combined(state), D * COMMON_DENOMINATOR)
            - gmpy2.mpq(
                (static[0] - static[1]) * int(density_numerators[state]), D * 45
            )
        )
        loads.append(load)
    residuals = [rho - load for load in loads]
    minimum_residual = min(residuals)
    argmax = residuals.index(minimum_residual)
    if minimum_residual < 0:
        raise RuntimeError(f"negative q-coordinate residual at state {argmax}")

    objective = (
        rho
        + gmpy2.mpq(static[0] * 3197, D * 10000)
        - gmpy2.mpq(static[1] * 2486, D * 10000)
        - gmpy2.mpq(2 * static[4], D * 25)
    )
    stored_objective = gmpy2.mpq(*certificate["objective"])
    if objective != stored_objective:
        raise RuntimeError("stored objective does not match exact replay")
    global_closure = objective <= 0
    if "global_closure" in certificate and bool(certificate["global_closure"]) != global_closure:
        raise RuntimeError("stored global-closure flag does not match exact objective")
    threshold = gmpy2.mpq(2, 25 * target_n * target_n)
    margin = threshold - objective
    if margin <= 0:
        raise RuntimeError("exact objective does not cross the requested finite-order threshold")
    integrality_error = gmpy2.mpq(25 * target_n * target_n, 2) * objective
    if integrality_error >= 1:
        raise RuntimeError("finite-graph integrality error is not strictly below one")
    if certificate.get("argmax_state") != argmax:
        raise RuntimeError("stored argmax coordinate does not match replay")

    print("EXACT_DUAL_REPLAY_OK", flush=True)
    print(f"descriptors={len(descriptors)} nnz={exact_matrix.nnz} digits={len(pieces)}", flush=True)
    print(f"digit_product_audit={len(audit_states)}/{len(audit_states)}", flush=True)
    print(f"min_q_residual={minimum_residual} at state={argmax}", flush=True)
    print(f"objective={objective} ~= {float(objective):.15e}", flush=True)
    print(f"global_closure={global_closure}", flush=True)
    print(f"threshold_n{target_n}={threshold} ~= {float(threshold):.15e}", flush=True)
    print(f"strict_margin={margin} ~= {float(margin):.15e}", flush=True)
    print(
        f"(25/2)*{target_n}^2*objective={integrality_error} ~= "
        f"{float(integrality_error):.15e} < 1",
        flush=True,
    )
    if global_closure:
        print(
            "GLOBAL_CONCLUSION: beta(G)<=N^2/25 for every finite triangle-free "
            "graph G on N vertices (using the cited density-tail and blow-up transfer)."
        )
    else:
        print(
            f"CONCLUSION: beta(G)<={target_n * target_n} for every triangle-free G on "
            f"{5 * target_n} vertices; C5[{target_n}] is sharp"
        )


if __name__ == "__main__":
    main()
