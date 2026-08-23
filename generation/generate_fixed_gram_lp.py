#!/usr/bin/env python3
"""Generate an order-10 certificate using the public fixed Gram moment row.

The public arXiv bundle exposes the final manifest-PSD moment functional but not
the 6,359 dynamically generated MaxCut/Horn rows.  This program regenerates a
new row pool from scratch.  It replaces the non-deterministic collection of
individual moment eigenvector rows by the already independently reconstructed
single Gram functional.  Every generated row is stored by a compact exact
descriptor, not as an opaque floating-point sparse vector.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import pickle
import sys
import time
from pathlib import Path

import highspy
import numpy as np
from joblib import Parallel, delayed
from scipy.sparse import coo_matrix, csr_matrix, vstack


INF = highspy.kHighsInf
LO = 0.2486
HI = 0.3197
K7_DENOMINATOR = 10 * 25 * 181440


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--public-anc", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--max-iterations", type=int, default=30)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--k7-restarts", type=int, default=96)
    parser.add_argument("--k7-keep", type=int, default=3)
    parser.add_argument("--horn-keep", type=int, default=2)
    parser.add_argument("--solver", choices=("ipm", "simplex"), default="ipm")
    parser.add_argument("--pool-cap", type=int, default=700)
    parser.add_argument("--slack-drop", type=float, default=1e-5)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--crossover",
        action="store_true",
        help="request a basic primal/dual solution; use for the final dual extraction",
    )
    return parser.parse_args()


def canonical_binary_rule(rule):
    rule = tuple(int(value) for value in rule)
    complement = tuple(1 - value for value in rule)
    return min(rule, complement)


def canonical_cycle(cycle):
    cycle = tuple(tuple(profile) for profile in cycle)
    variants = []
    for base in (cycle, tuple(reversed(cycle))):
        for shift in range(5):
            variants.append(base[shift:] + base[:shift])
    return min(variants)


class CertificateGenerator:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.flagsdp = args.flagsdp.resolve()
        self.public_anc = args.public_anc.resolve()
        sys.path.insert(0, str(self.flagsdp))

        import flag_cutgen as flag_cutgen
        from cutting_plane_u8 import maxcut_coloring
        from envelope_horn import horn_tuples_for_R
        from run_k7b import sep_multi

        self.fc = flag_cutgen
        self.maxcut_coloring = maxcut_coloring
        self.horn_tuples_for_R = horn_tuples_for_R
        self.sep_multi = sep_multi

        with (self.flagsdp / "cache_n9.pkl").open("rb") as handle:
            self.cache = pickle.load(handle)
        lift = np.load(self.flagsdp / "c5lift_cache.npz", allow_pickle=True)
        with (self.flagsdp / "u8_decomp.pkl").open("rb") as handle:
            edge_decomp = pickle.load(handle)
        with (self.flagsdp / "u8_decomp_all.pkl").open("rb") as handle:
            all_decomp = pickle.load(handle)
        self.k7_raw_counts = False
        compact_k7 = self.args.state.with_name("k7_compact_v1")
        compact_manifest = compact_k7 / "manifest.json"
        if compact_manifest.exists():
            manifest = json.loads(compact_manifest.read_text())
            if (
                manifest.get("format") != "erdos23-k7-exact-count-cache-v1"
                or int(manifest.get("denominator", -1)) != 181440
                or int(manifest.get("n_states", -1)) != len(self.cache["states"])
                or int(manifest.get("n_roots", -1)) != 107
            ):
                raise RuntimeError("invalid compact K7 cache manifest")
            types7 = self.fc.fe.enumerate_graphs(7, triangle_free=True)
            self.k7_types = []
            for root, (_k, adjacency) in enumerate(types7):
                with np.load(compact_k7 / f"root_{root:03d}.npz", allow_pickle=False) as data:
                    stored_adjacency = tuple(int(x) for x in data["adjacency"])
                    if stored_adjacency != tuple(int(x) for x in adjacency):
                        raise RuntimeError(f"compact K7 root order mismatch at {root}")
                    edge_raw = np.asarray(data["edge_raw"], dtype=np.uint32)
                    mass_raw = np.asarray(data["mass_raw"], dtype=np.uint32)
                classes = self.fc.profile_classes(7, adjacency)
                if edge_raw.shape != (len(self.cache["states"]), len(classes), len(classes)):
                    raise RuntimeError(f"compact K7 tensor shape mismatch at {root}")
                self.k7_types.append((7, adjacency, edge_raw, mass_raw, classes))
            self.k7_raw_counts = True
            print(f"loaded {len(self.k7_types)} exact-count K7 roots from {compact_k7}", flush=True)
        else:
            with (self.flagsdp / "k7_precompute.pkl").open("rb") as handle:
                self.k7_types = pickle.load(handle)
        with (self.public_anc / "mom_term_exact.pkl").open("rb") as handle:
            moment_raw = pickle.load(handle)
        with (self.public_anc / "horn_dual.pkl").open("rb") as handle:
            public_dual = pickle.load(handle)

        self.states9 = self.cache["states"]
        self.n9 = len(self.states9)
        self.n10 = int(lift["nJ"])
        self.n7 = len(self.k7_types)
        self.n8 = int(edge_decomp["nR"])
        if (self.n9, self.n10, self.n7, self.n8) != (1897, 12172, 107, 410):
            raise RuntimeError(
                f"unexpected dimensions {(self.n9, self.n10, self.n7, self.n8)}"
            )

        self.eta_col = self.n10
        self.u7_col = self.n10 + 1
        self.u8_col = self.u7_col + self.n7
        self.nvars = self.n10 + 1 + self.n7 + self.n8

        self.deletion = csr_matrix(
            (lift["Dval"], (lift["Drow"], lift["Dcol"])),
            shape=(self.n9, self.n10),
        )
        count_data = np.rint(np.asarray(lift["Dval"]) * 10).astype(np.int64)
        if not np.allclose(count_data / 10.0, lift["Dval"], atol=0, rtol=0):
            raise RuntimeError("9->10 deletion lift is not integral over denominator 10")
        self.deletion_count_t = csr_matrix(
            (count_data, (lift["Drow"], lift["Dcol"])),
            shape=(self.n9, self.n10),
            dtype=np.int64,
        ).T.tocsr()

        self.edge_decomp = edge_decomp["decomp"]
        self.all_decomp = all_decomp["decomp"]
        self.root_profiles = {
            int(root): tuple(
                sorted(
                    {tuple(profile) for profile in profiles},
                    key=lambda profile: (len(profile), profile),
                )
            )
            for root, profiles in all_decomp["Rprofiles"].items()
        }
        self.edge_profiles = {
            int(root): tuple(
                sorted(
                    {tuple(profile) for profile in profiles},
                    key=lambda profile: (len(profile), profile),
                )
            )
            for root, profiles in edge_decomp["Rprofiles"].items()
        }
        self.profile_indices = [
            {profile: index for index, profile in enumerate(self.root_profiles[root])}
            for root in range(self.n8)
        ]
        self.edge_profile_indices = [
            {profile: index for index, profile in enumerate(self.edge_profiles[root])}
            for root in range(self.n8)
        ]
        pair_cache_path = self.args.state.with_name("u8_pair_matrices_v1.pkl")
        if pair_cache_path.exists():
            with pair_cache_path.open("rb") as handle:
                pair_cache = pickle.load(handle)
            if (
                pair_cache.get("profiles") != self.root_profiles
                or pair_cache.get("edge_profiles") != self.edge_profiles
            ):
                raise RuntimeError("cached U8 profile order does not match regenerated data")
            self.edge_pair_matrices = pair_cache["edge"]
            self.all_pair_matrices = pair_cache["all"]
            print(f"loaded exact U8 pair matrices from {pair_cache_path}", flush=True)
        else:
            started = time.time()
            self.edge_pair_matrices = self.build_pair_matrices(
                self.edge_decomp, self.edge_profiles, self.edge_profile_indices
            )
            self.all_pair_matrices = self.build_pair_matrices(
                self.all_decomp, self.root_profiles, self.profile_indices
            )
            pair_cache = {
                "format": "erdos23-u8-pair-matrices-v1",
                "profiles": self.root_profiles,
                "edge_profiles": self.edge_profiles,
                "edge": self.edge_pair_matrices,
                "all": self.all_pair_matrices,
            }
            with pair_cache_path.open("wb") as handle:
                pickle.dump(pair_cache, handle, protocol=4)
            print(
                f"built and cached exact U8 pair matrices in {time.time() - started:.1f}s",
                flush=True,
            )

        self.moment_float = np.asarray(
            [numerator / denominator for numerator, denominator in moment_raw], dtype=float
        )
        self.moment_raw = moment_raw
        density9 = np.asarray(self.cache["dedge"], dtype=float)
        self.density10 = np.asarray(self.deletion.T @ density9).ravel()
        density_numerators = np.rint(self.density10 * 45).astype(np.int64)
        if not np.allclose(density_numerators / 45.0, self.density10, atol=2e-15):
            raise RuntimeError("order-10 edge densities are not integral over denominator 45")
        self.density_numerators = density_numerators
        self.mass10 = np.asarray(self.deletion.sum(axis=0)).ravel()
        if not np.allclose(self.mass10, 1.0, atol=2e-15):
            raise RuntimeError("deletion lift columns do not have mass one")

        public_x = np.asarray(public_dual["x"], dtype=float)
        self.public_q = np.maximum(public_x[: self.n10], 0.0)
        self.public_q /= self.public_q.sum()

        self.descriptors = []
        self.descriptor_keys = set()
        self.static_row_count = 6
        if args.resume and args.state.exists():
            with args.state.open("rb") as handle:
                checkpoint = pickle.load(handle)
            self.descriptors = checkpoint["descriptors"]
            if checkpoint.get("format") == "erdos23-fixed-gram-descriptors-v1":
                migrated = []
                for descriptor in self.descriptors:
                    if descriptor["kind"] != "k8":
                        migrated.append(descriptor)
                        continue
                    old_side_map = {
                        profile: int(side)
                        for profile, side in zip(
                            self.root_profiles[int(descriptor["root"])], descriptor["sides"]
                        )
                    }
                    migrated.append(
                        {
                            **descriptor,
                            "sides": tuple(
                                old_side_map.get(profile, 0)
                                for profile in self.edge_profiles[int(descriptor["root"])]
                            ),
                        }
                    )
                self.descriptors = migrated
                print("migrated legacy K8 descriptors to the exact edge-profile order", flush=True)
            elif checkpoint.get("edge_profiles") != self.edge_profiles:
                raise RuntimeError("checkpoint K8 profile order does not match regenerated data")
            self.descriptor_keys = {self.descriptor_key(desc) for desc in self.descriptors}
            print(
                f"resuming {len(self.descriptors)} exact row descriptors "
                f"from iteration {checkpoint.get('iteration')}",
                flush=True,
            )

        self.highs = None
        self.last_solve_optimal = False

    def build_pair_matrices(self, decomposition, profiles_by_root, profile_indices):
        rows = [[] for _ in range(self.n8)]
        columns = [[] for _ in range(self.n8)]
        for state, contributions in enumerate(decomposition):
            for root, profile_a, profile_b in contributions:
                profile_index = profile_indices[root]
                width = len(profile_index)
                pair_index = (
                    profile_index[tuple(profile_a)] * width
                    + profile_index[tuple(profile_b)]
                )
                rows[root].append(pair_index)
                columns[root].append(state)
        matrices = []
        for root in range(self.n8):
            width = len(profiles_by_root[root])
            values = np.ones(len(rows[root]), dtype=np.int16)
            matrix = coo_matrix(
                (values, (rows[root], columns[root])),
                shape=(width * width, self.n10),
                dtype=np.int16,
            ).tocsr()
            matrix.sum_duplicates()
            matrices.append(matrix)
        return matrices

    def sparse_row(self, entries):
        columns = np.fromiter(entries.keys(), dtype=np.int32, count=len(entries))
        values = np.fromiter(entries.values(), dtype=float, count=len(entries))
        return csr_matrix(
            (values, (np.zeros(len(columns), dtype=np.int32), columns)),
            shape=(1, self.nvars),
        )

    def descriptor_key(self, descriptor):
        kind = descriptor["kind"]
        if kind == "k7":
            return (kind, int(descriptor["root"]), tuple(descriptor["rule"]))
        if kind == "k8":
            return (kind, int(descriptor["root"]), tuple(descriptor["sides"]))
        return (kind, int(descriptor["root"]), tuple(descriptor["cycle"]))

    def k7_numerator(self, root, rule):
        _k, _adjacency, edge_tensor, root_mass, _classes = self.k7_types[root]
        rule = np.asarray(rule, dtype=np.int8)
        same_mask = np.equal.outer(rule, rule)
        same_mask = np.triu(same_mask)
        if self.k7_raw_counts:
            # Each state/root total is at most 9P7=181440, so uint32 is
            # exact here and avoids copying a whole root tensor to int64 for
            # every active rule.
            same_raw = np.tensordot(
                edge_tensor,
                same_mask.astype(np.uint32),
                axes=([1, 2], [0, 1]),
            ).astype(np.int64)
            root_raw = np.asarray(root_mass, dtype=np.int64)
        else:
            same_density = np.tensordot(
                edge_tensor, same_mask.astype(float), axes=([1, 2], [0, 1])
            )
            same_raw = np.rint(same_density * 181440).astype(np.int64)
            root_raw = np.rint(np.asarray(root_mass) * 181440).astype(np.int64)
            if not np.allclose(same_raw / 181440.0, same_density, atol=3e-13):
                raise RuntimeError(f"k7 edge counts failed exact recovery for root {root}")
            if not np.allclose(root_raw / 181440.0, root_mass, atol=3e-13):
                raise RuntimeError(f"k7 root counts failed exact recovery for root {root}")
        order9_numerator = 25 * same_raw - 2 * root_raw
        return np.asarray(self.deletion_count_t @ order9_numerator).ravel().astype(np.int64)

    def k8_numerator(self, root, sides):
        sides = np.asarray(sides, dtype=np.int8)
        coefficients = np.equal.outer(sides, sides).ravel().astype(np.int64)
        return np.asarray(
            self.edge_pair_matrices[root].T @ coefficients, dtype=np.int64
        ).ravel()

    def horn_numerator(self, root, cycle):
        cycle = tuple(tuple(profile) for profile in cycle)
        profile_index = self.profile_indices[root]
        width = len(profile_index)
        coefficients = np.zeros((width, width), dtype=np.int64)
        indices = [profile_index[profile] for profile in cycle]
        coefficients[np.ix_(indices, indices)] = 1
        for index in range(5):
            coefficients[indices[index], indices[(index + 1) % 5]] -= 4
        return np.asarray(
            self.all_pair_matrices[root].T @ coefficients.ravel(), dtype=np.int64
        ).ravel()

    def functional(self, descriptor):
        kind = descriptor["kind"]
        if kind == "k7":
            return self.k7_numerator(descriptor["root"], descriptor["rule"]), K7_DENOMINATOR
        if kind == "k8":
            return self.k8_numerator(descriptor["root"], descriptor["sides"]), 90
        return self.horn_numerator(descriptor["root"], descriptor["cycle"]), 90

    def descriptor_row(self, descriptor):
        numerator, denominator = self.functional(descriptor)
        nonzero = np.flatnonzero(numerator)
        entries = {int(index): -float(numerator[index]) / denominator for index in nonzero}
        if descriptor["kind"] == "k7":
            entries[self.u7_col + int(descriptor["root"])] = 1.0
        elif descriptor["kind"] == "k8":
            entries[self.u8_col + int(descriptor["root"])] = 1.0
        return self.sparse_row(entries)

    def add_descriptors(self, descriptors):
        new_descriptors = []
        new_rows = []
        for descriptor in descriptors:
            key = self.descriptor_key(descriptor)
            if key in self.descriptor_keys:
                continue
            self.descriptor_keys.add(key)
            self.descriptors.append(descriptor)
            new_descriptors.append(descriptor)
            if self.highs is not None:
                new_rows.append(self.descriptor_row(descriptor))
        if self.highs is not None and new_rows:
            self.add_rows(new_rows, [-INF] * len(new_rows), [0.0] * len(new_rows))
        return len(new_descriptors)

    def add_rows(self, rows, lower, upper):
        if not rows:
            return
        matrix = vstack(rows, format="csr")
        self.highs.addRows(
            matrix.shape[0],
            np.asarray(lower, dtype=float),
            np.asarray(upper, dtype=float),
            matrix.nnz,
            matrix.indptr[:-1].astype(np.int32),
            matrix.indices.astype(np.int32),
            matrix.data.astype(float),
        )

    def build_model(self):
        highs = highspy.Highs()
        highs.setOptionValue("output_flag", False)
        highs.setOptionValue("solver", self.args.solver)
        if self.args.solver == "ipm":
            highs.setOptionValue("presolve", "off")
            highs.setOptionValue("run_crossover", "on" if self.args.crossover else "off")
        else:
            highs.setOptionValue("simplex_strategy", 1)
        solve_tolerance = 1e-9
        highs.setOptionValue("primal_feasibility_tolerance", solve_tolerance)
        highs.setOptionValue("dual_feasibility_tolerance", solve_tolerance)
        highs.setOptionValue("ipm_optimality_tolerance", solve_tolerance)
        highs.setOptionValue("threads", max(1, self.args.threads))
        eta_lower = -INF if self.args.crossover else -1.0
        lower = np.asarray(
            [0.0] * self.n10 + [eta_lower] + [0.0] * (self.n7 + self.n8),
            dtype=float,
        )
        # In generation mode these redundant probability-scale bounds make a
        # heavily pruned IPM model much better conditioned.  They are removed
        # again for the final crossover/dual extraction.
        upper = np.asarray(
            [INF] * self.n10
            + [INF if self.args.crossover else 1.0]
            + [INF if self.args.crossover else 1.0] * (self.n7 + self.n8),
            dtype=float,
        )
        highs.addVars(self.nvars, lower, upper)
        objective = np.zeros(self.nvars)
        objective[self.eta_col] = -1.0
        highs.changeColsCost(
            self.nvars, np.arange(self.nvars, dtype=np.int32), objective
        )
        self.highs = highs

        static = [
            self.sparse_row(
                {int(index): value for index, value in enumerate(self.density10) if value}
            ),
            self.sparse_row(
                {int(index): -value for index, value in enumerate(self.density10) if value}
            ),
            self.sparse_row(
                {int(index): -value for index, value in enumerate(self.moment_float) if value}
            ),
            self.sparse_row(
                {self.eta_col: 1.0, **{self.u7_col + index: -1.0 for index in range(self.n7)}}
            ),
            self.sparse_row(
                {self.eta_col: 1.0, **{self.u8_col + index: -1.0 for index in range(self.n8)}}
            ),
            self.sparse_row(
                {int(index): value for index, value in enumerate(self.mass10) if value}
            ),
        ]
        self.add_rows(
            static,
            [-INF, -INF, -INF, -INF, -INF, 1.0],
            [HI, -LO, 0.0, 0.0, -2.0 / 25.0, 1.0],
        )
        if self.descriptors:
            rows = [self.descriptor_row(descriptor) for descriptor in self.descriptors]
            self.add_rows(rows, [-INF] * len(rows), [0.0] * len(rows))

    def solve(self):
        self.last_solve_optimal = False
        self.highs.run()
        status = self.highs.getModelStatus()

        def primal_payload(solution):
            lp = self.highs.getLp()
            column_values = np.asarray(solution.col_value, dtype=float)
            row_values = np.asarray(solution.row_value, dtype=float)
            column_lower = np.asarray(lp.col_lower_, dtype=float)
            column_upper = np.asarray(lp.col_upper_, dtype=float)
            row_lower = np.asarray(lp.row_lower_, dtype=float)
            row_upper = np.asarray(lp.row_upper_, dtype=float)
            violation = max(
                0.0,
                float(np.max(column_lower - column_values)),
                float(np.max(column_values - column_upper)),
                float(np.max(row_lower - row_values)),
                float(np.max(row_values - row_upper)),
            )
            return violation, (
                float(column_values[self.eta_col]),
                column_values[: self.n10],
                column_values[self.u7_col : self.u7_col + self.n7],
                column_values[self.u8_col : self.u8_col + self.n8],
            )

        if status not in (
            highspy.HighsModelStatus.kOptimal,
            highspy.HighsModelStatus.kObjectiveBound,
        ):
            use_crossover_fallback = self.args.solver == "ipm" and not self.args.crossover
            if use_crossover_fallback:
                interior_solution = self.highs.getSolution()
                primal_violation, payload = primal_payload(interior_solution)
                if interior_solution.value_valid and primal_violation <= 2e-7:
                    print(
                        f"accepted IPM primal despite status "
                        f"{self.highs.modelStatusToString(status)}; "
                        f"max violation {primal_violation:.2e}",
                        flush=True,
                    )
                    return payload
                print(
                    f"HiGHS IPM returned {self.highs.modelStatusToString(status)}; "
                    "crossing over its interior solution",
                    flush=True,
                )
                crossover_status = self.highs.crossover(interior_solution)
                status = self.highs.getModelStatus()
                crossed_solution = self.highs.getSolution()
                primal_violation, payload = primal_payload(crossed_solution)
                if (
                    crossed_solution.value_valid
                    and primal_violation <= 2e-7
                ):
                    print(
                        f"accepted crossed primal (call {crossover_status}) with "
                        f"max violation {primal_violation:.2e}",
                        flush=True,
                    )
                    return payload
                print(
                    f"direct crossover not accepted: call={crossover_status}, "
                    f"model={self.highs.modelStatusToString(status)}, "
                    f"value_valid={crossed_solution.value_valid}, "
                    f"primal_violation={primal_violation:.2e}",
                    flush=True,
                )
            action = "retrying with crossover" if use_crossover_fallback else "retrying"
            print(
                f"HiGHS first attempt returned {self.highs.modelStatusToString(status)}; "
                f"clearing solver state and {action}",
                flush=True,
            )
            self.highs.clearSolver()
            if use_crossover_fallback:
                self.highs.setOptionValue("run_crossover", "on")
            self.highs.run()
            status = self.highs.getModelStatus()
            if use_crossover_fallback:
                self.highs.setOptionValue("run_crossover", "off")
        if status not in (highspy.HighsModelStatus.kOptimal, highspy.HighsModelStatus.kObjectiveBound):
            raise RuntimeError(f"HiGHS status {self.highs.modelStatusToString(status)}")
        self.last_solve_optimal = True
        solution = self.highs.getSolution()
        columns = np.asarray(solution.col_value, dtype=float)
        eta = float(columns[self.eta_col])
        return (
            eta,
            columns[: self.n10],
            columns[self.u7_col : self.u7_col + self.n7],
            columns[self.u8_col : self.u8_col + self.n8],
        )

    def separate_k7(self, q, u7=None, force=False):
        order9_distribution = np.asarray(self.deletion @ q).ravel()

        def one(root, item):
            _k, _adjacency, edge_tensor, root_mass, _classes = item
            np.random.seed(1729 + root)
            rules, _values = self.sep_multi(
                edge_tensor,
                root_mass,
                order9_distribution,
                self.cache["t"],
                1e9,
                1e-10,
                restarts=self.args.k7_restarts,
                keep=self.args.k7_keep,
            )
            if not rules:
                return []
            descriptors = []
            for raw_rule in rules:
                rule = canonical_binary_rule(raw_rule)
                descriptor = {"kind": "k7", "root": root, "rule": rule}
                if force:
                    descriptors.append(descriptor)
                    continue
                numerator = self.k7_numerator(root, rule)
                value = float(numerator @ q) / K7_DENOMINATOR
                if u7[root] > value + 1e-9:
                    descriptors.append(descriptor)
            return descriptors

        results = Parallel(n_jobs=self.args.threads, prefer="threads")(
            delayed(one)(root, item) for root, item in enumerate(self.k7_types)
        )
        return [descriptor for descriptors in results for descriptor in descriptors]

    def k8_colorings(self, q):
        descriptors = []
        for root, matrix in enumerate(self.edge_pair_matrices):
            profiles = self.edge_profiles[root]
            width = len(profiles)
            pair_weights = np.asarray(matrix @ q).ravel().reshape(width, width) / 90.0
            off_diagonal = {}
            active = np.flatnonzero(pair_weights + pair_weights.T)
            for flat_index in active:
                index_a, index_b = divmod(int(flat_index), width)
                if index_a < index_b:
                    value = pair_weights[index_a, index_b] + pair_weights[index_b, index_a]
                    if value:
                        off_diagonal[(profiles[index_a], profiles[index_b])] = float(value)
            coloring = self.maxcut_coloring(list(profiles), off_diagonal)
            sides = tuple(int(coloring.get(profile, 0)) for profile in profiles)
            complement = tuple(1 - side for side in sides)
            sides = min(sides, complement)
            descriptors.append({"kind": "k8", "root": root, "sides": sides})
        return descriptors

    def separate_k8(self, q, u8=None, force=False):
        descriptors = self.k8_colorings(q)
        if force:
            return descriptors
        violated = []
        for descriptor in descriptors:
            numerator = self.k8_numerator(descriptor["root"], descriptor["sides"])
            value = float(numerator @ q) / 90.0
            if u8[descriptor["root"]] > value + 1e-9:
                violated.append(descriptor)
        return violated

    def separate_horn(self, q, keep=1):
        descriptors = []
        for root, pair_matrix in enumerate(self.all_pair_matrices):
            profiles = self.root_profiles[root]
            width = len(profiles)
            matrix = np.asarray(pair_matrix @ q).ravel().reshape(width, width) / 90.0
            if not np.any(matrix):
                continue
            matrix = 0.5 * (matrix + matrix.T)
            tuples = self.horn_tuples_for_R(matrix, KEEP=keep)
            for _value, indices in tuples:
                cycle = canonical_cycle(tuple(profiles[index] for index in indices))
                descriptors.append({"kind": "horn", "root": root, "cycle": cycle})
        return descriptors

    def checkpoint(self, iteration, eta):
        payload = {
            "format": "erdos23-fixed-gram-descriptors-v2",
            "iteration": iteration,
            "eta": eta,
            "solve_optimal": self.last_solve_optimal,
            "descriptors": self.descriptors,
            "edge_profiles": self.edge_profiles,
        }
        with self.args.state.open("wb") as handle:
            pickle.dump(payload, handle, protocol=4)

    def separation_distribution(self, q):
        q = np.maximum(np.asarray(q, dtype=float), 0.0)
        total = float(q.sum())
        if total <= 0:
            raise RuntimeError("LP returned a zero primal distribution")
        q /= total
        if np.count_nonzero(q > 1e-12) <= 800:
            return q
        order = np.argsort(q)[::-1]
        cumulative = np.cumsum(q[order])
        keep_count = int(np.searchsorted(cumulative, 1.0 - 1e-4)) + 1
        threshold = q[order[min(keep_count, self.n10 - 1)]]
        sparse_q = np.where(q >= threshold, q, 0.0)
        sparse_q /= sparse_q.sum()
        return sparse_q

    def prune_pool(self):
        if len(self.descriptors) <= self.args.pool_cap:
            return 0
        solution = self.highs.getSolution()
        row_values = np.asarray(solution.row_value, dtype=float)
        descriptor_values = row_values[
            self.static_row_count : self.static_row_count + len(self.descriptors)
        ]
        keep_mask = descriptor_values >= -self.args.slack_drop
        drop_positions = np.flatnonzero(~keep_mask)
        if not len(drop_positions):
            return 0
        self.descriptors = [
            descriptor
            for descriptor, keep in zip(self.descriptors, keep_mask)
            if keep
        ]
        self.descriptor_keys = {
            self.descriptor_key(descriptor) for descriptor in self.descriptors
        }
        # A fresh IPM model is more reliable than reusing its internal state after
        # row deletion followed by row addition.  Rebuilding costs only a few seconds.
        self.highs = None
        return len(drop_positions)

    def run(self):
        if not self.descriptors:
            seeds = [np.ones(self.n10) / self.n10, self.public_q]
            for seed_number, seed in enumerate(seeds, start=1):
                started = time.time()
                added7 = self.add_descriptors(self.separate_k7(seed, force=True))
                added8 = self.add_descriptors(self.separate_k8(seed, force=True))
                added_horn = self.add_descriptors(self.separate_horn(seed, keep=2))
                print(
                    f"seed {seed_number}: +{added7} k7 +{added8} k8 "
                    f"+{added_horn} Horn; pool={len(self.descriptors)}; "
                    f"elapsed={time.time() - started:.1f}s",
                    flush=True,
                )
        self.build_model()
        eta, q, u7, u8 = self.solve()
        initial_pruned = (
            self.prune_pool()
            if self.args.max_iterations > 0 and self.last_solve_optimal
            else 0
        )
        initial_quality = "optimal" if self.last_solve_optimal else "feasible-only"
        initial_bound = (
            math.floor(math.sqrt(2 / (25 * eta)))
            if self.last_solve_optimal and eta > 0
            else ("closed" if self.last_solve_optimal else "unavailable")
        )
        print(
            f"initial solve ({initial_quality}): eta={eta:+.10e}; "
            f"pool={len(self.descriptors)}; "
            f"-{initial_pruned} slack; "
            f"n-bound={initial_bound}",
            flush=True,
        )
        iteration = 0
        model_dirty = bool(initial_pruned)
        for iteration in range(1, self.args.max_iterations + 1):
            started = time.time()
            separation_q = self.separation_distribution(q)
            descriptors7 = self.separate_k7(separation_q, u7=u7)
            descriptors8 = self.separate_k8(separation_q, u8=u8)
            descriptors_horn = self.separate_horn(
                separation_q, keep=self.args.horn_keep
            )
            added7 = self.add_descriptors(descriptors7)
            added8 = self.add_descriptors(descriptors8)
            added_horn = self.add_descriptors(descriptors_horn)
            if added7 + added8 + added_horn == 0:
                print(f"converged at iteration {iteration}", flush=True)
                break
            if self.highs is None:
                self.build_model()
            eta, q, u7, u8 = self.solve()
            target_reached = self.last_solve_optimal and eta < 4.74e-5
            pruned = (
                0
                if target_reached or not self.last_solve_optimal
                else self.prune_pool()
            )
            model_dirty = bool(pruned)
            self.checkpoint(iteration, eta)
            solve_quality = "optimal" if self.last_solve_optimal else "feasible-only"
            print(
                f"iteration {iteration}: +{added7} k7 +{added8} k8 +{added_horn} Horn; "
                f"-{pruned} slack; pool={len(self.descriptors)}; "
                f"eta={eta:+.10e} ({solve_quality}); "
                f"elapsed={time.time() - started:.1f}s",
                flush=True,
            )
            if target_reached:
                print("target float margin reached; stopping generation", flush=True)
                break
        if model_dirty:
            if self.highs is None:
                self.build_model()
            eta, q, u7, u8 = self.solve()
            model_dirty = False
            print(
                f"post-prune solve: eta={eta:+.10e}; pool={len(self.descriptors)}",
                flush=True,
            )
        self.checkpoint(iteration, eta)
        row_duals = np.asarray(self.highs.getSolution().row_dual, dtype=float)
        dual_payload = {
            "format": "erdos23-fixed-gram-dual-v1",
            "eta": eta,
            "solve_optimal": self.last_solve_optimal,
            "row_duals": row_duals,
            "descriptors": self.descriptors,
            "static_row_count": self.static_row_count,
        }
        dual_path = self.args.state.with_name(self.args.state.stem + "_dual.pkl")
        with dual_path.open("wb") as handle:
            pickle.dump(dual_payload, handle, protocol=4)
        print(f"saved {self.args.state} and {dual_path}", flush=True)


def main() -> int:
    args = parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))
    os.environ.setdefault("OPENBLAS_NUM_THREADS", str(args.threads))
    generator = CertificateGenerator(args)
    generator.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
