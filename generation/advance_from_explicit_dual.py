#!/usr/bin/env python3
"""Add one exact cutting-plane batch from an optimal explicit-dual solution."""

from __future__ import annotations

import argparse
import importlib.util
import pickle
import time
from pathlib import Path

import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--public-anc", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--dual", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--k7-restarts", type=int, default=64)
    parser.add_argument("--k7-keep", type=int, default=1)
    parser.add_argument("--horn-keep", type=int, default=1)
    parser.add_argument("--slack-drop", type=float, default=1e-5)
    return parser.parse_args()


def main():
    args = parse_args()
    spec = importlib.util.spec_from_file_location("fixed_gram_generator", args.generator)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    namespace = argparse.Namespace(
        flagsdp=args.flagsdp,
        public_anc=args.public_anc,
        state=args.state,
        max_iterations=0,
        threads=args.threads,
        k7_restarts=args.k7_restarts,
        k7_keep=args.k7_keep,
        horn_keep=args.horn_keep,
        pool_cap=700,
        slack_drop=1e-5,
        resume=True,
        crossover=False,
        solver="ipm",
    )
    generator = module.CertificateGenerator(namespace)
    with args.dual.open("rb") as handle:
        dual = pickle.load(handle)
    dual_violation = max(
        float(dual.get("lower_violation", float("inf"))),
        float(dual.get("upper_violation", float("inf"))),
        float(dual.get("multiplier_violation", float("inf"))),
        float(dual.get("equality_error", float("inf"))),
    )
    if dual_violation > 1e-7:
        raise RuntimeError(
            f"explicit dual is not sufficiently feasible: violation={dual_violation}"
        )
    if len(dual["descriptors"]) != len(generator.descriptors):
        raise RuntimeError("dual and row-pool descriptor counts differ")
    if any(
        generator.descriptor_key(left) != generator.descriptor_key(right)
        for left, right in zip(dual["descriptors"], generator.descriptors)
    ):
        raise RuntimeError("dual and row-pool descriptor order differs")

    row_duals = np.asarray(dual["row_duals"], dtype=float)
    activities = np.asarray(dual["primal_row_activities"], dtype=float)
    primal_upper = np.concatenate(
        [np.asarray([0.3197, -0.2486, 0.0, 0.0, -0.08, 1.0]),
         np.zeros(len(activities) - generator.static_row_count)]
    )
    primal_violation = max(
        0.0,
        float(np.max(activities - primal_upper)),
        abs(float(activities[generator.static_row_count - 1]) - 1.0),
        float(-np.min(row_duals[: generator.eta_col])),
        float(-np.min(row_duals[generator.eta_col + 1 :])),
    )
    if primal_violation > 1e-7:
        raise RuntimeError(
            f"recovered primal is not sufficiently feasible: violation={primal_violation}"
        )
    q = np.maximum(row_duals[: generator.n10], 0.0)
    q /= q.sum()
    u7 = np.maximum(
        row_duals[generator.u7_col : generator.u7_col + generator.n7], 0.0
    )
    u8 = np.maximum(
        row_duals[generator.u8_col : generator.u8_col + generator.n8], 0.0
    )
    separation_q = generator.separation_distribution(q)

    pruned = 0
    if "primal_row_activities" in dual:
        descriptor_values = activities[
            generator.static_row_count : generator.static_row_count
            + len(generator.descriptors)
        ]
        keep_mask = descriptor_values >= -args.slack_drop
        pruned = int(np.count_nonzero(~keep_mask))
        if pruned:
            generator.descriptors = [
                descriptor
                for descriptor, keep in zip(generator.descriptors, keep_mask)
                if keep
            ]
            generator.descriptor_keys = {
                generator.descriptor_key(descriptor)
                for descriptor in generator.descriptors
            }

    started = time.time()
    descriptors7 = generator.separate_k7(separation_q, u7=u7)
    descriptors8 = generator.separate_k8(separation_q, u8=u8)
    descriptors_horn = generator.separate_horn(
        separation_q, keep=args.horn_keep
    )
    added7 = generator.add_descriptors(descriptors7)
    added8 = generator.add_descriptors(descriptors8)
    added_horn = generator.add_descriptors(descriptors_horn)

    with args.state.open("rb") as handle:
        previous = pickle.load(handle)
    generator.last_solve_optimal = False
    generator.checkpoint(int(previous.get("iteration", 0)) + 1, float(dual["objective"]))
    print(
        f"explicit-dual status={dual['model_status']} objective={dual['objective']:+.12e}; "
        f"dual/primal violations={dual_violation:.1e}/{primal_violation:.1e}; "
        f"q support={np.count_nonzero(q > 1e-12)}; "
        f"separation support={np.count_nonzero(separation_q)}; "
        f"+{added7} k7 +{added8} k8 +{added_horn} Horn; -{pruned} slack; "
        f"pool={len(generator.descriptors)}; elapsed={time.time() - started:.1f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
