#!/usr/bin/env python3
"""Solve the explicit dual of the fixed-Gram Erdős 23 cutting-plane LP.

The primal basis is expensive and unnecessary for certification.  This model
has one nonnegative multiplier for every upper row and one bounded free mass
multiplier rho.  A slightly inaccurate solution is still useful: the exact
verifier rationalizes it, repairs the two envelope legs, and pays the exact
minimum coordinate residual in the final objective.
"""

from __future__ import annotations

import argparse
import importlib.util
import pickle
import time
from pathlib import Path

import highspy
import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, hstack


INF = highspy.kHighsInf


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--public-anc", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--multiplier-upper", type=float, default=1000.0)
    parser.add_argument("--solver", choices=("ipm", "simplex"), default="ipm")
    parser.add_argument("--presolve", choices=("off", "on", "choose"), default="off")
    parser.add_argument("--crossover", action="store_true")
    parser.add_argument("--log", action="store_true")
    return parser.parse_args()


def load_generator(args):
    spec = importlib.util.spec_from_file_location("fixed_gram_generator", args.generator)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    namespace = argparse.Namespace(
        flagsdp=args.flagsdp,
        public_anc=args.public_anc,
        state=args.state,
        max_iterations=0,
        threads=args.threads,
        k7_restarts=1,
        k7_keep=1,
        horn_keep=1,
        pool_cap=700,
        slack_drop=1e-5,
        resume=True,
        crossover=True,
        solver="ipm",
    )
    generator = module.CertificateGenerator(namespace)
    generator.build_model()
    return generator


def highs_matrix(lp):
    matrix = lp.a_matrix_
    starts = np.asarray(matrix.start_, dtype=np.int64)
    indices = np.asarray(matrix.index_, dtype=np.int32)
    values = np.asarray(matrix.value_, dtype=float)
    shape = (lp.num_row_, lp.num_col_)
    if matrix.format_ == highspy.MatrixFormat.kRowwise:
        return csr_matrix((values, indices, starts), shape=shape)
    return csc_matrix((values, indices, starts), shape=shape).tocsr()


def main():
    args = parse_args()
    started = time.time()
    generator = load_generator(args)
    primal_lp = generator.highs.getLp()
    primal_matrix = highs_matrix(primal_lp)

    equality_row = generator.static_row_count - 1
    inequality_rows = np.concatenate(
        [np.arange(equality_row), np.arange(equality_row + 1, primal_lp.num_row_)]
    )
    upper_matrix = primal_matrix[inequality_rows]
    upper_rhs = np.asarray(primal_lp.row_upper_, dtype=float)[inequality_rows]
    mass_column = primal_matrix[equality_row].T.tocsr()

    # Variables are (all upper-row multipliers y >= 0, rho in [-1,1]).
    dual_matrix = hstack([upper_matrix.T, mass_column], format="csr")
    n_upper = upper_matrix.shape[0]
    n_dual = n_upper + 1
    rho_col = n_upper

    highs = highspy.Highs()
    highs.setOptionValue("output_flag", bool(args.log))
    highs.setOptionValue("log_to_console", bool(args.log))
    highs.setOptionValue("solver", args.solver)
    highs.setOptionValue("presolve", args.presolve)
    highs.setOptionValue("run_crossover", "on" if args.crossover else "off")
    if args.solver == "simplex":
        highs.setOptionValue("simplex_strategy", 1)
    highs.setOptionValue("threads", max(1, args.threads))
    highs.setOptionValue("primal_feasibility_tolerance", 1e-9)
    highs.setOptionValue("dual_feasibility_tolerance", 1e-9)
    highs.setOptionValue("ipm_optimality_tolerance", 1e-9)

    lower = np.asarray([0.0] * n_upper + [-1.0], dtype=float)
    upper = np.asarray(
        [args.multiplier_upper] * n_upper + [1.0], dtype=float
    )
    highs.addVars(n_dual, lower, upper)
    objective = np.concatenate([upper_rhs, [1.0]])
    highs.changeColsCost(
        n_dual, np.arange(n_dual, dtype=np.int32), objective
    )

    constraint_lower = np.zeros(primal_lp.num_col_, dtype=float)
    constraint_upper = np.full(primal_lp.num_col_, INF, dtype=float)
    constraint_lower[generator.eta_col] = 1.0
    constraint_upper[generator.eta_col] = 1.0
    highs.addRows(
        primal_lp.num_col_,
        constraint_lower,
        constraint_upper,
        dual_matrix.nnz,
        dual_matrix.indptr[:-1].astype(np.int32),
        dual_matrix.indices.astype(np.int32),
        dual_matrix.data.astype(float),
    )

    print(
        f"explicit dual: variables={n_dual}, constraints={primal_lp.num_col_}, "
        f"nnz={dual_matrix.nnz}; build={time.time() - started:.1f}s",
        flush=True,
    )
    solve_started = time.time()
    highs.run()
    status = highs.getModelStatus()
    solution = highs.getSolution()
    values = np.asarray(solution.col_value, dtype=float)
    row_duals = np.asarray(solution.row_dual, dtype=float)
    primal_row_activities = np.asarray(primal_matrix @ row_duals).ravel()
    activities = np.asarray(dual_matrix @ values).ravel()
    lower_violation = float(np.max(constraint_lower - activities))
    upper_violation = float(np.max(activities - constraint_upper))
    multiplier_violation = max(0.0, float(-np.min(values[:n_upper])))
    equality_error = abs(float(activities[generator.eta_col]) - 1.0)
    objective_value = float(objective @ values)
    info = highs.getInfo()
    print(
        f"status={highs.modelStatusToString(status)}; objective={objective_value:+.12e}; "
        f"lower_violation={max(0.0, lower_violation):.3e}; "
        f"upper_violation={max(0.0, upper_violation):.3e}; "
        f"multiplier_violation={multiplier_violation:.3e}; "
        f"eta_equality_error={equality_error:.3e}; solve={time.time() - solve_started:.1f}s",
        flush=True,
    )

    payload = {
        "format": "erdos23-fixed-gram-explicit-dual-v1",
        "model_status": highs.modelStatusToString(status),
        "objective": objective_value,
        "values": values,
        "row_duals": row_duals,
        "primal_row_activities": primal_row_activities,
        "n_upper": n_upper,
        "rho_col": rho_col,
        "inequality_rows": inequality_rows,
        "descriptors": generator.descriptors,
        "static_row_count": generator.static_row_count,
        "lower_violation": max(0.0, lower_violation),
        "upper_violation": max(0.0, upper_violation),
        "multiplier_violation": multiplier_violation,
        "equality_error": equality_error,
        "highs_info": {
            "primal_solution_status": int(info.primal_solution_status),
            "dual_solution_status": int(info.dual_solution_status),
            "ipm_iteration_count": int(info.ipm_iteration_count),
            "primal_dual_objective_error": float(
                getattr(info, "primal_dual_objective_error", float("nan"))
            ),
        },
    }
    with args.output.open("wb") as handle:
        pickle.dump(payload, handle, protocol=4)
    print(f"saved {args.output}", flush=True)


if __name__ == "__main__":
    main()
