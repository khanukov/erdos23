#!/usr/bin/env python3
"""The band SDP on the rebuilt 8-root rows, solved with SCS.

Same primal as `build_u8_lp.py` -- maximise `eta` over order-10 state vectors
`q` in the band with `eta <= -2/25 + sum_tau u_tau`, `u_tau >= <h_{tau,c}, q>`
for every K8 MaxCut rule in the pool -- but the positive-semidefiniteness of
the moment blocks is imposed exactly instead of through eigenvector cuts:

  * the 8-root pair matrices `M_tau(q)` (type size 8, widths up to 256),
  * the order-10 blocks of type size 0 and 2 (widths 14, 245, 135), read from
    the explicit triplets written by `blocks10 triplets`.

The order-10 blocks of type size 4 and 6 (widths up to 2,445) stay as
eigenvector cuts, as do the rooted-Horn rows and the rules, all separated at
the SDP optimum in a cutting-plane loop.  Every solve's dual is written in the
format `verify_u8_certificate.py` reads: the PSD duals are eigen-decomposed
into `rpsd` / `b10` rows with their eigenvalues as multipliers, so the exact
rational check needs nothing new.
"""

from __future__ import annotations

import argparse
import pickle
import random
import struct
import sys
import time
from pathlib import Path

import numpy as np
import scs
from scipy.sparse import coo_matrix, csc_matrix, csr_matrix, vstack

sys.path.insert(0, str(Path(__file__).resolve().parent))
from blocks10 import SCALE, Blocks10  # noqa: E402
from build_u8_lp import N9, N10, MomentBlocks, RootData, q10_of  # noqa: E402
from k7_leg_ceiling import cycle, grotzsch, petersen  # noqa: E402
from rebuild_u8 import Catalogue  # noqa: E402

HI, LO = 0.3197, 0.2486
SQ2 = np.sqrt(2.0)


def svec_index(n, i, j):
    """SCS packs the lower triangle column by column; entry (i >= j)."""
    return j * n - j * (j - 1) // 2 + (i - j)


def read_triplets(path):
    data = Path(path).read_bytes()
    pos = 0
    (n_levels,) = struct.unpack_from("i", data, pos)
    pos += 4
    blocks = []
    for _ in range(n_levels):
        s, n_types = struct.unpack_from("ii", data, pos)
        pos += 8
        for t in range(n_types):
            nF, = struct.unpack_from("i", data, pos)
            pos += 4
            n, = struct.unpack_from("q", data, pos)
            pos += 8
            arrays = []
            for _k in range(4):
                arrays.append(np.frombuffer(data, dtype=np.int32, count=n, offset=pos))
                pos += 4 * n
            blocks.append(((s, t), nF, arrays))
    return blocks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--public-anc", type=Path, required=True)
    parser.add_argument("--rebuilt", type=Path, required=True)
    parser.add_argument("--allpairs", type=Path, required=True)
    parser.add_argument("--triplets", type=Path, required=True, help="blocks10 triplets output for levels 02")
    parser.add_argument("--work", type=Path, default=Path(".work"))
    parser.add_argument("--pool", type=Path, default=None, help="LP checkpoint whose rows seed the pool")
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument("--starts", type=int, default=6)
    parser.add_argument("--horn-restarts", type=int, default=12)
    parser.add_argument("--horn-per-root", type=int, default=3)
    parser.add_argument("--max-rules", type=int, default=200)
    parser.add_argument("--max-horn", type=int, default=200)
    parser.add_argument("--max-b10", type=int, default=150)
    parser.add_argument("--b10-per-block", type=int, default=4)
    parser.add_argument("--eps", type=float, default=1e-6)
    parser.add_argument("--max-iters", type=int, default=50000)
    parser.add_argument("--indirect", action="store_true")
    parser.add_argument("--time-limit", type=float, default=20000.0)
    parser.add_argument("--no-root-blocks", action="store_true")
    parser.add_argument("--checkpoint", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    started = time.time()

    with args.rebuilt.open("rb") as handle:
        payload = pickle.load(handle)
    states10, roots8, edge_pairs = payload["states10"], payload["roots8"], payload["pairs"]
    with args.allpairs.open("rb") as handle:
        all_pairs = pickle.load(handle)["pairs"]
    with (args.flagsdp / "cache_n9.pkl").open("rb") as handle:
        cache9 = pickle.load(handle)
    lift = np.load(args.flagsdp / "c5lift_cache.npz", allow_pickle=True)
    moments = MomentBlocks(cache9, lift)
    cat9 = Catalogue([list(masks) for _o, masks in cache9["states"]])
    counts = np.rint(np.asarray(lift["Dval"]) * 10).astype(np.int64)
    columns = csr_matrix((counts, (lift["Drow"], lift["Dcol"])), shape=(N9, N10), dtype=np.int64).tocsc()
    profile_to_state = {}
    for j in range(N10):
        s0, e0 = columns.indptr[j], columns.indptr[j + 1]
        profile_to_state[tuple(sorted(zip(columns.indices[s0:e0].tolist(),
                                          columns.data[s0:e0].tolist())))] = j
    d_edge = np.asarray([sum(bin(m).count("1") for m in adj) / 2 / 45.0 for adj in states10])
    with (args.public_anc / "mom_term_exact.pkl").open("rb") as handle:
        moment = np.asarray([n / d for n, d in pickle.load(handle)], dtype=float)

    edge_roots = [RootData(tau, roots8[tau], edge_pairs[tau]) for tau in range(410) if edge_pairs[tau]]
    all_roots = [RootData(tau, roots8[tau], all_pairs[tau]) for tau in range(410) if all_pairs[tau]]
    edge_by_tau = {r.tau: r for r in edge_roots}
    all_by_tau = {r.tau: r for r in all_roots}
    blocks46 = Blocks10(states10, args.work, "46", verbose=False)
    print(f"root data: {len(edge_roots)} rule roots, {len(all_roots)} pair roots; "
          f"order-10 cut blocks (4,6): {len(blocks46.keys)}  [{time.time() - started:.0f}s]", flush=True)

    seeds = [("C5", cycle(5), np.full(5, 0.2)),
             ("Grotzsch", grotzsch(), np.asarray(
                 [10173, 9717, 8691, 10166, 7628, 8213, 7344, 9352, 12084, 10443, 6189], float)),
             ("Petersen", petersen(), np.full(10, 0.1))]
    seed_vectors = [(name, q10_of(cat9, profile_to_state, base, w / w.sum())) for name, base, w in seeds]
    rng = random.Random(2026)

    # ---- columns -------------------------------------------------------------
    col_eta = N10
    col_u = {r.tau: N10 + 1 + i for i, r in enumerate(edge_roots)}
    ncols = N10 + 1 + len(edge_roots)
    c = np.zeros(ncols)
    c[col_eta] = -1.0

    # ---- explicit PSD blocks (built once) ----------------------------------------
    psd_blocks = []          # (label, n, A_block (svec rows x ncols) as coo parts)

    def add_block(label, n, rows, cols, vals):
        psd_blocks.append((label, n, np.asarray(rows), np.asarray(cols), np.asarray(vals, dtype=float)))

    for key, nF, (st, f1, f2, cnt) in read_triplets(args.triplets):
        s = key[0]
        norm = float(blocks46.normaliser.get(s) or {0: 252.0, 2: 6300.0}[s])
        idx = svec_index(nF, f1.astype(np.int64), f2.astype(np.int64))
        scale = np.where(f1 > f2, SQ2, 1.0)
        add_block(("b10", key), nF, idx, st.astype(np.int64), -cnt.astype(float) / norm * scale)
    if not args.no_root_blocks:
        for r in all_roots:
            n = r.width
            rows, cols, vals = [], [], []
            N = r.N.tocsr()
            for o, m in enumerate(r.members):
                a, b = m // n, m % n
                keep = a >= b
                if not keep.any():
                    continue
                idx = svec_index(n, a[keep], b[keep])
                scale = np.where(a[keep] > b[keep], SQ2, 1.0)
                sl = slice(N.indptr[o], N.indptr[o + 1])
                st, v = N.indices[sl], N.data[sl] / (90.0 * r.sizes[o])
                rows.append(np.repeat(idx, len(st)))
                cols.append(np.tile(st, len(idx)))
                vals.append(-np.outer(scale, v).ravel())
            add_block(("rpsd", r.tau), n, np.concatenate(rows), np.concatenate(cols), np.concatenate(vals))
    psd_dims = [n for _l, n, _r, _c, _v in psd_blocks]
    psd_rows = sum(n * (n + 1) // 2 for n in psd_dims)
    print(f"explicit PSD blocks: {len(psd_blocks)}, svec rows {psd_rows}, "
          f"nonzeros {sum(len(v) for *_r, v in psd_blocks)}  [{time.time() - started:.0f}s]", flush=True)

    # ---- linear rows ---------------------------------------------------------------
    # each entry: (kind, a, b, indices, values, rhs); meaning  values . x <= rhs  (nonneg cone)
    lin_rows = []

    def add_lin(kind, a, b, indices, values, rhs):
        lin_rows.append((kind, a, b, np.asarray(indices, dtype=np.int64), np.asarray(values, dtype=float), float(rhs)))

    qcols = np.arange(N10)
    add_lin("band_hi", None, None, qcols, d_edge, HI)
    add_lin("band_lo", None, None, qcols, -d_edge, -LO)
    add_lin("gram", None, None, qcols, -moment, 0.0)
    add_lin("leg8", None, None, [col_eta] + [col_u[r.tau] for r in edge_roots],
            [1.0] + [-1.0] * len(edge_roots), -2.0 / 25.0)
    n_static = 4
    pool = {r.tau: [] for r in edge_roots}

    def add_rule(r, colouring):
        row = r.rule_row(colouring)
        nz = np.nonzero(row)[0]
        add_lin("rule", r.tau, np.asarray(colouring, dtype=np.int8),
                np.concatenate([nz, [col_u[r.tau]]]), np.concatenate([-row[nz], [1.0]]), 0.0)
        pool[r.tau].append(np.asarray(colouring, dtype=np.int8))

    def add_ge_row(kind, a, b, row):
        nz = np.nonzero(row)[0]
        add_lin(kind, a, b, nz, -row[nz], 0.0)

    def gate(row, label):
        for name, qs in seed_vectors:
            check = float(row @ qs)
            if check < -1e-9:
                raise RuntimeError(f"INVALID row {label}: {check:+.3e} at {name}")

    if args.pool:
        with args.pool.open("rb") as handle:
            saved = pickle.load(handle)
        pending = []
        for kind, a, b in saved["rows"]:
            if kind == "rule":
                add_rule(edge_by_tau[a], np.asarray(b, dtype=np.int8))
            elif kind == "horn":
                add_ge_row("horn", a, list(b), all_by_tau[a].form_row(all_by_tau[a].horn_coefficients(list(b))))
            elif kind == "psd":
                add_ge_row("psd", a, (b[0], np.asarray(b[1])), moments.cut_row(b[0], np.asarray(b[1])))
            elif kind == "b10" and tuple(a)[0] in (4, 6):
                pending.append((tuple(a), np.asarray(b, dtype=np.int64)))
        if pending:
            for (key, vec), row in zip(pending, blocks46.float_rows(pending)):
                add_ge_row("b10", key, vec, row)
        print(f"pool from {args.pool}: {len(lin_rows) - n_static} rows  [{time.time() - started:.0f}s]", flush=True)
    for r in edge_roots:
        if not pool[r.tau]:
            add_rule(r, np.zeros(r.width, dtype=np.int8))

    # ---- assemble and solve ----------------------------------------------------------
    def assemble():
        # zero cone: mass; nonneg: q >= 0, u >= 0, linear rows; PSD blocks
        rows_i, cols_i, vals_i = [np.zeros(N10, dtype=np.int64)], [qcols], [np.ones(N10)]
        b = [1.0]
        r0 = 1
        # q >= 0  ->  -q + s = 0
        rows_i.append(r0 + qcols)
        cols_i.append(qcols)
        vals_i.append(-np.ones(N10))
        b.extend([0.0] * N10)
        r0 += N10
        ucols = np.arange(N10 + 1, ncols)
        rows_i.append(r0 + np.arange(len(ucols)))
        cols_i.append(ucols)
        vals_i.append(-np.ones(len(ucols)))
        b.extend([0.0] * len(ucols))
        r0 += len(ucols)
        lin_start = r0
        for _k, _a, _b, idx, val, rhs in lin_rows:
            rows_i.append(np.full(len(idx), r0))
            cols_i.append(idx)
            vals_i.append(val)
            b.append(rhs)
            r0 += 1
        n_nonneg = r0 - 1
        psd_start = r0
        for _label, n, br, bc, bv in psd_blocks:
            rows_i.append(r0 + br)
            cols_i.append(bc)
            vals_i.append(bv)
            b.extend([0.0] * (n * (n + 1) // 2))
            r0 += n * (n + 1) // 2
        A = coo_matrix((np.concatenate(vals_i), (np.concatenate(rows_i), np.concatenate(cols_i))),
                       shape=(r0, ncols)).tocsc()
        cone = {"z": 1, "l": n_nonneg, "s": psd_dims}
        return A, np.asarray(b), cone, lin_start, psd_start

    def certificate(y, x, lin_start, psd_start, eta):
        """rows_log + duals in the verifier's (HiGHS) sign convention."""
        rows_log = [("band", None, None), ("gram", None, None), ("mass", None, None), ("leg8", None, None)]
        y_hi, y_lo = max(0.0, y[lin_start]), max(0.0, y[lin_start + 1])
        y_g, y_L = max(0.0, y[lin_start + 2]), max(0.0, y[lin_start + 3])
        duals = [-(y_hi - y_lo), y_g, -float(y[0]), -y_L]
        for k, (kind, a, b, _i, _v, _r) in enumerate(lin_rows[n_static:]):
            mult = max(0.0, float(y[lin_start + n_static + k]))
            rows_log.append((kind, a, b))
            duals.append(-mult if kind == "rule" else mult)
        pos = psd_start
        n_eig = 0
        for label, n, *_rest in psd_blocks:
            m = n * (n + 1) // 2
            Y = np.zeros((n, n))
            iu = np.tril_indices(n)
            # svec is column-major lower triangle: reorder
            order = svec_index(n, iu[0], iu[1])
            vals = np.asarray(y[pos:pos + m])[order]
            Y[iu] = vals
            Y = Y + Y.T - np.diag(np.diag(Y))
            off = ~np.eye(n, dtype=bool)
            Y[off] /= SQ2
            pos += m
            w, V = np.linalg.eigh((Y + Y.T) / 2)
            for j in range(n):
                if w[j] > 1e-12 * max(1.0, w[-1]):
                    v = V[:, j]
                    if label[0] == "rpsd":
                        rows_log.append(("rpsd", label[1], v.copy()))
                    else:
                        rows_log.append(("b10", label[1], np.rint(v * SCALE).astype(np.int64)))
                    duals.append(float(w[j]))
                    n_eig += 1
        return {"eta": eta, "q": x[:N10], "u": {r.tau: float(x[col_u[r.tau]]) for r in edge_roots},
                "duals": np.asarray(duals), "rows": rows_log, "col_u": col_u, "n_eig": n_eig}

    warm = None
    history = []
    for it in range(args.iterations):
        A, b, cone, lin_start, psd_start = assemble()
        t0 = time.time()
        solver = scs.SCS({"A": A, "b": b, "c": c}, cone, eps_abs=args.eps, eps_rel=args.eps,
                         max_iters=args.max_iters, verbose=True, use_indirect=args.indirect,
                         normalize=True, acceleration_lookback=10)
        if warm is not None and len(warm["y"]) == A.shape[0]:
            sol = solver.solve(warm_start=True, x=warm["x"], y=warm["y"], s=warm["s"])
        else:
            sol = solver.solve()
        info = sol["info"]
        x, y = np.asarray(sol["x"]), np.asarray(sol["y"])
        eta = float(x[col_eta])
        q = np.maximum(x[:N10], 0.0)
        bound = float(b @ y)      # dual objective: eta <= b.y up to the residuals
        print(f"iter {it:3d}: {info['status']} iters={info['iter']} eta={eta:+.6e} dual={bound:+.6e} "
              f"pres={info['res_pri']:.1e} dres={info['res_dual']:.1e} gap={info['gap']:.1e} "
              f"rows={A.shape[0]} band={float(d_edge @ q):.4f}  [{time.time() - t0:.0f}s solve]", flush=True)
        cert = certificate(y, x, lin_start, psd_start, eta)
        if args.checkpoint:
            with args.checkpoint.open("wb") as handle:
                pickle.dump(cert, handle)
        # warm start for the next round: extend y and s with zeros for new rows
        warm = {"x": x, "y": y, "s": np.asarray(sol["s"])}
        # ---- separation ----
        candidates = []
        for r in edge_roots:
            u_star = float(x[col_u[r.tau]])
            colouring, value = r.best_rule(q, rng, args.starts)
            if value < u_star - 1e-7:
                candidates.append((u_star - value, r, colouring))
        candidates.sort(key=lambda t: -t[0])
        added = 0
        for _viol, r, colouring in candidates[:args.max_rules]:
            add_rule(r, colouring)
            added += 1
        worst = candidates[0][0] if candidates else 0.0
        horn_candidates = []
        for r in all_roots:
            val, cyc = r.horn_search(q, rng, args.horn_restarts)
            if cyc is None:
                continue
            for val, cyc in r.last_cycles[:args.horn_per_root]:
                horn_candidates.append((val, r, cyc))
        horn_candidates.sort(key=lambda t: t[0])
        horns = 0
        for val, r, cyc in horn_candidates[:args.max_horn]:
            row = r.form_row(r.horn_coefficients(cyc))
            gate(row, f"horn {r.tau}")
            add_ge_row("horn", r.tau, list(cyc), row)
            horns += 1
        horn_worst = horn_candidates[0][0] if horn_candidates else 0.0
        b10_cuts = 0
        b10_worst = 0.0
        for key, lam, row, vec in sorted(blocks46.separate(q, per_block=args.b10_per_block),
                                         key=lambda t: t[1])[:args.max_b10]:
            gate(row, f"b10 {key}")
            add_ge_row("b10", key, vec, row)
            b10_cuts += 1
            b10_worst = min(b10_worst, lam)
        psd_cuts = 0
        psd_worst = 0.0
        for label, lam, row, vec in moments.separate(q):
            gate(row, f"psd {label}")
            add_ge_row("psd", label, vec, row)
            psd_cuts += 1
            psd_worst = min(psd_worst, lam)
        history.append((it, eta, bound, added, horns, b10_cuts, psd_cuts))
        print(f"          cuts: rules +{added} (viol {worst:.1e}) horn +{horns} ({horn_worst:+.1e}) "
              f"b10 +{b10_cuts} ({b10_worst:+.1e}) psd +{psd_cuts} ({psd_worst:+.1e})  "
              f"[{time.time() - started:.0f}s]", flush=True)
        if warm is not None:
            extra = added + horns + b10_cuts + psd_cuts
            if extra:
                # new nonneg rows sit between the old linear rows and the PSD part
                ins = psd_start
                warm["y"] = np.concatenate([warm["y"][:ins], np.zeros(extra), warm["y"][ins:]])
                warm["s"] = np.concatenate([warm["s"][:ins], np.zeros(extra), warm["s"][ins:]])
        if added == 0 and horns == 0 and b10_cuts == 0 and psd_cuts == 0:
            print("no violated rows found", flush=True)
            break
        if time.time() - started > args.time_limit:
            print("time limit reached", flush=True)
            break
    if args.out:
        with args.out.open("wb") as handle:
            pickle.dump({**cert, "history": history}, handle)
    print(f"total {time.time() - started:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
