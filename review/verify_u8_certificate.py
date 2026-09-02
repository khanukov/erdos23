#!/usr/bin/env python3
"""Exact rational verification of a dual certificate produced by build_u8_lp.py.

The LP solver's duals are floating point.  A proof needs none of that: any
nonnegative multipliers give a valid bound, and the bound they give is computed
here in exact arithmetic from integer data.

Primal rows, written as `<row, x> <= rhs`:

    (B+)  <d, q>            <= HI            y_hi >= 0
    (B-)  -<d, q>           <= -LO           y_lo >= 0
    (G)   -<m, q>           <= 0             y_g  >= 0
    (L)   eta - sum_tau u   <= -2/25         y_L  >= 0   (must equal 1)
    (R_i) u_tau - <h_i, q>  <= 0             z_i  >= 0
    (H_j) -<horn_j, q>      <= 0             w_j  >= 0
    (P_k) -<psd_k, q>       <= 0             v_k  >= 0
    (Q_l) -<rpsd_l, q>      <= 0             t_l  >= 0   (v^T M_tau(q) v >= 0)
    (T_m) -<b10_m, q>       <= 0             b_m  >= 0   (order-10 flag moment blocks)
    (W_n) -<agg_n, q>       <= 0             a_n  >= 0   (tr(U U^T M(q)) >= 0, aggregated)
    (M)   sum q              = 1             rho free

Summing with these multipliers, and using `u >= 0`, `q >= 0`, gives `eta <= delta`
provided every `u_tau` coefficient `sum_{i in tau} z_i - 1` is `>= 0` and every
`q_s` coefficient is `>= 0`; the latter fixes `rho` as the largest "load"

    load_s = sum_i z_i h_i[s] + sum_j w_j horn_j[s] + sum_k v_k psd_k[s] + y_g m_s - (y_hi - y_lo) d_s ,

and then  `delta = y_hi HI - y_lo LO - 2/25 + rho`.  The script rationalises the
solver's duals, scales each root's rule multipliers up to sum to at least one,
rebuilds every row from the integer decomposition, and prints `delta` exactly.
`delta <= 0` is an exact certificate of `d_mono <= 2/25` on the closed band,
relative to the validity of the row families.
"""

from __future__ import annotations

import argparse
import pickle
import sys
import time
from fractions import Fraction
from math import comb, prod
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_u8_lp import RootData  # noqa: E402

N9, N10 = 1897, 12172
HI, LO = Fraction(3197, 10000), Fraction(2486, 10000)


def rational(x, denominator=10 ** 12):
    return Fraction(int(round(float(x) * denominator)), denominator)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--public-anc", type=Path, required=True)
    parser.add_argument("--rebuilt", type=Path, required=True)
    parser.add_argument("--allpairs", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--work", type=Path, default=Path(".work"))
    args = parser.parse_args()
    started = time.time()

    with args.result.open("rb") as handle:
        result = pickle.load(handle)
    rows_log, duals = result["rows"], np.asarray(result["duals"], dtype=float)
    if len(rows_log) != len(duals):
        raise RuntimeError(f"{len(rows_log)} rows logged but {len(duals)} duals")
    with args.rebuilt.open("rb") as handle:
        payload = pickle.load(handle)
    states10, roots8, edge_pairs = payload["states10"], payload["roots8"], payload["pairs"]
    with args.allpairs.open("rb") as handle:
        all_pairs = pickle.load(handle)["pairs"]
    with (args.flagsdp / "cache_n9.pkl").open("rb") as handle:
        cache9 = pickle.load(handle)
    lift = np.load(args.flagsdp / "c5lift_cache.npz", allow_pickle=True)
    with (args.public_anc / "mom_term_exact.pkl").open("rb") as handle:
        gram = [Fraction(int(n), int(d)) for n, d in pickle.load(handle)]

    # exact static data
    d_edge = [Fraction(sum(bin(m).count("1") for m in adj) // 2, 45) for adj in states10]
    counts = np.rint(np.asarray(lift["Dval"]) * 10).astype(np.int64)
    deletion = csr_matrix((counts, (lift["Drow"], lift["Dcol"])), shape=(N9, N10)).tocsc()
    moment_blocks = []
    for label, _size, sigma, _flags, free_size, _flat, mats in cache9["moments"]:
        root_size = sigma[0]
        den = [prod(order - i for i in range(root_size)) * comb(order - root_size, free_size) ** 2
               if order - root_size >= free_size else 1 for order, _adj in cache9["states"]]
        moment_blocks.append((label, np.asarray(mats, dtype=object), den))
    labels = [b[0] for b in moment_blocks]

    # root data, exactly as the LP built it
    edge_roots, all_roots = {}, {}
    for tau in range(410):
        if edge_pairs[tau]:
            edge_roots[tau] = RootData(tau, roots8[tau], edge_pairs[tau])
        if all_pairs[tau]:
            all_roots[tau] = RootData(tau, roots8[tau], all_pairs[tau])
    print(f"root data rebuilt [{time.time() - started:.0f}s]", flush=True)

    def integer_incidence(root):
        """orbit x state integer counts, as Python ints (exact)."""
        N = root.N.tocoo()
        rows = [dict() for _ in range(len(root.members))]
        for o, s, v in zip(N.row.tolist(), N.col.tolist(), N.data.tolist()):
            rows[o][s] = rows[o].get(s, 0) + int(round(v))
        return rows

    incidence_cache = {}

    def orbit_rows(root, key):
        if key not in incidence_cache:
            incidence_cache[key] = integer_incidence(root)
        return incidence_cache[key]

    # multipliers
    y_hi = max(Fraction(0), -rational(duals[0]))    # HiGHS: negative dual on an active upper bound (min form)
    y_lo = max(Fraction(0), rational(duals[0]))
    y_g = max(Fraction(0), rational(duals[1]))
    y_L = max(Fraction(0), -rational(duals[3]))
    print(f"static multipliers: y_hi={float(y_hi):.6e} y_lo={float(y_lo):.6e} "
          f"y_g={float(y_g):.6e} y_L={float(y_L):.6e}")
    if y_L == 0:
        raise RuntimeError("the leg multiplier is zero; the dual does not bound eta")
    # normalise so the eta coefficient is exactly one
    scale = 1 / y_L
    y_hi, y_lo, y_g = y_hi * scale, y_lo * scale, y_g * scale
    # HiGHS (minimisation): a row at its upper bound carries a dual <= 0, a row at
    # its lower bound a dual >= 0.  Rule rows are `<= 0`, Horn and PSD rows `>= 0`.
    multipliers = []
    for index, (kind, _a, _b) in enumerate(rows_log):
        if kind in ("band", "gram", "mass", "leg8"):
            multipliers.append(Fraction(0))
        elif kind == "rule":
            multipliers.append(max(Fraction(0), -rational(duals[index])) * scale)
        else:
            multipliers.append(max(Fraction(0), rational(duals[index])) * scale)

    # every root's rule multipliers must sum to at least one
    per_root = {}
    for index, (kind, tau, _c) in enumerate(rows_log):
        if kind == "rule":
            per_root[tau] = per_root.get(tau, Fraction(0)) + multipliers[index]
    scaled_roots = 0
    for tau, total in per_root.items():
        if total < 1:
            factor = 1 / total if total > 0 else None
            if factor is None:
                raise RuntimeError(f"root {tau} carries no rule multiplier at all")
            for index, (kind, t2, _c) in enumerate(rows_log):
                if kind == "rule" and t2 == tau:
                    multipliers[index] *= factor
            scaled_roots += 1
    print(f"rule multipliers: {sum(1 for r in rows_log if r[0] == 'rule')} rows, "
          f"{scaled_roots} roots scaled up to reach a sum of one")

    # order-10 block rows: exact integer numerators from the C helper, all at once
    b10_rows = [(index, tuple(a), np.asarray(b, dtype=np.int64))
                for index, (kind, a, b) in enumerate(rows_log) if kind == "b10" and multipliers[index] != 0]
    b10_exact = {}
    if b10_rows:
        from blocks10 import SCALE, Blocks10
        levels = "".join(sorted({str(a[0]) for _i, a, _v in b10_rows}))
        blocks = Blocks10(states10, args.work, levels, verbose=False)
        numerators = blocks.rows([(a, v) for _i, a, v in b10_rows])
        for (index, a, _v), num in zip(b10_rows, numerators):
            b10_exact[index] = (num, SCALE ** 2 * blocks.normaliser[a[0]])
        print(f"order-10 block rows: {len(b10_rows)} exact numerators computed [{time.time() - started:.0f}s]",
              flush=True)

    b10w_rows = [(index, tuple(a), np.asarray(b, dtype=np.int64))
                 for index, (kind, a, b) in enumerate(rows_log) if kind == "b10w" and multipliers[index] != 0]
    b10w_exact = {}
    if b10w_rows:
        from blocks10 import AGG_SCALE, Blocks10
        levels = "".join(sorted({str(a[0]) for _i, a, _u in b10w_rows}))
        blocks = Blocks10(states10, args.work, levels, verbose=False)
        mats = [(a, U @ U.T) for _i, a, U in b10w_rows]          # exact integer Gram matrices
        numerators = blocks.form_rows(mats)
        for (index, a, _u), num in zip(b10w_rows, numerators):
            b10w_exact[index] = (num, AGG_SCALE ** 2 * blocks.normaliser[a[0]])
        print(f"aggregated order-10 rows: {len(b10w_rows)} exact numerators computed "
              f"[{time.time() - started:.0f}s]", flush=True)

    # loads
    load = [Fraction(0)] * N10
    for s in range(N10):
        load[s] = y_g * gram[s] - (y_hi - y_lo) * d_edge[s]
    n_rule = n_horn = n_psd = n_rpsd = n_b10 = n_agg = 0
    for index, (kind, a, b) in enumerate(rows_log):
        mult = multipliers[index]
        if mult == 0 or kind in ("band", "gram", "mass", "leg8"):
            continue
        if kind == "rule":
            root = edge_roots[a]
            colouring = np.asarray(b, dtype=np.int8)
            rows = orbit_rows(root, ("e", a))
            for o, (ma, mb) in enumerate(zip(root.member_a, root.member_b)):
                mono = int(np.sum(colouring[ma] == colouring[mb]))
                if mono == 0:
                    continue
                coeff = mult * Fraction(mono, 90 * len(ma))
                for s, cnt in rows[o].items():
                    load[s] += coeff * cnt
            n_rule += 1
        elif kind == "horn":
            root = all_roots[a]
            cyc = list(b)
            C = root.horn_coefficients(cyc)
            flat = np.rint(C.ravel()).astype(np.int64)
            rows = orbit_rows(root, ("a", a))
            for o, m in enumerate(root.members):
                total = int(flat[m].sum())
                if total == 0:
                    continue
                coeff = mult * Fraction(total, 90 * len(m))
                for s, cnt in rows[o].items():
                    load[s] += coeff * cnt
            n_horn += 1
        elif kind == "b10w":
            num, den = b10w_exact[index]
            for s, value in enumerate(num.tolist()):
                if value:
                    load[s] += mult * Fraction(value, den)
            n_agg += 1
        elif kind == "rpsdw":
            root = all_roots[a]
            U = np.asarray(b, dtype=np.int64)
            W = U @ U.T                                    # integer, entries over 1000^2
            rows = orbit_rows(root, ("a", a))
            for o, (ma, mb) in enumerate(zip(root.member_a, root.member_b)):
                total = int(W[ma, mb].sum())
                if total == 0:
                    continue
                coeff = mult * Fraction(total, 1000 ** 2 * 90 * len(ma))
                for s, cnt in rows[o].items():
                    load[s] += coeff * cnt
            n_agg += 1
        elif kind == "psdw":
            block_index, U = b
            label, mats, den = moment_blocks[block_index]
            W = np.asarray(U, dtype=np.int64) @ np.asarray(U, dtype=np.int64).T
            k = W.shape[0]
            Wl = [[int(W[i, j]) for j in range(k)] for i in range(k)]
            D = deletion.tocsr()
            for s9 in range(N9):
                M = mats[s9]
                acc = 0
                for i in range(k):
                    row_i = M[i]
                    acc += sum(Wl[i][j] * int(row_i[j]) for j in range(k))
                if acc == 0:
                    continue
                w = Fraction(acc, 1000 ** 2 * den[s9])
                for s10, cnt in zip(D.indices[D.indptr[s9]:D.indptr[s9 + 1]].tolist(),
                                    D.data[D.indptr[s9]:D.indptr[s9 + 1]].tolist()):
                    load[s10] += mult * w * Fraction(int(cnt), 10)
            n_agg += 1
        elif kind == "b10":
            num, den = b10_exact[index]
            for s, value in enumerate(num.tolist()):
                if value:
                    load[s] += mult * Fraction(value, den)
            n_b10 += 1
        elif kind == "rpsd":
            root = all_roots[a]
            v = [rational(x, 10 ** 6) for x in b]
            rows = orbit_rows(root, ("a", a))
            for o, (ma, mb) in enumerate(zip(root.member_a, root.member_b)):
                total = sum(v[i] * v[j] for i, j in zip(ma.tolist(), mb.tolist()))
                if total == 0:
                    continue
                coeff = mult * total / (90 * len(ma))
                for s, cnt in rows[o].items():
                    load[s] += coeff * cnt
            n_rpsd += 1
        elif kind == "psd":
            block_index, vec = b
            label, mats, den = moment_blocks[block_index]
            v = [rational(x, 10 ** 6) for x in vec]
            k = len(v)
            w = []
            for s9 in range(N9):
                M = mats[s9]
                acc = Fraction(0)
                for i in range(k):
                    if v[i] == 0:
                        continue
                    row_i = M[i]
                    inner = sum(int(row_i[j]) * v[j] for j in range(k) if v[j] != 0)
                    acc += v[i] * inner
                w.append(acc / den[s9])
            # lift through the deletion matrix (counts / 10)
            D = deletion.tocsr()
            for s9 in range(N9):
                if w[s9] == 0:
                    continue
                for s10, cnt in zip(D.indices[D.indptr[s9]:D.indptr[s9 + 1]].tolist(),
                                    D.data[D.indptr[s9]:D.indptr[s9 + 1]].tolist()):
                    load[s10] += mult * w[s9] * Fraction(int(cnt), 10)
            n_psd += 1
        if (index + 1) % 500 == 0:
            print(f"  accumulated {index + 1}/{len(rows_log)} rows [{time.time() - started:.0f}s]",
                  flush=True)
    rho = max(load)
    delta = y_hi * HI - y_lo * LO - Fraction(2, 25) + rho
    print(f"\nrows used: {n_rule} rules, {n_horn} Horn, {n_psd} PSD, {n_rpsd} root-PSD, "
          f"{n_b10} order-10 block rows, {n_agg} aggregated rows")
    print(f"rho (largest load) = {float(rho):+.12e}")
    print(f"delta exact        = {delta}")
    print(f"delta              = {float(delta):+.12e}")
    print(f"solver's eta       = {result['eta']:+.12e}")
    print("\nVERDICT:", "delta <= 0 -- exact certificate of d_mono <= 2/25 + delta on the band"
          if delta <= 0 else "delta > 0 -- no certificate")
    print(f"total {time.time() - started:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
