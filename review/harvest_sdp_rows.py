#!/usr/bin/env python3
"""Turn an `sdp_u8.py` result into a row pool the LP builder can resume from.

The SDP's PSD duals are eigen-decomposed into `rpsd` / `b10` rows by
`sdp_u8.py`; there can be tens of thousands.  This keeps, per block, the
eigen-rows with the largest multipliers (`--per-block`) and every linear row
with a positive multiplier, and writes a checkpoint-shaped pickle.
"""

from __future__ import annotations

import argparse
import pickle
from collections import defaultdict
from pathlib import Path

import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdp", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--per-block", type=int, default=12)
    parser.add_argument("--min-mult", type=float, default=1e-9)
    args = parser.parse_args()
    with args.sdp.open("rb") as handle:
        result = pickle.load(handle)
    rows, duals = result["rows"], np.asarray(result["duals"], dtype=float)
    keep = []
    by_block = defaultdict(list)
    for index, (kind, a, b) in enumerate(rows):
        if kind in ("band", "gram", "mass", "leg8"):
            keep.append(index)
        elif kind in ("rpsd", "b10"):
            by_block[(kind, a if kind == "rpsd" else tuple(a))].append((abs(float(duals[index])), index))
        elif abs(float(duals[index])) > args.min_mult or kind == "rule":
            keep.append(index)
    eig_rows = 0
    for _key, items in by_block.items():
        items.sort(reverse=True)
        for mult, index in items[:args.per_block]:
            if mult > args.min_mult:
                keep.append(index)
                eig_rows += 1
    keep.sort()
    out_rows = [rows[i] for i in keep]
    out_duals = duals[keep]
    with args.out.open("wb") as handle:
        pickle.dump({"eta": result["eta"], "q": result["q"], "u": result["u"],
                     "duals": out_duals, "rows": out_rows, "col_u": result["col_u"]}, handle)
    kinds = defaultdict(int)
    for kind, _a, _b in out_rows:
        kinds[kind] += 1
    print(f"kept {len(out_rows)} of {len(rows)} rows ({eig_rows} eigen-rows from {len(by_block)} blocks): "
          + ", ".join(f"{k} {v}" for k, v in sorted(kinds.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
