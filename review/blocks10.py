#!/usr/bin/env python3
"""Order-10 flag moment blocks (types of size 0, 2, 4, 6) via the C helper blocks10.c.

`Blocks10(states, work_dir, levels)` compiles the helper if needed, exports the
states, builds (or loads) the kept-tuple index, and offers

    matrices(q)            -> list of (key, M) with key = (s, type), M dense float64
    separate(q, tol, k)    -> eigenvector cuts: (key, lam, row, int_vector)
    rows(vectors)          -> exact int64 numerators per state for (key, int_vector)

A cut's row over states is  r_G = v^T M_sigma(G) v / (10^12 * normaliser)  with
`v` the eigenvector rounded to integers over 10^6, so the exact verifier can
reproduce it bit for bit from the same integer data.
"""

from __future__ import annotations

import struct
import subprocess
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SCALE = 10 ** 6


class Blocks10:
    def __init__(self, states, work_dir: Path, levels: str = "0246", verbose=True):
        self.levels = levels
        self.work = Path(work_dir)
        self.work.mkdir(parents=True, exist_ok=True)
        self.binary = self.work / "blocks10"
        src = HERE / "blocks10.c"
        if not self.binary.exists() or self.binary.stat().st_mtime < src.stat().st_mtime:
            subprocess.run(["gcc", "-O2", "-march=native", "-Wall", "-o", str(self.binary), str(src)],
                           check=True)
        self.states_bin = self.work / "states10.bin"
        arr = np.asarray(states, dtype=np.uint16)
        self.n_states = len(arr)
        payload = struct.pack("i", self.n_states) + arr.tobytes()
        if not self.states_bin.exists() or self.states_bin.read_bytes() != payload:
            self.states_bin.write_bytes(payload)
        self.index_bin = self.work / f"tuples10_{levels}.bin"
        if not self.index_bin.exists():
            t0 = time.time()
            subprocess.run([str(self.binary), "index", str(self.states_bin), str(self.index_bin), levels],
                           check=True, capture_output=not verbose)
            if verbose:
                print(f"blocks10: index built [{time.time() - t0:.0f}s]", flush=True)
        info = subprocess.run([str(self.binary), "info", str(self.states_bin), levels],
                              check=True, capture_output=True, text=True).stdout
        self.keys, self.sizes, self.normaliser = [], [], {}
        current = None
        for line in info.splitlines():
            if line.startswith("level"):
                parts = line.split()
                current = int(parts[1][2:])
                self.normaliser[current] = int(float(parts[-1]))
            elif line.strip().startswith("type"):
                parts = line.split()
                self.keys.append((current, int(parts[1])))
                self.sizes.append(int(parts[-1]))
        self._calls = 0

    def _tmp(self, name):
        self._calls += 1
        return self.work / f"b10_{name}_{self._calls}.bin"

    def matrices(self, q):
        qf = self._tmp("q")
        out = self._tmp("M")
        qf.write_bytes(np.ascontiguousarray(q, dtype=np.float64).tobytes())
        subprocess.run([str(self.binary), "moment", str(self.states_bin), str(self.index_bin),
                        str(qf), str(out), self.levels], check=True)
        data = out.read_bytes()
        qf.unlink()
        out.unlink()
        pos, result = 0, []
        for key, n in zip(self.keys, self.sizes):
            (nF,) = struct.unpack_from("i", data, pos)
            pos += 4
            if nF != n:
                raise RuntimeError(f"block {key}: size {nF} != {n}")
            M = np.frombuffer(data, dtype=np.float64, count=nF * nF, offset=pos).reshape(nF, nF)
            pos += 8 * nF * nF
            result.append((key, M))
        return result

    def rows(self, vectors):
        """vectors: list of (key, int64 array) -> int64 array (len(vectors), n_states)."""
        if not vectors:
            return np.zeros((0, self.n_states), dtype=np.int64)
        vf = self._tmp("v")
        out = self._tmp("r")
        chunks = [struct.pack("i", len(vectors))]
        for (s, t), v in vectors:
            v = np.ascontiguousarray(v, dtype=np.int64)
            chunks.append(struct.pack("iii", s, t, len(v)) + v.tobytes())
        vf.write_bytes(b"".join(chunks))
        subprocess.run([str(self.binary), "rows", str(self.states_bin), str(self.index_bin),
                        str(vf), str(out), self.levels], check=True)
        acc = np.frombuffer(out.read_bytes(), dtype=np.int64).reshape(len(vectors), self.n_states)
        vf.unlink()
        out.unlink()
        return acc

    def float_rows(self, vectors):
        acc = self.rows(vectors)
        rows = np.empty(acc.shape, dtype=np.float64)
        for i, ((s, _t), _v) in enumerate(vectors):
            rows[i] = acc[i] / (float(SCALE) ** 2 * self.normaliser[s])
        return rows

    def separate(self, q, tol=1e-9, per_block=3):
        from scipy.linalg import eigh as dense_eigh
        from scipy.sparse.linalg import ArpackNoConvergence, eigsh

        cuts = []
        for key, M in self.matrices(q):
            M = (M + M.T) / 2.0
            n = M.shape[0]
            if not M.any():
                continue
            if n <= 400:
                vals, vecs = np.linalg.eigh(M)
            else:
                k = min(per_block, n - 2)
                try:
                    vals, vecs = eigsh(M, k=k, which="SA", tol=1e-10, maxiter=20000)
                except ArpackNoConvergence:
                    vals, vecs = dense_eigh(M, subset_by_index=[0, k - 1])
                order = np.argsort(vals)
                vals, vecs = vals[order], vecs[:, order]
            for j in range(min(per_block, len(vals))):
                if vals[j] < -tol:
                    v = np.rint(vecs[:, j] * SCALE).astype(np.int64)
                    cuts.append((key, float(vals[j]), v))
        if not cuts:
            return []
        rows = self.float_rows([(key, v) for key, _lam, v in cuts])
        return [(key, lam, rows[i], v) for i, (key, lam, v) in enumerate(cuts)]


if __name__ == "__main__":
    import argparse
    import pickle
    import sys

    from scipy.sparse import csr_matrix

    sys.path.insert(0, str(HERE))
    from build_u8_lp import N9, N10, q10_of  # noqa: E402
    from k7_leg_ceiling import cycle, grotzsch, petersen  # noqa: E402
    from rebuild_u8 import Catalogue  # noqa: E402

    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuilt", type=Path, required=True)
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--work", type=Path, default=Path(".work"))
    parser.add_argument("--levels", default="0246")
    parser.add_argument("--q", type=Path, default=None, help="pickle with an LP 'q' to test")
    args = parser.parse_args()
    t0 = time.time()
    with args.rebuilt.open("rb") as handle:
        states10 = pickle.load(handle)["states10"]
    blocks = Blocks10(states10, args.work, args.levels)
    print(f"{len(blocks.keys)} blocks, sizes up to {max(blocks.sizes)}, "
          f"sum of squares {sum(n * n for n in blocks.sizes)}  [{time.time() - t0:.0f}s]")
    with (args.flagsdp / "cache_n9.pkl").open("rb") as handle:
        cache9 = pickle.load(handle)
    lift = np.load(args.flagsdp / "c5lift_cache.npz", allow_pickle=True)
    cat9 = Catalogue([list(masks) for _o, masks in cache9["states"]])
    counts = np.rint(np.asarray(lift["Dval"]) * 10).astype(np.int64)
    columns = csr_matrix((counts, (lift["Drow"], lift["Dcol"])), shape=(N9, N10), dtype=np.int64).tocsc()
    profile_to_state = {}
    for j in range(N10):
        s0, e0 = columns.indptr[j], columns.indptr[j + 1]
        profile_to_state[tuple(sorted(zip(columns.indices[s0:e0].tolist(),
                                          columns.data[s0:e0].tolist())))] = j
    seeds = [("C5", cycle(5), np.full(5, 0.2)),
             ("Grotzsch", grotzsch(), np.asarray(
                 [10173, 9717, 8691, 10166, 7628, 8213, 7344, 9352, 12084, 10443, 6189], float)),
             ("Petersen", petersen(), np.full(10, 0.1))]
    for name, base, w in seeds:
        q = q10_of(cat9, profile_to_state, base, w / w.sum())
        t1 = time.time()
        mats = blocks.matrices(q)
        worst = 0.0
        sums = {}
        for (s, t), M in mats:
            vals = np.linalg.eigvalsh((M + M.T) / 2)
            worst = min(worst, float(vals[0]))
            sums[(s, t)] = float(M.sum())
        print(f"{name:9s}: min eigenvalue over all blocks {worst:+.3e}; "
              f"sum of entries s=0: {sums[(0, 0)]:.6f}, s=2: "
              f"{sums.get((2, 0), 0):.6f} + {sums.get((2, 1), 0):.6f}  [{time.time() - t1:.1f}s]")
    if args.q:
        with args.q.open("rb") as handle:
            q = np.asarray(pickle.load(handle)["q"], dtype=float)
        t1 = time.time()
        cuts = blocks.separate(q)
        print(f"LP point: {len(cuts)} cuts, most negative eigenvalue "
              f"{min((c[1] for c in cuts), default=0.0):+.3e}  [{time.time() - t1:.1f}s]")
        for key, lam, row, v in sorted(cuts, key=lambda c: c[1])[:8]:
            print(f"  block {key}: lambda {lam:+.3e}  row.q = {float(row @ q):+.3e}")
    print(f"total {time.time() - t0:.0f}s")
