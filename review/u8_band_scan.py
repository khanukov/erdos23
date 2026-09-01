#!/usr/bin/env python3
"""sup of the rebuilt 8-root envelope over the band, on the points most likely
to be worst: for each base graph, the in-band weights that maximise d_mono
(closed form, cheap), then U8 by MaxCut per root on the rebuilt rows."""
from __future__ import annotations
import argparse, pickle, random, sys, time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_u8_lp import RootData, q10_of  # noqa: E402
from rebuild_u8 import Catalogue  # noqa: E402
from k7_leg_ceiling import cycle, edges_of, grotzsch, petersen, triangles  # noqa: E402
N9, N10 = 1897, 12172
LO, HI = 0.2486, 0.3197

def densities(base, w):
    E = edges_of(base); pr = [w[u] * w[v] for u, v in E]
    return 2 * sum(pr), 2 * min(sum(p for (u, v), p in zip(E, pr) if ((t >> u) & 1) == ((t >> v) & 1))
                                for t in range(1 << len(base)))

def maximise_dmono(base, seed, restarts=6, steps=300):
    h = len(base); rng = random.Random(seed); best = (-1, None)
    for t in range(restarts):
        w = np.full(h, 1.0 / h) if t == 0 else np.asarray([rng.random() + 1e-3 for _ in range(h)]); w /= w.sum()
        de, dm = densities(base, w); cur = dm if LO <= de <= HI else -1; sc = 0.3
        for st in range(steps):
            i, j = rng.randrange(h), rng.randrange(h)
            if i == j: continue
            T = w.copy(); d = sc * rng.random() * T[i]; T[i] -= d; T[j] += d
            de, dm = densities(base, T)
            if not LO <= de <= HI: continue
            if dm > cur: w, cur = T, dm
            if st % 50 == 49: sc *= 0.65
        if cur > best[0]: best = (cur, w.copy())
    return best

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--flagsdp", type=Path, required=True)
    ap.add_argument("--rebuilt", type=Path, required=True); ap.add_argument("--starts", type=int, default=6)
    ap.add_argument("--bases", type=int, default=12); args = ap.parse_args(); t0 = time.time()
    payload = pickle.load(args.rebuilt.open("rb")); states10, roots8, pairs = payload["states10"], payload["roots8"], payload["pairs"]
    cache = pickle.load((args.flagsdp / "cache_n9.pkl").open("rb")); states9 = [list(m) for _o, m in cache["states"]]; cat9 = Catalogue(states9)
    lift = np.load(args.flagsdp / "c5lift_cache.npz", allow_pickle=True)
    cnt = np.rint(np.asarray(lift["Dval"]) * 10).astype(np.int64)
    cols = csr_matrix((cnt, (lift["Drow"], lift["Dcol"])), shape=(N9, N10), dtype=np.int64).tocsc()
    p2s = {}
    for j in range(N10):
        s, e = cols.indptr[j], cols.indptr[j + 1]
        p2s[tuple(sorted(zip(cols.indices[s:e].tolist(), cols.data[s:e].tolist())))] = j
    roots = [RootData(t, roots8[t], pairs[t]) for t in range(410) if pairs[t]]
    print(f"root data ready [{time.time() - t0:.0f}s]", flush=True)
    rng = random.Random(3)
    def U8(q): return sum(r.best_rule(q, rng, args.starts)[1] for r in roots)
    cands = [("C5", cycle(5)), ("Petersen", petersen()), ("Grotzsch", grotzsch())]
    ranked = sorted(range(N9), key=lambda i: -sum(bin(m).count("1") for m in states9[i]))
    cands += [(f"tf9#{i}", states9[i]) for i in ranked[:args.bases] if edges_of(states9[i])]
    print(f"{'graphon':>12} {'d_edge':>8} {'d_mono':>9} {'U8':>10} {'U8-2/25':>10}  valid")
    worst = -1
    for name, base in cands:
        assert not triangles(base)
        dm, w = maximise_dmono(base, sum(map(ord, name)))
        if w is None: continue
        de, dm = densities(base, w); q = q10_of(cat9, p2s, base, w)
        u = U8(q); worst = max(worst, u)
        print(f"{name:>12} {de:>8.4f} {dm:>9.6f} {u:>10.6f} {u - 0.08:>+10.6f}  {u >= dm - 1e-9}", flush=True)
    print(f"\nlargest U8 in the band over {len(cands)} points: {worst:.6f}  (2/25 = 0.080000)  [{time.time() - t0:.0f}s]")
main()
