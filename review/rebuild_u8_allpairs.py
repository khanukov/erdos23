#!/usr/bin/env python3
"""All-pairs companion to rebuild_u8.py: every ordered pair of vertices, not
just edges, so that rooted-Horn (copositivity) rows can be generated on
trustworthy data.  Reuses the validated order-10 catalogue and root types."""
from __future__ import annotations
import argparse, pickle, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rebuild_u8 import Catalogue  # noqa: E402

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuilt", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    with args.rebuilt.open("rb") as handle:
        payload = pickle.load(handle)
    states10, roots8 = payload["states10"], payload["roots8"]
    cat8 = Catalogue(roots8)
    pairs = [dict() for _ in range(410)]
    for state, adjacency in enumerate(states10):
        for u in range(10):
            for v in range(u + 1, 10):
                rest = [x for x in range(10) if x != u and x != v]
                position = {x: i for i, x in enumerate(rest)}
                induced = [0] * 8
                for x in rest:
                    for y in rest:
                        if (adjacency[x] >> y) & 1:
                            induced[position[x]] |= 1 << position[y]
                tau, mapping = cat8.find(induced)
                prof_u = tuple(sorted(mapping[position[x]] for x in rest if (adjacency[u] >> x) & 1))
                prof_v = tuple(sorted(mapping[position[x]] for x in rest if (adjacency[v] >> x) & 1))
                for key in ((prof_u, prof_v), (prof_v, prof_u)):
                    bucket = pairs[tau].setdefault(key, {})
                    bucket[state] = bucket.get(state, 0) + 1
        if (state + 1) % 3000 == 0:
            print(f"  {state + 1}/12172 states [{time.time() - started:.0f}s]", flush=True)
    total = sum(sum(sum(b.values()) for b in root.values()) for root in pairs)
    print(f"all-pairs decomposition: {total} ordered incidences (expected {12172 * 90})")
    with args.out.open("wb") as handle:
        pickle.dump({"pairs": pairs}, handle)
    print(f"saved {args.out}  [{time.time() - started:.0f}s]")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
