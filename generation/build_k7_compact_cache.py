#!/usr/bin/env python3
"""Build a resumable exact-count cache for the 107 order-7 root types.

The historical ``k7_precompute.pkl`` stores normalized float64 tensors in one
multi-gigabyte pickle.  A killed writer leaves that file unusable.  Here every
root is written atomically and independently, and the normalized quantities
are stored as their exact integer numerators over 9P7 = 181440.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np


DENOMINATOR = 181440
_STATES = None
_OUTPUT = None
_CPP_THREADS = None


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--flagsdp", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--cpp-threads", type=int, default=2)
    return parser.parse_args()


def root_path(output: Path, root: int) -> Path:
    return output / f"root_{root:03d}.npz"


def valid_root(path: Path, root: int, adjacency) -> bool:
    if not path.exists():
        return False
    try:
        with np.load(path, allow_pickle=False) as payload:
            edge = payload["edge_raw"]
            mass = payload["mass_raw"]
            stored_root = int(payload["root"])
            stored_adjacency = tuple(int(x) for x in payload["adjacency"])
            return (
                stored_root == root
                and stored_adjacency == tuple(int(x) for x in adjacency)
                and edge.dtype == np.uint32
                and mass.dtype == np.uint32
                and edge.ndim == 3
                and mass.shape == (edge.shape[0],)
                and edge.shape[1] == edge.shape[2]
            )
    except Exception:
        return False


def build_one(task):
    root, adjacency = task
    import cpp_precompute as cpp

    edge, mass, classes = cpp.precompute_type_cpp(
        _STATES, 7, adjacency, nthreads=_CPP_THREADS
    )
    edge_scaled = edge * DENOMINATOR
    mass_scaled = mass * DENOMINATOR
    edge_raw = np.rint(edge_scaled).astype(np.uint32)
    mass_raw = np.rint(mass_scaled).astype(np.uint32)
    edge_error = float(np.max(np.abs(edge_scaled - edge_raw)))
    mass_error = float(np.max(np.abs(mass_scaled - mass_raw)))
    if edge_error > 5e-8 or mass_error > 5e-8:
        raise RuntimeError(
            f"root {root}: integer recovery failed: edge={edge_error}, mass={mass_error}"
        )

    target = root_path(_OUTPUT, root)
    temporary = target.with_suffix(f".tmp-{os.getpid()}")
    with temporary.open("wb") as handle:
        np.savez(
            handle,
            root=np.asarray(root, dtype=np.int32),
            adjacency=np.asarray(adjacency, dtype=np.uint32),
            edge_raw=edge_raw,
            mass_raw=mass_raw,
        )
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, target)
    return root, len(classes), tuple(edge_raw.shape), target.stat().st_size


def main():
    global _STATES, _OUTPUT, _CPP_THREADS
    args = parse_args()
    flagsdp = args.flagsdp.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(flagsdp))

    import cpp_precompute as cpp
    import flag_engine as fe

    cpp.compile_cpp()
    with (flagsdp / "cache_n9.pkl").open("rb") as handle:
        _STATES = pickle.load(handle)["states"]
    types = fe.enumerate_graphs(7, triangle_free=True)
    if len(types) != 107:
        raise RuntimeError(f"expected 107 root types, got {len(types)}")

    _OUTPUT = output
    _CPP_THREADS = max(1, args.cpp_threads)
    pending = [
        (root, adjacency)
        for root, (_k, adjacency) in enumerate(types)
        if not valid_root(root_path(output, root), root, adjacency)
    ]
    print(
        f"compact K7 cache: valid={len(types) - len(pending)}, pending={len(pending)}, "
        f"workers={args.workers}, cpp_threads={_CPP_THREADS}",
        flush=True,
    )
    started = time.time()
    if pending:
        with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
            futures = {executor.submit(build_one, task): task[0] for task in pending}
            completed = 0
            for future in as_completed(futures):
                root, nclasses, shape, size = future.result()
                completed += 1
                print(
                    f"root={root:03d} classes={nclasses:3d} shape={shape} "
                    f"bytes={size} ({completed}/{len(pending)})",
                    flush=True,
                )

    missing = [
        root
        for root, (_k, adjacency) in enumerate(types)
        if not valid_root(root_path(output, root), root, adjacency)
    ]
    if missing:
        raise RuntimeError(f"invalid or missing roots after build: {missing}")
    manifest = {
        "format": "erdos23-k7-exact-count-cache-v1",
        "denominator": DENOMINATOR,
        "n_states": len(_STATES),
        "n_roots": len(types),
        "adjacency": [[int(x) for x in adjacency] for _k, adjacency in types],
    }
    temporary_manifest = output / f"manifest.tmp-{os.getpid()}.json"
    temporary_manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    os.replace(temporary_manifest, output / "manifest.json")
    total_size = sum(root_path(output, root).stat().st_size for root in range(len(types)))
    print(
        f"complete: roots={len(types)} bytes={total_size} elapsed={time.time() - started:.1f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
