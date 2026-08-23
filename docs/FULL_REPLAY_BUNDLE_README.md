# Full replay bundle

The complete replay archive supplements the compact certificate package with
the four public flag-data files and the derived exact K7 count cache required
by the two independent verifiers.  Its uncompressed contents are large, but
the sparse count arrays compress well.

From the unpacked bundle root, run:

```bash
python verify_moment_vector.py \
  --flagsdp ./public_flagsdp_data \
  --public-anc ./public_anc --jobs 8

python verify_exact_fixed_gram_dual.py \
  --certificate ./erdos23_global_exact_dual.json \
  --flagsdp ./public_flagsdp_data \
  --public-anc ./public_anc \
  --k7-cache ./k7_compact_v1 --target-n 1
```

The four files in `public_flagsdp_data` are copied unchanged from the public
computational data associated with arXiv:2606.28041v1.  The `k7_compact_v1`
directory was generated deterministically by `build_k7_compact_cache.py`;
its manifest records the exact denominator, state count, root count, and root
adjacencies.  The top-level checksum manifest covers every bundled file.
