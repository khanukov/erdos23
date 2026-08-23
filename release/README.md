# Release artifacts

The two logical gzip archives are stored as ordered `.part-*` chunks:

- `erdos23_full_solution_exact_certificate_2026-08-23.tar.gz`;
- `erdos23_full_replay_bundle_2026-08-23.tar.gz`.

Reconstruct and verify them with:

```bash
./scripts/assemble_artifacts.sh
```

The full archive expands to about 1.8 GB and contains every data file required
by the independent exact replay. `ASSEMBLED_SHA256SUMS` fixes the digest of
each logical artifact.
