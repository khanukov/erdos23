# Certificate files

The primary transparent certificate is
`erdos23_global_exact_dual.json`. It is stored here as two numbered
`.part-*` files. Run:

```bash
./scripts/assemble_artifacts.sh
```

The reconstructed JSON appears at
`.replay-work/erdos23_global_exact_dual.json` and is checked against the
published SHA-256 digest.

The equivalent pickle certificate, floating source dual, descriptor state,
public moment files, and exact verifier are all contained in both logical
release archives. The complete replay archive additionally contains the four
public flag-data files and 107 derived K7 exact-count caches.
