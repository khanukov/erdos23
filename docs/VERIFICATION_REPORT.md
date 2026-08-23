# Verification report

Date: 2026-08-23

## Exact fixed-Gram reconstruction

Command: `verify_moment_vector.py --jobs 16`

```text
order-9 states: 1897
order-10 states: 12172
manifest-Gram atoms: 87
nonnegative rational weights: True
exact coordinate matches: 12172/12172
PASS: public exact moment vector was independently reconstructed
```

## Exact dual replay

Primary input: `erdos23_global_exact_dual.json`

```text
replayed exact matrix: shape=(10188, 12172) nnz=11560170
EXACT_DUAL_REPLAY_OK
descriptors=10188 nnz=11560170 digits=3
digit_product_audit=9/9
min_q_residual=0 at state=46
objective ~= -9.878886951679021e-04
global_closure=True
GLOBAL_CONCLUSION: beta(G)<=N^2/25 for every finite triangle-free graph G on N vertices
```

Descriptor counts:

```text
K7 MaxCut: 2790
K8 MaxCut: 2385
Horn:      5013
Total:    10188
```

The JSON and pickle certificates produced the same exact objective and both
passed the independent replay. Rationalization used denominator `10^12` and
required zero root repairs.

The self-contained full replay bundle was rechecked on 2026-08-23 using only
its bundled flag data and K7 cache. Its complete checksum manifest, moment
reconstruction, JSON replay, and pickle replay all passed.

## Search provenance

The previous exact certificate had positive objective
`3.4394777894128489e-05` and proved the conjectured equality through 240
vertices. Separating a later numerical witness added 218 K7, 352 K8, and 2045
Horn rows. The resulting 10,188-row dual had floating objective
`-9.878886947326e-04`, which rationalized to the negative exact objective above.

Fifty-eight adaptive Gram directions were also found as a diagnostic, but
they are not used by the final certificate.

## Trust boundary

The replay verifies the finite certificate, exact arithmetic, descriptor
syntax, root normalizations, fixed Gram vector, and objective. The deduction
to the graph theorem additionally cites the generic validity of the
per-root-MaxCut, rooted-Horn, and moment rows and the density-tail theorem, as
proved in the references listed in `SCIENTIFIC_NOTE.md`.

External peer review has not yet occurred.
