# Verification status

Status date: 2026-08-23.

## Completed

- Exact JSON certificate replay: **PASS**.
- Equivalent pickle certificate replay: **PASS**.
- Descriptor pool: 2,790 K7 rows, 2,385 K8 rows, 5,013 Horn rows.
- Reconstructed matrix: \(10188\times12172\), 11,560,170 nonzeros.
- Minimum exact dual residual: zero, at state 46.
- Exact objective:
  \(\delta_\star\approx-9.878886951679021\times10^{-4}<0\).
- Fixed Gram moment reconstruction: **12,172/12,172 exact matches**.
- Full-bundle checksum replay: **PASS**.
- Compact and full gzip archives: **PASS**.
- LaTeX compilation: two clean passes, no warnings or overfull boxes.
- PDF visual inspection: all five pages checked.

## Still required for an accepted solution

- Independent execution by a third party on a clean machine.
- Human audit of the generic K7/K8 MaxCut-row validity.
- Human audit of the rooted-Horn and fixed-Gram interfaces.
- Human audit of normalization, closed-band boundaries, density-tail transfer,
  and the final graphon-to-finite-graph deduction.
- Final authorship, attribution, license, and submission decisions.
- External peer review and public preprint submission.

The repository therefore records a complete, falsifiable, reproducible proof
candidate. It does not claim that external verification has already occurred.
