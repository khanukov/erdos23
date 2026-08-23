# Erdős Problem 23: exact order-10 certificate

[![Exact replay](https://github.com/khanukov/erdos23/actions/workflows/exact-replay.yml/badge.svg)](https://github.com/khanukov/erdos23/actions/workflows/exact-replay.yml)

This repository contains a reproducible computer-assisted proof candidate for
the Erdős-Győri-Simonovits triangle-free bipartization conjecture:

\[
  \operatorname{bip}(G)\le \frac{|V(G)|^2}{25}
\]

for every finite triangle-free graph \(G\), where
\(\operatorname{bip}(G)=|E(G)|-\operatorname{MaxCut}(G)\).

> **Review status.** The finite certificate has been replayed repeatedly in
> exact rational arithmetic, including from the self-contained archive. The
> argument has not yet received independent third-party peer review. Until
> that happens, describe it as an externally reviewable proof candidate, not
> as an accepted solution.

## Exact result

The order-10 dual certificate has

\[
  \delta_\star
  =-9.878886951679021\ldots\times10^{-4}
  <-\frac{9}{10000}<0.
\]

It contains 10,188 finitely described K7, K8, and rooted-Horn rows. The
independent verifier reconstructs an integer matrix of shape
\(10188\times12172\) with 11,560,170 nonzero entries and checks every dual
residual exactly. The fixed Gram moment vector is separately reconstructed
from 87 nonnegative rational rank-one atoms and matches all 12,172
coordinates.

The negative objective proves the middle-density band with no finite-order
integrality correction. Ferudun's exact blow-up identity and the
Balogh-Clemen-Lidický density-tail theorem then cover every finite order and
every edge density. The full deduction and its assumptions are in
[the review draft](paper/erdos23_full_solution_draft.pdf).

## One-command replay

The complete archive is only about 20 MB compressed but expands to about
1.8 GB. Allow at least 3 GB of free disk space and approximately 4 GB of RAM.

```bash
python3 -m pip install -r requirements-replay.txt
./scripts/check_artifacts.sh
./scripts/run_full_replay.sh
```

Successful output includes:

```text
PASS: public exact moment vector was independently reconstructed
EXACT_DUAL_REPLAY_OK
min_q_residual=0 at state=46
objective ~= -9.878886951679021e-04
global_closure=True
GLOBAL_CONCLUSION: beta(G)<=N^2/25 for every finite triangle-free graph G on N vertices
```

GitHub Actions runs the same checksum, JSON replay, pickle replay, and
moment-vector reconstruction on every push and pull request.

## Repository map

- `paper/`: LaTeX source and rendered five-page review draft.
- `certificate/`: numbered chunks of the transparent JSON certificate and
  exact assembly instructions. The equivalent pickle, source dual, descriptor
  state, and moment data are in the release archives.
- `verifier/`: independent exact dual replay and independent moment-vector
  reconstruction. The dual verifier does not import the LP generator.
- `generation/`: search, separation, solving, rationalization, and K7-cache
  regeneration scripts.
- `release/`: numbered chunks of the compact certificate archive and the
  self-contained full replay bundle.
- `docs/SCIENTIFIC_NOTE.md`: concise theorem and proof bridge.
- `docs/VERIFICATION_REPORT.md`: recorded checks and trust boundary.
- `REVIEW_GUIDE.md`: suggested external-review protocol.

## Trust boundary

The executable replay verifies:

1. certificate syntax and nonnegative multipliers;
2. all 517 per-root lower bounds and the exact envelope-leg equality;
3. exact reconstruction of every selected K7, K8, and Horn row;
4. all 12,172 coordinate residuals and the exact rational objective;
5. exact reconstruction of the public fixed Gram moment vector.

The deduction from this finite certificate to the graph theorem additionally
uses the generic validity of the per-root MaxCut, rooted-Horn, and positive
semidefinite moment rows, plus the cited density-tail theorem and blow-up
transfer. Those mathematical interfaces are the main targets for external
human review.

## References and provenance

- A. Ferudun, *The Erdős \(n^2/25\) max-cut conjecture for small multiples
  of five, via a per-root-MaxCut envelope and blow-up integrality*,
  [arXiv:2606.28041v1](https://arxiv.org/abs/2606.28041).
- J. Balogh, F. C. Clemen, and B. Lidický, *Max cuts in triangle-free
  graphs*, [arXiv:2103.14179](https://arxiv.org/abs/2103.14179).
- [Erdős Problem 23](https://www.erdosproblems.com/23).

The present certificate extends Ferudun's public framework. Manuscript
authorship and final attribution should be settled before formal submission.
