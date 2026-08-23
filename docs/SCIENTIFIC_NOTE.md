# An exact order-10 certificate for Erdos Problem 23

## Statement

For a finite simple graph `G`, let

\[
\operatorname{bip}(G)=|E(G)|-\operatorname{MaxCut}(G),
\]

the minimum number of edges whose deletion makes `G` bipartite.

**Theorem (computer-assisted).** Every triangle-free graph `G` on `N`
vertices satisfies

\[
\operatorname{bip}(G)\le \frac{N^2}{25}.
\]

The constant is sharp: if `5` divides `N`, the balanced blow-up of `C_5`
has `bip(G)=N^2/25`.

## New exact band certificate

The order-10 per-root-MaxCut/Horn/fixed-Gram relaxation of Ferudun's public
paper bounds the monochromatic-pair density on the closed edge-density band

\[
0.2486\le d_{\rm edge}\le0.3197
\]

by

\[
d_{\rm mono}(W)\le \frac{2}{25}+\delta.
\]

The certificate in this package has the exact rational objective stored in
`erdos23_global_exact_dual.json`. Its decimal value is

\[
\delta_\star=-9.878886951679021\ldots\times10^{-4}
             <-\frac{9}{10000}<0.
\]

It uses 10,188 exact row descriptors:

- 2,790 K7 per-root MaxCut rows;
- 2,385 K8 per-root MaxCut rows;
- 5,013 rooted Horn rows.

All non-mass multipliers have denominator `10^12`. The five static
numerators (high density, low density, fixed Gram, K7 leg, K8 leg) are

```text
             0
             0
53956267978680
             1
  999999999999
```

Thus the two envelope legs sum to `10^12` exactly. No root repair was
required. The exact q-coordinate residual is nonnegative on all 12,172
states and is zero at state 46.

## Deduction of the full theorem

Write

\[
d_{\rm mono}(G)=\frac{2\operatorname{bip}(G)}{N^2}.
\]

Ferudun proves the exact blow-up identity

\[
\operatorname{bip}(G[t])=t^2\operatorname{bip}(G),
\]

so the step graphon `W_G` satisfies

\[
d_{\rm mono}(W_G)=\frac{2\operatorname{bip}(G)}{N^2}
\]

with no limiting error.

If the edge density of `G` is in the closed middle band, the exact dual gives

\[
\frac{2\operatorname{bip}(G)}{N^2}
\le \frac{2}{25}+\delta_\star
<\frac{2}{25},
\]

and hence `bip(G)<N^2/25`.

For densities below `0.2486` or above `0.3197`, the exact-constant tail cases
of Balogh-Clemen-Lidicky, first stated for all sufficiently large orders,
transfer to every finite graph by applying them to `G[t]` and using the
blow-up identity. The two tail regions and the closed certificate band cover
all possible densities. Therefore

\[
\operatorname{bip}(G)\le \frac{N^2}{25}
\]

for every finite triangle-free graph.

## Exact verification

Two independent exact gates were run.

1. `verify_moment_vector.py` reconstructed the fixed moment vector from 87
   nonnegative rational rank-one Gram atoms and matched all 12,172
   coordinates.
2. `verify_exact_fixed_gram_dual.py` did not import the LP generator. It
   rebuilt every K7, K8, and Horn functional from finite combinatorial data,
   constructed an integer matrix with 11,560,170 nonzeros, audited its
   digit-split product on nine columns, and checked every rational dual
   residual and the exact objective.

The verifier printed:

```text
EXACT_DUAL_REPLAY_OK
descriptors=10188 nnz=11560170 digits=3
digit_product_audit=9/9
min_q_residual=0 at state=46
objective ~= -9.878886951679021e-04
global_closure=True
GLOBAL_CONCLUSION: beta(G)<=N^2/25 for every finite triangle-free graph G on N vertices
```

Both the JSON and pickle representations passed the replay.

## Scope and review status

This package contains an internally independently replayed exact certificate
and a complete deduction from the public envelope and density-tail results.
It has not yet received external peer review. Until that review is complete,
it should be described as a reproducible computer-assisted proof candidate,
not as an externally verified solution.

## References

- A. Ferudun, *The Erdos n^2/25 max-cut conjecture for small multiples of
  five, via a per-root-MaxCut envelope and blow-up integrality*,
  arXiv:2606.28041v1, <https://arxiv.org/abs/2606.28041>.
- J. Balogh, F. C. Clemen, and B. Lidicky, *Max Cuts in Triangle-Free
  Graphs*, arXiv:2103.14179, <https://arxiv.org/abs/2103.14179>.
- T. F. Bloom, Erdos Problem 23, <https://www.erdosproblems.com/23>.
