# Independent review of the order-10 certificate

Review date: 2026-08-23. Commit reviewed: `e35fca0`.

## Verdict

**The certificate does not prove Erdős Problem 23, and this repository should
not be described as a proof candidate for the full conjecture in its present
form.**

The machine replay is real: the stored multipliers are an exactly feasible dual
of the finite linear program that `verifier/verify_exact_fixed_gram_dual.py`
reconstructs, and its objective is exactly
`delta* = -9.878886951679021...e-04 < 0`. I reproduced this from a clean
checkout (`ARTIFACT_INTEGRITY_OK`, `PASS: public exact moment vector was
independently reconstructed`, `EXACT_DUAL_REPLAY_OK` for both the JSON and the
pickle certificate).

What fails is the step the replay cannot check: the passage from "this LP has a
negative dual objective" to "`bip(G) <= N^2/25`". That passage needs every LP
row to be a valid inequality for genuine triangle-free graphons, and needs the
LP optimum to dominate `d_mono(W) - 2/25`. Both requirements are false, and
they are false on the very graphon the conjecture is sharp for.

Concretely, evaluated at the exact order-10 state vector of the `C5` blow-up
graphon (`d_edge = 2/5`, `d_mono = 2/25`), using the verifier's *own* row
reconstruction and exact rational arithmetic:

| check | claim in the manuscript | measured |
|---|---|---|
| K8 per-root MaxCut envelope `U8` | `>= d_mono` | `0.055899136 < 0.08` |
| rooted-Horn rows | `>= 0` for every triangle-free graphon | 38 of 5,013 are **strictly negative** (min `-7.74e-03`), carrying 20.4% of the Horn dual weight |
| rooted pair-density matrices behind those rows | a nonnegative mixture of `p p^T`, hence PSD | 14 of 410 root types have a **negative eigenvalue** (min `-5.06e-03`), with an exact rational witness |
| dual weight on the two edge-density band rows | the certificate covers `[0.2486, 0.3197]` | both multipliers are exactly **zero** — the band is never used |

The same failures occur *inside* the certificate's own band: for the Petersen
graphon (`d_edge = 3/10`, `d_mono = 3/50`) 330 rooted-Horn rows are negative and
`U8 = 0.055884 < 0.06`; for an unbalanced `C5` blow-up with weights
`(3/5, 1/10, 1/10, 1/10, 1/10)` (`d_edge = 3/10`, `d_mono = 1/50`) 17 rows are
negative and `U8 = 0.013718 < 0.02`.

Everything below is reproducible with the two scripts in this directory.

---

## 1. What the replay actually verifies

Reading `generation/generate_fixed_gram_lp.py:420-484` together with
`generation/solve_fixed_gram_dual.py`, the primal LP is

```
maximise  eta
over      q in R^12172 (q >= 0),  eta free,  u7 in R^107 (>=0),  u8 in R^410 (>=0)
subject to
  (0)  <d_edge, q>  <=  HI = 3197/10000            [dual multiplier y0]
  (1)  <d_edge, q>  >=  LO = 2486/10000            [dual multiplier y1]
  (2)  <m, q>       >=  0        (fixed Gram)      [dual multiplier y2]
  (3)  eta          <=  sum_sigma u7_sigma         [dual multiplier y3, "K7 leg"]
  (4)  eta          <=  -2/25 + sum_tau u8_tau     [dual multiplier y4, "K8 leg"]
  (5)  <1, q>       =   1                          [dual multiplier rho]
       u7_sigma  <= <g_{sigma,c}, q>   for every selected K7 rule c
       u8_tau    <= <h_{tau,c},   q>   for every selected K8 rule c
       <horn_{tau,C}, q> >= 0          for every selected rooted-Horn cycle C
```

The verifier checks, exactly, that `(y0..y4, y_rows, rho)` is dual feasible:
nonnegativity, `y3 + y4 = 1` (`static[3] + static[4] = D`), the 517 per-root
multiplier lower bounds, all 12,172 coordinate residuals, and the objective
`y0*HI - y1*LO - (2/25)*y4 + rho = delta*`. All of that is correct and I
reproduced it.

Dual feasibility means exactly one thing: `eta <= delta*` for every point of the
LP's feasible region. It says nothing about `bip(G)` until one shows that the
true value `d_mono(W) - 2/25` is *attained inside that region*. That is the
claim in §3 of the manuscript ("Every genuine band step graphon induces a
feasible point"), and it is what fails.

## 2. Finding A — the certificate carries zero weight on the density band

```
static_numerators = [0, 0, 53956267978680, 1, 999999999999]
                     ^  ^
                     |  y1 (d_edge >= 0.2486)
                     y0 (d_edge <= 0.3197)
```

Both band multipliers are exactly zero. Deleting rows (0) and (1) from the LP
therefore leaves the same multipliers dual feasible with the same objective, so
the bound the certificate proves is *independent of the edge density*. The
middle band `[0.2486, 0.3197]` — the entire point of the construction, and the
only region the conjecture is still open in — plays no role at all.

This is already fatal on its own, before any computation. If the remaining rows
were valid and the legs dominated `d_mono - 2/25`, the certificate would prove
`d_mono(W) < 2/25` for **every** triangle-free graphon, including the balanced
`C5` blow-up, where `d_mono = 2/25` exactly (manuscript §1; `bip(C5[n]) = n^2`).
A correct band certificate cannot have this property: its bound has to become
false somewhere outside the band, so at least one band multiplier must be
positive.

For contrast, the referenced prior work (arXiv:2606.28041) reports a *positive*
`delta ~ +4.8558e-05` for the 7-root envelope on the same band, and needs
integrality of `bip` to reach a conclusion for `N <= 200`. A sign flip to
`-9.88e-04` together with zero band weight is not a strengthening of that
result; it is a symptom.

## 3. Finding B — the K8 per-root MaxCut envelope does not dominate `d_mono`

The manuscript's §3 justification is that a Boolean profile rule `c` "defines a
genuine global two-coloring and hence a valid per-root MaxCut row". If that
were so, then for *any* assignment of one rule `c_tau` to each root type,

```
sum_tau <h_{tau,c_tau}, q>  =  E_R[ mono-edge density of the colouring chi_{R,c_tau(R)} ]  >=  d_mono(W),
```

because each summand is the mono-edge density of an honest two-colouring, and
`d_mono` is the minimum over all two-colourings. In particular the per-root
minimum `U8 = sum_tau min_c <h_{tau,c}, q>` would satisfy `U8 >= d_mono`.

Measured on the balanced `C5` graphon, with the certificate's own pool:

```
sum over the 410 roots of the pool minimum : U8 = 109178/1953125 = 0.055899136
sum over the 410 roots of the pool maximum :      0.067798016
d_mono(C5)                                 :      0.08
```

Every single rule in the pool undershoots. Since the true `U8` (minimum over
*all* rules) is at most the pool minimum, `U8(C5) <= 0.0559 < 2/25`. The K8 leg
therefore permits `eta = -0.0241`, whereas the quantity to be bounded is
`d_mono - 2/25 = 0`. The K8 leg carries dual weight `1 - 10^-12`; it is
essentially the whole certificate.

Two controls confirm the measurement rather than the method:

* the same computation with the constant rule reproduces `d_edge = 2/5`
  exactly, and the all-pairs total reproduces `1` exactly;
* the **K7** leg, computed the same way from `k7_compact_v1`, gives exactly
  `2/25` at the `C5` graphon — i.e. the 7-root envelope is tight at the
  extremal point, exactly as arXiv:2606.28041 claims. So the framework is not
  wrong in principle; the 8-root data is.

The K7 leg passes the validity test at all three test graphons, the K8 leg
fails at all three.

## 4. Finding C — the rooted-Horn rows are not valid, and why

The copositivity inequality itself is correct: for `x >= 0`,
`(sum_i x_i)^2 - 4 sum_i x_i x_{i+1} >= 0` on a 5-cycle is Motzkin–Straus with
`omega(C5) = 2`. The defect is in its flag-algebra realisation.

For a genuine graphon, the rooted pair-density matrix

```
M^tau_{a,b} = E_R[ 1(type(R) = tau) * p_a(R) * p_b(R) ]
```

is a nonnegative mixture of rank-one matrices `p p^T` with `p >= 0` — hence
symmetric, positive semidefinite, and in fact *completely positive*, which is
precisely what makes `<C, M^tau> >= 0` hold for every copositive `C`.

Rebuilding `M^tau` from `public_flagsdp_data/u8_decomp_all.pkl` at the `C5`
graphon (via the verifier's own `exact_pair_matrices`): the matrices are
symmetric, but **14 of the 410 root types are indefinite**, the worst being
root 401 with smallest eigenvalue `-5.06e-03` against a root mass of `0.0645`.
An exact rational witness is printed by the audit script:

```
root 401 (width 33):  x^T M x = -30911132658357/6103515625000000 < 0,  x^T x ~ 1
```

So `u8_decomp_all.pkl` does not encode a rooted pair density: the profile of a
free vertex is not a function of the root alone. That single defect explains
both symptoms — it breaks complete positivity (hence the negative Horn rows,
whose bad roots are exactly `{0, 1, 7, 18, 43, 53, 123, 199, 230, 364, 367,
401, 405, 406}`) and it breaks the "one rule = one global colouring" reading of
the K8 rows.

The invalid Horn rows are not marginal: they carry 20.4% of the total Horn dual
weight at `C5` (23.7% at Petersen), and one of them (row 2305, root 43) has
multiplier `18.837`, among the largest in the certificate.

## 5. Finding D — the trust boundary is wider than advertised

`README.md` says "The dual verifier does not import the LP generator". True,
but the verifier takes the entire combinatorial semantics on trust from five
opaque binary inputs:

```
public_flagsdp_data/cache_n9.pkl        public_flagsdp_data/c5lift_cache.npz
public_flagsdp_data/u8_decomp.pkl       public_flagsdp_data/u8_decomp_all.pkl
k7_compact_v1/root_XXX.npz
```

Nothing in the verifier checks that the 12,172 columns are the triangle-free
order-10 graphs, that `u8_decomp*.pkl` are the 8-rooted pair decompositions, or
that the profiles are consistently labelled. The modules that produce them —
`flag_cutgen`, `cutting_plane_u8`, `envelope_horn`, `run_k7b` — are in neither
archive and neither is `public_anc/horn_dual.pkl`, so `generate_fixed_gram_lp.py`
cannot even be run from the released material. The files that turn out to be
wrong are precisely two of these unregenerable inputs.

Smaller points, in decreasing order of importance:

1. `verify_exact_fixed_gram_dual.py` prints `GLOBAL_CONCLUSION: beta(G)<=N^2/25
   for every finite triangle-free graph`. The script verifies dual feasibility;
   it does not and cannot verify that sentence. The line should be removed or
   demoted to "dual feasible with negative objective".
2. With `--target-n 1` (what `scripts/run_full_replay.sh` passes) the
   `threshold`, `strict_margin` and integrality checks reduce to
   `objective < 2/25` and `(25/2)*objective < 1`, which any negative objective
   satisfies. They report nothing.
3. The two envelope legs use different normalisations — the K7 row carries its
   `-2/25` internally and is scaled by `1/10`, the K8 leg subtracts `2/25`
   outside. Because the LP takes the *minimum* of the legs, the mismatch makes
   the K7 leg nearly vacuous whenever the true value is negative, which is why
   the optimiser could put weight `1 - 10^-12` on the K8 leg. This should be
   stated and justified in the manuscript, not left implicit in the code.
4. The manuscript's claim that "a later exact certificate extended this finite
   range through N=240" has no citation.

## 6. What does hold up

To be clear about the parts I checked and found correct:

* **Lemma 2 (blow-up identity)** `bip(G[t]) = t^2 bip(G)`. The proof is right:
  the cut value of `G[t]` divided by `t^2` is multilinear in the twin-class
  fractions, a multilinear function on a cube is maximised at a vertex, so
  `MaxCut(G[t]) = t^2 MaxCut(G)`, and `|E(G[t])| = t^2 |E(G)|`.
* **Normalisation** `d_edge = 2|E|/N^2`, `d_mono = 2 bip/N^2`, and the
  identity `d_mono(W_G) = 2 bip(G)/N^2` with no finite-size error. Consistent
  throughout, and consistent with the Balogh–Clemen–Lidický thresholds
  (`C5` sits at `d_edge = 0.4`, in their high-density tail).
* **The density-tail transfer argument.** Applying a large-`n` tail theorem to
  `G[t]` and dividing by `t^2` is a valid way to remove the "`n` sufficiently
  large" hypothesis, given the blow-up identity.
* **The fixed-Gram moment row.** `<m, q> >= 0` held at all three test
  graphons (`+1.73e-04`, `+2.76e-05`, `+1.48e-03`), consistent with the
  manifest-PSD construction, and `verify_moment_vector.py` reconstructs all
  12,172 coordinates from 87 nonnegative rank-one atoms.
* **The K7 per-root MaxCut leg.** Valid at all three test graphons and exactly
  tight (`2/25`) at the `C5` blow-up.
* **The exact arithmetic and the replay infrastructure** — checksums, chunked
  artefacts, digit-split sparse product with its `int64` bound, the nine-column
  audit, both certificate encodings. No defect found.

## 7. Does this close Erdős Problem 23?

No. Erdős Problem 23 asks whether every triangle-free graph on `N` vertices can
be made bipartite by deleting at most `N^2/25` edges. It is open; the tails
`d_edge <= 0.2486` and `d_edge >= 0.3197` are settled by Balogh, Clemen and
Lidický (arXiv:2103.14179), and arXiv:2606.28041 settles balanced `C5`
multiples up to `N = 200`. The middle band is the open part.

This certificate does not touch the middle band: its dual assigns weight zero
to both band constraints, and the leg that carries all of its strength rests on
an 8-rooted decomposition that provably fails a necessary structural property
of rooted pair densities. The negative objective `delta*` measures the slack of
a broken envelope, not any property of triangle-free graphs.

## 8. Reproducing this review

```bash
python3 -m pip install -r requirements-replay.txt
./scripts/check_artifacts.sh            # confirms artefact integrity
./scripts/run_full_replay.sh            # confirms the exact dual replay passes

mkdir -p .review-work
tar -xzf .replay-work/erdos23_full_replay_bundle_2026-08-23.tar.gz -C .review-work
B=.review-work/erdos23_full_replay_bundle

# exact order-10 state vectors of three genuine triangle-free graphons
python3 review/graphon_state_vector.py --flagsdp $B/public_flagsdp_data \
    --base c5 --out .review-work/c5_balanced.pkl                 # d_edge 2/5,  d_mono 2/25
python3 review/graphon_state_vector.py --flagsdp $B/public_flagsdp_data \
    --base petersen --out .review-work/petersen.pkl              # d_edge 3/10, d_mono 3/50
python3 review/graphon_state_vector.py --flagsdp $B/public_flagsdp_data \
    --base c5 --weights 6,1,1,1,1 --out .review-work/c5_inband.pkl  # d_edge 3/10, d_mono 1/50

python3 review/check_row_validity.py \
    --certificate .replay-work/erdos23_global_exact_dual.json \
    --flagsdp $B/public_flagsdp_data --public-anc $B/public_anc \
    --k7-cache $B/k7_compact_v1 \
    --graphon .review-work/c5_balanced.pkl \
    --graphon .review-work/petersen.pkl \
    --graphon .review-work/c5_inband.pkl
```

The two scans behind sections 9 and 10 run the same way:

```bash
python3 review/verify_aut_fix.py --flagsdp $B/public_flagsdp_data \
    --graphon .review-work/c5_balanced.pkl        # root cause and its repair

python3 review/band_ceiling.py --flagsdp $B/public_flagsdp_data   # ~4 min

python3 review/k7_leg_ceiling.py \
    --certificate .replay-work/erdos23_global_exact_dual.json \
    --flagsdp $B/public_flagsdp_data --k7-cache $B/k7_compact_v1 --sample 60
```

`graphon_state_vector.py` needs no order-10 catalogue: it builds each blow-up
`H[a]` directly, identifies it through its vertex-deletion profile against the
repository's own 9→10 deletion lift (verified to separate all 12,172 states),
and sums exact multinomial weights. The resulting vectors are validated three
ways — they sum to `1`, their edge density matches the closed form
`2 sum_{uv in E} w_u w_v` exactly, and `D q_10 = q_9` holds exactly against an
independently built order-9 vector.

`check_row_validity.py` reuses the verifier's row reconstruction verbatim, so a
disagreement is a disagreement with the certificate, not with a re-implementation.

## 9. Root cause, and a fix that is verified to work

`u8_decomp*.pkl` stores exactly **90 contributions per order-10 state** — 45
free pairs times 2 orderings, i.e. one labelling of each 8-element root subset.
The K7 cache is normalised by `9P7 = 181440`: it sums over *all* ordered
injections of the 7 root vertices. That is the whole difference between the
family that works and the family that does not.

With a single labelling per subset, the labelling a canonical form assigns
depends on how the subset sits inside the state, hence on the two free
vertices. Writing `lambda(R,u,v) = alpha . lambda_0` with `alpha in Aut(tau)`,
the recorded matrix is `E[e_{alpha a0} e_{alpha b0}^T]` with `alpha` correlated
with `(u,v)`, rather than the true rooted pair density
`M_0 = E_R[1(tau) p p^T]`. Two consequences, both observed:

* `M` need not be PSD, so copositive `C` no longer gives `<C, M> >= 0` — the
  negative rooted-Horn rows;
* a "rule" no longer defines one colouring per root, it defines a colouring
  that can adapt to the pair being scored — so the minimum over rules can slip
  *below* `d_mono`, which is exactly what `U8` does.

Averaging the recorded pairs over `Aut(tau)` cancels the arbitrary `alpha`:

```
(1/|Aut|) sum_beta e_{beta alpha a0} e_{beta alpha b0}^T
    = (1/|Aut|) sum_gamma e_{gamma a0} e_{gamma b0}^T,
```

so the symmetrised matrix is `(1/|Aut|) sum_gamma gamma M_0 gamma^T`, a
nonnegative mixture of PSD matrices; and the symmetrised MaxCut row is an
average over `gamma` of `<h_{tau, c . gamma}, q>`, each term being an average
over labelled roots of the mono-edge density of **one fixed** two-colouring,
hence `>= d_mono` times the root mass. Summing over root types restores
`U8 >= d_mono`.

`review/verify_aut_fix.py` checks the first half numerically. At the `C5`
graphon:

```
non-PSD before Aut-averaging : 14  [0, 1, 7, 18, 43, 53, 123, 199, 230, 364, 367, 401, 405, 406]
non-PSD after  Aut-averaging : 0   []
|Aut| of the failing types   : {0: 40320, 1: 1440, 7: 5040, 18: 192, 43: 72, 53: 1440,
                                123: 48, 199: 144, 230: 720, 364: 32, 367: 16, 401: 16,
                                405: 72, 406: 1152}
```

Every failing type has a nontrivial automorphism group — the defect bites
exactly where the labelling is ambiguous — and averaging over that group
removes all 14 failures. The production fix is to regenerate `u8_decomp.pkl`
and `u8_decomp_all.pkl` summing over all `10P8 = 1,814,400` ordered injections
per state, the convention the K7 cache already uses; Aut-averaging the existing
tables is the equivalent cheap route.

## 10. What would still be needed after the fix

Repairing the decomposition makes the rows valid. It does **not** by itself
produce a proof, and the remaining gap is the hard part.

**a. The certificate must actually use the band.** A valid band certificate has
to become false outside the band, because `U7` is exactly `2/25` at the `C5`
blow-up (`d_edge = 0.4`). So at least one of `y0, y1` must be strictly
positive. Any solve that returns `y0 = y1 = 0` is a build failure by
construction, not a result.

**b. The K7 leg goes positive inside the band, at the Grotzsch graphon.**
`eta <= leg7` is a constraint of the LP, and once the Horn and K8 rows are
repaired every in-band triangle-free graphon is a feasible point, so a
certificate in the published **7-root framework** — all dual weight on the K7
leg, no K8 leg — satisfies `delta >= max{ leg7(W) : W in-band }`.
`review/k7_leg_ceiling.py` hill-climbs the class weights of blow-up graphons
over `C5`, `C7`, `C9`, Petersen, Grotzsch and 60 random triangle-free
9-vertex bases:

```
          leg7       base   d_edge   d_mono
  3.649113e-05   Grotzsch   0.3196   0.0556      <-- positive
 -4.448000e-05   Petersen   0.3000   0.0600
 -2.415751e-04   tf9#1128   0.3197   0.0453
 -2.914771e-04    tf9#663   0.3197   0.0416
 ...
 -1.443000e-03   C5, weights 1/4,1/4,1/12,1/3,1/12
```

The Grotzsch blow-up sits inside the band, satisfies the conjecture with room
to spare (`d_mono = 0.0556`, well under `2/25`), and yet drives the K7 leg
**above zero**. With this row pool, a 7-root certificate therefore cannot reach
`delta <= 0` at all, let alone `-9.88e-04` — it is stuck at `>= +3.65e-05`.

Enlarging the pool does not obviously rescue it. `review/grotzsch_envelope.py`
re-optimises the Grotzsch weights and then replaces the pool minimum by a
multi-start local search over **all** Boolean profile rules — one weighted
MaxCut per root type, up to 128 profile classes wide:

```
Grotzsch blow-up, d_edge = 0.319385, d_mono = 0.060024
  leg7, certificate's K7 pool      : +5.503e-05
  leg7, local search over all rules: +3.158e-05     <-- still positive
control, Petersen
  leg7, certificate's K7 pool      : -4.448e-05
  leg7, local search over all rules: -1.280e-04     <-- the search does bite
```

The control matters: at Petersen the same search improves the pool value by a
factor of three, so it is finding genuinely better rules. At Grotzsch it
improves the value by 43% and still cannot cross zero. Hardening the search to
120 restarts with Kernighan-Lin sweeps moves it only to `+3.142e-05`, i.e. it
has converged (`review/grotzsch_envelope_hardened.py`).

A local search only upper-bounds the true minimum, so this is strong evidence,
not proof, that the 7-root envelope genuinely exceeds `2/25` at that point. The
obvious cheap rigorous bound is far too weak to settle it — the spectral bound
`min_c mono >= diag + S/2 + (n/4) lambda_min(B)` gives `leg7 >= -2.75e-02`, four
hundred times looser than needed. Closing the gap needs a real MaxCut solver or
an SDP bound; see **c**.

That number sits right next to the `+4.8558e-05` reported in
arXiv:2606.28041. Read together, the two say something important: that paper's
positive `delta` is very probably not slack in its relaxation but a real
obstruction, which is exactly why it needed integrality of `bip` and could only
reach `N <= 200`. Whether the obstruction survives an *unlimited* K7 rule pool
is the open question in **c** below — the pool minimum is only an upper
estimate of the true envelope.

Note also that weighted `C5` blow-ups are *not* the binding configurations:
Petersen is thirty times closer to zero than the best of them, Grotzsch clears
zero entirely. Any serious search has to sweep a far richer family than `C5`
blow-ups.

Third — and this is the trap — a repaired K8 leg does *not* restore the claimed
objective's plausibility, it relocates the entire burden of proof onto the
8-root envelope. Larger roots make the envelope *smaller* (more root types,
each minimised independently), so `U8 <= U7` and the K8 leg is the binding one;
`eta = min(leg7, leg8) = leg8`. A `delta` near `-1e-03` is then attainable only
because of the 8-root leg — exactly the leg whose validity fails today and
whose repaired form nobody has yet exhibited. Note also that the two legs are
not on the same scale: `leg7 = (mass7/10) * (U7 - 2/25)` carries a
graphon-dependent factor `mass7/10` (0.014 at Petersen, 0.019 at `C5`, 0.17 at
an unbalanced `C5`) that compresses it toward zero, while `leg8 = U8 - 2/25`
does not. Any future manuscript has to justify that mismatch explicitly rather
than let the optimiser exploit it.

**b'. There is a hard ceiling on how negative any valid `delta` can be.**
"The LP is a relaxation" means that at every in-band triangle-free graphon `W`
the point `(q(W), eta = d_mono(W) - 2/25, ...)` is feasible, so
`delta >= d_mono(W) - 2/25` for every such `W`, hence

```
delta  >=  sup{ d_mono(W) : W triangle-free, d_edge(W) in [0.2486, 0.3197] } - 2/25.
```

No extra rows, higher flag order or better Gram block can get past that.
`review/band_ceiling.py` lower-bounds the supremum by hill-climbing the class
weights of blow-up graphons over every triangle-free graph on 9 vertices (all
1,897, read out of `cache_n9.pkl`) plus `C5`, `C7`, `C9`, `C11`, Petersen and
Grotzsch, and by allowing dilutions `theta*W` — which stay triangle-free, since
`t(K3, theta W) = theta^3 t(K3, W) = 0`, and pull a graphon above the band back
into it along its own ray. Every base is asserted triangle-free first (an
earlier draft of this scan used a wrong Grotzsch edge list that contained
triangles and duly reported a "counterexample" to the conjecture).

```
   d_mono       base   h   m   d_edge   theta   bip/m
 0.063940         C5   5   5   0.3197  0.7992  0.2000
 0.063940   Grotzsch  11  20   0.3197  0.9671  0.2000
 0.063864   tf9#1866   9  14   0.3197  0.9929  0.1998
 ...
best in-band d_mono found : 0.063940 = 0.3197/5
best bip/m ratio seen     : exactly 1/5, never exceeded
```

So `delta >= -0.016060` for any valid certificate, and the entire slack of the
conjecture inside the band is at most `2/25 - 0.06394 = 0.01606`. This does
*not* exclude the claimed `-9.88e-04` — it is only 6% of the available room —
which is worth stating plainly: the defect in this certificate is the validity
of its rows, not the magnitude of its objective. It also calibrates the prior
work: `+4.8558e-05` misses the target by `0.3%` of the available room, so the
7-root route is not obviously capped, it is close.

**c. The open question this route turns on** is whether
`sup{ U7(W) : W triangle-free, d_edge(W) in [0.2486, 0.3197] } <= 2/25`, where
`U7` is the envelope with the minimum taken over *all* Boolean profile rules,
not just the pooled ones. The Grotzsch result in **b** shows the certificate's
pool fails that test; it does not yet settle whether an unlimited pool would.
That is now the single most valuable thing to compute, and the Grotzsch
blow-up is the point to compute it at.

Each root type turns into one weighted MaxCut instance over its profile
classes, and those are large — up to 128 classes for the empty 7-root — so
exhaustive enumeration is out and a real MaxCut routine (or an SDP bound) is
needed. Concretely: run an exact solver, or a Goemans-Williamson style SDP
dual, on the 107 instances that `review/grotzsch_envelope_hardened.py` already
builds, and see whether `leg7` at the Grotzsch point is provably positive.
Two possible outcomes:

* the true `U7` at Grotzsch still exceeds `2/25` — then **no** 7-root envelope
  certificate at any flag order can prove the conjecture on the band, the
  obstruction is in the envelope itself, and the route has to be replaced
  rather than strengthened;
* the true `U7` drops below `2/25` — then the route survives, the current row
  pool is simply too small, and the target is `delta <= 0` (which suffices;
  strict negativity was never required).

**d. If the envelope is capped, the sound alternative is the two-coloured flag
algebra**: carry the max-cut colouring as part of the structure, so `d_mono` is
a *linear* functional of coloured densities rather than an adaptive minimum,
and add the local-optimality inequalities of the cut (no vertex improves the
cut by flipping). That is the formulation Balogh, Clemen and Lidický use, it
has no adaptivity leak, and it is the natural place to spend order-11 effort.

**e. Independently of which route is taken, the following must accompany any
future claim:**

1. A written proof, with explicit quantifiers over labellings, that
   `U_k(W) >= d_mono(W)` for the rows actually used.
2. A written proof that the LP is a relaxation: exhibit the feasible primal
   point `(q(W), eta = d_mono(W) - 2/25, u7, u8)` for an arbitrary band
   graphon, and check that *both* legs admit that `eta`. The two legs currently
   use different normalisations (K7 carries its `-2/25` internally and is
   scaled by `1/10`; K8 subtracts `2/25` outside), which is what let the
   optimiser park weight `1 - 1e-12` on one leg.
3. A validity gate in the cutting-plane loop: no row enters the pool unless it
   is nonnegative at a battery of exact graphon state vectors. This is the
   single highest-value change — a separator that emits invalid rows will
   always find a spurious optimum, and this one did.
4. A precise citation of the Balogh–Clemen–Lidický theorem with matching
   normalisation, and an explicit check that the closed band plus the open
   tails cover every density.
5. Either publish the generators for the four `public_flagsdp_data` files, or
   add verifier-side structural checks — state catalogue, decomposition
   completeness, PSD of the rooted pair densities — so the trust boundary
   matches what `README.md` describes.
6. Drop `GLOBAL_CONCLUSION` from the verifier output, or demote it to what is
   actually checked.
