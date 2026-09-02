# Routes to an unconditional proof of Erdős Problem 23

Companion to `INDEPENDENT_REVIEW.md`. Erdős Problem 23 is **open**. Nothing
below is a proof or a guarantee that one exists along these lines — it is what
the measurements in the review say about which routes are still alive, what
each would cost, and which are provably dead.

## 1. What is already unconditional

Two of the three pieces are done and correct.

* **The tails.** Balogh, Clemen and Lidický prove `bip(G) <= n^2/25` for
  `d_edge <= 0.2486` and for `d_edge >= 0.3197`, for `n` large
  ([arXiv:2103.14179](https://arxiv.org/abs/2103.14179)).
* **The transfer to every finite order.** `bip(G[t]) = t^2 bip(G)` is correct
  (the cut value of `G[t]` divided by `t^2` is multilinear in the twin-class
  fractions, and a multilinear function on a cube is maximised at a vertex).
  Applying a large-`n` tail theorem to `G[t]` and dividing by `t^2` removes the
  "`n` sufficiently large" hypothesis. Densities are blow-up invariant, so this
  costs nothing.

What is missing is exactly one statement:

> For every triangle-free graphon `W` with `d_edge(W)` in the closed band
> `[0.2486, 0.3197]`, `d_mono(W) <= 2/25`.

Everything else is bookkeeping.

## 2. The number that reframes the problem

The band is **not** a knife-edge. Searching weighted blow-ups of every
triangle-free graph on 9 vertices plus `C5`, `C7`, `C9`, `C11`, Petersen and
Grotzsch, with dilutions (`review/band_ceiling.py`):

```
sup{ d_mono(W) : W triangle-free, d_edge(W) in the band }  >=  0.063940 = 0.3197/5
2/25                                                       =  0.080000
```

and the `bip/m` ratio never exceeded `1/5` anywhere in the family. So inside
the band the conjecture has roughly **20% of slack**. A certificate there does
not have to be sharp — it has to be a relaxation that loses less than about a
quarter. That is a far softer requirement than at the extremal `C5` blow-up,
which sits at `d_edge = 0.4`, outside the band, and is handled by the tails.

Two immediate corollaries, both useful as build-time assertions:

* No valid certificate can have `delta < -0.016060`. A more negative objective
  is proof of a bug, not of strength.
* A valid band certificate must give **positive dual weight to at least one
  band row**. With zero weight on both, the same multipliers certify the same
  bound with the band deleted, which would contradict `d_mono(C5) = 2/25`.

Note also that the target is `delta <= 0`, not `delta < 0`. `d_mono <= 2/25` is
the conjecture; strict negativity was never required, and chasing it is what
pushed the current certificate onto a broken row family.

## 3. Step 0 — SETTLED: the 7-root envelope is dead

`review/k7_leg_ceiling.py` and `review/grotzsch_envelope.py` found a Grotzsch
blow-up inside the band at which the 7-root per-root-MaxCut leg looks positive.
`review/envelope_maxcut_bound.py` turns that into a decidable question: writing
`min_c mono = diag + S - MaxCut`, a lower bound on the envelope needs a
*certified upper* bound on one weighted MaxCut per root type. Components up to
20 vertices are settled by enumeration; the rest need a bound.

A plain Goemans-Williamson dual was not enough — 4-6% loose against exact
enumeration, where the decision needs about 1%, leaving a gap of 30.9 envelope
units against a budget of 5.70. `review/maxcut_triangle_bound.py` closes it with
the triangle inequalities of the cut polytope, dualised: for `lam >= 0` and
`W' = W - sum_t lam_t A_t`,

```
min_x x^T W x  >=  min_{X psd, diag 1} <W', X> - sum_t lam_t
               >=  - sum_i mu_i - sum_t lam_t      whenever  W' + Diag(mu) >= 0,
```

so `(lam, mu)` is a certificate that is *checked*, never trusted. The inner SDP
is solved by the mixing method, `mu` is read off stationarity and repaired by a
uniform shift, and `lam` is driven by Polyak subgradient ascent with violated
triples separated from `X = R R^T` each round. Calibrated against exact
enumeration it closes 81-100% of the gap; on the four instances small enough to
check, three bounds are exactly the true MaxCut.

```
                       plain SDP        + triangle inequalities
total gap over 33 instances   30.8993  ->   1.7022
budget needed to prove leg7 > 0            5.7016
root 0 (128 vertices, 1093 edges)  18.37  ->   0.2695   (98.5% closed)

leg7 lower bound : +2.204235e-05     VERDICT: provably positive
```

At the graphon

```
Grotzsch blow-up, weights (10173, 9717, 8691, 10166, 7628, 8213, 7344, 9352,
                           12084, 10443, 6189) / 100000
d_edge = 0.319384   inside the closed band [0.2486, 0.3197]
d_mono = 0.060023   comfortably under 2/25 = 0.08, so the conjecture holds here
```

the 7-root envelope satisfies `leg7 >= +2.20e-05 > 0`. Since `eta <= leg7` is a
constraint of the LP and every in-band graphon is a feasible point once the Horn
and K8 rows are repaired, **every certificate in the published 7-root
per-root-MaxCut framework has `delta >= +2.2e-05 > 0`** — with any rule set, and
at any flag order, because `U7` is a property of the graphon and no relaxation
can certify a bound the quantity itself violates.

That closes the route, and it explains the `+4.8558e-05` of
[arXiv:2606.28041](https://arxiv.org/abs/2606.28041): that objective is
positive not because the relaxation is slack but because the envelope it bounds
genuinely exceeds `2/25` inside the band. Integrality of `bip` was not a
convenience there, it was unavoidable — and it is why that route stops at
`N <= 200` and cannot be pushed to all `N` by adding rows.

One caveat, stated precisely: the pipeline is float64 throughout, and the PSD
checks use a tolerance of `1e-9 * max|W'|`. The margin is `4.0` envelope units
against accumulated rounding of order `1e-12`, so the conclusion is numerically
overwhelming but is not yet an exact-rational proof. The upgrade is mechanical:
rational class weights make every `N^sigma` rational, scaling makes the MaxCut
instances integral, and the `W' + Diag(mu) >= 0` checks become rational
`LDL^T`. Worth doing before this appears in a manuscript.

## 3b. The 8-root envelope is the cheap next move, not the coloured algebra

A correction to the reasoning above, and it changes the order of the plan.

`U_k` — the average over root types of the best profile-rule cut — is squeezed
between `d_mono` below (validity: every rule is a genuine colouring) and
whatever the certificate needs above. Larger roots express **more** colourings,
so `U_k` decreases toward `d_mono` as `k` grows: an 8-root envelope is
**tighter** than a 7-root one, not looser. The 7-root envelope going above
`2/25` at Grotzsch therefore says nothing about the 8-root one. The repository's
instinct — add 8-root rows to push the objective down — was directionally
right; only the data was broken.

At the balanced `C5` blow-up this is forced: `d_mono = U7 = 2/25` exactly and
`d_mono <= U8 <= U7`, so a correctly generated `U8` is exactly `2/25` there.

`review/u8_envelope.py` measures the shipped `U8` and the Aut-averaged repair of
it, by local search over rules (an upper bound on each per-root minimum):

```
graphon                    d_edge    d_mono   U8 shipped  U8 repaired
Grotzsch(leg7-worst)       0.3194  0.060023    0.045205    0.064016
Grotzsch(dmono-max)        0.3196  0.061458    0.044136    0.061200
Petersen(dmono-max)        0.3065  0.060883    0.041809    0.056847
tf9#1329(dmono-max)        0.3186  0.063394    0.039917    0.056820
...
balanced C5 (control)      0.4000  0.080000    0.055756    0.072201
```

Two readings. The repair does most of the work — at `C5` it lifts `U8` from
`0.0558` to `0.0722`, and the total pair mass reproduces `d_edge = 0.400000000`
exactly, so the decomposition is complete as a count and the defect really is in
the profile labels. But it is not enough: at `C5` the patched value is still
`9.75%` below the `0.08` validity forces. The tables had to be **rebuilt**, not
patched.

`review/rebuild_u8.py` does that, trusting nothing that ships with the
certificate. It re-derives the order-10 catalogue by extension from the 1,897
order-9 graphs (**12,172**, as required), identifies it with the certificate's
state numbering through the vertex-deletion profiles (a **bijection**, as
required), rebuilds the 107 order-7 and 410 order-8 root types, and for every
order-10 state and every edge computes the root, its type, an explicit
isomorphism to that type, and the two endpoint profiles under it — then averages
over `Aut(tau)`, which is exactly the sum over all isomorphisms the shipped data
gets wrong. 320,476 ordered (edge, root) incidences.

The acid test is the balanced `C5` blow-up, where `d_mono = U7 = 2/25` forces
`U8 = 2/25` exactly:

```
rebuilt U8 = 0.080000000   target 2/25 = 0.080000000   error -5.551e-17
```

Machine epsilon. The rebuild is correct, and the shipped tables are definitively
wrong. With it:

```
     graphon    d_edge    d_mono          U8     U8 - 2/25   valid  useful
    Grotzsch  0.319384  0.060023    0.068542   -0.011458    True    True
    Petersen  0.300000  0.060000    0.060653   -0.019347    True    True
```

At the very graphon where the 7-root envelope is provably **above** `2/25` (§3),
the correctly rebuilt 8-root envelope clears it by `0.0115`; at Petersen it sits
`1.1%` above `d_mono`, i.e. nearly tight. **The envelope route is alive at eight
roots.** That is what the repository was reaching for, and the only thing that
was wrong with it was the data.

Order of operations, then:

1. ~~regenerate the decomposition~~ — done, `review/rebuild_u8.py`, validated at
   `C5` to machine epsilon;
2. scan `sup U8` over the band properly. Two points are not a band: hill-climb
   the weights to *maximise* `U8` over every base, as `band_ceiling.py` does for
   `d_mono`, and confirm the supremum stays under `2/25`;
3. build the order-10 LP on the rebuilt rows, with positive dual weight on a
   band row, and check the objective lands at or below zero;
4. rationalise the dual and verify it in exact arithmetic, with the row-validity
   gate of §6 wired into the generator;
5. write the interfaces up as proofs: validity of the 8-root envelope with
   explicit quantifiers over labellings, the primal-feasibility argument, and
   the citation of the density-tail theorem with matching normalisation;
6. only if step 2 or 3 fails, go to §4.

## 3c. The LP on the rebuilt rows — first runs

`review/build_u8_lp.py` builds the order-10 band LP on the rebuilt 8-root rows
and grows it by cutting planes: K8 rules by weighted MaxCut at each optimum,
rooted-Horn 5-cycles by local search on the all-pairs decomposition
(`review/rebuild_u8_allpairs.py`, 1,095,480 incidences), and PSD cuts from the
most negative eigenvectors of the order-9 flag moment matrices lifted through
the deletion map. Every Horn row and PSD cut is gated — it must be nonnegative
at the `C5`, Grotzsch and Petersen state vectors before it enters — so the
mistake that sank the shipped certificate cannot recur silently. The orbit
bookkeeping reproduces `d_edge` from the constant rule to `1e-16`, and the
moment matrices are PSD to `1e-18` at the genuine graphons and indefinite (to
`-9e-4`) at a weak LP optimum, so the cuts have something to bite on.

| run | families | iterations | `eta` |
|---|---|---:|---:|
| 1 | band, mass, fixed Gram row, K8 rules | 40 | `+0.2157` |
| 2 | + Horn rows | 60 | `+0.1818` |
| 3 | + PSD eigenvector cuts | 15 (wall clock) | `+0.0564`, still falling |
| 4 | same, with checkpoints / capped additions / row dropping | 27 | `+0.0567` at 4,158 rows, still falling |
| 5 | + root-PSD cuts + order-10 flag blocks (types 0, 2, 4, 6) | 2 | `+0.0535`; 20 min per simplex solve |
| 6 | same, interior-point solver | 17 | `+0.0293` at 8,213 rows, gaining 0.0007 per iteration |
| 7 | same, eight eigenvectors per order-10 block per iteration | running | resumed from run 6 |

Run 3's trajectory — `0.2397, 0.2283, 0.2048, 0.1805, 0.1529, 0.1342, 0.1142,
0.0958, 0.0792, 0.0659, 0.0564` — is the first time this LP has been seen to
head toward zero, and it says two things. The fixed Gram row on its own is
nearly useless (run 1: `q` simply piles onto the densest states); the strength
is in the Horn rows and, above all, in letting the LP choose its own PSD
combinations rather than inheriting one vector. And the obstacle is now
engineering, not mathematics: the rows are nearly dense over 12,172 columns,
run 3 had reached ~15,000 of them, and the fifteenth solve took a quarter of an
hour. The current run adds per-iteration checkpoints, capped additions, and
dropping of rows that stay slack, and resumes across wall-clock limits.

Run 4 reached run 3's value with a quarter of the rows (the LP never exceeded
4,200 rows), at the price of slower progress per iteration; its band row stayed
inactive throughout (dual `0`), i.e. the LP is still fighting states that no
graphon realises rather than the band itself. Two facts settled on the side:

* `review/verify_u8_certificate.py` is validated. On a 3,106-row checkpoint the
  exact rational `delta` is `+0.1257866` against the solver's `eta =
  +0.1257862`; the difference is the rationalisation of the duals and the
  seventeen roots whose rule multipliers had to be scaled up to sum to one.
  (A first version discarded the Horn and PSD multipliers because HiGHS
  reports non-negative duals on rows at their lower bound; fixed.)
* `review/u8_band_scan.py`: the rebuilt 8-root envelope, maximised over the
  band at fifteen blow-up points (C5, Petersen, Grotzsch, twelve dense order-9
  bases), peaks at `0.069219` (Grotzsch), below `2/25 = 0.08`.

Two stronger cut families are now in, both valid by the same argument as the
moment blocks (an integral of `phi phi^T` is positive semidefinite):

* **root-PSD** (`--root-psd`): the 8-root pair matrices `M_tau(q)` themselves
  are completely positive at every graphon, so their negative eigenvectors give
  rows `v^T M_tau(q) v >= 0` through the orbit machinery already built. At the
  run-4 point they are violated to `-3.2e-3`.
* **order-10 flag blocks** (`--blocks10 0246`, `review/blocks10.c`): the
  genuine order-10 moment blocks for the 48 types of size 0, 2, 4, 6 (widths up
  to 2,445), computed exactly as integer counts by a small C program from the
  12,172 states. They are positive semidefinite to `1e-15` at the C5, Grotzsch
  and Petersen blow-ups and reproduce their edge densities, and the type-0
  block is indefinite to `-3.4e-3` at the run-4 point: seventy times the
  violation the lifted order-9 blocks show. The repository's data only ever
  contained four order-9 blocks; this is the missing bulk of the flag SDP.
  Cut vectors are rounded to integers over `1e6`, so the exact verifier
  recomputes every row bit for bit from the same integer data.
  `review/check_blocks10.py` re-derives the blocks of twelve sampled states in
  plain Python (its own canonical forms and flag indexing) and compares
  indexing-free signatures: 394 nonzero blocks agree.

The order-10 rows are dense over the 12,172 states, and the dual simplex
slowed to twenty minutes per solve at 5,000 rows. HiGHS's interior-point
method solves the same LP in two to three minutes (`--solver ipm`), and any
nonnegative interior duals are as good as vertex duals for the exact
certificate, so run 6 uses it.

Run 6's dual says where the weight is. At `eta = +0.0293` the 1,520 cuts from
the size-6 order-10 blocks carry a dual mass of `1.4e5` against `3.1e3` for
the size-4 blocks, `42` for size 2, `1.3e3` for the four lifted order-9 blocks
and `243` for the 8-root pair matrices (scales differ, but the ranking is
robust: 1,477 of the 1,520 size-6 cuts are active). The band row is still
inactive, so the LP's worst case sits at `d_edge = 0.28`, inside the band,
on a `q` that no graphon realises. Two consequences:

* An SDP with the small blocks explicit (`review/sdp_u8.py`, SCS: the 410
  pair matrices and the size-0/2 blocks as cones, 869,230 cone rows, 17M
  nonzeros) is the wrong lever: SCS did not finish one solve in two hours,
  and the blocks that matter -- size 6, widths up to 2,445 -- cannot be made
  explicit at all without a symmetry reduction. Parked.
* The lever is more size-6 cuts per iteration. Run 7 takes eight eigenvectors
  per block instead of three.

And the reason Kelley cutting planes crawl here, measured at the run-6 point:
the 48 order-10 blocks have 11,804 negative eigenvalues between them (the
size-6 empty type alone has 1,084 of its 2,445, the size-2 non-edge block 116
of 245), all tiny (`-1e-7` to `-2e-4`), summing to `-3.1e-3`. Three to eight
eigenvector cuts per block per iteration against twelve thousand negative
directions is why each iteration buys `0.001`. Two remedies are in the builder
now: the eight-per-block separation of run 7, and in-out separation
(`--in-out`): cuts are separated at `lam * q_out + (1 - lam) * q_in`, `q_in` a
mixture of 43 graphon state vectors (the three seeds plus forty random
weighted blow-ups of order-7 triangle-free graphs) at which every valid row
holds, and kept only if they still cut the LP optimum; such cuts are deeper
and the literature reports order-of-magnitude fewer rounds.

What a negative `eta` would and would not mean, stated before the number is
in: it would be a floating-point LP certificate over row families each of
which is valid by construction, i.e. the strongest computer evidence for the
band so far, and `review/verify_u8_certificate.py` turns it into an exact
rational `delta` from integer data. It would still not be a theorem until the
row validities and the primal-feasibility argument are written as proofs and
the tail theorem is cited with matching normalisation (§7).

## 4. Fallback — the two-coloured flag algebra

The fallback if §3b fails, and the formulation Balogh–Clemen–Lidický already use.
Carry the cut as part of the structure instead of trying to recover it from the
root.

**Objects.** Triangle-free graphs with their vertices 2-coloured, up to
isomorphism and colour swap. The limit object is a pair `(W, S)` with `W` a
triangle-free graphon and `S` a measurable subset. Then

```
d_mono(W, S) = ∫∫ W(x,y) [ 1_S(x) 1_S(y) + 1_{S^c}(x) 1_{S^c}(y) ] dx dy
```

is a **linear** functional of coloured subgraph densities. That single change
removes the adaptivity leak that killed the K8 leg: there is no minimum over
rules any more, so nothing can adapt to the pair being scored.

**The constraint that does the work — local optimality.** `S` is a maximum cut,
so no vertex improves by flipping sides: for a.e. `x`,

```
∫_S W(x,y) dy  <=  ∫_{S^c} W(x,y) dy       for x in S,
∫_{S^c} W(x,y) dy  <=  ∫_S W(x,y) dy       for x in S^c.
```

These are rooted inequalities with a single coloured vertex as the type — the
natural, valid analogue of what the per-root envelope was trying to fake. They
generalise: no *pair* improves, no set of size `k` improves, each giving more
rooted rows at larger types. Cheap to add, and they are where the extremal
structure enters.

**The rest of the LP/SDP.** The band `LO <= d_edge <= HI` as two ordinary rows
(with the assertion of §2 that at least one carries positive weight); one PSD
block per type; objective `max d_mono`. If the optimum over the closed band
comes out `<= 2/25`, the band is closed and, with §1, the theorem follows for
every finite triangle-free graph.

## 5. What it costs

`review/count_coloured_flags.py` counts the objects by Burnside over
`Aut(G) x Z2`:

| order | triangle-free | 2-coloured | growth |
|---:|---:|---:|---:|
| 5 | 14 | 106 | |
| 6 | 38 | 538 | 5.1x |
| 7 | 107 | 3,111 | 5.8x |
| 8 | 410 | 25,143 | 8.1x |
| 9 | 1,897 | **265,407** | 10.6x |
| 10 | 12,172 | ~3,000,000 | — |

The uncoloured column reproduces 1, 2, 3, 7, 14, 38, 107, 410, 1897
(OEIS A006785), which is a self-check on the enumeration and the canonical form.

This is the real trade-off, and it is worth stating plainly: the envelope trick
exists precisely to *avoid* carrying the colouring, which is what let the
current work run at order 10 with only 12,172 states. Going coloured buys
tightness and costs about two orders of resolution — order 8 (25k variables) is
comparable in size to the present order-10 uncoloured LP, order 9 (265k) is
heavy but has been done in the flag-algebra literature, order 10 is out of
reach.

So the practical question is whether order 8, or order 9, plus local optimality
plus the band restriction, is strong enough. The calibration:

* BCL's global bound is `n^2/23.5`, i.e. `d_mono <= 0.08511`.
* The band needs `d_mono <= 0.08000` — a **6.4% relative** improvement.
* The truth in the band is about `0.06394`, so their bound is already about
  33% above the truth there.

A 6.4% tightening from a higher order plus band-restricted rows is a plausible
ask, not an obviously hopeless one. It is also not a guarantee: this is where
the actual research risk lives, and no amount of engineering removes it.

## 6. Step 2 — making it unconditional

A floating SDP optimum is not a proof. The repository's existing machinery is
genuinely good here and is worth reusing wholesale:

1. rationalise the dual (round, project back to feasibility, repair);
2. rebuild every row from exact combinatorial data rather than a stored float
   matrix;
3. check all residuals and the objective in exact rational arithmetic;
4. ship chunked artefacts with checksums and a one-command replay.

The one new ingredient over the current LP setting is PSD certification: an
SDP dual needs an exact `LDL^T` (or rational Cholesky) of each block to prove
positive semidefiniteness, not just residual checks. Standard, but it is the
step to plan for.

Add the assertions the current build lacked:

* no row enters the pool unless it is nonnegative at a battery of exact graphon
  state vectors (`C5`, Petersen, Grotzsch, unbalanced `C5`) —
  `review/graphon_state_vector.py` produces them;
* reject any dual with zero weight on both band rows;
* reject any `delta < -0.016060`;
* for the coloured setting, check that the rooted pair-density matrices are PSD
  at those same test points before solving anything, which is what would have
  caught the present defect on day one.

## 7. What must still be proved by hand

* validity of the local-optimality rows (short);
* that the coloured limit object is correct and `d_mono` of the limit equals
  `lim 2 bip(G_n)/N^2` (max cut is a cut-metric-continuous parameter);
* the blow-up identity (already correct in the repository);
* the BCL statement quoted with matching normalisation, and an explicit check
  that the closed band plus the open tails cover every density.

## 8. Routes that are provably closed

* **`bip(G) <= |E(G)|/5` in general.** False. Alon's constructions give
  triangle-free graphs whose max-cut surplus over `m/2` is only `Θ(m^{4/5})`,
  so `bip/m -> 1/2`. A density hypothesis is mandatory; the counterexamples are
  sparse and land harmlessly in BCL's low-density tail. Note also that the
  density-restricted version *is* tight inside the band — the diluted `C5`
  attains `bip/m = 1/5` there — so it is a strictly harder target than
  `d_mono <= 2/25`, which has 0.016 of slack. Aim at the weaker one.
* **Any certificate with `delta < -0.016060`.** Impossible, by §2.
* **Any certificate with zero dual weight on both band rows.** Impossible, by §2.
* **The 7-root per-root-MaxCut envelope.** Settled in §3: `U7 > 2/25` at an
  in-band graphon, so no rule set and no flag order can rescue it. This is the
  framework of arXiv:2606.28041, and it is why that paper needs integrality.

## 9. Honest summary

The realistic plan is: regenerate the 8-root decomposition properly and measure
`sup U8` over the band (§3b) — days of work, and it may well be enough, since
larger roots make the envelope tighter and the shipped one is only broken, not
hopeless. In parallel, hand the 33 exported MaxCut instances (§3) to a real
solver to learn whether the 7-root envelope is genuinely dead. Only if the
8-root envelope also lands above `2/25` in the band does this become a
two-coloured flag SDP project (§4), with order 8 first and order 9 as the
expensive fallback. Either way the dual has to be rationalised and shipped with
an exact replay.

The encouraging part is §2: the band has real slack, so this is a question of
relaxation strength rather than of hitting an exact extremal configuration. The
sobering part is §5: the people who built the current best bound are expert
flag-algebra practitioners, they reached `n^2/23.5`, and the band held. A 6.4%
improvement is the whole remaining distance, and there is no guarantee that two
more orders deliver it.
