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

## 3. Step 0 — settle whether the envelope route is dead (cheap, decisive)

`review/k7_leg_ceiling.py` and `review/grotzsch_envelope.py` found a Grotzsch
blow-up inside the band (`d_edge = 0.3194`, `d_mono = 0.0600`, comfortably
inside the conjecture) at which the 7-root per-root-MaxCut leg is **positive**:

```
leg7, certificate's rule pool                      : +5.503e-05
leg7, multi-start local search over ALL rules      : +3.158e-05
leg7, 120 restarts with Kernighan-Lin sweeps       : +3.142e-05   (converged)
control at Petersen: pool -4.448e-05 -> all rules  : -1.280e-04   (the search bites)
```

If that positivity is real, **no per-root-MaxCut envelope certificate at any
flag order can prove the conjecture on the band** — the obstruction is in the
envelope, not in the relaxation strength, and the route must be replaced rather
than strengthened. It would also explain the `+4.8558e-05` of
[arXiv:2606.28041](https://arxiv.org/abs/2606.28041) as a genuine obstruction
rather than slack, and hence why that paper needed integrality of `bip` and
stopped at `N <= 200`.

A local search only upper-bounds the minimum over rules, so this is evidence,
not proof. Settling it needs an exact MaxCut or an SDP dual bound on the 107
weighted instances (up to 128 profile classes each) that
`review/grotzsch_envelope_hardened.py` already builds. The cheap spectral bound
is four hundred times too loose. **Do this before writing any more certificate
code** — it decides everything downstream, and it is hours of work, not weeks.

## 4. Step 1 — the two-coloured flag algebra

The sound replacement, and the formulation Balogh–Clemen–Lidický already use.
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
* **The per-root-MaxCut envelope**, if Step 0 confirms the Grotzsch positivity.

## 9. Honest summary

The realistic plan is: settle Step 0; if the envelope is dead, build a
band-restricted two-coloured flag SDP at order 8, with local-optimality rows,
and see how far under `0.08` it lands; if order 8 is short, go to order 9 and
budget serious compute; rationalise the dual and ship an exact replay.

The encouraging part is §2: the band has real slack, so this is a question of
relaxation strength rather than of hitting an exact extremal configuration. The
sobering part is §5: the people who built the current best bound are expert
flag-algebra practitioners, they reached `n^2/23.5`, and the band held. A 6.4%
improvement is the whole remaining distance, and there is no guarantee that two
more orders deliver it.
