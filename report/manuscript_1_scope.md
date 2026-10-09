# Manuscript 1 Scope: ICC Frontier Geometry and Stochastic Dynamics

## Manuscript status

**Theory freeze reached after H12 verification.** The first manuscript should now move from discovery to consolidation, literature positioning, figure/table construction, and writing. New mathematical branches are out of scope unless drafting exposes a specific logical gap.

## Working title

**Same ICC, Different Geometry: Directional Asymmetry and Stochastic Dynamics on an Intraclass-Correlation Frontier**

Alternative title:

**The Geometry of a Fixed ICC: Directional Asymmetry, Neutrality, and Stochastic Frontier Dynamics**

## Central manuscript question

What observation-space geometry and sequential dynamics are compatible with a fixed intraclass correlation coefficient when new observations are added to balanced groups?

The manuscript is not primarily about proposing ICC as a new clustering criterion. ICC is the deliberately fixed structural constraint used to expose a broader phenomenon: a scalar structural summary can conceal strongly direction-dependent admissible heterogeneity and nontrivial sequential dynamics.

## Central story

The paper should tell one story:

1. A fixed ICC threshold defines an anisotropic frontier in observation space.
2. The amount of admissible novelty depends on alignment with the existing between-group structure.
3. Repeated threshold-hugging updates collapse to an exact low-dimensional recurrence.
4. The recurrence is directionally asymmetric: orthogonality and mean-zero alignment are not generally neutral.
5. Under random alignment, the long-run exponent is governed by the nonlinear quantity \(E[h(C)]\), not \(E[C]\).
6. At stochastic neutrality, short-range dependence does not destroy neutrality but changes the distribution of the limiting structural scale.

The final point is the intended endpoint of the theory for Manuscript 1.

## Headline results to emphasize

The manuscript should foreground only a small number of results.

### Result A: the ICC frontier is directionally anisotropic

At a fixed state and ICC threshold, the admissible frontier radius increases with alignment to the current group-offset direction. The same ICC therefore permits substantially different observation-space heterogeneity depending on direction.

### Result B: repeated frontier-following has an exact normalized recurrence

Once the process is on the ICC threshold, absolute scale and detailed path history drop out of the normalized frontier. The dynamics reduce to a scalar multiplicative recurrence for group-offset geometry.

### Result C: structural neutrality is nonlinear and asymmetric

Orthogonality is not exactly neutral at finite group size. Equal positive and negative alignment also does not generally cancel. Neutral deterministic and mixed policies require a small or systematic opposing bias.

### Result D: stochastic drift is governed by \(E[h(C)]\)

For IID alignment, the almost-sure polynomial growth exponent is \(E[h(C)]\), not the ordinary mean \(E[C]\). Strict convexity of \(h\) implies that every non-degenerate mean-zero IID alignment distribution is expansive in the current regime.

### Result E: dependence controls neutral-boundary variability

For a stationary two-state Markov alignment process whose marginal distribution is Gotham-neutral, the structural scale remains asymptotically finite and positive under short-range dependence, while persistence changes the variance of the random limiting scale. This is the final major theoretical result for Manuscript 1 and has now been numerically verified.

## Supporting rather than headline results

The following results are important but should function as lemmas, corollaries, or supporting propositions rather than separate manuscript stories:

- zero-radius ICC headroom;
- the sufficient condition for a unique positive frontier crossing;
- the finite-horizon H4/H5 distinction;
- exact one-step neutral alignment;
- deterministic greedy neutral controller;
- numerical package-verification identities.

They support the main argument but should not each become a separate narrative branch.

## Scope exclusions for Manuscript 1

Do **not** expand Manuscript 1 to include:

- structural criteria other than ICC;
- unbalanced groups;
- non-centered proposal theory;
- multilevel-model generalizations beyond the balanced one-way ICC formulation;
- Fourier/spatial Gotham loops;
- KPZ or interface-growth theory;
- multiple simultaneous structural constraints;
- a large downstream simulation study of mixed-model inference;
- additional stochastic processes unless they answer a specific manuscript need.

Those can become later papers or extensions.

## Final theory target

The final new theorem for this manuscript is the neutral-boundary dependent result labeled H12:

> For a stationary irreducible two-state Markov alignment process with Gotham-neutral stationary marginal distribution, the threshold-hugging group-offset geometry converges to a finite positive random limit. The mean log scale is independent of persistence when the stationary marginal is fixed, while the variance has an explicit persistence-dependent form.

H12 is now derived and numerically verified. **Stop adding new theory to Manuscript 1** unless manuscript drafting reveals a genuine logical gap.

## Planned manuscript structure

1. **Introduction**
   - scalar structural summaries versus admissible data geometry;
   - motivating question: what does holding ICC fixed actually permit?
   - contributions.

2. **ICC frontier geometry**
   - balanced one-way setup;
   - exact update identities;
   - frontier and directional anisotropy;
   - unique-crossing condition.

3. **Sequential threshold dynamics**
   - normalized recurrence;
   - contraction, finite-limit, and expansion regimes;
   - finite-\(n\) versus asymptotic neutrality.

4. **Directional asymmetry and neutral control**
   - orthogonality is not finite-\(n\) neutrality;
   - mean-zero deterministic alignment is not neutral;
   - explicit neutral deterministic/mixed rules.

5. **Stochastic alignment**
   - IID theorem;
   - convexity and the failure of \(E[C]\);
   - same-mean/different-regime demonstration;
   - dependent neutral boundary / H12.

6. **Discussion**
   - same ICC, different geometry;
   - implications for treating ICC as a complete description of clustered structure;
   - limitations and future extensions.

Technical derivations and some numerical-verification details can move to an appendix or supplement.

## Core figures/tables

Target a compact set:

1. Static ICC frontier radius versus alignment, with representative profiles.
2. Sequential fixed-alignment trajectories showing contraction / finite limit / expansion.
3. Nonlinear drift \(h(c)\), including the distinction between raw alignment balance and structural balance.
4. Same-\(E[C]\) stochastic distributions with different \(E[h(C)]\) and different resulting regimes.
5. Neutral Markov dependence figure/table showing common mean drift but persistence-dependent limiting variability.

## Stopping rule

The first manuscript is theoretically complete when:

- H12 is derived; **complete**
- H12 is numerically verified against the exact recurrence; **complete**
- the five core figures/tables can be produced reproducibly;
- the theorem/lemma hierarchy has been consolidated from H1-H12;
- relevant literature has been reviewed for positioning and terminology.

At that point, development should shift from discovering more phenomena to writing, simplifying, and validating the manuscript.
