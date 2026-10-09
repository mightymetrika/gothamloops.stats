# Manuscript 1 Theorem Map

## Purpose

This document converts the exploratory H1-H12 notebook into a compact manuscript-facing theorem hierarchy. It is not a new research agenda. Its purpose is to decide what belongs in the main text, what becomes a corollary or lemma, and what moves to an appendix or supplement.

The manuscript story is:

> A fixed ICC threshold defines an anisotropic observation-space frontier. Repeated movement along that frontier has an exact scale-free recurrence, producing asymmetric notions of neutrality. Under stochastic alignment, the nonlinear marginal drift determines the long-run regime, while dependence determines neutral-boundary variability.

## Main-text theorem hierarchy

### Proposition 1. Exact ICC update geometry

For balanced groups and a centered proposal with radius \(r\) and alignment \(c\),

\[
W' = W + \frac{n}{n+1}r^2,
\]

\[
B' = (n+1)A + 2\sqrt{A}\,rc + \frac{r^2}{n+1},
\]

and

\[
A' = A + \frac{2\sqrt{A}\,rc}{n+1}
       + \frac{r^2}{(n+1)^2}.
\]

**Role in paper:** foundation. This establishes the geometry from which every later result follows.

**Notebook sources:** update identities, H3.

---

### Theorem 1. The fixed-ICC frontier is directionally anisotropic

Under the threshold-hugging centered regime and the H2 sufficient condition:

1. there is exactly one positive frontier radius for every alignment \(c\in[-1,1]\);
2. that radius is strictly increasing in \(c\).

Therefore the same ICC threshold permits systematically different observation-space heterogeneity depending on whether new observations oppose, are orthogonal to, or reinforce existing between-group structure.

**Role in paper:** first headline result and primary geometric motivation.

**Supporting results:**
- H1 zero-radius headroom;
- H2 unique positive crossing;
- H2 alignment monotonicity.

**Likely placement:** Section 2, with detailed algebra for H1/H2 in appendix if space is tight.

---

### Theorem 2. Threshold-hugging dynamics reduce to a scale-free recurrence

Once the state lies on the ICC threshold,

\[
z_n(c)=\frac{r_n}{\sqrt{A_n}}
\]

depends only on \(n,G,\tau,c\), not on the absolute scale of \(A_n\) or the detailed history that produced the threshold state.

The exact dynamics are

\[
A_{n+1}
=
A_n
\left[
1+\frac{2c z_n(c)}{n+1}
+\frac{z_n(c)^2}{(n+1)^2}
\right].
\]

For fixed \(c\), this yields the asymptotic regimes

\[
c<0 \Rightarrow A_n,r_n\to0,
\]

\[
c=0 \Rightarrow A_n,r_n\to\text{finite positive limits},
\]

\[
c>0 \Rightarrow A_n,r_n\to\infty.
\]

**Role in paper:** second headline result. This is the reduction that turns the frontier into an analytically tractable dynamical system.

**Notebook source:** H6.

**Corollaries:**
- H4 finite-horizon versus asymptotic structural neutrality;
- H5 finite-horizon radius versus structural neutrality.

---

### Corollary 2A. Orthogonality is not exact finite-n structural neutrality

At \(c=0\),

\[
A' = A + \frac{r^2}{(n+1)^2}>A
\]

for \(r>0\).

The exact one-step structurally neutral alignment satisfies

\[
c_n^*<0
\]

and

\[
c_n^*
=
-\frac{z_\infty(0)}{2n}+O(n^{-2}).
\]

**Role in paper:** concise illustration of nonlinear finite-sample neutrality.

**Notebook sources:** H3, H9.

---

### Theorem 3. Directional balance is not structural balance

For a variable alignment sequence, define

\[
h(c)=2c z_\infty(c).
\]

Then

\[
\log A_N
=
\log A_{n_0}
+
\sum_n \frac{h(c_n)}{n}
+
O(1).
\]

Hence neutrality is governed by harmonic-time accumulation of the nonlinear drift \(h(c)\), not by the ordinary mean of the alignments.

For symmetric magnitudes,

\[
h(c)+h(-c)>0
\]

for \(c\ne0\), so equal positive and negative directional alignment is generically expansive.

**Role in paper:** third headline result and conceptual transition from deterministic to stochastic alignment.

**Notebook source:** H7.

**Corollary 3A. Explicit neutral mixed policy**

For two available directions \(c_-<0<c_+\), asymptotic neutrality requires

\[
p_+^*
=
\frac{-h(c_-)}
{h(c_+)-h(c_-)}.
\]

For symmetric support \(\{-c,+c\}\), the neutral policy assigns more mass to the opposing direction.

**Notebook source:** H8.

**Main-text treatment:** formula and interpretation only. The greedy deterministic controller can move to supplement.

---

### Theorem 4. IID stochastic alignment is governed by nonlinear marginal drift

Let \(C_n\) be IID and

\[
\mu_h=E[h(C)].
\]

Then

\[
\log A_N-\mu_h\log N
\]

converges almost surely to a finite random limit, so

\[
\frac{A_N}{N^{\mu_h}}
\to K\in(0,\infty)
\quad\text{almost surely}.
\]

Thus

\[
\mu_h<0,\quad \mu_h=0,\quad \mu_h>0
\]

correspond respectively to contraction, a finite positive random limiting scale, and expansion.

Moreover, \(h\) is strictly convex. Therefore

\[
E[h(C)]\ge h(E[C])
\]

with strict inequality for non-degenerate \(C\). In particular,

\[
E[C]=0,\ \operatorname{Var}(C)>0
\Rightarrow
E[h(C)]>0.
\]

**Role in paper:** fourth headline result. This gives the paper's sharpest stochastic statement: ordinary average alignment can be misleading.

**Notebook source:** H10.

**Primary demonstration:** distributions with the same \(E[C]\) but different \(E[h(C)]\) and therefore different Gotham regimes.

---

### Theorem 5. Dependence changes neutral-boundary variability, not the neutral regime

For a stationary irreducible two-state Markov alignment process whose stationary marginal satisfies

\[
E[h(C)]=0,
\]

the threshold-hugging geometry satisfies

\[
A_N\to K_\rho\in(0,\infty)
\]

almost surely under the short-range Markov dependence considered in H12.

For exact log multipliers with

\[
\delta_n
=
\log m_n(c_+)-
\log m_n(c_-),
\]

the finite-horizon mean log ratio depends only on the stationary marginal and is therefore independent of persistence \(\rho\). Its variance is

\[
p^*(1-p^*)
\left[
\sum_n \delta_n^2
+
2\sum_{i<j}
\delta_i\delta_j\rho^{j-i}
\right].
\]

Thus persistence changes the distribution of the limiting structural scale even when the ICC threshold, alignment support, stationary marginal distribution, and first-order Gotham drift are all identical.

**Role in paper:** final headline result and endpoint of Manuscript 1 theory.

**Notebook sources:** H11 provides the general stationary first-order result; H12 supplies the neutral-boundary Markov result.

**Numerical verification:** theoretical and empirical finite-horizon standard deviations agree closely across \(\rho=-0.5,0,0.5,0.8\), while the spread in final \(A\) grows sharply with positive persistence.

---

## Results that should not be standalone manuscript sections

These results should remain available but should not fragment the narrative:

- H1 as a separate headline hypothesis;
- H4 and H5 as separate hypotheses;
- the 5,000-step deterministic greedy controller;
- every individual package-validation discrepancy;
- the general 0/1/2-root topology outside the manuscript's threshold-hugging regime;
- long tables of all alignment values.

They can appear as lemmas, remarks, appendix derivations, or reproducibility checks.

## Proposed main-text results sequence

### Section 2: Geometry of a fixed ICC

- Proposition 1: exact update geometry.
- Theorem 1: unique anisotropic ICC frontier.
- Figure 1: frontier radius versus alignment plus representative observation profiles.

### Section 3: Sequential frontier dynamics

- Theorem 2: scale-free recurrence.
- Fixed-alignment asymptotic regimes.
- Corollary 2A: finite-n neutral direction.
- Figure 2: representative contracting, neutral, and expanding trajectories.

### Section 4: Nonlinear neutrality

- Theorem 3: nonlinear harmonic-time drift.
- Corollary 3A: asymmetric neutral mixture.
- Figure 3: \(h(c)\) with raw zero alignment and Gotham-neutral balance distinguished.

### Section 5: Stochastic alignment

- Theorem 4: IID drift classification and strict convexity.
- Figure/Table 4: same \(E[C]\), different \(E[h(C)]\), different regimes.
- Theorem 5: neutral Markov convergence and variance.
- Figure/Table 5: persistence-dependent spread under a common neutral marginal.

## Appendix / supplement candidates

- complete derivation of zero-radius headroom;
- proof details for the unique-root sufficient condition;
- algebra for the H6 normalized recurrence;
- finite-horizon H4/H5 identities;
- exact H8 finite-\(n\) log-neutral probability;
- deterministic greedy controller;
- detailed Markov Poisson-equation argument;
- package-versus-formula numerical checks.

## Theory freeze

H12 closes the planned mathematical expansion for Manuscript 1.

The next scientific work should be:

1. verify that the five main figures/tables tell the manuscript story cleanly;
2. conduct a focused literature review for positioning and terminology;
3. draft the introduction around the "same ICC, different geometry" problem;
4. draft the mathematical methods/results around the five-result hierarchy above;
5. move algebraic detail out of the narrative wherever it obscures the story.

Do not create H13 unless manuscript drafting identifies a specific missing result that is necessary for the argument.
