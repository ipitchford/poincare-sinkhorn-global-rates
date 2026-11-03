---
title: "Curvature, relative entropy, and global Sinkhorn convergence"
author: "Anonymous"
date: "3 November 2026 (revised version)"
lang: en-GB
fontsize: 11pt
geometry: margin=25mm
colorlinks: true
linkcolor: black
urlcolor: violet
toc: false
header-includes:
  - \usepackage{amsmath,amssymb,amsthm}
  - \usepackage{microtype}
  - \usepackage{booktabs}
  - \usepackage{float}
  - \floatplacement{figure}{H}
---

## Abstract

One Poincaré inequality on the source measure is enough for Sinkhorn's algorithm to converge geometrically from its first iteration, at a rate polynomial in the regularisation parameter $\lambda$. Suppose the source has Poincaré constant $C_P$ and, at each source point, the source gradients of the cost lie in a ball of radius $L$. Then, under integrability assumptions that hold automatically for bounded costs, every full Sinkhorn iteration multiplies the dual gap by at most $C_PL^2/(C_PL^2+\lambda^2)$: one explicit factor valid for all $\lambda$. For bounded costs this holds for every $\lambda>0$, every bounded initial potential and every target. Boundedness of the cost is not necessary: for quadratic cost the bound covers every Poincaré source on $\mathbb R^d$, including exponential-tailed sources that satisfy no Talagrand inequality, and every compactly supported target, including point clouds; here $L$ is the target's Chebyshev radius. The same curvature-to-entropy argument turns the curvature estimate of Chizat, Delalande and Vaškevičius into a global rate $C/(C+\lambda)$, linear in $\lambda$, in their bounded-density and log-concave settings. Within the Poincaré class with bounded source gradients, the exponent two cannot be lowered, even for one fixed source. In a two-atom model, order $\lambda^2$ is attained uniformly over a Poincaré class but by no single fixed absolutely continuous source. The closest prior result, by Chiarini, Conforti, Greco and Tamanini, gives the same rational rate under a Talagrand inequality and Hessian bounds on the cost. Relative to it, our gains are the weaker functional inequality, first-derivative hypotheses, validity for all $\lambda$, and the sharpness results.

## 1. The paper in brief

Sinkhorn's algorithm, also called iterative proportional fitting or matrix scaling, looks for a joint distribution (a *coupling*) with two prescribed marginal distributions. It alternately rescales the coupling so that its first marginal is correct, then so that its second marginal is correct, and repeats. It is the standard solver for entropy-regularised optimal transport. There, a parameter $\lambda>0$ blurs the transport plan. Small $\lambda$ approximates unregularised optimal transport more closely but slows the algorithm. The practical question is therefore how the number of iterations grows as $\lambda$ shrinks. Classical guarantees grow exponentially in $1/\lambda$. Recent work obtained polynomial growth, but only under strong assumptions on one of the marginals.

We show that one comparatively mild property of the source distribution is enough: a Poincaré inequality, which says that a function with small gradient has small variance. Under it, every full Sinkhorn step shrinks the error by a fixed factor, starting from the first step and for any target, provided the source gradients of the cost lie in balls of a fixed radius and mild integrability conditions hold (automatic for bounded costs). The resulting iteration count grows like $1/\lambda^2$ (Theorem 2). This includes quadratic cost with unbounded sources, such as Gaussian or exponential-tailed laws, whenever the target is bounded, for example a point cloud (Theorem 3 and Corollary 4). When, in addition, the source lives on a compact convex domain with a density whose ratio of maximum to minimum is bounded, and the cost is bounded and semiconcave in the source variable, the count improves to order $1/\lambda$ (Corollary 5). Explicit one-dimensional examples show that $1/\lambda^2$ cannot be improved in general (Section 7). All results concern the exact algorithm on continuous distributions; a source known only through samples needs a different argument.

## 2. Problem, results, and relation to prior work

### 2.1. Setting

Let $\mu$ and $\nu$ be probability measures on a source space $X$ and a target space $Y$. For $\lambda>0$, entropy-regularised optimal transport minimises

$$
\int c(x,y)\,d\gamma(x,y)+\lambda\operatorname{KL}(\gamma\Vert\mu\otimes\nu)
\quad\text{over }\gamma\in\Pi(\mu,\nu).
\tag{1}
$$

Our $\lambda$ is the $\lambda$ of Chizat, Delalande and Vaškevičius (2026; hereafter CDV) and the $\varepsilon$ of Chiarini, Conforti, Greco and Tamanini (2025; hereafter CCGT). All logarithms are natural. One *iteration* means one source normalisation followed by one target normalisation. We measure progress by the semidual objective gap $\Delta$, defined in (3). When an optimal potential exists, this gap is $\lambda$ times the relative entropy between the optimal coupling and the current source-normalised coupling (equation (34)). We do not study the bias between (1) and unregularised transport, or the arithmetic cost of evaluating continuous integrals.

### 2.2. Main results at a glance

| Result | Assumptions (source $\mu$; cost $c$; target $\nu$) | Per-iteration factor on the dual gap |
|:------|:-----------|:------|
| Theorem 1 | curvature bound (8), with variance in the *current* target marginal | $1-\lambda/A$ |
| Theorem 2 | Poincaré constant $C_P$; source gradients of $c$ in balls of radius $L$; integrability (S1)--(S3) | $C_PL^2/(C_PL^2+\lambda^2)$ |
| Theorem 3, Corollary 4 | Poincaré source on $\mathbb R^d$; cost Lipschitz in $y$ with linear growth in $|x|$; compact target. Corollary 4: quadratic or bilinear cost, target Chebyshev radius $R$ | $C_PL^2/(C_PL^2+\lambda^2)$, every $\lambda$; $C_PR^2/(C_PR^2+\lambda^2)$ in Corollary 4 |
| Corollary 5 | CDV's bounded-density or log-concave settings | $C/(C+\lambda)$ |
| Theorem 6, Propositions 7--8, Corollary 9 | one-dimensional two-atom models | the exponent $2$ is sharp, also for Talagrand sources |

All upper contraction bounds hold from iteration zero and for every bounded initial potential; the sharpness results are local and asymptotic statements.

### 2.3. Relation to prior work

*Transport inequalities (CCGT).* The closest prior result is CCGT's. Their Propositions 1.8(i) and 1.9(i) assume that one marginal satisfies Talagrand's inequality $\mathrm{TI}(\tau)$: $W_2^2(\cdot,\nu)\le2\tau\operatorname{KL}(\cdot\Vert\nu)$. They also assume second-derivative bounds on the cost. For small regularisation, Proposition 1.8(i) (two-sided Hessian bounds, cost Lipschitz in the Talagrand variable) gives the factor $1-\varepsilon^2/(\varepsilon^2+\tau\varepsilon(H-h)+\tau\operatorname{Lip}^2)$. Proposition 1.9(i) (other marginal compactly supported, Hessian in the Talagrand variable constant) gives $1-\varepsilon^2/(\varepsilon^2+\tau\|g\|_\infty^2)$. In our orientation this is $\tau L^2/(\tau L^2+\lambda^2)$, which is Theorem 2's factor with $\tau$ in place of $C_P$. Section 9.1 quotes their statements and compares hypotheses in detail. Because $\mathrm{TI}(\tau)$ implies a Poincaré inequality with $C_P\le\tau$ (Gozlan 2010, Proposition 4.9, applied with $\min(|u|,u^2)\le u^2$), our contribution here is a weakening of hypotheses at an equal or better constant, not a new rate. The weakening is threefold: Poincaré in place of Talagrand, first source derivatives of the cost in place of Hessian bounds, and validity for every $\lambda$ and every bounded initial potential. Theorem 3 extends this to unbounded costs with compactly supported targets. Talagrand sources, such as Gaussians, with compact targets are already covered by CCGT's Proposition 1.9(i). The genuinely new cases are Poincaré sources without a Talagrand inequality, such as the two-sided exponential law.

*Curvature (CDV).* CDV's Proposition 5.1 bounds the curvature of the semidual objective below by a variance in the current target marginal. CDV's global rate passes through a bound on $\operatorname{KL}(\nu\Vert\nu[\psi_t+sv])$ in which, by their Remark 5.2, a factor $\lambda/c_\infty$ is lost. The remark adds that improving that bound to $\lambda^{-1}(\delta_t-\delta_{t+1})$ would give the rate $(1-\Theta(\lambda/c_\infty))^t$. Theorem 1 supplies that estimate, and Corollary 5 realises the anticipated rate globally. CDV's Remark 3.1 conjectures that a Poincaré inequality on the source "should be a sufficient condition for our results to hold". Theorem 2 confirms the conjecture for the global $\lambda^2$ rate, for costs satisfying (7). Theorem 6 shows that the linear rate does not extend to the whole Poincaré class (Section 7.6).

*Mirror descent (AFKL).* Aubin-Frankowski, Korba and Léger (2022; hereafter AFKL) write Sinkhorn's algorithm as mirror descent. In their framework, a reverse entropy inequality such as (10) is relative strong convexity. Their Theorem 4 then converts it into a linear rate. Their Proposition 7 supplies such an inequality with coefficient $1+4e^{3D_c/\varepsilon}$, which is exponential in $1/\varepsilon$. Our new step is (8) $\Rightarrow$ (9)--(10), which yields a polynomial coefficient. Once (10) holds, the contraction (11) also follows from AFKL's mirror-descent analysis (Remark 4.2). Section 4.3 gives a self-contained dual proof.

*Probabilistic reading.* Hypothesis (8) says that every source-normalised coupling has maximal correlation at most $\sqrt{1-\lambda/A}$ (Remark 4.1). Inequalities (25)--(27) are the classical decomposition of variance into a mean conditional variance and the variance of a conditional mean, which is used for functional inequalities of Gibbs measures and mixtures (Otto and Reznikoff 2007; Grunewald, Otto, Villani and Westdickenberg 2009; Chafaï and Malrieu 2010).

*Priority.* CDV report an announcement by S. Di Marino at a 2022 Pisa workshop of a rate $1-\Theta(\lambda^2/c_\infty^2)$ for bounded Lipschitz costs when both marginals satisfy Poincaré inequalities. A 2021 Simons Institute seminar abstract by Di Marino also mentions counterexamples. We could not inspect those proofs. The exponent $\lambda^2$ under Poincaré assumptions was therefore announced earlier. The gain claimed here is that one marginal suffices, together with the explicit constant, the generality in the target, and the unbounded costs. Priority for the exponent and for counterexamples remains uncertain.

### 2.4. Roadmap

Section 3 sets out definitions and integrability assumptions that allow unbounded costs. Section 4 proves the curvature-to-entropy transfer (Theorem 1) by finite-dimensional convex duality. Section 5 derives the Poincaré rate (Theorem 2), extends it to compactly supported targets with unbounded costs (Theorem 3 and Corollary 4), and records practical consequences: an iteration count, a stopping certificate and an inexact-iteration bound. Section 6 applies Theorem 1 to CDV's curvature estimate. Section 7 proves the sharpness results. Section 8 summarises the numerical checks, and Section 9 compares the results in detail with the literature and lists open questions.

## 3. Definitions and standing assumptions

Let $X$ and $Y$ be standard Borel spaces with probability measures $\mu$ and $\nu$, and let $c:X\times Y\to\mathbb R$ be jointly measurable. Fix $\lambda>0$. We use the following standing assumptions.

* **(S1)** $c\in L^1(\mu\otimes\nu)$.
* **(S2)** There is a measurable $\omega\in L^1(\mu)$ with $|c(x,y)-c(x,y')|\le\omega(x)$ for all $x\in X$ and $y,y'\in Y$.
* **(S3$_\lambda$)** There is a measurable $G_\lambda\in L^1(\mu)$ such that, for all $x\in X$ and $y\in Y$,
  $e^{-c(x,y)/\lambda}\le G_\lambda(x)\int_Ye^{-c(x,y')/\lambda}\,d\nu(y')$.

A bounded cost satisfies all three for every $\lambda$, with $\omega=\operatorname{osc}c$ and $G_\lambda=e^{\operatorname{osc}c/\lambda}$. More generally, (S2) implies (S3$_\lambda$) whenever $e^{\omega/\lambda}\in L^1(\mu)$. When $\nu$ is finitely supported, (S3$_\lambda$) holds with $G_\lambda=1/\min_j\nu(\{y_j\})$. Theorem 3 verifies (S1)--(S3$_\lambda$) for every $\lambda$ in the main unbounded examples. Adding to $c$ a term $a(x)+b(y)$ with $a\in L^1(\mu)$ and $b$ bounded changes none of the quantities below, apart from additive constants and a bounded shift of potentials. A source-only term $a\notin L^1(\mu)$ does not change the oscillation in (S2), but it makes (1) infinite, so (S1) cannot be dropped. In that case one can still remove $a$ and study the renormalised problem.

For a bounded measurable target potential $\psi:Y\to\mathbb R$, set

$$
\begin{aligned}
Z_\psi(x)&=\int_Y e^{(\psi(y)-c(x,y))/\lambda}\,d\nu(y),\\
\Phi(\psi)&=\lambda\int_X\log Z_\psi(x)\,d\mu(x),\\
q_x^\psi(dy)&=Z_\psi(x)^{-1}e^{(\psi(y)-c(x,y))/\lambda}\,d\nu(y),\\
r_\psi(dy)&=\int_X q_x^\psi(dy)\,d\mu(x),\qquad
\gamma_\psi(dx,dy)=\mu(dx)q_x^\psi(dy).
\end{aligned}
\tag{2}
$$

By Jensen's inequality and (S2), $\int\psi\,d\nu-\int c(x,y)\,d\nu(y)\le\lambda\log Z_\psi(x)\le\sup\psi-\int c(x,y)\,d\nu(y)+\omega(x)$. By (S1)--(S2), $\Phi(\psi)$ is therefore finite. The coupling $\gamma_\psi$ always has source marginal $\mu$; its target marginal $r_\psi$ need not equal $\nu$. For every $x$ and $\nu$-almost every $y$,

$$
e^{-(\operatorname{osc}\psi+\omega(x))/\lambda}\le
\frac{dq_x^\psi}{d\nu}(y)\le e^{\operatorname{osc}\psi/\lambda}G_\lambda(x).
$$

Integrating in $x$, and applying Jensen's inequality to the left side, gives

$$
e^{-(\operatorname{osc}\psi+\int\omega\,d\mu)/\lambda}\le\frac{dr_\psi}{d\nu}
\le e^{\operatorname{osc}\psi/\lambda}\|G_\lambda\|_{L^1(\mu)}.
$$

Thus $r_\psi$ and $\nu$ are mutually absolutely continuous with log-density bounded above and below. Define the semidual objective and its gap by

$$
E(\psi)=\int_Y\psi\,d\nu-\Phi(\psi),\qquad
E^*=\sup_{\eta\in L^\infty(\nu)}E(\eta),\qquad
\Delta(\psi)=E^*-E(\psi).
\tag{3}
$$

The supremum is finite since, by Jensen's inequality, $E(\psi)\le\int c\,d(\mu\otimes\nu)$ for every $\psi$. CDV's entropic Kantorovich functional is $K=-\Phi$, and their semidual objective $\langle\psi^{c,\lambda},\mu\rangle+\langle\psi,\nu\rangle$ coincides with $E$. The exact full Sinkhorn update is

$$
S\psi=\psi-\lambda\log\frac{dr_\psi}{d\nu}.
\tag{4}
$$

Here the source potential is $-\lambda\log Z_\psi$, and transforming it back to the target gives (4). By the bounds above, $S$ maps bounded potentials to bounded potentials; this is the only use of (S3$_\lambda$). Adding a constant to a potential changes none of $q_x^\psi$, $r_\psi$, $E$ and $\Delta$. Differentiation in bounded directions gives

$$
D\Phi(\psi)[h]=\int h\,dr_\psi,\qquad
D^2\Phi(\psi)[h,h]=\frac1\lambda\int_X
\operatorname{Var}_{q_x^\psi}(h)\,d\mu(x).
\tag{5}
$$

The $s$-derivatives of $\lambda\log Z_{\psi+sh}(x)$ are a mean and a variance under $q_x^{\psi+sh}$. They are bounded by $\|h\|_\infty$ and $\|h\|_\infty^2/\lambda$, uniformly in $x$, so differentiation under $\int d\mu$ needs no boundedness of $c$.

For the Poincaré results, $X$ is a connected open subset of $\mathbb R^d$ and

$$
\operatorname{Var}_\mu(f)\le C_P\int_X|\nabla f|^2\,d\mu
\tag{6}
$$

for bounded $C^1$ functions $f$ with bounded gradient. We assume that each $c(\cdot,y)$ is $C^1$ on $X$, that $\nabla_xc$ is jointly measurable, and that, for every $x\in X$, the set of source gradients lies in a closed ball of radius $L$:

$$
\inf_{b\in\mathbb R^d}\sup_{y\in Y}|\nabla_xc(x,y)-b|\le L.
\tag{7}
$$

If $|\nabla_xc|\le L$, then (7) holds with $b=0$. For $c(x,y)=|x-y|^2/2$ or $c(x,y)=-\langle x,y\rangle$ with $Y\subset\mathbb R^d$ bounded, (7) holds with $L$ the Chebyshev radius of $Y$ (the radius of its smallest enclosing ball). CCGT use the same device: their $g(x)=\nabla_2c(x,y)-\Sigma y$. No Lebesgue density is assumed for $\nu$, and the source domain may be nonconvex or unbounded.

## 4. From uniform curvature to a global entropy bound

**Theorem 1 (curvature transfer).** Fix $\lambda>0$ and assume (S1)--(S3$_\lambda$). Suppose there is a constant $A\ge\lambda$ such that, for *every* bounded $\psi$ and $h$,

$$
D^2\Phi(\psi)[h,h]\ge A^{-1}\operatorname{Var}_{r_\psi}(h).
\tag{8}
$$

Then, for every bounded $\psi$,

$$
\boxed{\Delta(\psi)\le A\operatorname{KL}(\nu\Vert r_\psi).}
\tag{9}
$$

For any two bounded potentials $\eta$ and $\psi$,

$$
\boxed{\operatorname{KL}(\gamma_\eta\Vert\gamma_\psi)
\le \frac A\lambda\operatorname{KL}(r_\eta\Vert r_\psi).}
\tag{10}
$$

Consequently, exact full Sinkhorn iterations satisfy

$$
\boxed{\Delta(S\psi)\le(1-\lambda/A)\Delta(\psi).}
\tag{11}
$$

The variance in (8) is taken in the target marginal $r_\psi$, which moves with the potential. The proof uses (8) at potentials $\psi+u$ near $\psi$, not only at $\psi$ (see (13)). A curvature bound known only along the Sinkhorn trajectory would therefore not suffice. Replacing $r_\psi$ by $\nu$ in (8) would also be unjustified.

### 4.1. Finite partitions and a genuine inverse

Fix a bounded $\psi$. Partition $Y$ into finitely many measurable cells $B_1,\ldots,B_m$ of positive $\nu$-measure, merging any null cells into a neighbour. Put

$$
Q_i(x)=q_x^\psi(B_i),\qquad
F(u)=\lambda\int_X\log\left(\sum_{i=1}^m Q_i(x)e^{u_i/\lambda}\right)d\mu(x),\qquad u\in\mathbb R^m.
\tag{12}
$$

The function $F$ is the increment $\Phi(\psi+h_u)-\Phi(\psi)$ for the cell-constant potential $h_u=\sum_iu_i\mathbf 1_{B_i}$. It satisfies $F(0)=0$, $F(u+a\mathbf1)=F(u)+a$ and $|F(u)|\le\max_i|u_i|$. Its $u$-derivatives are softmax moments, bounded uniformly in $x$, so $F$ is smooth. Write $r(u)=\nabla F(u)$, so that $r(u)_i=r_{\psi+h_u}(B_i)>0$. By (5), $z^\top\nabla^2F(u)z=D^2\Phi(\psi+h_u)[h_z,h_z]$. Applying (8) at the potential $\psi+h_u$ gives

$$
\nabla^2F(u)\succeq A^{-1}
\left(\operatorname{diag}r(u)-r(u)r(u)^\top\right).
\tag{13}
$$

Since $r(u)$ is a strictly positive probability vector, the right side of (13) has kernel exactly $\operatorname{span}(\mathbf1)$. The Hessian $\nabla^2F(u)$ annihilates $\mathbf1$ by translation invariance. Hence $\ker\nabla^2F(u)=\operatorname{span}(\mathbf1)$, and $F$ is strictly convex on the quotient of $\mathbb R^m$ by constants.

Let $p$ be a strictly positive probability vector, and use the gauge $\max_iu_i=0$. Choose $i^*$ with $u_{i^*}=0$. The lower bound for $dq_x^\psi/d\nu$ in Section 3 gives $Q_{i^*}(x)\ge\nu(B_{i^*})e^{-(\operatorname{osc}\psi+\omega(x))/\lambda}$, hence

$$
F(u)-p\cdot u\ge\lambda\min_i\log\nu(B_i)-\operatorname{osc}\psi-\int\omega\,d\mu+
(\min_i p_i)\operatorname{osc}u.
\tag{14}
$$

Thus $F(u)-p\cdot u$ is coercive modulo constants, and it has a unique minimiser modulo constants. Its first-order condition is $r(u)=p$. The gradient therefore maps the quotient bijectively onto the open probability simplex. Because the Hessian is nonsingular on the quotient, the implicit function theorem gives a smooth inverse.

Let $G=F^*$ on the open simplex. For a tangent vector $v$ with $\sum_iv_i=0$, we claim

$$
D^2G(p)[v,v]\le A\sum_{i=1}^m\frac{v_i^2}{p_i}.
\tag{15}
$$

Let $H=\nabla^2F(u)$ at $r(u)=p$. On tangent vectors, $D^2G(p)$ is the inverse quadratic form of $H$ on the quotient:

$$
\begin{aligned}
v^\top H^{-1}v
&=\sup_z\{2v\cdot z-z^\top Hz\}\\
&\le\sup_z\{2v\cdot z-A^{-1}\operatorname{Var}_p(z)\}
=A\sum_i v_i^2/p_i.
\end{aligned}
\tag{16}
$$

Here $H^{-1}$ denotes only this inverse form on tangent vectors, not an inverse of the singular matrix. The final equality follows by centring $z$ under $p$ and completing the square; the maximiser is $z_i=Av_i/p_i$. This proves (15).

### 4.2. Entropy integration and the infinite target

Write $r=r(0)$, so that $r_i=r_\psi(B_i)$. Then $G(r)=0$, and $DG(r)$ is represented by $u=0$ modulo constants. Since $\sum_i(p_i-r_i)=0$, the first-order term vanishes. The segment $r+t(p-r)$, $t\in[0,1]$, stays in the open simplex because $r$ and $p$ are strictly positive. Integrating (15) along it gives

$$
\begin{aligned}
G(p)
&\le A\int_0^1(1-t)\sum_i
\frac{(p_i-r_i)^2}{r_i+t(p_i-r_i)}\,dt
=A\sum_i p_i\log(p_i/r_i).
\end{aligned}
\tag{17}
$$

The equality is the integral form of Taylor's theorem for $p\mapsto\sum_ip_i\log p_i$. Hence

$$
\sup_u\{p\cdot u-F(u)\}\le A\operatorname{KL}(p\Vert r).
\tag{18}
$$

Take $p_i=\nu(B_i)$. Then $p\cdot u-F(u)=E(\psi+h_u)-E(\psi)$ is the semidual improvement due to the cell-constant increment $h_u$. Relative entropy cannot increase when each point is mapped to its cell, so $\operatorname{KL}(p\Vert r)\le\operatorname{KL}(\nu\Vert r_\psi)$. This entropy is finite by the density bounds in Section 3.

Every bounded measurable increment is a uniform limit of simple functions, and the level sets of a simple function form a finite partition. If two potentials differ by at most $s$ in the supremum norm, their values of $\lambda\log Z$ differ by at most $s$ at every $x$. Hence $\Phi$ is 1-Lipschitz and $E$ is 2-Lipschitz in the supremum norm. Passing to uniform limits and then taking the supremum over all increments proves (9).

For (10), take $p_i=r_\eta(B_i)$ and the increment $h=\eta-\psi$. The quantity being bounded is

$$
\int h\,dr_\eta-\Phi(\psi+h)+\Phi(\psi)
=\lambda\operatorname{KL}(\gamma_\eta\Vert\gamma_\psi).
\tag{19}
$$

To see this, integrate $\log(d\gamma_\eta/d\gamma_\psi)=h/\lambda-\log Z_\eta+\log Z_\psi$ against $\gamma_\eta$; both logarithms are $\mu$-integrable by Section 3. Uniform simple approximation of $h$ and data processing finish the proof of (10). Every inverse in this argument is finite-dimensional. In particular, nothing is assumed about the range of an infinite-dimensional gradient map.

### 4.3. The Sinkhorn ascent identity

Let $\rho_\psi=dr_\psi/d\nu$ and $\sigma_\psi(x)=\int\rho_\psi^{-1}\,dq_x^\psi$. The function $\sigma_\psi$ is the density, with respect to $\mu$, of the source marginal of the coupling obtained from $\gamma_\psi$ by one target normalisation; note that $\int\sigma_\psi\,d\mu=\int\rho_\psi^{-1}\,dr_\psi=1$. From (4), $Z_{S\psi}(x)=Z_\psi(x)\sigma_\psi(x)$ and $\int(S\psi-\psi)\,d\nu=\lambda\operatorname{KL}(\nu\Vert r_\psi)$. Therefore

$$
E(S\psi)-E(\psi)=\lambda\operatorname{KL}(\nu\Vert r_\psi)+\lambda\operatorname{KL}(\mu\Vert\sigma_\psi\mu)
\ge\lambda\operatorname{KL}(\nu\Vert r_\psi).
\tag{20}
$$

The two terms are the gains of the target and source half-steps. Combining (9) with (20) gives $\Delta(S\psi)\le\Delta(\psi)-\lambda\operatorname{KL}(\nu\Vert r_\psi)\le(1-\lambda/A)\Delta(\psi)$, which proves (11). No bound on the oscillation of the initial potential is used, and the argument applies from iteration zero. $\square$

### 4.4. Two readings of Theorem 1

**Remark 4.1 (maximal correlation).** Under $\gamma_\psi$, let $(X,Y)$ have source marginal $\mu$ and conditional law $q_X^\psi$ of $Y$ given $X$. By (5) and the law of total variance, $\lambda D^2\Phi(\psi)[h,h]=\operatorname{Var}h(Y)-\operatorname{Var}\mathbb E[h(Y)\mid X]$. Hence (8) holds for all bounded $h$ if and only if

$$
\varrho(\gamma_\psi)^2:=\sup_h\frac{\operatorname{Var}\mathbb E[h(Y)\mid X]}{\operatorname{Var}h(Y)}\le1-\frac\lambda A,
\tag{21}
$$

where $\varrho$ is the Hirschfeld--Gebelein--Rényi maximal correlation (Rényi 1959). Theorem 1 therefore says that if every source-normalised coupling has maximal correlation at most $\bar\varrho$, then every full iteration multiplies the dual gap by at most $\bar\varrho^2$. Inequality (27) below is the bound $\varrho(\gamma_\psi)^2\le B/(B+\lambda^2)$.

**Remark 4.2 (the mirror-descent route).** In AFKL's notation, take the constraint set to consist of couplings with source marginal $\mu$, take the objective to be the marginal entropy $\operatorname{KL}(r_\gamma\Vert\nu)$, and take relative entropy as Bregman divergence. Data processing gives relative smoothness with constant 1 (their Lemma 6). Inequality (10) gives relative strong convexity with constant $\lambda/A$, but only on the family $\mathcal R=\{\gamma_\psi:\psi\in L^\infty(\nu)\}$, not on all couplings with source marginal $\mu$; this mirrors AFKL's own restriction to $\Pi_c\cap C$. The family $\mathcal R$ is invariant under full Sinkhorn steps, since $S$ maps bounded potentials to bounded potentials. Suppose an optimal bounded potential $\psi^*$ exists, so that $\gamma^*=\gamma_{\psi^*}\in\mathcal R$. AFKL's Theorem 4, applied on $\mathcal R$ with reference entropy $\operatorname{KL}(\cdot\Vert\gamma_0)$, then bounds the objective after $n$ steps by $l\operatorname{KL}(\gamma^*\Vert\gamma_0)/((1+l/(L-l))^n-1)$, with $L=1$ and $l=\lambda/A$; this decays at rate $1-\lambda/A$. Under the same condition, the three-point inequality underlying that theorem gives $\operatorname{KL}(\gamma^*\Vert\gamma_{t+1})\le(1-\lambda/A)\operatorname{KL}(\gamma^*\Vert\gamma_t)$, which by (34) is (11). The proof in Section 4.3 does not presuppose an optimal potential.

## 5. One Poincaré inequality on the source

**Theorem 2 (explicit Poincaré rate).** Assume (6), (7) and the differentiability hypotheses of Section 3, and assume (S1)--(S3$_\lambda$) for a given $\lambda>0$. Put $B=C_PL^2$. For every bounded $\psi$,

$$
\Delta(\psi)\le\left(\lambda+\frac B\lambda\right)
\operatorname{KL}(\nu\Vert r_\psi).
\tag{22}
$$

For every bounded initial potential and every integer $t\ge0$,

$$
\boxed{\Delta(\psi_t)\le
\left(\frac B{B+\lambda^2}\right)^t\Delta(\psi_0),\qquad
\psi_{t+1}=S\psi_t.}
\tag{23}
$$

When $B=0$, one full iteration reaches zero gap. For bounded $\eta$ and $\psi$,

$$
\operatorname{KL}(\gamma_\eta\Vert\gamma_\psi)
\le\left(1+\frac B{\lambda^2}\right)
\operatorname{KL}(r_\eta\Vert r_\psi).
\tag{24}
$$

In particular, for a bounded cost satisfying (6)--(7), all three conclusions hold for every $\lambda>0$.

*Proof.* Fix a bounded $h$ and let $m(x)=\mathbb E_{q_x^\psi}h$. We first show that $m\in C^1(X)$ with

$$
\nabla m(x)=-\lambda^{-1}
\operatorname{Cov}_{q_x^\psi}(h,\nabla_xc(x,\cdot)).
\tag{25}
$$

Fix $x_0\in X$, a ball $U\subset X$ of radius $s$ around $x_0$, and a point $y_0\in Y$. Replacing $c(x,y)$ by $c'(x,y)=c(x,y)-c(x,y_0)$ leaves $q_x^\psi$ unchanged. By (7), $|\nabla_xc'(x,y)|\le2L$. By the mean value theorem and (S2), $|c'(x,y)|\le\omega(x_0)+2Ls$ on $U\times Y$. The integrands defining $m$, and their $x$-gradients, are therefore bounded on $U\times Y$. Dominated convergence justifies differentiation under the integral and continuity of the derivative, which gives (25). The term $\nabla_xc(x,y_0)$ cancels from the covariance.

Let $b(x)$ be a centre in (7). For a scalar-vector covariance, Cauchy--Schwarz gives

$$
\begin{aligned}
|\nabla m(x)|^2
&\le\lambda^{-2}\operatorname{Var}_{q_x^\psi}(h)\,
\mathbb E_{q_x^\psi}|\nabla_xc-\mathbb E_{q_x^\psi}\nabla_xc|^2\\
&\le\lambda^{-2}\operatorname{Var}_{q_x^\psi}(h)\,
\mathbb E_{q_x^\psi}|\nabla_xc-b(x)|^2
\le\lambda^{-2}L^2\operatorname{Var}_{q_x^\psi}(h).
\end{aligned}
\tag{26}
$$

Thus $m$ is bounded by $\|h\|_\infty$, is $C^1$, and has gradient bounded by $L\|h\|_\infty/\lambda$, so (6) applies to it. Let $D=\int\operatorname{Var}_{q_x^\psi}(h)\,d\mu$. The law of total variance and (6) give

$$
\operatorname{Var}_{r_\psi}(h)
=D+\operatorname{Var}_\mu(m)
\le\left(1+\frac B{\lambda^2}\right)D.
\tag{27}
$$

By (5), condition (8) holds with $A=\lambda+B/\lambda$, uniformly in $\psi$. Theorem 1 then gives (22)--(24). $\square$

### 5.1. Unbounded costs with compactly supported targets

The next two lemmas remove the bounded-cost hypothesis when the target is compact. The first is a form of the classical exponential concentration under a Poincaré inequality (Gromov and Milman 1983; Bobkov and Ledoux 1997). We include a proof because our test-function class in (6) is restricted.

**Lemma 1 (exponential integrability).** If $\mu$ satisfies (6) on a connected open $X\subset\mathbb R^d$, then $\int e^{t|x|}\,d\mu<\infty$ for every $t\ge0$ with $C_Pt^2<4$.

*Proof.* Let $f$ be bounded and $C^1$ with $|\nabla f|\le1$, and set $H(t)=\int e^{tf}d\mu$. Applying (6) to $e^{tf/2}$ gives $H(t)-H(t/2)^2\le C_P(t^2/4)H(t)$. Iterating $k$ times and letting $k\to\infty$, using $H(t/2^k)^{2^k}\to e^{t\int f\,d\mu}$, yields

$$
\int e^{t(f-\int f\,d\mu)}\,d\mu\le K_P(t):=\prod_{j\ge0}\left(1-\frac{C_Pt^2}{4^{j+1}}\right)^{-2^j}<\infty.
\tag{28}
$$

The same bound holds for $-f$. Let $\varrho_0(x)=(1+|x|^2)^{1/2}$, which is smooth and 1-Lipschitz. Let $\theta(s)=s$ for $s\le0$ and $\theta(s)=1-e^{-s}$ for $s>0$, and put $f_n=n+\theta(\varrho_0-n)$. Each $f_n$ is bounded and $C^1$, with $|\nabla f_n|\le1$, and $f_n\uparrow\varrho_0$ pointwise. Choose $R_0$ with $\mu(\varrho_0\le R_0)\ge3/4$, and choose $r_0$ with $K_P(t)e^{-tr_0}<3/4$. If $\int f_n\,d\mu>R_0+r_0$, then $\{f_n\le R_0\}\subset\{\int f_n-f_n\ge r_0\}$, whose measure is less than $3/4$ by Chebyshev's inequality and (28) for $-f_n$. This is a contradiction, so $\int f_n\,d\mu\le R_0+r_0$ for every $n$. Hence $\int e^{tf_n}d\mu\le K_P(t)e^{t(R_0+r_0)}$, and monotone convergence gives $\int e^{t\varrho_0}d\mu<\infty$. Finally, $|x|\le\varrho_0(x)$. $\square$

**Lemma 2 (cap bound).** Let $(Y,d_Y)$ be a compact metric space and $\nu$ a Borel probability measure with $\operatorname{supp}\nu=Y$. Then $m_\delta:=\inf_{y\in Y}\nu(B(y,\delta))>0$ for every $\delta>0$. If $|c(x,y)-c(x,y')|\le\ell(x)d_Y(y,y')$ for all $x,y,y'$, then for all $x$, $y$ and $\delta>0$,

$$
e^{-c(x,y)/\lambda}\le m_\delta^{-1}e^{\delta\ell(x)/\lambda}\int_Ye^{-c(x,y')/\lambda}\,d\nu(y').
\tag{29}
$$

*Proof.* Cover $Y$ by finitely many balls $B(y_j,\delta/2)$. Each has positive mass because $y_j\in\operatorname{supp}\nu$. Every $B(y,\delta)$ contains one of them, so $m_\delta\ge\min_j\nu(B(y_j,\delta/2))>0$. On $B(y,\delta)$ we have $c(x,y')\le c(x,y)+\delta\ell(x)$. Restricting the integral in (29) to this ball proves (29). $\square$

**Theorem 3 (Poincaré sources with compactly supported targets).** Let $X\subset\mathbb R^d$ be connected and open, let $\mu$ satisfy (6), and let $c$ satisfy the differentiability hypotheses of Section 3 and (7). Let $(Y,d_Y)$ be a compact metric space and $\nu$ a Borel probability measure with $\operatorname{supp}\nu=Y$. Assume that there are constants $\ell_0,\ell_1\ge0$ with

$$
|c(x,y)-c(x,y')|\le(\ell_0+\ell_1|x|)\,d_Y(y,y')\qquad(x\in X,\ y,y'\in Y),
\tag{30}
$$

and that $c(\cdot,y_0)\in L^1(\mu)$ for some $y_0\in Y$. Then (S1)--(S3$_\lambda$) hold for every $\lambda>0$. Consequently (22)--(24) hold for every $\lambda>0$ and every bounded initial potential, with $B=C_PL^2$.

*Proof.* The function $c$ is $C^1$ in $x$ and Lipschitz in $y$, hence jointly measurable. By Lemma 1, $\int|x|\,d\mu<\infty$. Since $|c(x,y)|\le|c(x,y_0)|+(\ell_0+\ell_1|x|)\operatorname{diam}Y$, (S1) holds, and (S2) holds with $\omega(x)=(\ell_0+\ell_1|x|)\operatorname{diam}Y$. Fix $t>0$ with $C_Pt^2<4$. Lemma 2 with $\delta=\lambda t/\ell_1$ (any $\delta>0$ if $\ell_1=0$) gives (S3$_\lambda$) with $G_\lambda(x)=m_\delta^{-1}e^{\delta\ell_0/\lambda}e^{t|x|}$. This function is $\mu$-integrable by Lemma 1. Theorem 2 applies. $\square$

**Corollary 4 (quadratic cost and semi-discrete transport).** Let $\mu$ satisfy (6) on a connected open $X\subset\mathbb R^d$. Let $\nu$ be any probability measure on $\mathbb R^d$ with compact support $K$, and let $R$ be the Chebyshev radius of $K$. For $c(x,y)=|x-y|^2/2$ or $c(x,y)=-\langle x,y\rangle$, every $\lambda>0$ and every bounded initial potential on $K$,

$$
\boxed{\Delta(\psi_t)\le\left(\frac{C_PR^2}{C_PR^2+\lambda^2}\right)^t\Delta(\psi_0).}
\tag{31}
$$

*Proof.* Take $Y=K$ with the Euclidean distance. The source gradients $x-y$ (or $-y$), $y\in K$, lie in a ball of radius $R$. Condition (30) holds with $\ell_0=\max_K|y|$ and $\ell_1=1$, because $\tfrac12\bigl||x-y|^2-|x-y'|^2\bigr|\le(|x|+\max_K|y|)|y-y'|$. Finally, $|x-y_0|^2\in L^1(\mu)$ by Lemma 1. $\square$

*Examples.* Two examples show the range of Corollary 4.

* A Gaussian source $N(m,\Sigma)$ has $C_P=\|\Sigma\|_{\mathrm{op}}$ and gives the factor $\|\Sigma\|R^2/(\|\Sigma\|R^2+\lambda^2)$. Since Gaussian measures also satisfy Talagrand's inequality $\mathrm{TI}(\|\Sigma\|)$, this factor is CCGT's Proposition 1.9(i). Corollary 4 adds only validity for every $\lambda$ and every bounded initial potential.
* The two-sided exponential law $\mu(dx)=\frac12e^{-|x|}dx$ on $\mathbb R$ has $C_P\le4$. For bounded $C^1$ $f$ with bounded derivative, let $g=f-f(0)$. Integration by parts on $(0,\infty)$ gives $\int_0^\infty g^2e^{-x}dx=2\int_0^\infty gg'e^{-x}dx$. Cauchy--Schwarz then gives $\int_0^\infty g^2e^{-x}\le4\int_0^\infty g'^2e^{-x}$, and the same holds on $(-\infty,0)$. This law satisfies no Talagrand inequality. Conditioning $\mu$ on $[n,\infty)$ gives a law $\mu_n$ with $\operatorname{KL}(\mu_n\Vert\mu)=n+\log2$. Its mean is $n+1$, whereas $\mu$ has mean $0$, so $W_2^2(\mu_n,\mu)\ge(n+1)^2$. This grows faster than $2\tau(n+\log2)$ for any $\tau$. Every target in an interval of half-length $R$ therefore gives the factor $4R^2/(4R^2+\lambda^2)$, a case outside CCGT's Proposition 1.9.

CCGT's generalised Theorem 1.4, which uses transport inequalities $\mathrm{TI}_\omega$ with general costs, reaches beyond $\mathrm{TI}(\tau)$. We therefore do not claim that convergence for such sources is new; the explicit factor in (31), valid for all $\lambda$ and all bounded initial potentials, is.

*What does not follow.* For quadratic cost with an unbounded target, for example Gaussian to Gaussian, the source gradients $x-y$ are not contained in any ball, so (7) fails. That case needs a different mechanism, such as the semiconcavity control used by CCGT and by Conforti, Durmus and Greco (2023). For general costs, (S3$_\lambda$) follows from $\int e^{\omega/\lambda}d\mu<\infty$. For oscillation of linear growth, Lemma 1 guarantees this only when $\lambda$ is not too small. Theorem 3 avoids that restriction through the compactness of the target and Lemma 2.

### 5.2. Geometry, iteration counts and stopping

The bound has no explicit dimension factor and no target-size factor beyond $L$. The constant $C_P$ may nevertheless depend on the source geometry and dimension. The radius $L$ may depend on the extent of the target. For $B>0$, (23) implies

$$
\Delta(\psi_t)\le
\exp\left[-\frac{t\lambda^2}{B+\lambda^2}\right]\Delta(\psi_0).
\tag{32}
$$

To reduce the gap by a factor $\tau_{\mathrm{tol}}\in(0,1)$, it therefore suffices to run $\lceil(1+B/\lambda^2)\log(1/\tau_{\mathrm{tol}})\rceil$ iterations. This is polynomial in $1/\lambda$. Hilbert-metric analyses instead guarantee a contraction whose deficit is exponentially small in $1/\lambda$ for bounded costs (Franklin and Lorenz 1989; see CDV, Section 1.1).

**Remark 5.1 (a computable stopping certificate).** Combining (9) with (20), under the hypotheses of Theorem 2,

$$
\Delta(\psi_t)\le\left(1+\frac B{\lambda^2}\right)\bigl[E(\psi_{t+1})-E(\psi_t)\bigr].
\tag{33}
$$

The observed increase of the dual objective, multiplied by $1+B/\lambda^2$, therefore bounds the current gap. Using it requires a valid $C_P$ and accurate evaluation of the population integrals.

If an optimal bounded potential $\psi^*$ exists, then (19) and $r_{\psi^*}=\nu$ give

$$
\Delta(\psi)=\lambda\operatorname{KL}(\gamma_{\psi^*}\Vert\gamma_\psi).
\tag{34}
$$

With total variation defined as the supremum over measurable events, Pinsker's inequality then gives $\|\gamma_{\psi^*}-\gamma_\psi\|_{\mathrm{TV}}\le\sqrt{\Delta(\psi)/(2\lambda)}$. The coupling $\gamma_\psi$ is source-normalised, but its target marginal is generally wrong; (34) does not make it feasible.

**Remark 5.2 (where the bound is pessimistic).** In the two-atom examples of Section 7, potentials modulo constants are scalars. A second-order expansion of (20) at the optimum shows that

$$
\frac{E(S\psi)-E(\psi)}{\lambda\operatorname{KL}(\nu\Vert r_\psi)}\to2-D_\lambda
\quad\text{as }\psi\to\psi^*,
$$

where $D_\lambda=D_\lambda(\theta^*_\lambda)$ is as in (41). Near the optimum, the source half-step, which (20) discards, contributes $1-D_\lambda$ times the target gain; this ratio tends to one as $\lambda\to0$. In maximal-correlation terms, the optimal coupling has $\varrho^2=1-D_\lambda$, and the true local factor $(1-D_\lambda)^2=\varrho^4$ is the square of the one-step quantity used by Theorem 1. Linearising Sinkhorn's map at an optimal potential suggests an eventual factor $(B/(B+\lambda^2))^2$ under the hypotheses of Theorem 2. We do not prove this. Nor can (20) be improved globally by a constant factor. For the uniform two-atom source, the increase $E(S\psi_\theta)-E(\psi_\theta)$ is at most $\Delta(\psi_\theta)=\theta+O(1)$, while $\lambda\operatorname{KL}(\nu\Vert r_{\psi_\theta})=\theta+O(1)$ as $\theta\to\infty$. Their ratio therefore tends to $1$. Section 7.4 quantifies the total gap between guarantee and truth on explicit examples.

### 5.3. A deterministic approximation statement

Suppose an implementation produces bounded potentials satisfying

$$
\inf_{a\in\mathbb R}\|\widehat\psi_{t+1}-S\widehat\psi_t-a\|_\infty\le\eta_t.
$$

Write $q=B/(B+\lambda^2)$. Since $E$ is 2-Lipschitz in the supremum norm and invariant under constants,

$$
\Delta(\widehat\psi_t)
\le q^t\Delta(\widehat\psi_0)+2\sum_{j=0}^{t-1}q^{t-1-j}\eta_j.
\tag{35}
$$

If $\eta_j\le\eta$ for all $j$, the error floor is at most $2\eta(1+B/\lambda^2)$. In a quadrature or minibatch implementation, $\eta_t$ is the supremum-norm error of the computed target update. CDV's Appendix C gives a statement of the same form for their rates. Inequality (35) is conditional. It supplies no quadrature or sampling error estimate by itself, and replacing the source by a finite empirical measure does not preserve (6).

## 6. A global rate linear in the regularisation for bounded densities

**Corollary 5 (CDV settings).** Assume CDV's condition (A). That is, $X\subset\mathbb R^d$ is convex, each $c(\cdot,y)$ is $\xi$-semiconcave for some $\xi\ge0$, and $c_\infty=\operatorname{osc}_{X\times Y}c<\infty$. Define

* $C=e\,(c_\infty+\xi R_X^2/2)\,\kappa$ under CDV's (A1): $X\subset\{|x|\le R_X\}$ is compact and $\mu$ has a density $f$ with $\sup f/\inf f=\kappa<\infty$;
* $C=e\max\{c_\infty,1\}$ under CDV's (A2): $\mu$ has density $e^{-V}$ with $V$ $\xi$-strongly convex; or under CDV's (A3): each $c(\cdot,y)$ is also $\zeta$-semiconvex for some $\zeta\ge0$, and $\mu$ has density $e^{-V}$ with $V$ $\max\{\xi,(\xi+\zeta)/\lambda\}$-strongly convex. When $\xi=0$ or $c_\infty\ge1$ this is CDV's constant $e\,c_\infty$.^[CDV state $C=e\,c_\infty$ under (A2)/(A3). Their proof rescales the problem to $(\alpha c,\alpha\lambda)$ with $\alpha=1/c_\infty$ and needs $V$ to be $\alpha\xi$-strongly convex, which (A2) guarantees only when $\alpha\le1$ or $\xi=0$. Taking $\alpha=\min\{1,1/c_\infty\}$ in the same argument gives $C=\alpha^{-1}e^{\alpha c_\infty}\le e\max\{c_\infty,1\}$, which is what we use. The same choice works under (A3).]

Then every bounded initial potential satisfies

$$
\boxed{\Delta(\psi_t)\le
\left(\frac C{C+\lambda}\right)^t\Delta(\psi_0).}
\tag{36}
$$

When $C=0$, one full iteration reaches zero gap.

*Proof.* CDV's Proposition 5.1 states that, under (A), $\frac{d^2}{ds^2}K(\psi+sv)\le-(C+\lambda)^{-1}\operatorname{Var}_{\nu[\psi+sv]}(v)$ for any $\psi$, $v$ and $s$, with the values of $C$ above. Here $K=-\Phi$, and $\nu[\psi]$ is the second marginal of the current coupling, which is our $r_\psi$. At $s=0$ this is (8) with $A=C+\lambda$, for every bounded $\psi$ and $h=v$. The cost is bounded, so (S1)--(S3$_\lambda$) hold. Theorem 1 gives (36). $\square$

For bilinear cost, $\xi=0$ and $C=e\kappa c_\infty$. The whole curvature estimate is CDV's; the new step is its conversion into a global marginal-entropy bound by Theorem 1, which is the improvement anticipated in CDV's Remark 5.2. Under (A1), CDV's Theorem 3.1 has global coefficient $\alpha=176\{1+(c_\infty+\xi R_X^2/2)\kappa\lambda^{-1}+c_\infty^2\lambda^{-2}\}$, a deficit of order $\lambda^2$. Their Corollary 3.2 improves this to $176\{1+(c_\infty+\xi R_X^2/2)\kappa\lambda^{-1}\}$ after an instance-dependent number of iterations, and needs either a Lipschitz cost with bounded $Y$ or a finitely supported $\nu$. Under (A2), the corresponding coefficients are $176\{1+c_\infty\lambda^{-1}+c_\infty^2\lambda^{-2}\}$ globally and $176\{1+c_\infty\lambda^{-1}\}$ eventually. Corollary 5 obtains linear dependence on $\lambda$ from iteration zero in both cases, without CDV's extra conditions and without the factor $176$. Under (A3), which includes CDV's Theorem 1.1 (log-concave source and bilinear cost), CDV's global coefficient is already $176\{1+c_\infty\lambda^{-1}\}$, so Corollary 5 improves only the constant. No lower bound on target atom masses is needed.

*The linear order is optimal.* Take the uniform source on $[-1,1]$, the target $(\delta_{-1}+\delta_1)/2$ and $c=-xy$. This satisfies (A1) with $\kappa=1$, $\xi=0$ and $c_\infty=2$. Section 7.2 shows that its local dual-gap contraction factor is $(1-D_\lambda)^2$ with

$$
D_\lambda=\lambda\tanh(1/\lambda),\qquad 1-(1-D_\lambda)^2\sim2\lambda.
\tag{37}
$$

No guarantee of the form $1-K\lambda^\beta$ with $\beta<1$ can therefore hold in the bounded-density class, so the order $\lambda$ in (36) cannot be improved. CDV's Theorem 1.3 gives a Gaussian lower bound with a deficit linear in $\lambda$, in a non-compact setting.

The table compares guarantees as $\lambda\downarrow0$, with the other parameters fixed. A deficit $d$ means a contraction factor $1-d$. "Worst case over class" means the order that the guarantee implies uniformly over its hypothesis class; individual instances may converge faster.

| Guarantee | Functional-inequality marginal | Cost | Small-$\lambda$ deficit | Validity |
|:----|:-------|:------|:--:|:-----|
| Hilbert metric (e.g. Franklin and Lorenz 1989) | none | bounded | $e^{-\Theta(c_\infty/\lambda)}$ | all iterations |
| CDV Thm 3.1, (A1)/(A2) | bounded density / $\xi$-strongly log-concave | bounded, semiconcave, convex $X$ | $\Theta(\lambda^2)$ | all iterations |
| CDV Cor 3.2, (A1)/(A2) | same | same, plus Lipschitz or finite $\nu$ | $\Theta(\lambda)$ | after an instance-dependent index |
| CDV Thm 3.1, (A3) (incl. Thm 1.1) | $\max\{\xi,(\xi+\zeta)/\lambda\}$-strongly log-concave | bounded, $\xi$-semiconcave and $\zeta$-semiconvex, convex $X$ | $\Theta(\lambda)$ | all iterations |
| CCGT Prop 1.8(i) | $\mathrm{TI}(\tau)$ | Lipschitz in Talagrand variable, Hessian bounds; unbounded allowed | $\Theta(\lambda^2)$ | all $n\ge1$; explicit factor for $\varepsilon\le\varepsilon_m$ |
| CCGT Prop 1.9(i) | $\mathrm{TI}(\tau)$ | constant Hessian (e.g. quadratic); other marginal compact | $\Theta(\lambda^2)$ | all $n\ge1$; explicit factor for $\varepsilon\le\sqrt\tau\|g\|$ |
| Theorem 2 here | Poincaré $C_P$ | first derivative, (7); (S1)--(S3) | $\Theta(\lambda^2)$ | all iterations, all $\lambda$ |
| Theorem 3, Corollary 4 here | Poincaré $C_P$ | e.g. quadratic, unbounded; target compact | $\Theta(\lambda^2)$ | all iterations, all $\lambda$ |
| Corollary 5 here | as CDV (A1)/(A2)/(A3) | as CDV (A) | $\Theta(\lambda)$ | all iterations |

For the Poincaré and Talagrand classes (CCGT's Proposition 1.9(i), Theorems 2--3 and Corollary 4) the order $\lambda^2$ is sharp as a worst-case statement over the class (Section 7); for CDV's bounded-density and log-concave classes, Corollary 5 shows that order $\lambda$ holds globally. Within the Poincaré class, the obstruction comes from sources with little mass where the conditional law of the target switches. Proposition 7 shows that density ratios $\kappa_\lambda=1/\lambda$ reconcile the two regimes.

## 7. The exponent two is sharp

All examples in this section use $X=(-1,1)$, $\nu=(\delta_{-1}+\delta_1)/2$ and $c(x,y)=-xy$, so that (7) holds with $L=1$. The source $\mu$ has a density $w$.

### 7.1. A fixed source that prevents every smaller exponent

**Theorem 6 (fixed-source obstruction).** Let

$$
w(x)=\frac25|x|\bigl[\log(e/|x|)\bigr]^2\quad(0<|x|<1),
\qquad w(0)=0.
\tag{38}
$$

1. The source satisfies (6) with $C_P\le3/4$.
2. For each $\lambda>0$, define the *local dual-gap contraction factor* as $\lim\Delta(S\psi)/\Delta(\psi)$ as $\psi$ tends to an optimal potential, with $\psi$ not optimal. This limit exists and equals $(1-D_\lambda)^2$, where $D_\lambda=\int_{-1}^1\operatorname{sech}^2(x/\lambda)w(x)\,dx$. Every Sinkhorn trajectory started from a non-optimal bounded potential satisfies $\Delta(\psi_{t+1})/\Delta(\psi_t)\to(1-D_\lambda)^2$.
3. As $\lambda\downarrow0$,
$$
1-(1-D_\lambda)^2\sim\frac85\log2\,\lambda^2\log^2(1/\lambda).
\tag{39}
$$

Consequently, there are no $\beta<2$, $K>0$ and $\lambda_0>0$ such that, for every $0<\lambda<\lambda_0$, some non-optimal trajectory obeys a bound $\Delta(\psi_t)\le K'(1-K\lambda^\beta)^t\Delta(\psi_0)$ for all $t$, with $K'$ allowed to depend on $\lambda$ and $\psi_0$. In particular, no eventual rate and no one-step bound of the form $1-K\lambda^\beta$ with $\beta<2$ holds for this fixed problem.

The last statement follows from part 2: along a non-optimal trajectory, $\Delta(\psi_t)^{1/t}\to(1-D_\lambda)^2$. Any bound of the stated form therefore forces $K\lambda^\beta\le1-(1-D_\lambda)^2$, which contradicts (39) for small $\lambda$.

### 7.2. Proof of Theorem 6

*Poincaré certificate.* The density is normalised because $\int_0^1x(1-\log x)^2\,dx=5/4$. For $0<x<1$,

$$
R(x):=\int_0^x\frac{dt}{w(t)}=\frac{5}{2\log(e/x)}.
$$

Weighted Cauchy--Schwarz gives $(f(x)-f(0))^2\le R(x)\int_0^xf'(t)^2w(t)\,dt$. Integrating and exchanging the nonnegative integrals yields

$$
\int_0^1(f-f(0))^2w
\le\int_0^1 f'(t)^2w(t)
\left(\int_t^1 R(x)w(x)\,dx\right)dt
\le\frac34\int_0^1 f'(t)^2w(t)\,dt,
$$

since $R(x)w(x)=x\log(e/x)$ and its integral on $(0,1)$ is $3/4$. The same bound holds on the negative half, and $\operatorname{Var}_\mu(f)\le\int(f-f(0))^2\,d\mu$. The vanishing of the density at zero therefore does not destroy the Poincaré inequality.

*An exact scalar Sinkhorn map.* The argument in this paragraph applies to any density $w$ on $(-1,1)$. Every target potential is, modulo constants, $\psi_\theta(y)=\theta y$. Let $m_\lambda(\theta)=\mathbb E_\mu\tanh((X+\theta)/\lambda)$. The target probabilities of $r_{\psi_\theta}$ are $(1\pm m_\lambda(\theta))/2$, so (4) gives

$$
T_\lambda(\theta)=\theta-\lambda\operatorname{artanh}m_\lambda(\theta).
\tag{40}
$$

The objective is $E(\psi_\theta)=-\lambda\,\mathbb E_\mu\log\cosh((X+\theta)/\lambda)$. It is maximised at the unique zero $\theta^*_\lambda$ of the strictly increasing function $m_\lambda$. Write $D_\lambda(\theta)=\mathbb E_\mu\operatorname{sech}^2((X+\theta)/\lambda)$. Then

$$
\Delta(\psi_\theta)=\frac{D_\lambda(\theta^*_\lambda)}{2\lambda}(\theta-\theta^*_\lambda)^2+O(|\theta-\theta^*_\lambda|^3),\qquad
T_\lambda'(\theta^*_\lambda)=1-D_\lambda(\theta^*_\lambda).
\tag{41}
$$

Since $T_\lambda(\theta^*_\lambda)=\theta^*_\lambda$, it follows that

$$
\lim_{\theta\to\theta^*_\lambda,\ \theta\ne\theta^*_\lambda}
\frac{\Delta(\psi_{T_\lambda(\theta)})}{\Delta(\psi_\theta)}
=(1-D_\lambda(\theta^*_\lambda))^2.
\tag{42}
$$

For the symmetric density (38), $\theta^*_\lambda=0$ and $D_\lambda(0)=D_\lambda$. For $\theta>0$, pairing $x$ with $-x$ shows that $0<m_\lambda(\theta)<\tanh(\theta/\lambda)$, since the paired average is

$$
\frac{\sinh(2\theta/\lambda)}{\cosh(2\theta/\lambda)+\cosh(2x/\lambda)}.
$$

Thus $0<T_\lambda(\theta)<\theta$. The iterates decrease to a fixed point of the continuous map $T_\lambda$, which can only be $0$, and they never reach it. Oddness handles negative initial values. By (42), ratios of successive gaps along every non-optimal trajectory converge to $(1-D_\lambda)^2$, which proves part 2.

*The small-regularisation limit.* Substituting $x=\lambda u$ gives

$$
D_\lambda=\frac45\lambda^2
\int_0^{1/\lambda}u\bigl[\log(e/\lambda)-\log u\bigr]^2
\operatorname{sech}^2u\,du.
\tag{43}
$$

After division by $\lambda^2\log^2(1/\lambda)$, the integrand converges to $(4/5)u\operatorname{sech}^2u$. For $\lambda\le e^{-1}$ it is dominated by a constant times $u(2+|\log u|)^2\operatorname{sech}^2u$, which is integrable. Since $\int_0^\infty u\operatorname{sech}^2u\,du=\log2$, dominated convergence gives $D_\lambda\sim(4/5)\log2\,\lambda^2\log^2(1/\lambda)$. As $1-(1-D)^2=D(2-D)$, (39) follows. $\square$

*Refined expansion.* The leading term in (39) is accurate only slowly: the ratio of $D_\lambda$ to $(4/5)\log2\,\lambda^2\log^2(1/\lambda)$ is $2.20$, $1.53$, $1.33$ and $1.16$ at $\lambda=10^{-1},10^{-2},10^{-3},10^{-6}$. Expanding the square in (43) and extending the integrals to $(0,\infty)$ gives, with $\Lambda=\log(e/\lambda)$,

$$
D_\lambda=\frac45\lambda^2\bigl[\log2\,\Lambda^2-2a\Lambda+b\bigr]+\mathcal R_\lambda,
\qquad |\mathcal R_\lambda|\le69\,\lambda^{-1}e^{-2/\lambda}\quad(0<\lambda\le1).
\tag{44}
$$

The constants are

$$
a=\int_0^\infty u\log u\operatorname{sech}^2u\,du=\log2\left(1-\tfrac32\log2\right)=-0.0275323403\ldots,
$$

$$
b=\int_0^\infty u\log^2u\operatorname{sech}^2u\,du=\log2\left(\tfrac{\pi^2}6-\gamma_E^2-2\gamma_1-3\log2+\tfrac73\log^22\right)=0.3458827476\ldots,
$$

where $\gamma_E$ is Euler's constant and $\gamma_1$ the first Stieltjes constant. These follow from $\int_0^\infty u^s\operatorname{sech}^2u\,du=2^{1-s}\Gamma(s+1)\eta(s)$, where $\eta(s)=(1-2^{1-s})\zeta(s)$, by differentiating at $s=1$. The expansion follows from $\operatorname{sech}^2u\le4e^{-2u}$ and $|\Lambda-\log u|\le1+2\log u\le3u$ for $u\ge1/\lambda\ge1$, and the integral $\int_M^\infty u^3e^{-2u}du\le\frac{19}8M^3e^{-2M}$ for $M\ge1$.

### 7.3. Exact quadratic order uniformly over a Poincaré class

The logarithmic factor disappears when the source may vary with the regularisation while its Poincaré constant stays bounded.

**Proposition 7 (uniform-class sharpness).** For $0<\lambda\le1/2$, let

$$
\rho_\lambda(x)=\min\{\lambda\cosh(x/\lambda),1\},\qquad
w_\lambda(x)=\frac{\rho_\lambda(x)}{2I_\lambda},\quad |x|<1,
\qquad
I_\lambda=1-\lambda\operatorname{arcosh}(1/\lambda)+\lambda\sqrt{1-\lambda^2}.
\tag{45}
$$

Every $w_\lambda$ satisfies (6) with $C_P\le1+\pi/2$. Its local dual-gap contraction deficit satisfies

$$
1-(1-D_\lambda)^2\sim\pi\lambda^2.
\tag{46}
$$

*Proof.* The transition point is $x_*=\lambda\operatorname{arcosh}(1/\lambda)<1$. Integrating on either side gives $\int_0^1\rho_\lambda=I_\lambda\le1$. The anchored Cauchy--Schwarz argument of Section 7.2, with $\int_t^1R\,w\le R(1)\int_0^1w$, gives

$$
C_P\le\left(\int_0^1\rho_\lambda\right)\left(\int_0^1\rho_\lambda^{-1}\right)\le1+\pi/2,
\qquad
\int_0^1\rho_\lambda(x)^{-1}\,dx=\arccos\lambda+1-\lambda\operatorname{arcosh}(1/\lambda).
$$

The exact conditional-variance integral is

$$
D_\lambda=
\frac{\lambda^2\arccos\lambda+
\lambda\bigl[\tanh(1/\lambda)-\sqrt{1-\lambda^2}\bigr]}{I_\lambda}.
\tag{47}
$$

Here $I_\lambda\to1$, $\arccos\lambda\to\pi/2$, and the second term of the numerator is $O(\lambda^3)$. Hence $D_\lambda\sim(\pi/2)\lambda^2$. The density is symmetric, so $\theta^*_\lambda=0$, and (42) proves (46). $\square$

Theorem 2 and Proposition 7 together show that the worst local deficit over this uniformly bounded Poincaré class has exact order $\lambda^2$. For each fixed $\lambda$, the density in Proposition 7 is bounded above and below with ratio $\kappa_\lambda=1/\lambda$. Corollary 5 then gives $C=2e/\lambda$ and a deficit of order $\lambda^2$, consistent with (46). Its linear conclusion applies when the density ratio is held fixed.

### 7.4. No fixed source attains order $\lambda^2$

Proposition 7 attains order $\lambda^2$ only with sources that vary with $\lambda$. This is unavoidable in the two-atom model.

**Lemma 3 (Poincaré forces integrable reciprocal density).** Let $\mu=w\,dx$ satisfy (6) on $(-1,1)$. Then (6) extends to Lipschitz $f$, with $f'$ the almost-everywhere derivative. Moreover, for $-1<a<b<1$ with $\mu((-1,a))>0$ and $\mu((b,1))>0$,

$$
\int_a^b\frac{dx}{w(x)}\le\frac{C_P}{\mu((-1,a))\,\mu((b,1))}.
\tag{48}
$$

*Proof.* For Lipschitz $f$, mollify $f'$ and integrate. This gives bounded smooth $f_k\to f$ uniformly, with $f_k'\to f'$ almost everywhere and $|f_k'|\le\|f'\|_\infty$. Dominated convergence passes (6) to the limit. Next, let $g_n=\min(1/w,n)\mathbf 1_{(a,b)}$, with $1/w=\infty$ where $w=0$, and let $f_n(x)=\int_{-1}^xg_n$. Then $f_n=0$ on $(-1,a]$ and $f_n=F_n:=\int_a^bg_n$ on $[b,1)$, so $\operatorname{Var}_\mu(f_n)\ge\mu((-1,a))\mu((b,1))F_n^2$. Since $g_nw\le1$, we have $\int f_n'^2w\le F_n$. Hence $F_n\le C_P/(\mu((-1,a))\mu((b,1)))$, and monotone convergence proves (48). $\square$

Lemma 3 is in the spirit of the necessity part of Muckenhoupt's (1972) criterion for weighted Hardy inequalities.

**Proposition 8 (no fixed source attains order $\lambda^2$).** Let $\mu=w\,dx$ be any probability measure on $(-1,1)$ satisfying (6). In the two-atom model, its local dual-gap contraction deficit $\delta_\lambda=1-(1-D_\lambda(\theta^*_\lambda))^2$, given by (42), satisfies

$$
\delta_\lambda/\lambda^2\to\infty\qquad(\lambda\downarrow0).
$$

*Proof.* Let $s_1<s_2$ be the endpoints of the convex hull of $\operatorname{supp}\mu$. For $s_1<a<b<s_2$ both side masses in (48) are positive, so $1/w\in L^1_{\mathrm{loc}}((s_1,s_2))$. In particular, $w>0$ almost everywhere there, and the distribution function $F$ of $\mu$ is continuous and strictly increasing on $[s_1,s_2]$, with a unique median $x_m\in(s_1,s_2)$. Put $x_\lambda=-\theta^*_\lambda$, so that $\mathbb E_\mu\tanh((X-x_\lambda)/\lambda)=0$. Choose $\varkappa>0$ with $[x_m-3\varkappa,x_m+3\varkappa]\subset(s_1,s_2)$. If $x_\lambda\ge x_m+2\varkappa$, then $\tanh((X-x_\lambda)/\lambda)\le-\tanh(\varkappa/\lambda)$ on $\{X\le x_m+\varkappa\}$, so

$$
0=\mathbb E_\mu\tanh((X-x_\lambda)/\lambda)\le1-F(x_m+\varkappa)-F(x_m+\varkappa)\tanh(\varkappa/\lambda).
$$

The right side tends to $1-2F(x_m+\varkappa)<0$ as $\lambda\downarrow0$. Hence $x_\lambda<x_m+2\varkappa$ for small $\lambda$, and symmetrically $x_\lambda>x_m-2\varkappa$. For $\lambda<\varkappa$, the window $J_\lambda=[x_\lambda-\lambda,x_\lambda+\lambda]$ lies in $[x_m-3\varkappa,x_m+3\varkappa]$. Using $\operatorname{sech}^2\ge\operatorname{sech}^21$ on $J_\lambda$ and Cauchy--Schwarz,

$$
D_\lambda(\theta^*_\lambda)\ge\operatorname{sech}^2(1)\int_{J_\lambda}w\ge\operatorname{sech}^2(1)\frac{(2\lambda)^2}{\int_{J_\lambda}w^{-1}}.
$$

The function $1/w$ is integrable on the compact interval $[x_m-3\varkappa,x_m+3\varkappa]$, so $\sup_{J}\int_Jw^{-1}\to0$ over subintervals $J$ of length $2\lambda$. Hence $D_\lambda(\theta^*_\lambda)/\lambda^2\to\infty$. Since $\delta_\lambda=D(2-D)\ge D$ with $D=D_\lambda(\theta^*_\lambda)\le1$, the proof is complete. $\square$

In the two-atom model, order $\lambda^2$ is therefore the worst case uniformly over a bounded Poincaré class (Theorem 2 and Proposition 7), but no single source attains it, even eventually. Theorem 6 shows that a fixed source can come within a factor $\log^2(1/\lambda)$.

**Remark 7.1 (power densities).** The elementary obstruction uses $w_a(x)=\frac{a+1}2|x|^a$ with $0\le a<1$. The anchored argument of Section 7.2 gives $R(x)w_a(x)=x/(1-a)$, hence $C_P\le1/(2(1-a))$. Substituting $x=\lambda u$ gives $D_\lambda=(a+1)\lambda^{1+a}\int_0^{1/\lambda}u^a\operatorname{sech}^2u\,du$. The local deficit is therefore asymptotic to $2(a+1)c_a\lambda^{1+a}$, with $c_a=2^{1-a}\Gamma(a+1)\eta(a)$; for example, $c_{1/2}=0.7581\ldots$. For any $\beta<2$, choosing $a\in(\beta-1,1)$ defeats a deficit $K\lambda^\beta$. The Poincaré constants of this family blow up as $a\uparrow1$. The logarithmic cusp (38) is the borderline weight that defeats every $\beta<2$ with one source.

**Remark 7.2 (how far Theorem 2 is from the truth).** We estimated the true Poincaré constants numerically by Sturm--Liouville shooting in natural-scale coordinates, with an independent finite-element check for the cusp (Section 8). These are numerical estimates, not certificates. The cusp (38) gives $C_P\approx0.634$, against the certificate $3/4$. The power sources with $a=1/2$ and $a=0.9$ give $C_P\approx0.893$ and $4.879$, against the certificates $1$ and $5$. The family of Proposition 7 gives $C_P\approx0.604,\,1.29,\,1.80,\,1.85$ at $\lambda=0.5,\,0.1,\,0.01,\,0.005$, against the certificate $1+\pi/2\approx2.571$. For that family, the ratio of the exact local deficit to Theorem 2's guaranteed deficit $\lambda^2/(C_P+\lambda^2)$ is about $8.3$ at $\lambda=0.005$ with the certified constant, tending to $\pi(1+\pi/2)\approx8.1$. With the numerical constant it is about $6.0$. A factor of about two is the squaring discussed in Remark 5.2. The rest comes from the Cauchy--Schwarz and Poincaré steps (26)--(27) and, with the certified constant, from slack in $C_P$. Theorem 2 is pessimistic by a bounded factor on this class, not by a power of $\lambda$.

### 7.5. The exponent is also sharp for Talagrand sources

**Corollary 9.** Let $\kappa_G=\sqrt{18e^{\sqrt5}}$. Every probability density on $(-1,1)$ satisfying (6) satisfies $\mathrm{TI}(\tau)$, that is $W_2^2(\cdot,\mu)\le2\tau\operatorname{KL}(\cdot\Vert\mu)$, with

$$
\tau=\tfrac12\max\{\kappa_G^2C_P,\ 2\kappa_G\sqrt{C_P}\}.
\tag{49}
$$

Hence the sources of Proposition 7 satisfy $\mathrm{TI}(\bar\tau)$ with a constant $\bar\tau$ independent of $\lambda$, and the fixed source (38) satisfies $\mathrm{TI}(\tau)$ for some finite $\tau$. In the setting of CCGT's Proposition 1.9(i), the exponent $2$ in their guaranteed deficit $\lambda^2/(\lambda^2+\tau\|g\|^2)$ can therefore not be lowered. Over the class $\mathrm{TI}(\bar\tau)$ the worst local deficit is $\Theta(\lambda^2)$, and for the fixed Talagrand source (38) no bound with deficit $K\lambda^\beta$, $\beta<2$, holds even eventually. The setting is: compact second marginal (our two-atom target), Talagrand marginal (our source), and $\nabla^2$ constant in the Talagrand variable (here $0$), with $\|g\|=1$.

*Proof.* By Lemma 3, (6) holds for Lipschitz functions on $(-1,1)$. A locally Lipschitz function on $\mathbb R$ restricts to such a function, and its almost-everywhere derivative is bounded by its local slope. Hence $\mu$, viewed as an absolutely continuous measure on $\mathbb R$, satisfies Gozlan's Poincaré inequality $\mathbb{SG}(C_P)$. Gozlan's Theorem 4.7 (Gozlan 2010), which restates Bobkov, Gentil and Ledoux (2001, Corollary 5.1), then gives $\inf_\pi\int\alpha(|x-y|/(\sqrt{C_P}\kappa_G))\,d\pi\le\operatorname{KL}(\cdot\Vert\mu)$, where $\alpha(u)=\min(|u|,u^2)$. On $[-1,1]$ we have $|x-y|\le2$, and with $s=\sqrt{C_P}\kappa_G$, $\alpha(r/s)\ge r^2/\max\{s^2,2s\}$ for $0\le r\le2$. This gives (49). CCGT's error measure $\operatorname{KL}(\pi^\varepsilon\Vert\pi^{n,n})$ is our $\Delta/\lambda$, by (34), along the same sequence of couplings. The final statements therefore follow from Proposition 7 and from Theorem 6, part 2. $\square$

The absolute constant $\kappa_G$ is large, so (49) is not a good numerical estimate of $\tau$. Only its uniformity in $\lambda$ is used. Three clarifications. First, what is sharp is the exponent $2$, not the constant in CCGT's factor. Second, in these symmetric examples CCGT's default initialisation $\varphi^0=0$, which is our $\psi_0=0$, is already optimal; the sharpness statements concern non-optimal starting potentials, which CCGT's results allow. Third, the pairing and trajectory argument of Section 7.2 uses only the symmetry of the density, so it applies verbatim to the family of Proposition 7: every non-optimal trajectory has gap ratios converging to $(1-D_\lambda)^2$.

### 7.6. What the examples say about CDV's conjecture

CDV's Remark 3.1 discusses the convexity assumption on the support and conjectures that a Poincaré inequality on the source "should be a sufficient condition for our results to hold". For costs satisfying (7), Theorem 2 confirms this for a global rate of order $\lambda^2$, with no convexity, density or Hessian assumption. The source (38) satisfies CDV's assumption (A): it lives on a convex interval and the cost is bilinear and bounded. Yet by Theorem 6 no eventual rate $1-\Theta(\lambda^\beta)$ with $\beta<2$ holds for it. CDV's eventual linear rate (Corollary 3.2) therefore does not extend to the Poincaré class once the density bound is dropped. A different reading of the remark keeps the bounded density and replaces only convexity of the support by a Poincaré inequality, as for a dumbbell-shaped domain. Under that reading the question of a global linear rate remains open; our examples do not address it.

## 8. Numerical checks and reproducibility

The proofs above do not use numerical evidence. Three Python programmes check the displayed formulas by floating-point computation; they are not interval certificates.

* `code/sharpness_check.py` treats the uniform source and the cusp (38) with two target atoms. It evaluates exact gaps, the error bound (22), the curvature condition and the ascent inequality over 180 combinations of $\lambda$ and potential. It also follows ten nonzero trajectories, whose gap ratios approach $(1-D_\lambda)^2$ to within $3.5\times10^{-10}$. The cusp quadrature gives $D_{0.1}=0.0647035371183$, $D_{0.01}=0.00179454249623$ and $D_{0.001}=0.0000353004893444$.
* `code/multitarget_gaussian_check.py` uses a Gaussian source, the bounded cost $\sin(x)y$ and five target atoms. It checks the same inequalities in a non-compact, multi-atom setting.
* `code/revision_checks.py` (new in this version) checks several quantities: the refined expansion (44), including both closed-form constants to 25 digits and the remainder bound; the constants $c_a$ of Remark 7.1; and the numerical Poincaré constants of Remark 7.2 at two solver tolerances, together with the uniform-source value $4/\pi^2$ as a sanity check. It also illustrates Corollary 4 with the two-sided exponential source, quadratic cost and three target atoms at $\lambda\in\{1,0.5,0.3\}$. There, the observed one-step gap ratios ($0.136$, $0.363$, $0.552$) stay below the guaranteed $4/(4+\lambda^2)$, and the certificate (33) holds. Finally, it checks the local limit $2-D_\lambda$ of Remark 5.2. An independent finite-element computation in the same natural-scale coordinates reproduces the cusp constant $0.634118$.

Discretising a source for these checks does not make the discrete measure satisfy (6). The family of Proposition 7 is verified analytically; its numerical constants enter only Remark 7.2.

![Asymptotic one-step improvement $1-r_\lambda^2$ of the dual gap near the optimum, plotted against $\lambda$ on logarithmic axes, where $r_\lambda=1-D_\lambda$ is the derivative of the scalar Sinkhorn map (40) at the optimum. Upper curve: uniform source, deficit about $2\lambda$. Middle curve: the fixed cusp (38), deficit of order $\lambda^2\log^2(1/\lambda)$. Dashed line: the guaranteed deficit $\lambda^2/(C_P+\lambda^2)$ of Theorem 2 with $C_P=3/4$.](figures/sharpness_plot.png){width=82%}

## 9. Detailed comparison with prior work, and open questions

### 9.1. Chiarini, Conforti, Greco and Tamanini

CCGT work in the orientation where the Talagrand marginal is called $\nu$ and the other marginal is called $\rho$, with regularisation $\varepsilon$. Their Theorem 1.3 assumes that $y\mapsto c(x,y)+\psi^n_\varepsilon(y)$ is $\Lambda$-semiconcave uniformly in $x$ and $n\ge N-1$. Under $\mathrm{TI}(\tau)$ it gives contraction of $\operatorname{KL}(\pi^\varepsilon\Vert\pi^{n+1,n+1})$ by $1-\min\{\tau\Lambda,\varepsilon\}/(\min\{\tau\Lambda,\varepsilon\}+\tau\Lambda)$. Proposition 1.8(i) treats open connected, possibly unbounded, domains in $\mathbb R^d$ with $h(c)\preceq\nabla_2^2c\preceq H(c)$ and $\nabla_1^2c\preceq H(c)$. It also assumes the Lipschitz bound $\sup_x|c(x,y)-c(x,\hat y)|\le\operatorname{Lip}_{\infty,2}|y-\hat y|$ in the Talagrand variable. For $\varepsilon$ below an explicit threshold, the rate is $1-\varepsilon^2/(\varepsilon^2+\tau\varepsilon(H-h)+\tau\operatorname{Lip}^2_{\infty,2})$. Proposition 1.9(i) assumes that $\rho$ has compact support and that $\nabla_2^2c=\Sigma$ is constant. For $\varepsilon\le\sqrt\tau\|g\|_{L^\infty(\rho)}$, with $g(x)=\nabla_2c(x,y)-\Sigma y$, the rate is $1-\varepsilon^2/(\varepsilon^2+\tau\|g\|^2_{L^\infty(\rho)})$. Proposition 1.9(ii), with both marginals compact, gives $1-\varepsilon^2/(\varepsilon^2+\tau\varepsilon(H-h)+\tau\|\nabla_2c\|^2_\infty)$. CCGT's convention (1.10) means that these rates hold for all $n\ge1$. For larger $\varepsilon$ their rate is still explicit; for example, in Proposition 1.9(i) with $\varepsilon>\sqrt\tau\|g\|$, Theorem 1.3 gives the factor $1/2$.

| Aspect | CCGT Prop. 1.9(i) | CCGT Prop. 1.8(i) | Theorems 2--3 here |
|:------|:------|:------|:------|
| Inequality on one marginal | $\mathrm{TI}(\tau)$ | $\mathrm{TI}(\tau)$ | Poincaré, $C_P\le\tau$ |
| Other marginal | compact support | open domain in $\mathbb R^d$ | arbitrary (Thm 2); compact (Thm 3) |
| Cost regularity | $\nabla^2c$ constant in Talagrand variable | two-sided Hessian bounds; Lipschitz | $C^1$ in source variable; (7) |
| Bounded cost needed | no | no | Thm 2: (S1)--(S3); Thm 3: no |
| Range of $\lambda$ for the stated factor | $\lambda\le\sqrt\tau\|g\|$ | below an explicit threshold | every $\lambda>0$ |
| Validity in $n$ | bound (1.10) for all $n\ge1$ | same | every bounded initial potential, $t\ge0$ |
| Error measure | $\operatorname{KL}(\pi^\varepsilon\Vert\pi^{n+1,n+1})$ | same | $\Delta=\lambda\operatorname{KL}(\gamma^*\Vert\gamma_\psi)$, same sequence |
| Factor | $\tau\|g\|^2/(\tau\|g\|^2+\lambda^2)$ | includes $\tau\lambda(H-h)$ | $C_PL^2/(C_PL^2+\lambda^2)$ |

On compact sources, Poincaré and Talagrand inequalities are equivalent up to constants that depend on the diameter (Bobkov, Gentil and Ledoux 2001; Gozlan 2010, Theorem 4.7). There, our gain lies in constants, cost regularity and the range of $\lambda$. Poincaré is strictly weaker on unbounded supports, as the two-sided exponential law shows. CCGT treat several settings that our results do not cover, such as log-concave or light-tailed second marginals, which may be unbounded (their Propositions 1.5 and 1.7), and Riemannian manifolds. They also obtain rates linear in $\varepsilon$ under log-concavity. Finally, CCGT note, citing Di Marino and Gerolin (2020, Proposition 2.4), that the gradient of a Sinkhorn potential is a conditional mean of the cost gradient. Identity (25) is the same mechanism, written for the conditional mean of a test function.

### 9.2. Other work

* **Relative geometry.** Léger (2021) gave a gradient-descent perspective on Sinkhorn with sublinear rates. AFKL (2022) proved the mirror-descent framework discussed in Remark 4.2. Dual-space preconditioning (Maddison et al. 2021) is the finite-dimensional counterpart of the conjugacy argument in Section 4. Chen and Xu (2026) use dual mirror geometry to accelerate finite Sinkhorn problems. Theorem 1 contributes a polynomial reverse-entropy coefficient obtained from tilted curvature. It does not claim to invent the geometry.
* **Hilbert-metric and exponential rates.** The matrix-scaling iteration goes back to Sinkhorn (1964) and Sinkhorn and Knopp (1967). Linear convergence for positive kernels follows from contraction in Hilbert's projective metric (Franklin and Lorenz 1989), with a guarantee that deteriorates exponentially as the regularisation decreases (CDV, Section 1.1). Carlier (2022) gave an elementary proof of linear convergence for multimarginal problems. Eckstein (2025) extended Hilbert-metric contraction to unbounded settings in which the marginals have sufficiently light tails relative to the growth of the cost.
* **Unbounded and smooth settings.** Conforti, Durmus and Greco (2023) proved exponential convergence for quadratic cost on $\mathbb R^d$ when both marginals have asymptotically positive log-concavity profiles. Greco, Noble, Conforti and Durmus (2023) analysed the torus by a coupling approach. Ghosal and Nutz (2026) proved non-asymptotic sublinear rates for a large class of costs and marginals, including quadratic cost with sub-Gaussian marginals. Delalande (2022) quantified semi-discrete entropic transport as a function of the regularisation. Entropic transport and its applications are surveyed by Peyré and Cuturi (2019). Cuturi (2013) popularised Sinkhorn's algorithm for computing transport distances.
* **Continuous-time structure.** Srinivasan and Slotine (2026) formulate conditional-expectation and spectral-gap structures for a continuous-time Sinkhorn flow. Equation (27) gives a uniform quantitative bound for the corresponding quantity, and Theorem 1 turns it into a discrete-iteration guarantee.

We searched the literature on 1 November 2026 and again during this revision (web search, arXiv full text, publisher pages and citation links). We found no earlier statement of the explicit bounds or the sharpness constructions proved here. A negative search is not a proof of priority, and the unpublished details of Di Marino's announcements remain a material uncertainty.

### 9.3. Open questions

1. Is the eventual factor $(B/(B+\lambda^2))^2$ suggested in Remark 5.2 correct, and what is the sharp universal constant?
2. Does a bounded density on a non-convex Poincaré domain give a global rate linear in $\lambda$ (Section 7.6)?
3. Can unbounded targets with quadratic cost, such as Gaussian to Gaussian, be treated under one Poincaré inequality?
4. Empirical sources do not satisfy (6). A discrete analogue would need a graph Poincaré inequality, together with a discrete version of the covariance identity (25).
5. Does Theorem 1 extend to multimarginal block-coordinate Sinkhorn with one Poincaré marginal per block? Di Marino's announced result extends to several marginals.

## Declarations

This is version 0.1.0-candidate, archived at https://doi.org/10.5281/zenodo.23929617 (all versions: https://doi.org/10.5281/zenodo.23929616); author Anonymous; AI-assisted and unrefereed; text and figures CC0-1.0, code MIT.

## References

Aubin-Frankowski, P.-C., Korba, A., & Léger, F. (2022). Mirror descent with relative smoothness in measure spaces, with application to Sinkhorn and EM. In *Advances in Neural Information Processing Systems 35* (pp. 17263--17275). <https://doi.org/10.52202/068431-1255>

Bobkov, S., & Ledoux, M. (1997). Poincaré's inequalities and Talagrand's concentration phenomenon for the exponential distribution. *Probability Theory and Related Fields, 107*(3), 383--400. <https://doi.org/10.1007/s004400050090>

Bobkov, S. G., Gentil, I., & Ledoux, M. (2001). Hypercontractivity of Hamilton--Jacobi equations. *Journal de Mathématiques Pures et Appliquées, 80*(7), 669--696. <https://doi.org/10.1016/S0021-7824(01)01208-9>

Carlier, G. (2022). On the linear convergence of the multimarginal Sinkhorn algorithm. *SIAM Journal on Optimization, 32*(2), 786--794. <https://doi.org/10.1137/21M1410634>

Chafaï, D., & Malrieu, F. (2010). On fine properties of mixtures with respect to concentration of measure and Sobolev type inequalities. *Annales de l'Institut Henri Poincaré, Probabilités et Statistiques, 46*(1), 72--96. <https://doi.org/10.1214/08-AIHP309>

Chen, L., & Xu, Z. (2026). *Accelerating Sinkhorn for entropy-regularized optimal transport* (arXiv:2605.30267). <https://arxiv.org/abs/2605.30267>

Chiarini, A., Conforti, G., Greco, G., & Tamanini, L. (2025). *A semiconcavity approach to stability of entropic plans and exponential convergence of Sinkhorn's algorithm* (arXiv:2412.09235, version 2). <https://arxiv.org/abs/2412.09235v2>

Chizat, L., Delalande, A., & Vaškevičius, T. (2026). Sharper exponential convergence rates for Sinkhorn's algorithm in continuous settings. *Mathematical Programming, 215*(1--2), 809--858. Published online 1 July 2025. <https://doi.org/10.1007/s10107-025-02242-z> (also arXiv:2407.01202, version 2)

Conforti, G., Durmus, A., & Greco, G. (2023). *Quantitative contraction rates for Sinkhorn's algorithm: beyond bounded costs and compact marginals* (arXiv:2304.04451, version 4, revised 2026). <https://arxiv.org/abs/2304.04451>

Cuturi, M. (2013). Sinkhorn distances: lightspeed computation of optimal transportation distances. *Advances in Neural Information Processing Systems 26*, 2292--2300 (arXiv:1306.0895). <https://arxiv.org/abs/1306.0895>

Delalande, A. (2022). Nearly tight convergence bounds for semi-discrete entropic optimal transport. *Proceedings of the 25th International Conference on Artificial Intelligence and Statistics*, PMLR 151, 1619--1642 (arXiv:2110.12678). <https://arxiv.org/abs/2110.12678>

Di Marino, S., & Gerolin, A. (2020). An optimal transport approach for the Schrödinger bridge problem and convergence of Sinkhorn algorithm. *Journal of Scientific Computing, 85*(2), Article 27. <https://doi.org/10.1007/s10915-020-01325-7>

Eckstein, S. (2025). Hilbert's projective metric for functions of bounded growth and exponential convergence of Sinkhorn's algorithm. *Probability Theory and Related Fields, 193*(1--2), 585--621. <https://doi.org/10.1007/s00440-025-01366-9>

Franklin, J., & Lorenz, J. (1989). On the scaling of multidimensional matrices. *Linear Algebra and its Applications, 114--115*, 717--735. <https://doi.org/10.1016/0024-3795(89)90490-4>

Ghosal, P., & Nutz, M. (2026). On the convergence rate of Sinkhorn's algorithm. *Mathematics of Operations Research, 51*(2), 1080--1096. <https://doi.org/10.1287/moor.2024.0427>

Gozlan, N. (2010). Poincaré inequalities and dimension free concentration of measure. *Annales de l'Institut Henri Poincaré, Probabilités et Statistiques, 46*(3), 708--739. <https://doi.org/10.1214/09-AIHP209>

Greco, G., Noble, M., Conforti, G., & Durmus, A. (2023). Non-asymptotic convergence bounds for Sinkhorn iterates and their gradients: a coupling approach. *Proceedings of Thirty Sixth Conference on Learning Theory*, PMLR 195, 716--746 (arXiv:2304.06549). <https://arxiv.org/abs/2304.06549>

Gromov, M., & Milman, V. D. (1983). A topological application of the isoperimetric inequality. *American Journal of Mathematics, 105*(4), 843. <https://doi.org/10.2307/2374298>

Grunewald, N., Otto, F., Villani, C., & Westdickenberg, M. G. (2009). A two-scale approach to logarithmic Sobolev inequalities and the hydrodynamic limit. *Annales de l'Institut Henri Poincaré, Probabilités et Statistiques, 45*(2), 302--351. <https://doi.org/10.1214/07-AIHP200>

Léger, F. (2021). A gradient descent perspective on Sinkhorn. *Applied Mathematics & Optimization, 84*(2), 1843--1855. <https://doi.org/10.1007/s00245-020-09697-w>

Maddison, C. J., Paulin, D., Teh, Y. W., & Doucet, A. (2021). Dual space preconditioning for gradient descent. *SIAM Journal on Optimization, 31*(1), 991--1016. <https://doi.org/10.1137/19M130858X>

Muckenhoupt, B. (1972). Hardy's inequality with weights. *Studia Mathematica, 44*(1), 31--38. <https://doi.org/10.4064/sm-44-1-31-38>

Otto, F., & Reznikoff, M. G. (2007). A new criterion for the logarithmic Sobolev inequality and two applications. *Journal of Functional Analysis, 243*(1), 121--157. <https://doi.org/10.1016/j.jfa.2006.10.002>

Peyré, G., & Cuturi, M. (2019). Computational optimal transport with applications to data sciences. *Foundations and Trends in Machine Learning, 11*(5--6), 355--607. <https://doi.org/10.1561/2200000073>

Rényi, A. (1959). On measures of dependence. *Acta Mathematica Academiae Scientiarum Hungaricae, 10*(3--4), 441--451. <https://doi.org/10.1007/BF02024507>

Sinkhorn, R. (1964). A relationship between arbitrary positive matrices and doubly stochastic matrices. *The Annals of Mathematical Statistics, 35*(2), 876--879. <https://doi.org/10.1214/aoms/1177703591>

Sinkhorn, R., & Knopp, P. (1967). Concerning nonnegative matrices and doubly stochastic matrices. *Pacific Journal of Mathematics, 21*(2), 343--348. <https://doi.org/10.2140/pjm.1967.21.343>

Srinivasan, A., & Slotine, J.-J. (2026). *Thermodynamic structure of the Sinkhorn flow* (arXiv:2510.12639, version 2). <https://arxiv.org/abs/2510.12639v2>

### Public announcements relevant to priority

Simons Institute for the Theory of Computing. (2021, November 10). *GMOS weekly seminar: The speed of Sinkhorn algorithm meets Poincaré inequality* (S. Di Marino) [Seminar abstract]. <https://simons.berkeley.edu/events/gmos-weekly-seminar-speed-sinkhorn-algorithm-meets-poincare-inequality>

The 2022 Pisa workshop announcement is known to us only through CDV's account (CDV, Section 1.1); we found no stable public URL for it.
