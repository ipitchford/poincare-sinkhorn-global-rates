---
title: "Curvature, relative entropy, and global Sinkhorn convergence"
date: "1 November 2026"
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
---

## Abstract

We prove explicit global convergence bounds for exact Sinkhorn iterations by converting a uniform curvature estimate, measured in the current marginal, into a relative-entropy error bound. If the source satisfies a Poincaré inequality with constant $C_P$ and the bounded cost has source gradient bounded by $L$, the dual gap contracts in each full iteration by at most $C_PL^2/(C_PL^2+\lambda^2)$, for every bounded initial potential and every target probability measure. Neither convexity of the source support nor a density assumption on the target is required. The argument uses finite-dimensional convex conjugacy on measurable partitions; it does not invert an infinite-dimensional Hessian. Applying the same transfer to an existing curvature theorem of Chizat, Delalande and Vaškevičius gives a global contraction $C/(C+\lambda)$ in their bounded-density setting, improving their global quadratic dependence on regularisation to linear dependence. One fixed Poincaré source has local contraction deficit asymptotic to $(8/5)\log(2)\lambda^2\log^2(1/\lambda)$, while a family with uniformly bounded Poincaré constants attains order $\lambda^2$. The polynomial exponent is therefore sharp. The literature already supplies reverse-entropy methods and, indirectly, qualitative polynomial convergence in overlapping smooth settings. Our claims concern the explicit bounds, their stated scope, the global bounded-density improvement, and the sharpness constructions.

## 1. Problem and contribution

Let $\mu$ and $\nu$ be probability measures on source and target spaces $X$ and $Y$. For $\lambda>0$, entropy-regularised optimal transport minimises

$$
\int c(x,y)\,d\gamma(x,y)+\lambda\operatorname{KL}(\gamma\Vert\mu\otimes\nu)
\quad\text{over }\gamma\in\Pi(\mu,\nu).
\tag{1}
$$

All logarithms in this paper are natural. One iteration means one source normalisation followed by one target normalisation. We study the dual objective gap, not the bias between (1) and unregularised transport, and not the arithmetic cost of evaluating a continuous integral.

Chizat, Delalande and Vaškevičius (2026, published online in 2025; hereafter CDV) ask in Remark 3.1 whether a Poincaré inequality on one marginal suffices for their convergence results. Their Theorem 1.2 gives a global contraction deficit of order $\lambda^2$ for a bounded positive source density and bilinear cost, and improves this to order $\lambda$ eventually. Their Proposition 5.1 gives a uniform curvature estimate that is stronger than the global rate extracted from it. We use that estimate directly through convex conjugacy.

Three results are established here.

1. A source Poincaré inequality yields the explicit global rate in Theorem 2, with no condition on the target and only a first derivative of the cost.
2. An abstract transfer from curvature to relative entropy gives a global linear-in-$\lambda$ rate in the CDV bounded-density setting, without an initial phase or a special initialisation.
3. A single fixed source proves that the power two cannot be improved throughout the Poincaré class, even for a bounded bilinear cost and two target atoms. A family with uniformly bounded Poincaré constants attains exact quadratic order.

The transfer uses established convex-duality and relative-smoothness ideas. Reverse relative-entropy comparisons for entropic plans appear in Aubin-Frankowski, Korba and Léger (2022), with a coefficient exponential in inverse regularisation. We replace that coefficient by explicit geometric quantities. Qualitative polynomial convergence in overlapping smooth settings also follows by combining the quadratic-linear transport inequality of Bobkov, Gentil and Ledoux (2001) with the generalised semiconcavity theorem of Chiarini, Conforti, Greco and Tamanini (2025). This overlap is not confined to compact supports. We therefore do not claim the first qualitative one-Poincaré polynomial convergence theorem. Section 8 and the accompanying novelty report separate these precedents from the conclusions proved here.

## 2. Definitions and assumptions

Initially, $X$ and $Y$ are standard Borel spaces, $\mu$ and $\nu$ are probability measures, and $c$ is a bounded measurable real-valued cost. For a bounded measurable potential $\psi:Y\to\mathbb R$, set

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

Thus $\gamma_\psi$ always has source marginal $\mu$; its target marginal is $r_\psi$, which need not equal $\nu$. Boundedness ensures that $r_\psi/\nu$ is bounded above and away from zero. Define the semidual objective and its gap by

$$
E(\psi)=\int_Y\psi\,d\nu-\Phi(\psi),\qquad
E^*=\sup_{\eta\in L^\infty(\nu)}E(\eta),\qquad
\Delta(\psi)=E^*-E(\psi).
\tag{3}
$$

The supremum is finite. For example, Jensen's inequality bounds $E(\psi)$ above by $\int c\,d(\mu\otimes\nu)$. Standard bounded-kernel Sinkhorn duality also supplies an optimal potential; the contraction proof below only needs the supremum in (3).

The exact full update is

$$
S\psi=\psi-\lambda\log\frac{dr_\psi}{d\nu}.
\tag{4}
$$

This is the usual pair of entropic transforms: the source potential is $-\lambda\log Z_\psi$, and transforming it back gives (4). Bounded potentials remain bounded. Constants in a potential do not change $q_x^\psi$, $r_\psi$, $E$, or $\Delta$.

Differentiation in bounded directions gives

$$
D\Phi(\psi)[h]=\int h\,dr_\psi,\qquad
D^2\Phi(\psi)[h,h]=\frac1\lambda\int_X
\operatorname{Var}_{q_x^\psi}(h)\,d\mu(x).
\tag{5}
$$

For the Poincaré application, let $X$ be a connected open subset of $\mathbb R^d$ carrying $\mu$. Assume

$$
\operatorname{Var}_\mu(f)\le C_P\int_X|\nabla f|^2\,d\mu
\tag{6}
$$

for bounded $C^1$ functions with bounded gradient. Assume $c(\cdot,y)$ is $C^1$ on $X$, the cost and its source gradient are jointly measurable, and

$$
\sup_{x\in X,\,y\in Y}|\nabla_xc(x,y)|\le L<\infty.
\tag{7}
$$

Essential suprema may replace suprema with the corresponding almost-everywhere differentiability and domination hypotheses. We use the classical $C^1$ formulation to keep those analytic domain issues explicit. No Lebesgue density is assumed for $\nu$. The source domain may be nonconvex or unbounded; the cost itself remains bounded. Boundary points of an interval with zero measure do not affect the statements.

## 3. From uniform curvature to a global entropy bound

**Theorem 1 (curvature transfer).** Suppose there is a constant $A\ge\lambda$ such that, for every bounded $\psi,h$,

$$
D^2\Phi(\psi)[h,h]\ge A^{-1}\operatorname{Var}_{r_\psi}(h).
\tag{8}
$$

Then

$$
\boxed{\Delta(\psi)\le A\operatorname{KL}(\nu\Vert r_\psi)}
\tag{9}
$$

for every bounded $\psi$. More generally, for any two bounded potentials $\eta,\psi$,

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

The marginal in (8) moves with the potential. The proof must preserve that dependence. Replacing it by $\nu$ before optimisation would be unjustified.

### 3.1. Finite partitions and a genuine inverse

Fix $\psi$. Partition $Y$ into finitely many measurable sets $B_1,\ldots,B_m$ of positive $\nu$-measure, discarding zero-measure sets. Put

$$
Q_i(x)=q_x^\psi(B_i),\qquad
F(u)=\lambda\int_X\log\left(\sum_{i=1}^m Q_i(x)e^{u_i/\lambda}\right)d\mu(x).
\tag{12}
$$

The function $F$ is the increment of $\Phi$ for a potential constant on each cell. It satisfies $F(0)=0$ and $F(u+a\mathbf1)=F(u)+a$. Boundedness of $c$ and $\psi$ gives a common $\epsilon>0$ such that $Q_i(x)\ge\epsilon$ for every cell and almost every $x$. One possible lower bound for cell $i$ is

$$
Q_i(x)\ge\nu(B_i)
\exp[-(\operatorname{osc}\psi+\operatorname{osc}c)/\lambda].
$$

Let $r(u)=\nabla F(u)$. Direct differentiation and (8) give

$$
\nabla^2F(u)\succeq A^{-1}
\left(\operatorname{diag}r(u)-r(u)r(u)^\top\right).
\tag{13}
$$

The Hessian is positive definite modulo constants. For every strictly positive probability vector $p$, $F(u)-p\cdot u$ is coercive modulo constants: using the gauge $\max_i u_i=0$,

$$
F(u)-p\cdot u\ge\lambda\log\epsilon+
(\min_i p_i)\operatorname{osc}u.
\tag{14}
$$

It therefore has a unique minimum modulo constants. Its first-order condition is $r(u)=p$. The gradient maps the quotient of $\mathbb R^m$ by constants onto the open probability simplex. The implicit function theorem supplies a smooth inverse there.

Let $G=F^*$ on that simplex. For a tangent vector $v$ with $\sum_i v_i=0$, (13) implies

$$
D^2G(p)[v,v]\le A\sum_{i=1}^m\frac{v_i^2}{p_i}.
\tag{15}
$$

To check the singular-matrix issue explicitly, let $H=\nabla^2F(u)$ at $r(u)=p$. Its inverse quadratic form on the quotient is

$$
\begin{aligned}
v^\top H^{-1}v
&=\sup_z\{2v\cdot z-z^\top Hz\}\\
&\le\sup_z\{2v\cdot z-A^{-1}\operatorname{Var}_p(z)\}
=A\sum_i v_i^2/p_i.
\end{aligned}
\tag{16}
$$

In this formula $H^{-1}$ means only the inverse quadratic form on tangent vectors, not an ordinary inverse of the singular full matrix. The final equality follows by centring $z$ under $p$ and completing the square. This proves (15).

### 3.2. Entropy integration and the infinite target

Write $r=r(0)$. Because $G(r)=0$ and $DG(r)$ is represented by $u=0$, integration of (15) along $r+t(p-r)$ yields

$$
\begin{aligned}
G(p)
&\le A\int_0^1(1-t)\sum_i
\frac{(p_i-r_i)^2}{r_i+t(p_i-r_i)}\,dt\\
&=A\sum_i p_i\log(p_i/r_i).
\end{aligned}
\tag{17}
$$

The last line is the integral Hessian formula for $p\mapsto\sum_i p_i\log p_i$. Hence

$$
\sup_u\{p\cdot u-F(u)\}\le A\operatorname{KL}(p\Vert r).
\tag{18}
$$

Take $p_i=\nu(B_i)$. The left side bounds every semidual improvement obtained by a partition-constant increment. Relative entropy decreases under the map sending a point to its cell, so the right side is at most $A\operatorname{KL}(\nu\Vert r_\psi)$.

Every bounded measurable increment is uniformly approximable by finite simple functions. The log-partition functional $\Phi$ is 1-Lipschitz in the supremum norm: if two potentials differ by at most $\varepsilon$, their log-partitions differ by at most $\varepsilon$. Thus $E$ is 2-Lipschitz. Passing to uniform limits and then taking the supremum of all improvements proves (9).

For (10), use $p_i=r_\eta(B_i)$ in the same argument. At increment $h=\eta-\psi$, the expression being bounded is

$$
\begin{aligned}
\int h\,dr_\eta-\Phi(\psi+h)+\Phi(\psi)
&=\lambda\operatorname{KL}(\gamma_\eta\Vert\gamma_\psi).
\end{aligned}
\tag{19}
$$

Uniform simple approximation and data processing finish the proof. All Hessian inverses used in this argument were finite-dimensional. In particular, no unproved assertion about the image of an infinite-dimensional gradient map is needed.

### 3.3. The Sinkhorn ascent identity

For completeness, write the two-potential dual functional as

$$
J(\phi,\psi)=\int\phi\,d\mu+\int\psi\,d\nu+
\lambda\left(1-\int e^{(\phi+\psi-c)/\lambda}\,d\mu\,d\nu\right).
$$

At $\phi_\psi=-\lambda\log Z_\psi$, $J(\phi_\psi,\psi)=E(\psi)$ and the exponential density has total mass one. Replacing $\psi$ by $S\psi$ makes its target marginal equal to $\nu$, so its total mass remains one. Equation (4) therefore gives

$$
J(\phi_\psi,S\psi)-J(\phi_\psi,\psi)
=\lambda\operatorname{KL}(\nu\Vert r_\psi).
$$

Optimising the source potential again can only increase this value. Consequently,

$$
E(S\psi)-E(\psi)\ge\lambda\operatorname{KL}(\nu\Vert r_\psi).
\tag{20}
$$

Combining (9) with (20) proves (11). This argument uses no assumption that the initial potential has oscillation bounded by the cost. It applies from iteration zero. $\square$

## 4. A Poincaré inequality on one marginal

**Theorem 2 (explicit Poincaré rate).** Under (6)--(7) and the bounded-cost assumptions of Section 2, put $B=C_PL^2$. For every bounded $\psi$,

$$
\Delta(\psi)\le\left(\lambda+\frac B\lambda\right)
\operatorname{KL}(\nu\Vert r_\psi).
\tag{21}
$$

For any bounded initial potential and every integer $t\ge0$,

$$
\boxed{\Delta(\psi_t)\le
\left(\frac B{B+\lambda^2}\right)^t\Delta(\psi_0),\qquad
\psi_{t+1}=S\psi_t.}
\tag{22}
$$

When $B=0$, the statement means that one full iteration reaches zero gap; the iteration-zero bound is equality.

The corresponding pairwise comparison is

$$
\operatorname{KL}(\gamma_\eta\Vert\gamma_\psi)
\le\left(1+\frac B{\lambda^2}\right)
\operatorname{KL}(r_\eta\Vert r_\psi).
\tag{23}
$$

*Proof.* For bounded $h$, define $m(x)=\mathbb E_{q_x^\psi}h$. Differentiating the normalised kernel gives the covariance identity

$$
\nabla m(x)=-\lambda^{-1}
\operatorname{Cov}_{q_x^\psi}(h,\nabla_xc).
\tag{24}
$$

The bounded gradient supplies domination, and $m$ belongs to the class in (6). For a scalar-vector covariance, Cauchy--Schwarz gives

$$
\begin{aligned}
|\nabla m(x)|^2
&\le\lambda^{-2}\operatorname{Var}_{q_x^\psi}(h)
\mathbb E_{q_x^\psi}|\nabla_xc-\mathbb E_{q_x^\psi}\nabla_xc|^2\\
&\le\lambda^{-2}L^2\operatorname{Var}_{q_x^\psi}(h).
\end{aligned}
\tag{25}
$$

Let $D=\int\operatorname{Var}_{q_x^\psi}(h)\,d\mu$. The Poincaré inequality and total variance imply

$$
\operatorname{Var}_{r_\psi}(h)
=D+\operatorname{Var}_\mu(m)
\le\left(1+\frac B{\lambda^2}\right)D.
\tag{26}
$$

By (5), condition (8) holds with $A=\lambda+B/\lambda$. Theorem 1 gives all three conclusions. $\square$

### 4.1. Geometry, invariance, and stopping

The bound has no explicit target-size or dimension factor. The constant $C_P$ may nevertheless depend on source geometry and dimension, and $L$ may depend on target extent. Thus the theorem does not assert dimension-independent computational complexity.

An additive source-only cost cancels from (24). More generally, the same proof works if a source vector field $b(x)$ satisfies $|\nabla_xc(x,y)-b(x)|\le L$, with differentiation justified as before. For $c(x,y)=|x-y|^2/2$ on bounded supports, choose $b(x)=x-y_0$. Then $L=\sup_y|y-y_0|$, an enclosing radius of the target. This is a refinement of the constant, not a theorem for an arbitrary unbounded quadratic cost.

For $B>0$, (22) implies

$$
\Delta(\psi_t)\le
\exp\left[-\frac{t\lambda^2}{B+\lambda^2}\right]\Delta(\psi_0).
\tag{27}
$$

To reduce the gap by a factor $\varepsilon\in(0,1)$, it suffices to use $(1+B/\lambda^2)\log(1/\varepsilon)$ iterations, rounded upward. Equation (21) also gives a stopping certificate from the current marginal relative entropy, if that population integral and a valid $C_P$ are available.

If $\psi^*$ is an optimal potential, then

$$
\Delta(\psi)=\lambda\operatorname{KL}(\gamma_{\psi^*}\Vert\gamma_\psi).
\tag{28}
$$

This follows directly from (19) and $r_{\psi^*}=\nu$. With total variation defined as the supremum over measurable events, Pinsker's inequality gives

$$
\|\gamma_{\psi^*}-\gamma_\psi\|_{\mathrm{TV}}
\le\sqrt{\Delta(\psi)/(2\lambda)}.
$$

The iterated plan $\gamma_\psi$ is source-normalised but generally has a target-marginal error; (28) does not turn it into an exactly feasible coupling.

### 4.2. A deterministic approximation statement

Suppose an implementation produces bounded potentials satisfying

$$
\inf_{a\in\mathbb R}\|\widehat\psi_{t+1}-S\widehat\psi_t-a\|_\infty\le\eta_t.
$$

The 2-Lipschitz property of $E$ and its invariance under constants give, for $q=B/(B+\lambda^2)$,

$$
\Delta(\widehat\psi_t)
\le q^t\Delta(\widehat\psi_0)+2\sum_{j=0}^{t-1}q^{t-1-j}\eta_j.
\tag{29}
$$

If $\eta_j\le\eta$, the error floor is at most $2\eta(1+B/\lambda^2)$. This is a conditional approximation guarantee. It does not supply a quadrature or sampling error estimate by itself. In particular, replacing the source by a finite empirical measure does not preserve the Euclidean-gradient Poincaré inequality (6).

## 5. A global linear regularisation rate for bounded densities

**Corollary 3 (bounded-density improvement).** Let $X,Y\subset\mathbb R^d$ be compact, with $X$ convex, let $c(x,y)=-\langle x,y\rangle$, and suppose $\mu$ has Lebesgue density $f$ on $X$ satisfying $0<m\le f\le M<\infty$. Set

$$
\kappa=M/m,\qquad c_\infty=\operatorname{osc}_{X\times Y}c,
\qquad C=e\kappa c_\infty.
$$

Then every bounded initial potential satisfies

$$
\boxed{\Delta(\psi_t)\le
\left(\frac C{C+\lambda}\right)^t\Delta(\psi_0).}
\tag{30}
$$

*Proof.* CDV Proposition 5.1 gives (8) with $A=C+\lambda$. Applying Theorem 1 proves (30). $\square$

This corollary credits the entire curvature estimate to CDV. The additional step is its conversion into a global marginal-entropy error bound. Their Theorem 1.2 has a global coefficient of the form $176(1+\kappa c_\infty/\lambda+c_\infty^2/\lambda^2)$ and an eventual coefficient $176(1+\kappa c_\infty/\lambda)$. Equation (30) obtains a linear dependence from the initial iteration. No lower bound on target atom masses is required.

For any other setting covered by CDV Proposition 5.1, Theorem 1 likewise gives the rate $C/(C+\lambda)$ with their stated value of $C$ and all their hypotheses retained. We state the bilinear case separately because it avoids importing additional semiconcavity notation and covers the standard quadratic-cost problem after removal of marginal-only terms.

## 6. A fixed source that prevents every smaller exponent

The quadratic power in Theorem 2 cannot be replaced by $\lambda^\beta$ for any $\beta<2$, even if the positive constant is allowed to depend on a fixed source and cost.

**Theorem 4 (polynomial exponent obstruction).** Take $X=(-1,1)$, $\nu=(\delta_{-1}+\delta_1)/2$, $c(x,y)=-xy$, and the fixed source density

$$
w(x)=\frac25|x|\bigl[\log(e/|x|)\bigr]^2\quad(0<|x|<1),
\qquad w(0)=0.
\tag{31}
$$

This source satisfies (6) with $C_P\le3/4$, and $L=1$. Its local dual-gap contraction factor is $(1-D_\lambda)^2$, where

$$
D_\lambda=\int_{-1}^1\operatorname{sech}^2(x/\lambda)w(x)\,dx,
$$

and

$$
1-(1-D_\lambda)^2
\sim\frac85\log2\,\lambda^2\log^2(1/\lambda)
\quad\text{as }\lambda\downarrow0.
\tag{32}
$$

Consequently, for this one fixed problem, there are no $\beta<2$, $K>0$ and $\lambda_0>0$ such that

$$
\Delta(S\psi)\le(1-K\lambda^\beta)\Delta(\psi)
$$

holds for every bounded $\psi$ and every $0<\lambda<\lambda_0$.

### 6.1. Normalisation and an explicit Poincaré certificate

The density is normalised because $\int_0^1x(1-\log x)^2\,dx=5/4$. For $0<x<1$,

$$
R(x):=\int_0^x\frac{dt}{w(t)}=
\frac{5}{2\log(e/x)}.
$$

Weighted Cauchy--Schwarz gives

$$
(f(x)-f(0))^2\le R(x)\int_0^x f'(t)^2w(t)\,dt.
$$

Integrating and exchanging the nonnegative integrals yields

$$
\begin{aligned}
\int_0^1(f-f(0))^2w
&\le\int_0^1 f'(t)^2w(t)
\left(\int_t^1 R(x)w(x)\,dx\right)dt\\
&\le\frac34\int_0^1 f'(t)^2w(t)\,dt,
\end{aligned}
$$

since $R(x)w(x)=x\log(e/x)$ and its integral on $(0,1)$ is $3/4$. The same bound holds on the negative half. Finally,
$\operatorname{Var}_\mu(f)\le\int(f-f(0))^2\,d\mu$ proves the claimed constant. The vanishing density at zero therefore does not destroy the Poincaré inequality.

### 6.2. An exact scalar Sinkhorn map

Represent potentials modulo constants as $\psi_\theta(y)=\theta y$. Symmetry makes $\theta=0$ optimal. Let

$$
m_\lambda(\theta)=\mathbb E_\mu\tanh((X+\theta)/\lambda).
$$

The target probabilities are $(1\pm m_\lambda(\theta))/2$, so (4) gives

$$
T_\lambda(\theta)=\theta-\lambda\operatorname{arctanh}m_\lambda(\theta).
\tag{33}
$$

Differentiation at zero yields

$$
m_\lambda'(0)=D_\lambda/\lambda,\qquad
T_\lambda'(0)=1-D_\lambda.
\tag{34}
$$

The exact gap is

$$
\Delta_\lambda(\theta)=\lambda\mathbb E_\mu
\left[\log\cosh((X+\theta)/\lambda)-\log\cosh(X/\lambda)\right]
=\frac{D_\lambda}{2\lambda}\theta^2+o(\theta^2).
\tag{35}
$$

Hence

$$
\lim_{\theta\to0,\,\theta\ne0}
\frac{\Delta_\lambda(T_\lambda(\theta))}{\Delta_\lambda(\theta)}
=(1-D_\lambda)^2.
\tag{36}
$$

This factor also describes actual nonzero trajectories near their limit. For $\theta>0$, pairing $x$ with $-x$ shows $0<m_\lambda(\theta)<\tanh(\theta/\lambda)$, since the paired average is

$$
\frac{\sinh(2\theta/\lambda)}
{\cosh(2\theta/\lambda)+\cosh(2x/\lambda)}.
$$

Thus $0<T_\lambda(\theta)<\theta$. Iterates converge to zero, and oddness handles negative initial values. In particular, the slow factor is not an artefact of starting exactly at the optimum.

### 6.3. The small-regularisation limit

Substitution $x=\lambda u$ gives

$$
D_\lambda=\frac45\lambda^2
\int_0^{1/\lambda}u\bigl[\log(e/\lambda)-\log u\bigr]^2
\operatorname{sech}^2u\,du.
\tag{37}
$$

After division by $\lambda^2\log^2(1/\lambda)$, the integrand converges to $(4/5)u\operatorname{sech}^2u$. For $\lambda\le e^{-1}$ it is dominated by a constant times
$u(2+|\log u|)^2\operatorname{sech}^2u$, an integrable function. Since

$$
\int_0^\infty u\operatorname{sech}^2u\,du=\log2,
$$

dominated convergence proves
$D_\lambda\sim(4/5)\log2\,\lambda^2\log^2(1/\lambda)$, and hence (32). Dividing (32) by $\lambda^\beta$ gives zero for every $\beta<2$. Taking the limit (36) contradicts any proposed uniform one-step bound with deficit $K\lambda^\beta$. $\square$

Theorem 4 proves optimality of the polynomial exponent. It does not assert that the constant in (22) is sharp, or that this fixed source has an exact $\Theta(\lambda^2)$ deficit. The logarithmic factor is stated explicitly. The density in (31) has no positive lower bound, so the example does not contradict Corollary 3.

### 6.4. Exact quadratic order uniformly over a Poincaré class

The logarithmic factor can be removed when the source is allowed to vary with regularisation, while its Poincaré constant stays uniformly bounded.

**Proposition 5 (uniform-class sharpness).** Keep the two-atom target and cost $c(x,y)=-xy$. For $0<\lambda\le1/2$, let

$$
\rho_\lambda(x)=\min\{\lambda\cosh(x/\lambda),1\},\qquad
w_\lambda(x)=\frac{\rho_\lambda(x)}{2I_\lambda},\quad |x|<1,
$$

where

$$
I_\lambda=1-\lambda\operatorname{arcosh}(1/\lambda)
+\lambda\sqrt{1-\lambda^2}.
\tag{38}
$$

Every $w_\lambda$ has $C_P\le1+\pi/2$. Its local dual-gap contraction deficit satisfies

$$
1-(1-D_\lambda)^2\sim\pi\lambda^2.
\tag{39}
$$

*Proof.* The transition point is $x_*=\lambda\operatorname{arcosh}(1/\lambda)<1$. Integrating on either side gives $\int_0^1\rho_\lambda=I_\lambda\le1$. The same anchored Cauchy--Schwarz proof used above gives

$$
C_P\le\left(\int_0^1\rho_\lambda\right)
\left(\int_0^1\rho_\lambda^{-1}\right)
\le1+\pi/2,
$$

because

$$
\int_0^1\rho_\lambda(x)^{-1}\,dx
=\arccos\lambda+1-\lambda\operatorname{arcosh}(1/\lambda)
\le1+\pi/2.
$$

The exact conditional-variance integral is

$$
D_\lambda=
\frac{\lambda^2\arccos\lambda+
\lambda\bigl[\tanh(1/\lambda)-\sqrt{1-\lambda^2}\bigr]}{I_\lambda}.
\tag{40}
$$

Here $I_\lambda\to1$, $\arccos\lambda\to\pi/2$, and the second numerator term is $O(\lambda^3)$. Hence $D_\lambda\sim(\pi/2)\lambda^2$. The scalar calculation (33)--(36) proves (39). $\square$

Theorem 2 and Proposition 5 establish the quadratic order of the worst local deficit over this uniformly bounded-Poincaré class. The sources in Proposition 5 depend on $\lambda$; Theorem 4 supplies the separate obstruction for a single fixed source. Neither result identifies the sharp universal multiplicative constant.

For each fixed $\lambda$, the density in Proposition 5 is positive and bounded, with ratio $\kappa_\lambda=1/\lambda$. Thus Corollary 3 gives $C=2e/\lambda$ and a contraction deficit of order $\lambda^2$ for this family, consistently with (39). Its linear-in-$\lambda$ conclusion applies when the density ratio is held fixed.

## 7. Numerical checks and reproducibility

The analytic arguments above do not use numerical evidence. The accompanying Python code checks their identities and asymptotics independently by quadrature. It treats the uniform source on $[-1,1]$ and the density (31), with a two-atom target, over 180 combinations of regularisation and nonzero potential. It also follows ten actual nonzero Sinkhorn trajectories. Stable paired formulas for the gap, substitutions $x=\lambda u$, and normalisation by the potential or its square avoid cancellation at small gaps.

For the uniform source, a separate exact check is available:

$$
D_\lambda=\lambda\tanh(1/\lambda),\qquad
\text{local gap factor}=(1-\lambda\tanh(1/\lambda))^2.
\tag{41}
$$

For the cusp source, quadrature gives $D_{0.1}=0.0647035371183$, $D_{0.01}=0.00179454249623$, and $D_{0.001}=0.0000353004893444$. These numbers are floating-point checks, not interval certificates. The source normalisation agreed with one to within $3\times10^{-16}$ in the recorded run. The checked global error-bound ratios remained below one; the largest were approximately $0.79150$ for the uniform source and $0.52720$ for the cusp source.

The initial numerical implementation exposed a cancellation-sensitive tiny-gap discrepancy. The final script evaluates scaled integrands before quadrature. The verification record preserves that correction and the final environment, tolerances, counts and output. A complementary test with a Gaussian source, a bounded sine cost and five target atoms checks the noncompact-source and multiple-target setting; its exact scope is recorded separately. Numerical discretisation of a source is only a check of the continuous formulas, not a claim that the discrete measure satisfies (6).

![Local contraction deficits. The uniform source has a linear deficit; the fixed cusp approaches a quadratic power with a logarithmic correction. The figure evaluates the exact scalar derivative, rather than fitting a rate from a short trajectory.](figures/sharpness_plot.png){width=92%}

## 8. Prior work, scope, and unresolved questions

The ingredients and the new quantitative deductions should be distinguished.

* **Existing geometry.** Convex-conjugate Hessian inversion and relative smoothness are established optimisation tools; see Maddison et al. (2021). Aubin-Frankowski et al. (2022, Proposition 7) already compare joint and marginal relative entropy for entropic plans. Xu and Chen (2026) use dual-relative geometry for finite Sinkhorn problems. Theorem 1 gives a self-contained finite-partition implementation suitable for the constants used here; it does not claim to invent that geometry.
* **Existing qualitative overlap.** The transport-inequality route through Bobkov et al. (2001) and Chiarini et al. (2025) overlaps sufficiently smooth one-Poincaré settings, including some unbounded supports. For globally Lipschitz costs with bounded Hessians, a quadratic semiconcavity estimate at short distances and a Lipschitz estimate at long distances match the quadratic-linear transport cost. This is a deduction from those papers, not a theorem located there verbatim. It does not provide the direct $C_PL^2$ estimate under only the first-source-derivative hypothesis used here.
* **Existing curvature, stronger consequence.** Corollary 3 uses CDV's curvature estimate unchanged. The global $C/(C+\lambda)$ consequence is the claimed strengthening. It avoids the interpolation loss discussed in their Remark 5.2.
* **Conditional expectations.** Srinivasan and Slotine (2026) formulate related conditional-expectation and spectral-gap structures for a continuous-time flow. Equation (26) supplies a uniform geometric bound, while Theorem 1 turns it into a discrete-iteration guarantee.

The search was conducted on 1 November 2026 using web search, arXiv full text, publisher pages, proceedings and citation-linked sources. No earlier statement of the exact global bounds or the explicit sharpness constructions was found within that recorded boundary. This is not a proof of historical priority. CDV reports an earlier two-marginal Poincaré result announced by Simone Di Marino. A public Simons Institute seminar abstract dated 10 November 2021 also announces polynomial rates and counterexamples, without giving the formulas. Those proof details were not available in the inspected record and remain a material priority uncertainty. The separate novelty report records queries, inspected formulas, indirect overlaps and unavailable evidence.

The theorem leaves several substantive questions open. The sharp multiplicative constant and exact worst fixed-source rate are not established. Unbounded costs require additional domain and integrability arguments. A finite empirical source does not inherit the required gradient inequality; an appropriate graph inequality would require a separate discrete proof. Finally, a computable stopping rule needs certified integration error and a defensible source Poincaré constant. These limitations do not alter the exact bounded-cost statements proved above.

## References

### Sinkhorn convergence and relative geometry

Aubin-Frankowski, P.-C., Korba, A., & Léger, F. (2022). Mirror descent with relative smoothness in measure spaces, with application to Sinkhorn and EM. *Advances in Neural Information Processing Systems, 35*. https://papers.nips.cc/paper_files/paper/2022/file/6e3daaeca6be8579573f69082b2dd58b-Paper-Conference.pdf

Chiarini, A., Conforti, G., Greco, G., & Tamanini, L. (2025). *A semiconcavity approach to stability of entropic plans and exponential convergence of Sinkhorn's algorithm* (arXiv:2412.09235, version 2). https://arxiv.org/abs/2412.09235v2

Chizat, L., Delalande, A., & Vaškevičius, T. (2026). Sharper exponential convergence rates for Sinkhorn's algorithm in continuous settings. *Mathematical Programming, 215*, 809--858. Advance online publication 1 July 2025. https://doi.org/10.1007/s10107-025-02242-z

Maddison, C. J., Paulin, D., Teh, Y. W., & Doucet, A. (2021). Dual space preconditioning for gradient descent. *SIAM Journal on Optimization, 31*(1), 991--1016. https://doi.org/10.1137/19M130858X

Srinivasan, A., & Slotine, J.-J. (2026). *Thermodynamic structure of the Sinkhorn flow* (arXiv:2510.12639, version 2). https://arxiv.org/abs/2510.12639v2

Xu, Z., & Chen, L. (2026). *Accelerating Sinkhorn for entropy-regularized optimal transport* (arXiv:2605.30267, version 1). https://arxiv.org/abs/2605.30267v1

### Functional inequalities

Bobkov, S. G., Gentil, I., & Ledoux, M. (2001). Hypercontractivity of Hamilton--Jacobi equations. *Journal de Mathématiques Pures et Appliquées, 80*(7), 669--696. https://doi.org/10.1016/S0021-7824(01)01208-9

### Public announcement relevant to priority

Simons Institute for the Theory of Computing. (2021, November 10). *GMOS weekly seminar: The speed of Sinkhorn algorithm meets Poincaré inequality* [Seminar abstract]. https://simons.berkeley.edu/events/gmos-weekly-seminar-speed-sinkhorn-algorithm-meets-poincare-inequality
