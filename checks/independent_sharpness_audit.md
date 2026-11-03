# Independent spectral and sharpness audit

Date: 1 November 2026. This note supplies independent mathematical derivations; it is not an external peer review or a claim of priority.

## 1. Conditional variance argument

Let \(\pi(dx,dy)=\mu(dx)p_x(dy)\), where
\[
p_x(dy)=\frac{\exp((v(y)-c(x,y))/\lambda)\nu(dy)}{\int\exp((v(z)-c(x,z))/\lambda)\nu(dz)}.
\]
Write \(q\) for its actual second marginal. Suppose \(\mu\) has Poincaré constant \(C_P\), and \(c(\cdot,y)\) is uniformly \(L\)-Lipschitz. For bounded \(h\), put \(g(x)=E[h(Y)\mid X=x]\). The differentiation identity is
\[
\nabla g=-\lambda^{-1}\operatorname{Cov}(h(Y),\nabla_xc(x,Y)\mid X=x).
\]
Cauchy–Schwarz gives \(\|\nabla g\|^2\le L^2\lambda^{-2}\operatorname{Var}(h(Y)\mid X=x)\). Hence, with \(A=C_PL^2/\lambda^2\),
\[
V:=\operatorname{Var}_\mu(g)\le A D,\qquad
D:=E_\mu\operatorname{Var}(h(Y)\mid X),\qquad
\operatorname{Var}_q(h)=V+D.
\]
Consequently
\[
\operatorname{Var}_\mu(E[h\mid X])\le\frac{A}{1+A}\operatorname{Var}_q(h).
\]
The use of \(q\), rather than \(\nu\), is essential away from the optimal potential. At the optimum, \(q=\nu\), and the reversible two-step conditional expectation operator on \(L^2(\nu)\) has spectral gap at least \(1/(1+A)\).

## 2. Exact two-point Sinkhorn dynamics

Let \(\mu\) be any symmetric nondegenerate law on \([-1,1]\), \(\nu=(\delta_{-1}+\delta_1)/2\), and \(c(x,y)=-xy\). Normalise the target potential as \(v_t(y)=ty\); adding a constant changes neither coupling nor semidual gap. Conditional probabilities satisfy
\[
P_t(Y=y\mid X=x)=\frac{e^{y(x+t)/\lambda}}{2\cosh((x+t)/\lambda)}.
\]
Put \(m_\lambda(t)=E_\mu\tanh((X+t)/\lambda)\). A full Sinkhorn update, with the source potential eliminated, is exactly
\[
T_\lambda(t)=t-\lambda\operatorname{arctanh}(m_\lambda(t)).
\]
Symmetry gives \(m_\lambda(0)=0\). Differentiation yields
\[
T_\lambda'(0)=1-E_\mu\operatorname{sech}^2(X/\lambda)
=E_\mu\tanh^2(X/\lambda)=:r_\lambda.
\]
More generally, writing \(Z=\tanh((X+t)/\lambda)\),
\[
T_\lambda'(t)=1-\frac{E(1-Z^2)}{1-(EZ)^2}
=\frac{\operatorname{Var}(Z)}{1-(EZ)^2}\in(0,1).
\]
Thus for every \(t>0\), \(0<T_\lambda(t)<t\): positivity follows from monotonicity and \(T_\lambda(0)=0\), and the upper inequality follows from \(m_\lambda(t)>0\). Iterates started at any \(t_0>0\) decrease to zero; continuity and the strict inequality exclude a positive limit. Negative initial conditions follow by oddness.

The semidual objective, up to an irrelevant constant, is
\[
F_\lambda(t)=-\lambda E_\mu\log\cosh((X+t)/\lambda).
\]
Its unique maximiser is zero. The gap obeys
\[
\Delta_\lambda(t):=F_\lambda(0)-F_\lambda(t)
=\frac{1-r_\lambda}{2\lambda}t^2+o(t^2).
\]
For every nonzero initial value, no iterate hits zero and
\[
\lim_{n\to\infty}\frac{\Delta_\lambda(t_{n+1})}{\Delta_\lambda(t_n)}=r_\lambda^2.
\]
This is the actual asymptotic gap ratio of Sinkhorn, not merely an eigenvalue of a related operator.

## 3. Power-cusp family

For \(0<a<1\), let
\[
\mu_a(dx)=\frac{a+1}{2}|x|^a\mathbf 1_{[-1,1]}(x)\,dx.
\]
Writing \(w(x)=(a+1)x^a/2\) on the positive half interval, weighted Cauchy–Schwarz gives
\[
|f(x)-f(0)|^2\le\left(\int_0^x f'(s)^2w(s)\,ds\right)
\left(\int_0^x\frac{ds}{w(s)}\right).
\]
After integration, Tonelli's theorem bounds the positive-half contribution by
\[
\int_0^1f'(s)^2w(s)\left(\int_s^1\frac{x}{1-a}\,dx\right)ds
\le\frac1{2(1-a)}\int_0^1f'(s)^2w(s)\,ds.
\]
Apply the same bound on the negative half interval and use \(\operatorname{Var}_\mu f\le E(f-f(0))^2\). Therefore
\[
C_P(\mu_a)\le\frac1{2(1-a)}.
\]
The explicit spectral gap is
\[
1-r_\lambda=(a+1)\lambda^{a+1}\int_0^{1/\lambda}u^a\operatorname{sech}^2u\,du
\sim(a+1)\lambda^{a+1}\int_0^\infty u^a\operatorname{sech}^2u\,du.
\]
The integral is finite and strictly positive. In particular the exponent can approach two while every fixed member of the family satisfies a Poincaré inequality.

## 4. A single fixed law rules out every exponent below two

Define
\[
\mu_\log(dx)=\frac25 |x|\log^2\!\left(\frac e{|x|}\right)\mathbf 1_{0<|x|\le1}\,dx.
\]
Set its density to zero at \(x=0\). Normalisation follows from
\[
\int_0^1x(1-\log x)^2\,dx=\frac54.
\]
On the positive half interval, \(w(x)=\tfrac25x\log^2(e/x)\), and
\[
\int_0^x\frac{ds}{w(s)}=\frac5{2\log(e/x)}.
\]
The preceding anchored argument yields
\[
C_P(\mu_\log)\le\int_0^1x\log(e/x)\,dx=\frac34.
\]
For \(L_\lambda=\log(e/\lambda)\),
\[
1-r_\lambda=\frac45\lambda^2\int_0^{1/\lambda}u(L_\lambda-\log u)^2\operatorname{sech}^2u\,du.
\]
Let \(A_j=\int_0^\infty u(\log u)^j\operatorname{sech}^2u\,du\), for \(j=0,1,2\). All three integrals are absolutely convergent. The missing tail after expansion is exponentially small, so
\[
1-r_\lambda=\frac45\lambda^2(A_0L_\lambda^2-2A_1L_\lambda+A_2)+o(\lambda^2).
\]
Integration by parts, or the antiderivative \(u\tanh u-\log\cosh u\), gives \(A_0=\log2\). It follows that
\[
1-r_\lambda\sim\frac45\log2\,\lambda^2\log^2(1/\lambda),\qquad
1-r_\lambda^2\sim\frac85\log2\,\lambda^2\log^2(1/\lambda).
\]
For every \(\beta<2\), \((1-r_\lambda^2)/\lambda^\beta\to0\). Thus this one fixed source, this one fixed two-point target, and this one fixed cost preclude any universally valid per-iteration dual-gap contraction whose small-\(\lambda\) improvement is bounded below by a positive constant times \(\lambda^\beta\), even when that constant is allowed to depend on these fixed data. The conclusion concerns polynomial exponent: the example leaves a logarithmic gap from an \(\lambda^2\) guarantee.

## 5. Material qualifications

* These sources are absolutely continuous but vanish at the origin; they do not have log-density bounded below. The logarithmic density is continuous at the origin, but its derivative there is unbounded.
* The target is atomic and does not satisfy a Euclidean-gradient Poincaré inequality. That is deliberate: it tests the one-marginal hypothesis.
* Starting the symmetric example at the zero target potential starts at the optimum. A lower-bound experiment must start at any nonzero potential, for example \(v_0(y)=y\), or use the local perturbation argument explicitly.
* No finite set of numerical iterates is required for the asymptotic ratio: the exact scalar map establishes it analytically.
* No claim is made here that the main global convergence theorem follows solely from the conditional-variance estimate. Any global proof must handle the moving second marginal \(q_v\).

## 6. Initial prior-art lead, superseded by the final novelty report

The initial search located Srinivasan and Slotine (2025), *Contraction and entropy production in continuous-time Sinkhorn dynamics*, arXiv:2510.12639, which studies conditional expectation coercivity and positive spectral gaps for the continuous-time Sinkhorn flow. At that stage this auditor had inspected its abstract, not its full text. The paper was subsequently retitled *Thermodynamic structure of the Sinkhorn flow*, version 2 (2026), https://arxiv.org/abs/2510.12639v2. Other members of the research team inspected the later full text. Their final comparison is recorded in `../prior_art/novelty_report.md` and supersedes this historical search warning. This note does not claim that the sharpness auditor independently repeated that full-text inspection.

## 7. Exact quadratic order uniformly over a bounded Poincaré class

The following strengthens the polynomial-exponent statement when the source may vary with the regularisation parameter. It is a separate quantifier statement from the fixed-source logarithmic example.

For \(0<\lambda\le1/2\), define on \([-1,1]\)
\[
\rho_\lambda(x)=\min\{\lambda\cosh(x/\lambda),1\},\qquad
\mu_\lambda(dx)=\frac{\rho_\lambda(x)}{2I_\lambda}\,dx,
\]
where
\[
I_\lambda=1-\lambda\operatorname{arcosh}(1/\lambda)
             +\lambda\sqrt{1-\lambda^2}.
\]
Then all these sources satisfy a Poincaré inequality with the common bound
\[
C_P(\mu_\lambda)\le 1+\frac\pi2.
\]
For the same fixed two-point target and cost as above, the exact scalar spectral gap is
\[
q_\lambda:=1-r_\lambda
=\frac{\lambda^2\arccos\lambda
+\lambda\{\tanh(1/\lambda)-\sqrt{1-\lambda^2}\}}{I_\lambda},
\]
and therefore
\[
q_\lambda\sim\frac\pi2\lambda^2,\qquad
1-r_\lambda^2\sim\pi\lambda^2.
\]

**Proof.** Put \(u_\lambda=\operatorname{arcosh}(1/\lambda)\) and \(a_\lambda=\lambda u_\lambda\). For the given range, \(0<a_\lambda<1\), and the two pieces of \(\rho_\lambda\) meet at \(a_\lambda\). Direct integration on the positive half interval gives
\[
\int_0^1\rho_\lambda(x)dx
=\lambda^2\sinh u_\lambda+1-a_\lambda=I_\lambda.
\]
The reciprocal-weight integral is
\[
\int_0^1\frac{dx}{\rho_\lambda(x)}
=\int_0^{u_\lambda}\operatorname{sech}u\,du+1-a_\lambda
=\arccos\lambda+1-a_\lambda\le1+\frac\pi2.
\]
The anchored Cauchy–Schwarz proof in Section 3, applied to the unnormalised weight \(\rho_\lambda\), bounds its Poincaré constant by
\[
\int_0^1\rho_\lambda(x)\left(\int_0^x\rho_\lambda(s)^{-1}ds\right)dx
\le\left(\int_0^1\rho_\lambda\right)
\left(\int_0^1\rho_\lambda^{-1}\right)\le1+\frac\pi2,
\]
because \(\rho_\lambda\le1\). The normalising constant cancels from the Poincaré ratio. Splitting the integral for \(q_\lambda=I_\lambda^{-1}\int_0^1\rho_\lambda(x)\operatorname{sech}^2(x/\lambda)dx\) at \(a_\lambda\) gives
\[
q_\lambda=I_\lambda^{-1}\left\{
\lambda^2\int_0^{u_\lambda}\operatorname{sech}u\,du
+\lambda\int_{u_\lambda}^{1/\lambda}\operatorname{sech}^2u\,du\right\},
\]
which is the displayed formula since \(\tanh u_\lambda=\sqrt{1-\lambda^2}\). Finally \(I_\lambda\to1\), \(\arccos\lambda\to\pi/2\), and the second term in braces is \(O(\lambda^3)\), proving the asymptotics.

**Implication and limit.** The \(\lambda^2\) order of a guarantee uniform over sources with a common Poincaré bound is optimal: the asymptotic improvement cannot be replaced by a function that grows faster than a constant times \(\lambda^2\) as \(\lambda\downarrow0\). This does not identify the best numerical constant. The source here changes with \(\lambda\); Section 4 supplies the distinct fixed-source statement.
