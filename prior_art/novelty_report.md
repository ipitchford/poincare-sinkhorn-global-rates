# Prior-art and object-identity audit: geometric Sinkhorn rates

> **Revision note (3 November 2026).** This report records the literature audit of 1 November 2026 and is kept as written, apart from this note, one deleted phrase, and the addendum at the end. Its description of the CCGT overlap as indirect or qualitative (Sections 1, 4.2, 4.3, 6 and 9 below) is superseded: CCGT's Propositions 1.8(i) and 1.9(i) state the same rational factor, with the Talagrand constant tau in place of C_P, under a Talagrand inequality and Hessian bounds. The paper now presents its Poincaré result as a weakening of CCGT's hypotheses at an equal or better constant, not as a new rate. The current positioning, including AFKL's Theorem 4, CDV's Remark 5.2 and the Di Marino announcements, is in Sections 2.3 and 9.1 of the revised manuscript (with Section 9.2 for other work). Result numbers in this report refer to the original manuscript: its Theorem 2 is still Theorem 2, its Corollary 3 is now Corollary 5, its Theorem 4 is now Theorem 6, and its Proposition 5 is now Proposition 7. The recommended claim language in Section 9 below predates the revision and is not the paper's current wording.

**Audit date:** 1 November 2026.  
**Status:** bounded literature audit supporting further expert review; not a certificate of historical priority.  
**Target:** the exact statements and normalisations in `target.yaml`.

## 1. Decision in brief

The most defensible contribution is the **global contraction factor `C/(C+lambda)` obtained from the published Chizat–Delalande–Vaškevičius (CDV) curvature estimate**, including their bounded-log-density case. Their stated global guarantee in that case has quadratic dependence on the regularisation parameter, and their linear dependence is stated only eventually. The inspected later literature did not supply this global improvement under the same assumptions.

The direct **one-source Poincaré estimate with the explicit factor `C_P L^2/(C_P L^2+lambda^2)`** is a second plausible contribution. Its formula, first-derivative assumptions, arbitrary measurable target, absence of a source diameter, and reverse-entropy estimate are more precise than the closest results located. However, a broad claim to the first qualitative polynomial convergence theorem under a single Poincaré marginal is **not justified**. There is both an earlier two-marginal announcement and an indirect route through published transportation inequalities and the generalized semiconcavity results of Chiarini–Conforti–Greco–Tamanini (CCGT). The latter overlap is explained below and must remain in the manuscript.

The conversion from Hessian inequalities to a relative-entropy error bound uses established convex-conjugate and dual-relative geometry. It should be presented as a useful continuous-space implementation of that machinery, **not as a new general principle of convex duality**. The finite-partition passage and the sharp analytic constants are the substantive work in this implementation.

The fixed logarithmic-cusp example appears to be a useful sharpness result: the single source has a finite Poincaré constant, while the local dual-gap contraction deficit has order `lambda^2 log^2(1/lambda)`. This rules out every power smaller than two over the full Poincaré class. Proposition 5 adds a separate source family, allowed to depend on regularisation, whose Poincaré constants stay bounded and whose local dual-gap contraction deficit is asymptotic to `pi lambda^2`. Together these establish fixed-source exponent sharpness and exact quadratic order uniformly over a class with bounded Poincaré constants. They do not identify the sharp universal coefficient or establish an exact order `lambda^2` for a single fixed source.

## 2. Exact objects and claim units

Let `lambda>0` be the coefficient of relative entropy in the primal problem. Use

\[
q_x^\psi(dy)=\frac{\exp((\psi(y)-c(x,y))/\lambda)\nu(dy)}
 {\int\exp((\psi(z)-c(x,z))/\lambda)\nu(dz)},\qquad
\pi_\psi(dx,dy)=\mu(dx)q_x^\psi(dy),\qquad r_\psi=\pi_{\psi,Y}.
\]

Write `Phi=-K` for the convex negative entropic Kantorovich functional, and `delta(psi)=E(psi*)-E(psi)` for the semidual gap. A full Sinkhorn iteration means one target-to-source transform followed by one source-to-target transform. This convention matters when comparing contraction exponents or local eigenvalues.

| Unit | Exact proposed statement | Classification after audit |
| --- | --- | --- |
| U1: source Poincaré curvature | If `Var_mu f <= C_P int |grad f|^2 dmu` and `|grad_x c|<=L`, then `Phi''(psi)[h,h] >= Var_(r_psi)(h)/(lambda+C_P L^2/lambda)` uniformly over bounded tilts. | Plausibly new quantitative application; covariance differentiation and Poincaré are standard. No identical formula located. |
| U2: curvature to reverse entropy | Uniform `Phi'' >= k Var_(r_psi)` gives `KL(pi_eta || pi_psi) <= (k lambda)^(-1) KL(r_eta || r_psi)` in the bounded entropic family. | A continuous-space extension/application of known dual-relative curvature geometry. Do not claim the underlying Hessian inversion or entropy-relative framework as new. |
| U3: global one-Poincaré rate | `delta(psi) <= (lambda+C_P L^2/lambda) KL(nu || r_psi)` and `delta_(t+1) <= C_P L^2/(C_P L^2+lambda^2) delta_t` for every bounded initial potential. | Plausible quantitative theorem. Narrow priority to the explicit constant, weak derivative assumptions, measurable target, and all-iteration guarantee. Qualitative one-Poincaré polynomial convergence has substantial indirect prior overlap. |
| U4: global CDV improvement | If CDV Proposition 5.1 supplies `Phi'' >= Var_(r_psi)/(C+lambda)`, then `delta_(t+1) <= C/(C+lambda) delta_t` from every bounded initial potential. | Strongest apparent new consequence in this audit. It applies known curvature through known convex geometry, with a rigorous continuous-space bridge. |
| U5: fixed-source exponent obstruction | For `rho(x)=(2/5)|x| log^2(e/|x|)` on `(-1,1)`, target `(delta_-1+delta_1)/2`, and `c=-xy`, the local gap factor `r_lambda^2` satisfies `1-r_lambda^2 ~ (8/5)log(2) lambda^2 log^2(1/lambda)`. | No identical example located. It is an exponent obstruction, not an exact matching lower bound at order `lambda^2`. |
| U6: binary spectrum/formulas | For symmetric source and the binary target, `r_lambda=int tanh^2(x/lambda) dmu(x)`; the nonconstant local potential mode has eigenvalue `r_lambda`. | Elementary specialization of the known conditional-expectation/Sinkhorn linearisation. Use as a proof tool, not a stand-alone novelty claim. |
| U7: uniform-class quadratic sharpness (Proposition 5) | For `rho_lambda(x)=min{lambda cosh(x/lambda),1}` on `(-1,1)` and `w_lambda=rho_lambda/(2I_lambda)`, with the same binary target and bilinear cost, `C_P<=1+pi/2` uniformly and the local dual-gap contraction deficit is asymptotic to `pi lambda^2`. | No identical family or formula located in the bounded rerun. Plausible explicit uniform-class sharpness construction; source dependence on `lambda` and unpublished counterexample priority uncertainty must remain explicit. |

The project has no natural integer sequence whose OEIS identity would resolve the main theorem. Its diagnostic matrices are covariance matrices on finite partitions and the two-state alternating-conditional operator. The latter has eigenvalues `1,r_lambda`; changing the iteration convention squares the local dual-gap factor. These identities and the quotient by constants have been kept explicit.

## 3. The originating open question and the exact improvement

### CDV: primary source inspected

Lénaïc Chizat, Alex Delalande, Tomas Vaškevičius, **Sharper exponential convergence rates for Sinkhorn’s algorithm in continuous settings**, *Mathematical Programming* **215** (2026), 809–858; published online 1 July 2025. DOI: [10.1007/s10107-025-02242-z](https://doi.org/10.1007/s10107-025-02242-z). Full text: <https://link.springer.com/article/10.1007/s10107-025-02242-z>. Preprint: <https://arxiv.org/abs/2407.01202>.

Inspected: introduction and related work; Theorems 1.1, 1.2 and 3.1; Corollary 3.2; Remarks 3.1–3.5; Proposition 5.1 and its proof setting; the variance-transfer step; Gaussian sharpness discussion.

Remark 3.1 conjectures that a Poincaré inequality on the source should suffice. Its wording is not a fully parameterised conjecture and should not be enlarged to a claim that every rate in the paper follows from Poincaré alone. Proposition 5.1 gives

\[
\Phi''(\psi)[h,h]\geq\frac{\operatorname{Var}_{r_\psi}(h)}{C+\lambda},\qquad
C=e\kappa\left(c_\infty+\frac\xi2 R_\mathcal X^2\right)
\]

under their compact convex bounded-density assumptions, where `kappa=sup f/inf f`, `c_infty=osc(c)`, and `xi` is source semiconcavity. The bilinear-cost Theorem 1.2 states a global quadratic-parameter rate and an eventual linear-parameter rate. U4 removes that global-to-eventual loss. The curvature estimate itself is prior work and is attributed as such.

**Caveat on the conjecture:** the result developed here establishes the robust polynomial conclusion for bounded costs with a bounded source derivative. It does not prove every conceivable reading of “our results” in Remark 3.1. In particular U5 prevents a linear-in-`lambda` bound for the entire Poincaré class.

## 4. Closest overlaps that constrain the novelty claim

### 4.1 An earlier Poincaré/Sinkhorn announcement

CDV's related-work section reports an announcement by Simone Di Marino at the 2022 Pisa workshop *Optimal Transportation and Applications*: polynomial dual convergence for bounded Lipschitz costs when **both marginals** satisfy Poincaré inequalities, also extending to multiple marginals. CDV attributes the details to personal communication; that proof was not inspected.

A primary event page predating the cited workshop also exists: [The Speed of Sinkhorn Algorithm Meets Poincaré Inequality](https://simons.berkeley.edu/events/gmos-weekly-seminar-speed-sinkhorn-algorithm-meets-poincare-inequality), Simons Institute GMOS seminar, 10 November 2021. Its abstract describes algebraic dependence on the regularisation under Poincaré assumptions on the marginals and mentions counterexamples. No complete proof or slides were retrieved in this audit. Therefore neither a first Poincaré connection nor a first polynomial-rate claim is appropriate.

### 4.2 CCGT's semiconcavity and transportation route

Alberto Chiarini, Giovanni Conforti, Giacomo Greco, Luca Tamanini, **A semiconcavity approach to stability of entropic plans and exponential convergence of Sinkhorn’s algorithm**, [arXiv:2412.09235v2](https://arxiv.org/html/2412.09235v2), 3 October 2025.

Inspected: Theorems 1.1–1.4; Proposition 1.8, Proposition 1.9; their proof subsections and the definition of the asserted rate. Their Lipschitz-cost result assumes a `T2(tau)` inequality on one marginal and bounded cost Hessians and yields, for small `epsilon`, the factor

\[
1-\frac{\epsilon^2}{\epsilon^2+\tau\epsilon(H(c)-h(c))+\tau L^2}.
\]

Theorem 1.4 extends their method to a general transport cost `omega` and a matching generalized semiconcavity assumption. This is a genuine nearby prior mechanism, including a reverse entropy stability estimate for entropic plans. Its assumptions and parameter differ from U1–U3; it does not state the `C_P L^2` theorem or U4.

### 4.3 Why compact one-Poincaré qualitative priority cannot be claimed

Sergey G. Bobkov, Ivan Gentil, Michel Ledoux, **Hypercontractivity of Hamilton–Jacobi equations**, *Journal de Mathématiques Pures et Appliquées* **80** (2001), 669–696. Primary author PDF: <https://www.math.univ-toulouse.fr/~ledoux/HamiltonJacobi.pdf>. Inspected Corollary 5.1 and its cost definition (5.10).

That corollary converts a Poincaré inequality to a transportation inequality with a quadratic-linear cost. For a source of finite diameter, this cost is bounded below by a positive multiple of squared distance; consequently the source satisfies a `T2` inequality. Swapping the marginals in CCGT then yields qualitative polynomial Sinkhorn convergence for overlapping smooth compact cases. This is an **inference obtained by combining the two sources**, rather than a one-Poincaré theorem found verbatim in CCGT.

For an explicit normalization, Nathael Gozlan, **Poincaré inequalities and dimension free concentration of measure**, *Annales de l'Institut Henri Poincaré, Probabilités et Statistiques* **46** (2010), 708–739, DOI [10.1214/09-AIHP209](https://doi.org/10.1214/09-AIHP209), restates BGL in Theorem 4.7 and Corollary 4.8. Primary PDF: <https://www.numdam.org/article/AIHPB_2010__46_3_708_0.pdf>. It supplies a transport cost `alpha(|x-y|/sqrt(kappa C_P))`, where `alpha(u)=min(u^2,|u|)` and `kappa` is an absolute constant. If the support diameter is `D`, then

\[
\alpha\left(\frac r{\sqrt{\kappa C_P}}\right)
\geq\frac{r^2}{\max\{\kappa C_P,D\sqrt{\kappa C_P}\}},\qquad 0\leq r\leq D.
\]

There is also broader overlap on unbounded supports under additional global smoothness. This follows from the exact assumptions of CCGT Theorem 1.4; their Section 2.1 permits the opposite space to be Polish. Suppose the cost is globally twice differentiable in the Poincaré variable, with `|grad c|<=L` and `hI<=Hess c<=HI`. In CCGT's sign convention, each normalized potential gives `f=c+psi` satisfying `|grad f|<=2L` and `Hess f<=A I`, where `A=H-h+L^2/lambda`. Therefore its semiconcavity remainder is at most

\[
\min\{A r^2/2,4Lr\}.
\]

For `s^2=kappa C_P`, BGL supplies transport with `omega(r)=min(r^2/s^2,r/s)` and `W_omega<=KL`. The choice `Lambda=max(A s^2,8Ls)` gives the generalized semiconcavity bound required by CCGT, with `tau=1/2`. Their factor is then `(1+min{1,2lambda/Lambda})^(-1)`, which has a quadratic contraction deficit as `lambda` tends to zero. This is an **inferred corollary obtained here from the two published results**, not a verbatim statement found in either source. It confirms qualitative one-Poincaré overlap for that smoother class without a compactness requirement.

What U1–U3 add relative to these combinations is a direct all-tilt estimate using the source Poincaré constant itself, requiring only the source cost derivative and allowing a completely measurable target. Their precise reverse-entropy coefficient has not been located in the inspected sources.

### 4.4 Reverse entropy for entropic families was already known

Pierre-Cyril Aubin-Frankowski, Anna Korba, Flavien Léger, **Mirror Descent with Relative Smoothness in Measure Spaces, with application to Sinkhorn and EM**, *NeurIPS 2022*. Primary paper: <https://papers.nips.cc/paper_files/paper/2022/file/6e3daaeca6be8579573f69082b2dd58b-Paper-Conference.pdf>.

Inspected Proposition 7, equation (20), its proof sketch, and Proposition 8. For two entropic plans with a common marginal, Proposition 7 proves a reverse entropy inequality with coefficient `1+4 exp(3D_c/epsilon)`, using a projective-metric stability estimate and Pinsker. The associated framework is relative strong convexity of the marginal objective with respect to joint relative entropy.

Therefore U2 is not the first reverse data-processing inequality for this family. The issue is obtaining polynomial geometric coefficients and retaining the best coefficient furnished by a uniform tilted curvature estimate.

### 4.5 General dual-relative geometry is established machinery

Chris J. Maddison, Daniel Paulin, Yee Whye Teh, Arnaud Doucet, **Dual Space Preconditioning for Gradient Descent**, *SIAM Journal on Optimization* **31** (2021), 991–1016, DOI [10.1137/19M130858X](https://doi.org/10.1137/19M130858X). Primary author PDF: <https://www.stats.ox.ac.uk/~doucet/MaddisonPaulinTehDoucet_dualspacepreconditioningforgradientdescent.pdf>.

Inspected equation (1.5), Proposition 3.4, and Theorem 3.9. A lower Hessian bound relative to the inverse Hessian of a dual reference function is converted by convex conjugacy into a Bregman inequality; this gives global linear convergence with the corresponding condition number. That finite-dimensional mechanism encompasses the algebra used in U2. The paper does not derive the geometric Sinkhorn constants in U1 or U4.

Emanuel Laude and Panagiotis Patrinos, **Anisotropic proximal gradient**, *Mathematical Programming* **214** (2025), 801–845; DOI [10.1007/s10107-025-02204-5](https://doi.org/10.1007/s10107-025-02204-5), [arXiv:2210.15531v4](https://arxiv.org/html/2210.15531v4). Inspected Propositions 5.10–5.11 and the conjugate characterization of anisotropic strong convexity. This further confirms that conjugate-relative smoothness and its error-bound implications must be treated as prior machinery.

### 4.6 A 2026 finite-Sinkhorn application of that geometry

Zeyi Xu and Long Chen, **Accelerating Sinkhorn for Entropy-Regularized Optimal Transport**, [arXiv:2605.30267v1](https://arxiv.org/html/2605.30267v1), 28 May 2026.

Inspected Section 3, inverse Hessian and Schur-complement formulas, and Appendix A.3. The paper explicitly writes finite-dimensional Sinkhorn as a dual mirror method. Its linear-convergence argument obtains a positive but unspecified Hessian minimum on a compact sublevel set, then a dual PL inequality. It neither obtains a one-source Poincaré coefficient nor gives the global CDV linear-in-regularisation improvement. It is nevertheless directly relevant to the interpretation of U2.

## 5. Other inspected contemporary corpus

The following records were used to test whether the contribution had already appeared under other terminology. They are not all equally close, and the depth of inspection is stated to avoid representing an abstract skim as a complete proof audit.

| Primary work and inspected version | Inspection | Relation to target |
| --- | --- | --- |
| O. Deniz Akyildiz, Pierre Del Moral, Joaquín Miguez, **On the contraction properties of Sinkhorn semigroups**, [arXiv:2503.09887v3](https://arxiv.org/html/2503.09887v3), 27 January 2026; *Foundations of Data Science* **11** (2026), 150–206, [DOI 10.3934/fods.2026006](https://www.aimsciences.org/article/doi/10.3934/fods.2026006) | Introduction, stated contributions, bounded-cost contraction discussion, whole-text keyword searches | Lyapunov/minorisation and weighted contraction methods. No Poincaré term or matching `C_P L^2` theorem located; not every proof was read. |
| Pierre Del Moral, **A Contraction Theory for Sinkhorn and Schrödinger Bridges via Log-Sobolev Inequalities**, [arXiv:2503.15963v6](https://arxiv.org/html/2503.15963v6), 18 January 2026; [DOI 10.1080/07362994.2026.2619443](https://doi.org/10.1080/07362994.2026.2619443) | Abstract, introduction, model/assumption overview, whole-text searches | Uses logarithmic Sobolev, curvature and Riccati tools. No apparent mere one-source Poincaré theorem. Earlier versions have a different title and were not conflated with the latest version. |
| Anand Srinivasan and Jean-Jacques Slotine, **Thermodynamic structure of the Sinkhorn flow**, [arXiv:2510.12639v2](https://arxiv.org/html/2510.12639v2), 25 January 2026 | Latest Section 3.3, conditional-expectation Dirichlet form and Lemma 1, introduction, searches | Continuous-time operator spectral-gap and LSI viewpoint. Does not supply the source-Poincaré bounded-cost derivative estimate. Version 1 had a different title. |
| Stephan Eckstein and Aziz Lakhal, **Exponential convergence of general iterative proportional fitting procedures**, [arXiv:2502.20264v2](https://arxiv.org/html/2502.20264v2); *SIAM Journal on Optimization* **36** (2026), 912–937, [DOI 10.1137/25M1752584](https://doi.org/10.1137/25M1752584) | Main theorem overview, subspace-angle formulation, explicit bounded-cost constants, whole-text searches | General IPFP and Friedrichs angles; the inspected OT specialization carries exponential cost constants. No matching polynomial geometric coefficient. |
| Giacomo Greco and Luca Tamanini, **Hessian stability and convergence rates for entropic and Sinkhorn potentials via semiconcavity**, primary [author-repository PDF](https://cvgmt.sns.it/media/doc/paper/7087/stab_hessian.pdf), dated 16 April 2025; arXiv:2504.11133 | Introduction, assumptions, Theorems 1.1–1.2, searches | Gradient/Hessian convergence uses semiconcavity and transport inequalities. Different derivative-level stability target; no source-Poincaré-only coefficient. |
| Cyril Letrouit, **Quantitative stability in optimal transport**, September 2025 lecture notes, [primary PDF](https://www.imo.universite-paris-saclay.fr/~cyril.letrouit/teaching/Peccotfinal.pdf) | Chapter 5.3, Theorem 5.11 and Propositions 5.12–5.13; relevant Poincaré discussion | Exposition of the CDV global quadratic-parameter argument. Does not contain U4. |
| Bo'az Klartag and Or Ordentlich, **The strong data processing inequality under the heat flow**, [primary PDF, arXiv:2406.03427v1](https://www.weizmann.ac.il/math/klartag/sites/math.klartag/files/uploads/2406.03427v1.pdf); *IEEE Transactions on Information Theory* **71** (2025), 3317–3333, DOI 10.1109/TIT.2025.3548961 | Theorem 1.1, maximal-correlation identities, comparison of Poincaré and LSI divergence contractions | Gaussian-channel chi-square contraction is bounded using Poincaré. It is an important conceptual neighbor, but not the all-tilt entropic family or U3. |
| Cyril Letrouit and Tudor Manole, **Sharp Stability of Brenier Potentials and its Link to Poincaré Inequalities**, [arXiv:2610.37659v1](https://arxiv.org/html/2610.37659v1), 27 October 2026 | Abstract, introduction, full-text Sinkhorn search | Recent and relevant terminology, but unregularised Brenier stability. No Sinkhorn result located. |
| Promit Ghosal and Marcel Nutz, **On the Convergence Rate of Sinkhorn’s Algorithm**, [DOI 10.1287/moor.2024.0427](https://doi.org/10.1287/moor.2024.0427); *Mathematics of Operations Research* **51** (2026), 1080–1096 | Primary abstract and bibliographic metadata only | Polynomial-in-iteration entropy/dual bounds for unbounded costs, not the claimed geometric contraction. Full paper not proof-audited. |
| Giacomo Greco, Maxence Noble, Giovanni Conforti, Alain Durmus, **Non-asymptotic convergence bounds for Sinkhorn iterates and their gradients: a coupling approach**, [COLT 2023 primary page](https://proceedings.mlr.press/v195/greco23a.html) | Primary abstract and hypotheses overview | Torus densities, stochastic-control and reflection-coupling approach. Not a matched theorem; full paper not proof-audited. |

Searches also surfaced 2026 papers on martingale Sinkhorn, asymptotically scalable matrices, operator scaling, and the Sinkhorn Riemannian metric. Their retrieved statements concern different constraints, possibly zero kernels, noncommutative scaling, or derivatives of the divergence. They were not relied on to establish absence of prior art for the present theorem. In particular an attempted primary fetch of arXiv:2609.18850 was blocked, so that item remains uninspected beyond its retrieved abstract.

## 6. Structural reruns of the gate

1. **Original target:** CDV Remark 3.1, a Poincaré replacement for geometric density assumptions. Searched the explicit conjecture and Poincaré/Sinkhorn variants; located the two-marginal announcement and modern semiconcavity papers.
2. **First simplification:** source covariance differentiation plus total variance gives U1 directly. Searched conditional expectation, maximal correlation, Fisher information, spectral-gap and heat-flow aliases. The Gaussian-channel result is a conceptual precedent; no identical bounded-score all-tilt estimate was found.
3. **Second simplification:** conjugate Hessian inversion gives U2, hence U3 without a variance-transfer loss. Reran reverse/inverse data processing, entropy-relative smoothness, entropic PL, generalized Brascamp–Lieb, dual mirror, anisotropic convexity and all-tilt aliases. Located AFKL's exact reverse-entropy form and the MPTD/Xu–Chen geometric machinery; narrowed the novelty claim accordingly.
4. **Third simplification:** apply U2 to CDV's existing curvature. Reran bounded-log-density, global versus eventual, linear-in-regularisation and later-2026 searches. No exact global improvement under that source class was located.
5. **Fourth simplification:** binary-target logistic model and a logarithmic cusp give a fixed-source exponent obstruction. Searched two-point, vanishing-density, bottleneck, sharp Poincaré/Sinkhorn rate aliases. No identical example was located. The algebraic two-state spectrum is elementary and is not claimed as new.
6. **Unexpected overlap:** BGL's Poincaré-to-transport theorem plus CCGT proves a significant part of the qualitative claim indirectly. The gate was rerun after recognizing this relation; the broad “first one-Poincaré polynomial result” description was rejected.
7. **Additional sharpness construction:** Proposition 5 replaces the fixed logarithmic cusp with a regularisation-dependent capped hyperbolic-cosine density. A bounded formula/alias rerun checked `cosh`, `arcosh`, `sech`, spectral-gap, quadratic-counterexample and Di Marino variants. No exact source family was located. The nearest new hyperbolic-cosine hit was inspected and concerned the geometry of Sinkhorn divergence rather than iteration-rate sharpness; details are recorded in Section 10.

## 7. Exact focused query log

The following are the **exact query strings retained from the final focused audit**. Searches were run on 1 November 2026 using the two provided search engines; `S1` denotes the stronger engine and `S2` the routine engine. Earlier exploratory searches selected the target and populated the source ledger above; their full exact strings are not reconstructed here. The log is intentionally not represented as an exhaustive machine export of every search in the conversation.

### Poincaré/transport overlap

`S1`:

- `"Poincare" "bounded diameter" "T2"`
- `"Poincaré" "bounded support" "transportation inequality"`
- `"Poincare" "compact" "Talagrand" "implies"`

`S2`:

- `Poincare inequality equivalent transport quadratic linear cost Bobkov Gentil Ledoux`
- `"Poincaré" "T2" "bounded" "implies"`
- `"Poincaré" "transportation" "quadratic-linear"`

Outcome: BGL 2001, Gozlan 2010, and the quantitative transport-variance literature. BGL Corollary 5.1 and Gozlan Theorem 4.7 were opened and inspected; the crucial overlap follows from their formulas, not from a search snippet.

### Direct statement and newest-paper searches

`S1`:

- `"Sinkhorn" "relative strong convexity" Poincare`
- `"Sinkhorn" "Poincaré" "2026"`
- `"Sharper exponential convergence rates" "2026" "linear"`
- `"Sinkhorn" "reverse" "entropy" "Poincaré"`
- `"Sinkhorn" "bounded" "density" "linear in" convergence` — recency 730 days
- `"Sinkhorn" "global" "lambda" "Chizat" convergence` — recency 730 days
- `"Sinkhorn" "global" "epsilon" "relative strong convexity"` — recency 730 days
- `"Sinkhorn" "Poincaré" "2026" convergence` — recency 366 days

`S2`:

- `"Exponential Convergence of General Iterative Proportional Fitting Procedures"`
- `"Sinkhorn" "conjugate" "relative" "curvature"`
- `"strong data processing" "exponential family" "reverse"`
- `"Sinkhorn" "bounded log-density" "2026"`

Outcome: CDV, AFKL, Eckstein–Lakhal, Akyildiz–Del Moral–Miguez, Greco–Tamanini, Xu–Chen, and constrained variants. The relevant primary statements were inspected as recorded above. Repeated searches returning CDV are bounded negative evidence only.

### Mechanism and identity aliases

`S2`:

- `"reverse data processing" "Poincare"`
- `"strong data processing" "all" "tilts"`
- `"entropy" "conjugate" "Hessian" "Sinkhorn"`
- `"Sinkhorn" "Brascamp" "reverse"`
- `"Poincare" "maximal correlation" "Fisher"`
- `"Poincaré" "conditional expectation" "covariance" "Sinkhorn"`
- `"Sinkhorn" "log-density" "2026" convergence linear rate`
- `"Poincare" "channel" "score" "contraction"`
- `"Sinkhorn" "dual relative strong convexity"`
- `"Sinkhorn" "anisotropic" "convergence"`
- `"Sinkhorn" "relative strong convexity" "Poincaré"`
- `"Sinkhorn" "entropy" "inverse Hessian" convergence`
- `"Dual Space Preconditioning for Gradient Descent" arxiv`
- `"Anisotropic Proximal Gradient" Sinkhorn relative strong convexity`

`S1`:

- `"Sinkhorn" "Poincare" "variance"`
- `"Sinkhorn" "inverse data processing"`
- `"strong data processing" "exponential tilt"`
- `"reverse" "relative entropy" "conditional" "Poincaré"`
- `"Sinkhorn" "spectral independence"`
- `"Sinkhorn" "entropic independence"`
- `"reverse data processing" "Hessian"`
- `"relative smoothness" "conjugate" "entropy" "Hessian"`

Outcome: MPTD, Laude–Patrinos, Xu–Chen, AFKL, and heat-flow SDPI. The generic finite-dimensional bridge is established machinery, so U2's novelty description was narrowed.

### Sharpness aliases

`S1`:

- `"Sinkhorn" "vanishing density" convergence`
- `"Sinkhorn" "sharp" "Poincaré"`
- `"Sinkhorn" "bottleneck" "Poincare"`
- `"Sinkhorn" "log" "two-point" convergence rate`

Outcome: the 2021 seminar abstract mentions counterexamples without enough detail to compare; CDV has Gaussian lower bounds; no logarithmic-cusp fixed-source example was located.

## 8. Coverage limits and remaining missingness

- This was a focused primary-source web/arXiv/proceedings/publisher audit. It did not include a complete MathSciNet or zbMATH query, a full citation-graph export, author correspondence, or unpublished seminar notes.
- Di Marino's announced proof and counterexamples remain a material priority uncertainty. Their existence is explicitly acknowledged; their scope cannot be reconstructed from the abstracts alone.
- The indirect BGL–CCGT overlap was checked against their exact generalized-cost assumptions and covers the globally smooth class described above, including unbounded supports. It is a corollary inferred in this audit, not a separately published statement or externally reviewed theorem. It must not be silently omitted because it weakens a broad novelty slogan.
- Some 2026 arXiv/publisher records had version-date discrepancies. This report identifies the actual inspected version; it does not infer that every published version is identical to it.
- For several neighboring papers only the main statements and relevant sections were read. A negative keyword search is used to locate passages, not as proof that no equivalent result exists.
- The search log above records exact strings from the final focused phase, with honest missingness for earlier exploratory queries. It is not a reproducible search-engine snapshot: indexes and rankings can change.
- No journal acceptance or independent research priority is established by this audit. A prospective assessment should distinguish theorem correctness, importance of the rate improvement, and priority risk.

## 9. Recommended claim language

A defensible abstract-level description is:

> We give a direct quantitative Poincaré bound for the tilted entropic conditional kernels and use entropy-relative convex duality to obtain global Sinkhorn contraction. Applying the same argument to the curvature estimate of Chizat, Delalande and Vaškevičius upgrades their bounded-log-density convergence bound to a global rate linear in the regularisation parameter. A fixed one-dimensional Poincaré example rules out every smaller polynomial exponent. A separate regularisation-dependent source family with uniformly bounded Poincaré constants attains an exactly quadratic local contraction deficit, establishing the quadratic order uniformly over that class.

Follow immediately in the related-work section with the earlier Di Marino announcement, AFKL reverse entropy, MPTD/Xu–Chen dual-relative geometry, and the BGL–CCGT indirect overlap. Say that an identical theorem was **not located in the inspected corpus**, not that no such theorem exists. Describe the originating conjecture as resolved in the explicitly stated bounded-cost, source-derivative setting, with its qualitative overlap and limits made clear.

## 10. Proposition 5 addendum: bounded formula and alias rerun

**Date:** 1 November 2026. **Scope:** only the new uniform-class sharpness construction in manuscript Section 6.4. This addendum does not reopen the whole search or change the classification of U1–U6.

### Exact identity to distinguish

For `0<lambda<=1/2`, the raw density, normalizing half-integral and transition point are

\[
\rho_\lambda(x)=\min\{\lambda\cosh(x/\lambda),1\},\qquad
x_*=\lambda\operatorname{arcosh}(1/\lambda),
\]

\[
I_\lambda=\int_0^1\rho_\lambda(x)\,dx
=1-\lambda\operatorname{arcosh}(1/\lambda)
 +\lambda\sqrt{1-\lambda^2},\qquad
w_\lambda(x)=\frac{\rho_\lambda(x)}{2I_\lambda}.
\]

The diagnostic integral from the binary-target linearisation is

\[
D_\lambda=
\frac{\lambda^2\arccos\lambda+
\lambda\left[\tanh(1/\lambda)-\sqrt{1-\lambda^2}\right]}{I_\lambda}.
\]

The normalization tends to one, `D_lambda~(pi/2)lambda^2`, and the local dual-gap deficit is `1-(1-D_lambda)^2~pi lambda^2`. The source satisfies `C_P<=1+pi/2` uniformly, while the support, target, dimension, bounded cost and `L=1` all stay fixed. The source **does change with regularisation**. Each density is positive, but its lower bound deteriorates as `lambda` tends to zero; therefore this example does not contradict the global linear-rate theorem with a fixed bounded density ratio.

The object is a one-dimensional continuous probability density, not an integer sequence or a new matrix-scaling polynomial. The distinguishing formula is the cap of `lambda cosh(x/lambda)` and its elementary `arcosh/arccos` integrals; the two-state spectrum itself is U6 and is not claimed as new.

### Exact queries and inspected result

The root agent reported these exact `S2` queries for the new family:

- `"Sinkhorn" "cosh" "Poincare"`
- `"Sinkhorn" "arcosh" convergence`
- `"Sinkhorn" "counterexample" "Poincaré" "Di Marino"`

They returned unrelated records or already-known seminar material, without the family formula. The present bounded stronger-engine rerun used exactly three `S1` queries:

- `"Sinkhorn" "cosh" "Poincaré"`
- `"Sinkhorn" "sech" "spectral gap"`
- `"Sinkhorn" "quadratic" "counterexample" "Poincare"`

These returned CDV, the previously inspected continuous-time spectral-gap paper, transport-inequality papers concerning a different quantity, and one potentially relevant new formula hit. That hit was followed to its primary source:

Hugo Lavenant, Jonas Luckhardt, Bernhard Schmitzer, **Approximation and computation of the geodesic Sinkhorn distance**, [arXiv:2610.02007v1](https://arxiv.org/html/2610.02007v1), 1 October 2026. Inspected abstract and Example 4.13, including its hyperbolic-cosine formula. The example compares a Dirac mass with two separating Dirac masses and computes `S_epsilon(delta_0,(delta_-r+delta_r)/2)=(epsilon/2)log cosh(2r^2/epsilon)`. It concerns the divergence and its geodesic metric, not a capped source density, a Poincaré constant, or the Sinkhorn iteration-rate obstruction here. The shared binary/cosh structure does not identify the objects.

### Decision and remaining uncertainty

**Retain U7 as a plausible explicit new construction within the audited boundary.** The mathematical conclusion is exact quadratic order for the worst local deficit uniformly over a class whose Poincaré constants are bounded. It is separate from U5's obstruction for one fixed source, and neither conclusion identifies the sharp universal multiplicative constant.

The unpublished counterexamples mentioned by the 2021 Simons seminar abstract, and the Di Marino work reported by CDV, remain a material priority uncertainty. The search did not retrieve their formulas; absence of a matching indexed formula cannot establish independence or first priority. No further broad search was undertaken after this three-query rerun because no new close candidate warranted it.

## Addendum (3 November 2026): GitHub and corpus search

A further bounded search was run on 3 November 2026, after the revision. It had two parts.

**Evidence Press corpus.** The live Evidence Press feed (224 items on that date) was scanned for the keywords Sinkhorn, entropic, optimal transport, Poincaré, Schrödinger, transport, Wasserstein, matrix scaling and IPF. The hits concern different objects: Gaussian Gromov–Wasserstein exactness, non-integer sliced Wasserstein distances, sensitivity certificates for planar transport maps, and entropy contraction for a different dynamics, together with keyword-only matches. No earlier release concerns Sinkhorn iterations or their convergence rates. This release is therefore not a duplicate or an extension of an earlier release.

**GitHub.** Issue searches (`gh search issues`) used these queries: sinkhorn poincare; sinkhorn poincaré inequality; sinkhorn convergence rate poincare; sinkhorn global convergence rate; sinkhorn exponential convergence rate; sinkhorn convergence curvature; sinkhorn linear convergence lambda; sinkhorn dual gap contraction; sinkhorn relative entropy convergence; chizat delalande vaskevicius; sinkhorn spectral gap; sinkhorn convergence rate epsilon dependence; IPFP convergence rate poincare; schrodinger bridge sinkhorn rate poincare; sinkhorn semidual gap; sinkhorn convergence; sinkhorn rate; sinkhorn algorithm convergence proof; entropic OT convergence rate; sinkhorn iterations contraction; poincare inequality entropic. Discussion searches (GitHub GraphQL) used: sinkhorn poincare; sinkhorn convergence rate; sinkhorn curvature; entropic optimal transport convergence. One code search used `gh search code "sinkhorn poincare"`.

The hits were software tickets in optimal-transport and related libraries (for example POT, ott, moscot, sbi, otdd, bipca and dualsimplex) about non-convergence in code, and one unrelated issue in a hyperbolic recommender project, whose body was checked: there "Poincaré" refers to the Poincaré ball model. None states a Sinkhorn rate derived from a Poincaré inequality, a curvature-to-entropy transfer, or the sharpness constructions.

**Conclusion.** No public prior statement of the paper's results was found in these sources. This was a bounded search of named sources with the listed queries. It is not proof of novelty, and it does not change the priority uncertainty about Di Marino's announced results recorded above and in the paper.
