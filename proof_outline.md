# Proof outline: curvature, relative entropy, and global Sinkhorn convergence

Informal mathematical proofs, revised 3 November 2026. This outline is not an external human review or a proof-assistant certificate. Equation and result numbers refer to the revised `manuscript.md`.

## Statement

For the exact population Sinkhorn algorithm, measured by the semidual objective gap Δ, the paper proves:

* **Theorem 1.** Assume the integrability conditions (S1)–(S3_λ) of Section 3, which allow unbounded costs. If the semidual curvature dominates A⁻¹ times the variance in the *current* target marginal r_ψ, for every bounded potential, then Δ(ψ) ≤ A KL(ν‖r_ψ) and Δ(Sψ) ≤ (1 − λ/A)Δ(ψ).
* **Theorem 2.** Suppose the source satisfies a Poincaré inequality with constant C_P, and for each x the source gradients ∇ₓc(x,·) lie in a ball of radius L. Then A = λ + C_P L²/λ, so every iteration contracts by C_P L²/(C_P L² + λ²).
* **Theorem 3 and Corollary 4.** If the target is a compact metric space and the cost is Lipschitz in y with constant of linear growth in |x|, then (S1)–(S3_λ) hold for every λ. Quadratic cost, any Poincaré source on R^d and any compactly supported target give C_P R²/(C_P R² + λ²), where R is the Chebyshev radius of the target.
* **Corollary 5.** CDV's Proposition 5.1, fed into Theorem 1, gives the global factor C/(C + λ) under CDV's (A1), (A2) or (A3).
* **Theorem 6, Proposition 7, Proposition 8 and Corollary 9.** In a one-dimensional two-atom model, one fixed source has local deficit (8/5) log 2 λ² log²(1/λ), so no smaller exponent is possible even eventually. A λ-dependent family with uniformly bounded C_P has local deficit ~ πλ². No single fixed absolutely continuous source attains order λ². All of these sources satisfy a Talagrand inequality with uniform constant, so CCGT's λ² exponent is also sharp.

## Idea

The source Poincaré inequality bounds the variance of a conditional mean by the mean conditional variance. That is exactly the comparison between the semidual Hessian and the variance under the moving target marginal (a maximal-correlation bound). Convex conjugacy reverses the Hessian comparison into a relative-entropy bound. The inversion is done on each finite partition of the target, where it is finite-dimensional, and is then extended by uniform simple approximation. One exact Sinkhorn step gains at least λ KL(ν‖r_ψ), which closes the contraction argument. Unbounded costs enter only through integrability: lower density bounds by Jensen's inequality with an integrable oscillation, and upper density bounds from the compactness of the target (a cap argument) together with the exponential integrability that a Poincaré inequality provides.

## Outline and dependencies

1. **Differential quantities (standard; Section 3, eq. (5)).** The gradient of Φ is the current target marginal, and its Hessian is the mean conditional variance divided by λ. The derivatives are bounded uniformly in x, so no boundedness of c is needed.
2. **Integrability (new; Section 3).** (S1)–(S2) make Φ finite and bound dr_ψ/dν below; (S3_λ) bounds it above, so S maps bounded potentials to bounded potentials.
3. **Finite-partition conjugacy (Section 4.1, eqs. (12)–(16); coercivity rewritten).** Coercivity (14) now uses the integrated lower bound λ min log ν(B_i) − osc ψ − ∫ω dμ. Its kernel is span(1), and the inverse quadratic form on tangent vectors gives (15).
4. **Entropy integration and limits (Section 4.2, eqs. (17)–(19)).** Data processing and the 2-Lipschitz property of E in the supremum norm give (9) and (10).
5. **Ascent identity (Section 4.3, eq. (20); rewritten).** E(Sψ) − E(ψ) = λ KL(ν‖r_ψ) + λ KL(μ‖σ_ψ μ). This identity gives (11) and explains the factor-two slack discussed in Remark 5.2.
6. **Poincaré curvature (Section 5, eqs. (25)–(27); domination argument rewritten).** Local domination after subtracting c(x, y₀), the covariance identity, Cauchy–Schwarz with the centre of the gradient ball, and total variance.
7. **Unbounded costs (Section 5.1; new).** Lemma 1 (exponential integrability under (6), proved for the restricted test class), Lemma 2 (cap bound) and Theorem 3. Corollary 4 then follows for quadratic and bilinear costs.
8. **CDV settings (Section 6; extended).** CDV's Proposition 5.1 is quoted with its constants. The uniform source shows that the linear order is optimal.
9. **Sharpness (Section 7).** These steps use the exact scalar map (40)–(42) for any density, and the anchored Cauchy–Schwarz Poincaré certificates.
   * Theorem 6: dominated convergence, plus the trajectory statement.
   * Refined expansion (44): closed-form constants with an explicit remainder bound.
   * Proposition 7: closed form (47).
   * Lemma 3 (new): reciprocal density integrable.
   * Proposition 8 (new): no fixed source attains λ², via localisation of the moving optimum near the median and Cauchy–Schwarz on a window of length 2λ.
   * Corollary 9 (new): Gozlan's Theorem 4.7, restating Bobkov–Gentil–Ledoux, applied on a bounded interval.

## Prior work that the argument relies on

* CDV (Mathematical Programming 215, 2026), Proposition 5.1, Remarks 3.1 and 5.2, Theorems 1.2 and 3.1, and Corollary 3.2, read from arXiv:2407.01202v2.
* CCGT (arXiv:2412.09235v2), Theorems 1.3 and 1.4 and Propositions 1.8 and 1.9, including the proofs in Sections 2.5–2.6.
* AFKL (NeurIPS 2022), Lemma 6, Proposition 7, Theorem 4 and Proposition 8.
* Gozlan (2010), Theorem 4.7 and Proposition 4.9.

Convex-conjugate Hessian inversion, entropy Hessians, data processing, the variance decomposition and exponential concentration under Poincaré inequalities are established ingredients and are credited as such in the paper.

## Verification status

No proof-assistant verification was performed. The new and rewritten proofs are listed in the round-1 revision response for a separate proof audit. Floating-point checks are summarised in `VERIFICATION.md`. They support the displayed formulas and constants within their tested ranges; they do not prove the theorems.
