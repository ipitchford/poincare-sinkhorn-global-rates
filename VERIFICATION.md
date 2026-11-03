# Verification and provenance

Original record: 1 November 2026. Revised: 3 November 2026, after one review round. This record separates analytic proof inspection, floating-point checks, source comparison and typesetting inspection.

## 1. What was established and how it was checked

| Result (revised numbering) | Verification basis | Important boundary |
| --- | --- | --- |
| Theorem 1: curvature in the moving marginal implies a global entropy bound and contraction | Finite-partition proof. In the first version it was reconstructed independently and checked line by line in review. In this revision it was re-checked under the integrability conditions (S1)–(S3), with coercivity rewritten as (14) | Informal mathematics; no proof-assistant certificate |
| Theorem 2: explicit Poincaré rate | Covariance identity, Cauchy–Schwarz, Poincaré inequality and total variance. The domination argument was rewritten for unbounded costs, and the constant uses the radius in (7). Binary-target, five-target and Laplace-source numerical checks | Requires (S1)–(S3) for the given λ |
| Theorem 3 and Corollary 4: Poincaré sources, compactly supported targets, unbounded (e.g. quadratic) costs | New proofs: Lemma 1 (exponential integrability), Lemma 2 (cap bound). Numerical illustration with a Laplace source, quadratic cost and three atoms | Compact targets only; Gaussian to Gaussian is not covered |
| Corollary 5: global rate C/(C + λ) in CDV's settings | Constants and hypotheses compared directly with CDV's Proposition 5.1 (arXiv:2407.01202v2); then Theorem 1 | The curvature input is CDV's published result |
| Theorem 6: fixed-source exponent obstruction | Normalisation, Poincaré certificate, exact scalar map, trajectory limit, dominated convergence. The refined expansion (44) has closed-form constants and a remainder bound, checked numerically | Local and eventual statement; it has a logarithmic correction |
| Proposition 7: uniform-class order λ² | Exact normaliser, reciprocal integral, Poincaré certificate, closed form (47) | Sources vary with λ |
| Proposition 8 and Lemma 3: no fixed source attains order λ² | New proofs (Section 7.4) | Two-atom bilinear model; local deficit |
| Corollary 9: sharpness for Talagrand sources | Gozlan (2010), Theorem 4.7 (Bobkov–Gentil–Ledoux), read from the primary text, applied on a bounded interval | Constant κ_G is large; only uniformity is used |

The new and rewritten proofs are listed, with their locations, in the round-1 revision response prepared for a separate proof audit. The earlier derivation in `checks/independent_sharpness_audit.md` uses the original numbering. (An earlier reconstruction of the curvature-to-entropy proof for bounded costs, in the original numbering, is not included in this archive; the paper's proofs supersede it.)

## 2. Numerical environment

The original runs used Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0, with Matplotlib 3.10.8 for the figure. The revision run of `code/revision_checks.py` used Python 3.11.14, NumPy 2.3.5, SciPy 1.17.0 and mpmath 1.4.1. It ran via `uv run ... python3 -I -B`, and again with `-O` into a scratch directory, with identical results. All programmes are local computations with no network access. Each exits nonzero if any check fails; none uses `assert`.

### Original programmes (unchanged)

* `code/sharpness_check.py` covers the uniform source and the cusp (38) at nine values of λ from 0.5 to 0.001 and ten potential amplitudes, which is 180 cases, plus ten trajectories. Its outputs are `checks/sharpness_checks.json` and `.log`, with 0 failures.
  * The largest error-bound ratios were 0.791499226 (uniform source) and 0.527191967 (cusp).
  * The trajectory gap ratios match the exact limit (1 − D_λ)² to within 3.46 × 10⁻¹⁰.
* `code/multitarget_gaussian_check.py` uses a source N(0,1), cost sin(x)y, five target atoms and 36 potential cases at 4096 and 16384 Gauss–Hermite nodes, plus six trajectories. Its outputs are `checks/multitarget_gaussian_checks.json` and `.log`, with 0 failures.

In this round, the review lead replayed both programmes from a fresh copy, with and without `-O`. Every reported value matched apart from roundoff-level fields. `checks/reviewer_replay_provenance.json` keeps the earlier replay record.

### New programme: `code/revision_checks.py`

Outputs are `checks/revision_checks.json` and `.log`, with 0 failures. The run took 13.5 s.

* **A. Refined cusp expansion (44).**
  * Quadrature gives a = −0.027532340317356827583 and b = 0.34588274767781242384, agreeing with the closed forms to 25 digits.
  * Direct quadrature of D_λ agrees with (44) within the stated remainder bound at λ = 10⁻¹, 5·10⁻², 10⁻², 10⁻³, 10⁻⁴ and 10⁻⁶. The difference at λ = 0.1 is 3.1 × 10⁻¹⁰, against the bound 1.4 × 10⁻⁶.
  * The ratio of D_λ to the leading term is 2.20, 1.87, 1.53, 1.33, 1.24 and 1.16 at the same values of λ.
* **B. Power-family constants** c_a = 2^{1−a} Γ(a+1) η(a) agree with quadrature at a = 0, 0.25, 0.5, 0.9 and 1.
* **C. Numerical Poincaré constants.** These come from Sturm–Liouville shooting in natural-scale coordinates (DOP853), in both parity sectors, at two tolerances. The two tolerances agree to within 10⁻⁷.
  * The uniform-source sanity check gives 4/π² = 0.4052847346.
  * Cusp (38): 0.634118, against the certificate 3/4. An independent P1 finite-element computation, not shipped, gave 0.6341182, 0.6341183 and 0.6341183 with 1000, 2000 and 4000 elements.
  * Power densities: 0.892507 (certificate 1) and 4.879012 (certificate 5).
  * Proposition 7 family at λ = 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005: 0.60350, 0.99835, 1.29493, 1.52314, 1.71564, 1.79975, 1.84999. Each is below the product certificate (∫ρ)(∫ρ⁻¹) ≤ 1 + π/2.
  * Earlier finite-element estimates in x-coordinates, which reviewers reported as 0.52–0.60 for the cusp, converge slowly because of the logarithmic weight; the natural-scale computation is regular.
* **D. Guarantee versus exact deficit** for the Proposition 7 family, using the certified constant: 7.24, 9.99, 9.72, 9.17, 8.65, 8.41 and 8.27 at the same values of λ. The limit is π(1 + π/2) ≈ 8.08. With the numerical constant the ratios are 2.19, 3.97, 4.91, 5.44, 5.77, 5.89 and 5.95.
* **E. Illustration of Corollary 4.** The source is the two-sided exponential law (quadrature on [−60, 60] with 20-point Gauss–Legendre panels of width 0.25, total mass 1 to 10⁻¹³). The cost is quadratic, and the target has atoms (−1, 0.3, 1) with masses (0.2, 0.5, 0.3), so R = 1 and C_P = 4.
  * The worst observed one-step gap ratios are 0.136, 0.363 and 0.552 at λ = 1, 0.5 and 0.3, below the guaranteed 0.800, 0.941 and 0.978.
  * The stopping certificate (33) held, with worst ratios 0.231, 0.092 and 0.049.
* **F. Local ascent ratio of Remark 5.2.** The ratio (E(Sψ) − E(ψ))/(λ KL) equals 2 − D_λ to 6 digits at λ = 0.5 and 0.1 (uniform source, θ = 10⁻⁶).

All values are floating-point approximations, not interval certificates. Discretising a source for a check does not make the discrete measure satisfy the Euclidean Poincaré inequality.

## 3. Source comparison

The positioning statements were checked against primary texts downloaded during the revision:

* CDV, arXiv:2407.01202v2: Theorems 1.1–1.3 and 3.1, Corollary 3.2, Remarks 3.1 and 5.2, Proposition 5.1, and Section 1.1 on the Di Marino announcement;
* CCGT, arXiv:2412.09235v2: Theorems 1.3–1.4, Propositions 1.5, 1.8 and 1.9, Sections 2.5–2.6, and convention (1.10);
* AFKL, NeurIPS 2022 proceedings PDF: Theorem 4, Lemma 6, Propositions 5, 7 and 8;
* Gozlan (2010): Theorem 4.7 and Proposition 4.9;
* Chafaï and Malrieu (2010), arXiv version: Section 4.1, decomposition of variance;
* the Simons Institute seminar page of 10 November 2021.

The CRM Pisa talk page returned HTTP 404, so the 2022 Pisa announcement is cited only through CDV's account. Every bibliographic entry was copied from a Crossref or arXiv record; the raw records are kept with the revision materials outside this bundle.

## 4. Version boundaries

* `versions/manuscript_round1_input.md` is the manuscript as it entered review on 3 November 2026.
* `versions/round0_to_round1.diff` records all changes in this revision.
* Earlier snapshots and diffs from 1 November are retained unchanged.
* Final byte hashes are in `MANIFEST.sha256`.

## 5. Artifact inspection

* Pandoc and pdfLaTeX produced `Sinkhorn_Convergence.pdf`.
* Every page was rendered to an image and inspected for visible TeX, overflow and broken mathematics.
* The Evidence Press PDF TeX preflight reported no raw TeX findings.

## 6. Remaining limitations

* No formal proof assistant, interval arithmetic, external human review, exhaustive bibliographic search or unpublished-proof comparison is supplied.
* The following are not covered:
  * unbounded targets with quadratic cost;
  * gradient inequalities for empirical sources;
  * sampling complexity;
  * the sharp universal constant;
  * a proof of the eventual factor (B/(B + λ²))²;
  * the bias relative to unregularised transport.
* Priority for the λ² exponent relative to Di Marino's announcements remains uncertain.
