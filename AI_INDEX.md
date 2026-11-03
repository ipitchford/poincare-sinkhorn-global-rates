# AI index: curvature, relative entropy, and global Sinkhorn convergence (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with
[README.md](README.md) and the paper [Sinkhorn_Convergence.pdf](Sinkhorn_Convergence.pdf).
DOI: https://doi.org/10.5281/zenodo.23929617 (concept DOI 10.5281/zenodo.23929616).

## Normalisation

Entropic optimal transport between probability measures mu (source, on X) and nu (target, on Y) with cost c and regularisation
lambda > 0, objective int c dgamma + lambda KL(gamma | mu x nu). One iteration = one source normalisation then one target
normalisation (exact continuous transforms). Error measure: the semidual objective gap Delta(psi) (equation (3)); when an optimal
potential exists, Delta = lambda KL(gamma* | gamma_psi) (equation (34)). Poincare inequality (6): Var_mu(f) <= C_P int |grad f|^2 dmu.
Gradient condition (7): source gradients of c lie in balls of radius L. B = C_P L^2. Exact claims and fingerprints:
[prior_art/target.yaml](prior_art/target.yaml). Standing integrability assumptions (S1)-(S3): Section 3 of [manuscript.md](manuscript.md).

## Claims and where they are proved

| ID | Claim | Location in [manuscript.md](manuscript.md) | Evidence type |
|---|---|---|---|
| C1 | Theorem 1 (curvature transfer): a uniform curvature bound (8) with constant A >= lambda gives the reverse-entropy inequality (9)-(10) and the per-iteration contraction 1 - lambda/A of the dual gap (11) | Section 4 (4.1-4.3) | proof |
| C2 | Theorem 2 (Poincare rate): under (6), (7), (S1)-(S3), Delta(S psi) <= B/(B + lambda^2) Delta(psi) for every bounded psi; iteration count ceil((1 + B/lambda^2) log(1/tol)) (32) | Section 5, eqs (22)-(23), (32) | proof; floating-point checks of the inequalities in examples |
| C3 | Lemma 1 (exponential integrability of Poincare sources), Lemma 2 (cap bound for compactly supported targets) | Section 5.1 | proof |
| C4 | Theorem 3: Poincare sources on connected open X in R^d, cost Lipschitz in y with linear growth in x, compact target: factor C_P L^2/(C_P L^2 + lambda^2) for every lambda, without bounded cost | Section 5.1 | proof |
| C5 | Corollary 4: quadratic or bilinear cost, any Poincare source on R^d, any compactly supported target (including point clouds): factor C_P R^2/(C_P R^2 + lambda^2), R the target Chebyshev radius | Section 5.1 | proof |
| C6 | Laplace example: the two-sided exponential law (C_P <= 4, no Talagrand inequality) is covered by Corollary 4; Gaussian example | Section 5.1, Examples | proof; floating-point illustration (revision_checks.py part E: observed one-step ratios 0.136, 0.363, 0.552 below 4/(4 + lambda^2) at lambda = 1, 0.5, 0.3; certificate (33) holds) |
| C7 | Remark 5.1 (stopping certificate): Delta(psi_t) <= (1 + B/lambda^2)[E(psi_(t+1)) - E(psi_t)] | Section 5.2, eq (33) | proof; floating-point check |
| C8 | Inexact iterations: error floor 2 eta (1 + B/lambda^2) (35); conditional, no quadrature or sampling estimate | Section 5.3 | proof |
| C9 | Corollary 5 (CDV settings A1/A2/A3): CDV's curvature estimate plus Theorem 1 gives the global factor C/(C + lambda), about 1/lambda iterations, from iteration zero | Section 6 | proof (curvature input is CDV's Proposition 5.1) |
| C10 | Theorem 6 (fixed-source obstruction): the cusp source w(x) = (2/5)\|x\| log^2(e/\|x\|) on (-1,1) with two target atoms has local deficit ~ (8/5) log 2 lambda^2 log^2(1/lambda); no exponent below 2 is possible for a fixed source, even eventually | Section 7.1-7.2 | proof; floating-point checks (sharpness_check.py) |
| C11 | Proposition 7: order lambda^2 attained uniformly over a bounded Poincare class (local deficit ~ pi lambda^2) | Section 7.3 | proof; numerical constants only in Remark 7.2 |
| C12 | Proposition 8 (with Lemma 3): no single fixed absolutely continuous Poincare source attains order lambda^2 in the two-atom model | Section 7.4 | proof |
| C13 | Corollary 9: the sources above satisfy a Talagrand inequality with uniform constant, so the lambda^2 exponent of the Talagrand-inequality rates (CCGT) is also sharp | Section 7.5 | proof (uses Gozlan 2010, Theorem 4.7) |
| C14 | Remarks 4.1-4.2, 5.1, 7.1-7.2: maximal-correlation reading; mirror-descent route (AFKL); stopping certificate; power densities; numerical Poincare constants (how far Theorem 2 is from the truth) | Sections 4.4, 5.2, 7.4 | proof (4.1, 4.2, 5.1, 7.1); NUMERICAL estimate (7.2) |
| C15 | Remark 5.2: local limit 2 - D_lambda of the one-step ratio; the eventual factor (B/(B + lambda^2))^2 is a CONJECTURE (Open question 1) | Section 5.2, Section 9.3 | local limit: proof and floating-point check (revision_checks.py part F); eventual factor: CONJECTURE, not proved |

Practical reading: about (B/lambda^2) log(1/eps) iterations; about 1/lambda in the bounded-density settings of Corollary 5.

## Verification artefacts

Scripts are floating-point checks of selected formulas, constants and trajectories only; they are not proofs and not interval certificates.

- [code/sharpness_check.py](code/sharpness_check.py): uniform and cusp sources with two target atoms; gaps, error bound (22), curvature condition,
  ascent inequality over 180 (lambda, potential) combinations and ten trajectories. Output: [checks/sharpness_checks.json](checks/sharpness_checks.json).
- [code/multitarget_gaussian_check.py](code/multitarget_gaussian_check.py): Gaussian source, bounded cost sin(x) y, five target atoms.
  Output: [checks/multitarget_gaussian_checks.json](checks/multitarget_gaussian_checks.json).
- [code/revision_checks.py](code/revision_checks.py): refined expansion (44) (constants to 25 digits with mpmath), constants of Remark 7.1,
  numerical Poincare constants of Remark 7.2, the Laplace illustration of Corollary 4, the local limit of Remark 5.2.
  Output: [checks/revision_checks.json](checks/revision_checks.json).
- Each programme exits non-zero on any failed check, also under `python -O` (no assert statements). Outputs reproduce to
  roundoff, not bit for bit across Python/BLAS builds. Proof dependencies: [proof_outline.md](proof_outline.md);
  check scope and limits: [VERIFICATION.md](VERIFICATION.md).

## Exclusions

Exact algorithm on population (continuous) measures only: empirical sources do not satisfy the Euclidean Poincare inequality, and no
sampling complexity or bias relative to unregularised transport is bounded. Unbounded targets with quadratic cost (e.g. Gaussian to
Gaussian) are not covered. Constants are not sharp; the sharp universal constant and the eventual factor of Remark 5.2 are open.
Multimarginal extensions are open.

## Review and provenance

Internal AI-assisted review only: see [review/REVIEW_SUMMARY.md](review/REVIEW_SUMMARY.md). No human peer review,
independent human reproduction, proof assistant or interval arithmetic. Prior-art audit: [prior_art/novelty_report.md](prior_art/novelty_report.md).
Closest prior work: CCGT (same rate shape under a Talagrand inequality and Hessian bounds); AFKL (mirror-descent route);
CDV (curvature estimate used in Corollary 5); S. Di Marino's announced two-marginal result (not inspectable; priority for the
exponent uncertain). Bibliography: [refs.bib](refs.bib). Earlier versions: [versions/](versions/). Licences: [LICENSE.md](LICENSE.md).
