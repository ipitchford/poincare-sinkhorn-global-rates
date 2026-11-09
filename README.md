# Curvature, relative entropy, and global Sinkhorn convergence

Research bundle. First version 1 November 2026; revised 3 November 2026 after one round of review; version 0.1.1-candidate, 9 November 2026 (wording correction).

## Main results and status

The paper gives global, explicit convergence rates for the exact (continuous) Sinkhorn algorithm, measured by the entropic dual gap.

* **One Poincaré inequality on the source suffices.** Suppose the source has Poincaré constant C_P and, at each source point, the source gradients of the cost lie in a ball of radius L. Then every full iteration contracts the gap by C_P L²/(C_P L² + λ²), from any bounded initial potential (Theorem 2).
* **Unbounded costs with compact targets.** For quadratic cost, any Poincaré source on R^d and any compactly supported target, including point clouds, the factor is C_P R²/(C_P R² + λ²) for every λ > 0. Here R is the target's Chebyshev radius (Theorem 3 and Corollary 4). The two-sided exponential source is a covered case that satisfies no Talagrand inequality.
* **Global linear dependence on λ for bounded densities.** The curvature estimate of Chizat, Delalande and Vaškevičius (CDV) gives the global factor C/(C + λ) in their bounded-density and log-concave settings (Corollary 5).
* **Sharpness.** No exponent below 2 is possible for a fixed source, even eventually (Theorem 6). Order λ² is attained uniformly over a bounded Poincaré class (Proposition 7) but by no single fixed source (Proposition 8). The same examples show that the λ² exponent in the Talagrand-inequality rates of Chiarini, Conforti, Greco and Tamanini (CCGT) is sharp (Corollary 9).

Closest prior work: CCGT's Propositions 1.8(i) and 1.9(i) give the same rational factor, with a Talagrand constant τ ≥ C_P, under Hessian bounds on the cost. Here the Talagrand inequality is replaced by a Poincaré inequality and the Hessian bounds by a first-derivative condition, the factor holds for every λ, and the sharpness results are new. Section 9 of the paper gives a side-by-side comparison. CDV (Section 1.1) report an earlier announcement by S. Di Marino of λ² rates under Poincaré inequalities on both marginals; no priority is claimed for the exponent.

The mathematical statements are informal proofs. No proof assistant or interval arithmetic was used. The supplied floating-point programmes check selected formulas, constants and trajectories; they are not the proof of the theorems.

## Reading order

1. `Sinkhorn_Convergence.pdf`: the typeset paper. Section 1 is a plain-English summary and Section 2.4 a roadmap.
2. `proof_outline.md`: the mechanism and dependencies of the proofs.
3. `VERIFICATION.md`: what was checked, how, and the limits of those checks.
4. `prior_art/novelty_report.md`: the original (1 November) literature audit, with a dated note and a 3 November addendum. Its positioning against CCGT is superseded by Sections 2.3 and 9.1 of the paper.
5. `review/REVIEW_SUMMARY.md`: summary of the internal AI-assisted review, its objections and their resolution.
6. `AI_INDEX.md`: machine-oriented map of claims, evidence and limits.

## Bundle contents

| File or directory | Purpose |
| --- | --- |
| `manuscript.md` | Editable mathematical source for the paper |
| `prior_art/target.yaml` | Precise claims, normalisation, exclusions and formula fingerprints |
| `prior_art/novelty_report.md` | Dated literature audit (1 November 2026) with a 3 November note and addendum |
| `review/REVIEW_SUMMARY.md` | Rating-free summary of the internal AI-assisted review |
| `AI_INDEX.md`, `RESEARCH_GATES.json` | Machine-oriented claims map; contribution map, searches and limitation triage |
| `refs.bib` | The reference list as BibTeX, from Crossref and arXiv records |
| `LICENSE.md`, `CITATION.cff`, `.zenodo.json`, `environment.json` | Licences, citation and archive metadata, software versions |
| `.github/workflows/replay.yml` | CI replay of the three check programmes (normal and `-O`) |
| `checks/independent_sharpness_audit.md` | Earlier independent derivation of the sharpness examples (original numbering) |
| `code/` | Three standalone numerical programmes |
| `checks/*.json`, `checks/*.log` | Numerical outputs and the earlier replay provenance |
| `figures/` | PNG and SVG of the scalar contraction calculation |
| `versions/` | Earlier manuscripts, the round-1 input, and diffs between versions |
| `requirements.txt` | Numerical dependency versions |
| `MANIFEST.sha256` | SHA-256 of every file in the archive, excluding the manifest itself |

## Reproduce the numerical checks

The programmes need NumPy and SciPy. `revision_checks.py` also needs mpmath. Matplotlib is needed only to regenerate the plot. The recorded runs used NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8 and mpmath 1.4.1. From this directory:

```sh
uv run --with numpy==2.3.5 --with scipy==1.17.0 --with matplotlib==3.10.8 --with mpmath \
  python3 -I -B code/revision_checks.py --output-dir reproduction/revision
python code/sharpness_check.py --output-dir reproduction/binary
python code/multitarget_gaussian_check.py --output-dir reproduction/gaussian
```

Each command returns zero only when all its checks pass, including under `python -O`. Outputs are JSON and a text log in the selected directory. Without `--output-dir`, the programmes write into `checks/` and overwrite the delivered outputs. `sharpness_check.py` also writes `sharpness_plot.png` and `.svg` when Matplotlib is available.

## Rebuild the paper

With Pandoc and a LaTeX installation containing `pdflatex`, `amsmath`, `amssymb`, `amsthm`, `microtype`, `booktabs` and `float`:

```sh
pandoc manuscript.md --from=markdown+tex_math_dollars --standalone --shift-heading-level-by=-1 --pdf-engine=pdflatex --resource-path=. --output=Sinkhorn_Convergence.pdf
```

A rebuild may change PDF metadata and typesetting across versions, so the PDF's byte hash need not be stable even when the source is unchanged.

## Changes in 0.1.1-candidate

The wording about S. Di Marino's announced results (manuscript Sections 2.3 and 9, README, AI index, research gates, verification notes, review summary and prior-art notes) now cites CDV, Section 1.1, as the source of the report and states that no priority is claimed for the exponent. No theorem, proof, figure or code changed.

## How to cite

Anonymous (2026). *Curvature, relative entropy, and global Sinkhorn convergence* (version 0.1.1-candidate). Evidence Press. https://doi.org/10.5281/zenodo.24091432

Concept DOI (all versions): https://doi.org/10.5281/zenodo.23929616. Release page: https://evidencepress.org/releases/poincare-sinkhorn-global-rates/. Repository: https://github.com/ipitchford/poincare-sinkhorn-global-rates.

Licences: text and figures CC0-1.0; code MIT (see `LICENSE.md`).

## Scope that must accompany reuse

* All results concern the exact continuous Sinkhorn transforms and the regularised dual gap, from bounded initial potentials.
* Theorem 2 needs the integrability conditions (S1)–(S3) of Section 3, which hold automatically for bounded costs. Theorem 3 removes boundedness of the cost for compactly supported targets.
* Unbounded targets with quadratic cost, for example Gaussian to Gaussian, are not covered.
* Empirical sources do not satisfy the Euclidean Poincaré inequality, and nothing here bounds sampling complexity or the bias relative to unregularised transport.
* The sharp universal constant and the conjectured eventual factor (B/(B + λ²))² are not established.
* The curvature input to Corollary 5 is CDV's. The mirror-descent route from the reverse-entropy inequality to a linear rate is AFKL's. The variance decomposition is classical.

Text is released under CC0 and code under the MIT licence. The originating question is Remark 3.1 of Chizat, L., Delalande, A., and Vaškevičius, T. (2026), *Sharper exponential convergence rates for Sinkhorn's algorithm in continuous settings*, Mathematical Programming 215, 809–858, https://doi.org/10.1007/s10107-025-02242-z.
