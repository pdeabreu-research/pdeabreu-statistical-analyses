# Statistical Analysis Portfolio

Selected, reproducible **R and Stata** workflows demonstrating advanced quantitative methods used in applied behavioral and organizational research.

The public files are deliberately **sanitized and generalized** from prior research workflows. They contain no confidential datasets, participant identifiers, collaborator information, private file paths, proprietary stimuli, or publication-restricted materials. Each workflow uses synthetic data by default so that the methodological logic can be inspected without exposing underlying research data.

## Featured workflows

### 1. Longitudinal panel models with robust inference — R
`workflows/01_longitudinal_robust_bootstrap.R`

Demonstrates:
- baseline-adjusted prospective panel models
- reciprocal-direction tests
- HC3 heteroskedasticity-robust standard errors
- parallel indirect-effect estimation
- participant-level nonparametric bootstrap inference (5,000 resamples)
- comparison of indirect pathways
- Benjamini-Hochberg false-discovery-rate control for exploratory screens
- explicit separation of prospective association from causal mediation claims

### 2. Multilevel count models for repeated engagement data — Stata
`workflows/02_multilevel_count_models.do`

Demonstrates:
- repeated-observation / panel data engineering
- network-normalized predictors
- temporal lag and rolling-history features
- mixed-effects Poisson regression with random intercepts
- cross-level interaction testing
- marginal effects and predicted-count visualization
- mixed-effects negative-binomial sensitivity analysis
- reproducible export of model results

### 3. Moderated mediation and serial indirect effects — Stata
`workflows/03_moderated_mediation_sem.do`

Demonstrates:
- factorial interaction models with robust inference
- structural equation modeling (SEM)
- parallel multiple-mediator models
- conditional indirect effects
- 5,000-resample nonparametric bootstrap inference
- serial moderated mediation
- marginal-effects visualization

## Reproducibility and privacy

These scripts are **methodological portfolio artifacts**, not publication replication packages. Project-specific variable names and substantive hypotheses have been generalized, and synthetic data are generated in-script. Real research datasets are intentionally excluded.

The workflows preserve the statistical architecture of analyses used in prior research while prioritizing auditability, privacy, and reproducibility.

## Software

- R: `sandwich`, `lmtest`, `boot`
- Stata: written for Stata 18; core models use built-in `mepoisson`, `menbreg`, `sem`, `bootstrap`, `margins`, and `collect`

## Repository structure

```text
.
├── README.md
├── LICENSE
├── .gitignore
├── workflows/
│   ├── 01_longitudinal_robust_bootstrap.R
│   ├── 02_multilevel_count_models.do
│   └── 03_moderated_mediation_sem.do
└── outputs/
    └── README.md
```

## Notes

The scripts are intentionally readable: assumptions, estimands, data-processing decisions, and inferential choices are documented in comments rather than hidden behind automated wrappers.
