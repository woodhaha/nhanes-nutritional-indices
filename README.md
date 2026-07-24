# Nutritional Indices and Mortality Across Cancer Types — NHANES Analysis

> **Manuscript**: Composite Nutritional Indices Predict Mortality Across Cancer Types: A Pooled NHANES Analysis With In-Depth Gastrointestinal Validation
> **Target Journal**: Clinical Nutrition
> **Analysis code for**: PNI, CONUT, and GNRI as predictors of all-cause mortality in 2,942 cancer patients (NHANES III 1988–1994 + Continuous NHANES 2005–2016)

## Project Structure

```
NHANES/
├── Scripts/                       # Main analysis pipeline (current)
│   ├── 01_prepare.R               # Data preparation and merging
│   ├── 02_analyze.R               # Main GI analysis (Cox, time-dependent, landmark, competing risks, RCS)
│   ├── 03_advanced.R              # PNI decomposition, NLR/SII, CRP, CALLY
│   ├── 04_revision_sensitivity.R  # Cohort-stratified, survey-weighted, extended covariates
│   ├── 05_extra_confounders.R     # Additional covariate adjustments
│   ├── 06_pancancer_prep.R        # Pancancer data preparation
│   ├── 07_pancancer_analysis.R    # Cross-cancer Cox models + interaction tests
│   ├── 08_survey_weighted_comparison.R  # Side-by-side unweighted vs survey-weighted
│   ├── flow_diagram.R             # Study flow diagram
│   ├── figures_merge.R            # Figure composition
│   ├── run_all.R                  # Pipeline orchestrator
│   └── audit.py / audit_final.py  # Numerical consistency checks
├── R_scripts/                     # Alternative pipeline version
│   ├── 04_gi_tumor_analysis.R     # GI tumor classification + analysis
│   └── 05_gi_figures.R            # GI-specific figures
├── results/
│   └── gi_analysis/               # All analysis output CSVs
├── manuscript.tex                 # Main manuscript (LaTeX)
├── supplementary.tex              # Supplementary materials
├── findings.md                    # Key results summary
└── research-state.yaml            # Project state tracking
```

## Reproducibility

1. **Data**: Publicly available from CDC at `https://wwwn.cdc.gov/nchs/nhanes/`
2. **R version**: 4.6.0
3. **Key packages**: `survival`, `survey`, `dplyr`, `ggplot2`, `splines`, `survminer`
4. **Run order**: `Scripts/01_prepare.R` → `02_analyze.R` → `03_advanced.R` → subsequent scripts

Note: NHANES III data uses the public-use mortality linkage file; Continuous NHANES data downloads via `nhanesA` package (or direct `haven::read_xpt` with the updated CDC URL format: `/Nchs/Data/Nhanes/Public/{year}/DataFiles/{table}.xpt`).

## Key Results

### Cross-Cancer (n=2,942, 1,135 deaths)
| Index | Pooled HR (95% CI) | p |
|-------|-------------------|---|
| PNI | 0.78 (0.73–0.84) | <0.001 |
| CONUT | 0.78 (0.73–0.82) | <0.001 |
| GNRI | 0.80 (0.75–0.86) | <0.001 |

### GI Cancer (n=353, 189 deaths) — Survey-Weighted
| Index | HR (95% CI) | p |
|-------|------------|---|
| PNI | 0.87 (0.61–1.23) | 0.432 |
| CONUT | 0.76 (0.63–0.91) | 0.003 |
| GNRI | 0.80 (0.61–1.05) | 0.103 |

## Manuscript

- `manuscript.tex` — Full manuscript formatted for Clinical Nutrition (elsarticle)
- `supplementary.tex` — Supplementary materials (Tables S1–S8, Figure S1)
- `references.bib` — Reference database

## Data Availability

All NHANES data are publicly available from the CDC. Analysis code is provided in this repository. No patient-level data is included.

## Contact

Zhuha Zhou — zhouzhuha@wmu.edu.cn
