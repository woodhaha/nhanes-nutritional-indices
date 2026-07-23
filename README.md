# NHANES Nutrition + Cognition + Depression Study (老年人)

> **Research Question**: In US adults ≥60, what is the association between multi-index nutritional status (DII/PNI/CONUT/GNRI) and cognitive function, and does depression mediate this relationship?

**Data**: NHANES 2011-2014 (only cycles with CERAD/AFT/DSST cognitive assessments)
**Design**: Cross-sectional, survey-weighted
**Status**: Analysis pipeline ready. Data not yet downloaded.

## File Structure

```
NHANES_cognition_nutrition/
├── run_all.R                    # Master runner — source this to run everything
├── README.md                    # This file
├── R_scripts/
│   ├── 00_config.R              # Package loading, paths, NHANES variable lists
│   ├── 01_load_and_derive.R     # Data download (nhanesA), merge, DII/PNI/CONUT
│   ├── 02_analysis.R            # Survey design, Table 1, models, RCS, subgroup, SEM
│   └── 03_figures.R             # 6 main + 5 supplementary figures
├── data/                        # Raw and derived data (auto-created)
├── results/                     # CSV output tables
└── figures/                     # PDF + PNG figures
```

## Quick Start

```r
# In R or RStudio:
setwd("C:/Users/woodh/Documents")
source("NHANES_cognition_nutrition/run_all.R")
```

## What the Pipeline Produces

### Results Tables
| File | Description |
|------|-------------|
| `table1_baseline.csv` | Weighted baseline by PNI tertile |
| `main_results.csv` | All models × all exposures (4 models × 4 indices) |
| `head2head_comparison.csv` | Standardized β comparison (Model 3) |
| `r2_increment.csv` | Incremental R² beyond base model |
| `rcs_predictions.csv` | RCS predicted values + 95% CI |
| `subgroup_results.csv` | Stratified β by 8 subgroups |
| `interaction_results.csv` | Interaction p-values |
| `mediation_sem.csv` | SEM parameter estimates (1000 bootstrap) |
| `sensitivity_analysis.csv` | 7 sensitivity checks |
| `depression_outcome.csv` | PHQ-9 as secondary outcome |

### Figures
| Figure | Content |
|--------|---------|
| Fig 2 | Head-to-head forest: 4 nutritional indices |
| Fig 3 | RCS dose-response: DII ~ cognitive function |
| Fig 4 | Subgroup forest: 8 strata |
| Fig 5 | Sensitivity tornado plot |
| Fig 6 | Mediation path diagram (DII → PHQ-9 → cognition) |
| S1 | Raw scatter + LOESS |
| S2 | Cognitive domain distributions |
| S3 | PNI tertile boxplots |
| S4 | Incremental R² bar chart |
| S5 | Sequential model stability |

## Key Variables

### Exposures (4 nutritional indices)
| Index | Formula | Data Source |
|-------|---------|-------------|
| **E-DII** | Energy-adjusted Dietary Inflammatory Index (28 components) | Day 1 24h recall |
| **PNI** | 10 × Alb(g/dL) + 0.005 × Lymph(/μL) | BIOPRO + CBC |
| **CONUT** | Alb + Chol + Lymph score (0-12) | BIOPRO + TCHOL + CBC |
| **GNRI** | 14.89 × Alb(g/dL) + 41.7 × (BMI/22) | BIOPRO + BMX |

### Outcomes
- **Cognitive composite Z-score**: mean(z_CERAD_imm, z_CERAD_del, z_AFT, z_DSST)
- **PHQ-9 depression**: continuous (0-27) + binary (≥10)
- **Probable MCI**: composite Z < -1 SD

### Covariates
Model 1: age + sex
Model 2: + race + education + income (PIR)
Model 3: + BMI + smoking + alcohol + physical activity + comorbidity + CRP
Model 4: + PHQ-9 (for cognition outcome)

## Novelty vs. Published Literature

| Published | This Study |
|-----------|-----------|
| DII → cognition (Zhang 2024, Du 2024) | **4 indices head-to-head** comparison |
| Single-index analysis | **PNI/CONUT blood-based** (more objective than dietary recall) |
| No food security | **Food security as upstream determinant** |
| DII → depression mediation (Sun 2022) | **Bidirectional mediation** (cog ↔ depression) |
| No kidney function consideration | **CKD stratification** (major elderly confounder) |

## Notes

- **First run** downloads ~200MB from CDC; cached to `data/nhanes_*_raw.RData`
- **NHANES 2011-2014 is the ONLY period** with cognitive assessment data
- **Survey weights**: Pooled 4-year MEC weight (WTMEC2YR/2)
- **Strata**: Combined cycle+original strata to ensure uniqueness
- DII calculation requires ≥20 of 28 components to be non-missing
- SEM uses MLR estimator with 1000 bootstrap (survey weights not natively supported by lavaan)

## References

- Zhang et al. (2024) DII + cognitive impairment, NHANES 2011-2014. *Front Aging Neurosci*. doi:10.3389/fnagi.2024.1371873
- Du et al. (2024) DII + PA + cognitive function. *Gen Hosp Psychiatry*. doi:10.1016/j.genhosppsych.2024.09.003
- Sun et al. (2022) DII → cognition → depression mediation. *Nutrients*
- Shivappa et al. (2014) DII development & validation. *Public Health Nutr*. doi:10.1017/S1368980013002115
- Onodera et al. (1984) PNI formula. *Jpn J Surg*
- Ignacio de Ulíbarri et al. (2005) CONUT score. *Nutr Hosp*
