# NHANES GI Tumor × Nutrition × Survival — Results Summary
# Generated: 2026-07-23
# Data: NHANES 2005-2016 (6 cycles) + NCHS Mortality 2019
# Nutrition indices: PNI (Prognostic Nutritional Index), CONUT, GNRI

## Data Sources
- `data/gi_analysis/nhanes_gi_nutrition_raw.rds` — 9,835 participants ≥60yo with complete nutrition
- `data/gi_analysis/mort_2019.rds` — 70,190 mortality records (CDC FTP)
- `results/gi_analysis/*.csv` — all analysis outputs

## Key Findings

### 1. Sample
- Total N=9,835 (≥60yo, complete nutrition biomarkers)
- GI tumor patients: n=199 (156 colon, 11 esophageal, 11 stomach, 8 rectal, 7 liver, 5 pancreas, 1 gallbladder)
- Deaths: 99/199 (49.7%) with 30 GI cancer-specific deaths (UCOD=002 w/ baseline GI tumor)

### 2. Cross-sectional: Nutrition in GI tumor vs non-cancer
| Index | GI Tumor | Non-Cancer | Adjusted β | p |
|-------|----------|------------|-----------|---|
| PNI | 419.1 | 428.3 | −7.07 | 0.002 |
| CONUT | 1.3 | 0.9 | +0.29 | <0.001 |
| GNRI | 664.6 | 677.7 | −9.28 | 0.005 |

### 3. All-cause survival (Cox PH)
Per 1-SD increase, adjusted for age/sex/race/education:
| Index | HR | 95% CI | p |
|-------|-----|--------|---|
| PNI | 0.655 | 0.525–0.818 | <0.001 |
| CONUT (reversed) | 0.759 | 0.628–0.917 | 0.004 |
| GNRI | 0.638 | 0.514–0.791 | <0.001 |

### 4. Median survival by PNI tertile
| PNI Tertile | Median | 5-yr OS | 10-yr OS |
|-------------|--------|---------|----------|
| Low (T1) | 5.83 yr | ~67% | ~38% |
| Mid (T2) | 8.67 yr | — | — |
| High (T3) | 11.33 yr | ~76% | ~50% |

### 5. Time-dependent effects
PNI/GNRI show strongest protection in first 2 years post-diagnosis (HR≈0.36), with attenuation after (interaction p<0.05). Effect persists across landmark analyses at 1/3/5 years. No era interaction (p=0.99).

### 6. Competing risks (GI cancer death vs other death)
| Cause | PNI HR | p |
|-------|--------|---|
| GI cancer death | 0.84 | 0.36 |
| Non-cancer death | 0.59 | <0.001 |
Nutrition primarily protects against non-cancer mortality (infection, CVD) — power limited for GI cancer-specific due to n=30 events.

## Output Files
```
results/gi_analysis/
├── table1_descriptive.csv        — Baseline characteristics by GI status
├── cross_sectional.csv           — Nutrition index differences
├── cox_allcause.csv              — All-cause mortality Cox results
├── cox_causespecific.csv         — Competing risks (GI cancer vs other)
├── cox_timedependent.csv         — Time-split Cox (0-2, 2-5, 5+ yr)
├── cox_landmark.csv              — Landmark analysis (1/3/5 yr)
├── rmst_results.csv              — RMST at 5 and 10 years
└── km_fit.rds                    — KM by PNI tertile (R object)
```
