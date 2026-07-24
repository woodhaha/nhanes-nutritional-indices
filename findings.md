# Findings — NHANES GI Cancer Nutritional Indices

## Current Understanding

We have comprehensively validated three composite nutritional indices (PNI, CONUT, GNRI) as prognostic factors for all-cause mortality in GI cancer patients across two independent NHANES cohorts spanning 30 years (1988-2019).

### Key Results

| Index | HR (95% CI) | p-value | Survives Bonferroni? |
|-------|-------------|---------|---------------------|
| PNI | 0.78 (0.63-0.96) | 0.020 | No (α=0.0083) |
| CONUT | 0.71 (0.61-0.83) | <0.001 | Yes |
| GNRI | 0.67 (0.57-0.80) | <0.001 | Yes |
| HEI-2015 | 0.97 (0.70-1.35) | 0.850 | N/A (null) |

### Time Dynamics
- PNI effect is strongly time-dependent: HR 0.34 in first 2 years → attenuation after
- CONUT/GNRI effects are more stable over time
- Landmark analysis: CONUT/GNRI remain significant at 1, 3, 5 years; PNI does not in conditional survival

### Mechanism
- Non-cancer death pathway (HR 0.63) >> GI cancer death (HR 1.06)
- Albumin drives the PNI effect; lymphocyte count contributes nothing independently
- No PNI × CRP interaction (p=0.780)
- PNI × NLR: joint exposure shows highest risk, sub-multiplicative interaction
- Dietary quality (HEI-2015) is not prognostic → biochemical > dietary pathway

### Optimal Threshold
- PNI cutpoint 48.5 (cross-validated mean 48.3)
- Below threshold: median survival 6.6 vs 12.8 years
- RMST difference: 1.97 years (1.19-2.76)
- PAF: 13.4% of deaths potentially attributable to low PNI

### Robustness
- E-values: PNI 1.66, CONUT 1.85, GNRI 1.97
- Consistent across NHANES generations (era interaction p=0.480)
- Consistent across GI sites (CRC vs non-CRC)
- Survives extended confounder adjustment

## Patterns and Insights
1. CONUT and GNRI are more robust than PNI (survive multiple testing correction)
2. The mechanism is physiological resilience, not anti-tumor effect
3. Simple lab values (albumin) outperform complex dietary assessments
4. Inflammatory and nutritional pathways share prognostic information

## Open Questions
- Can these findings replicate in a larger GI cancer cohort (e.g., UK Biobank)?
- Would longitudinal PNI measurements improve prediction?
- Does early nutritional intervention in low-PNI patients improve outcomes?
