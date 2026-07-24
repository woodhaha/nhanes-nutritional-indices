# Findings — NHANES Nutrition + Cognition + Depression

**Last Updated**: 2026-07-24
**Status**: Inner loop 1 complete — baseline analysis done

## Current Understanding

### H1: All 4 indices predict cognitive function — PARTIALLY SUPPORTED
Only **CONUT** (β=0.052, p=0.001, ΔR²=0.0041) and **E-DII** (β=0.040, p=0.034, ΔR²=0.0032) are significantly associated with cognitive function in fully adjusted models (Model 3). PNI (β=0.041, p=0.066) and GNRI (β=0.040, p=0.059) show trends but do not reach significance.

**Ranking**: CONUT > E-DII > PNI ≈ GNRI

**CONUT** (albumin + cholesterol + lymphocyte composite) objectively outperforms the dietary recall measure, which is clinically important.

### H2: Depression mediates nutrition-cognition — REFUTED (for DII)
- DII → PHQ-9 path: a = -0.132, p = 0.127 (NS)
- PHQ-9 → cognition: b = -0.017, p < 0.001
- Indirect effect: 0.002 (NS)
- Proportion mediated: -6.8%

Depression does NOT mediate the DII→cognition relationship. The DII→depression path is non-significant and goes in the opposite direction from predicted. This contrasts with published AHEI and vitamin D mediation studies.

### H3: Blood-based indices vs dietary recall — INCONCLUSIVE
- CONUT (blood-based) strongest: β=0.052, p=0.001
- PNI/GNRI (blood-based): trends only, p>0.05
- E-DII (dietary): significant β=0.040, p=0.034
- CONUT objectively outperforms, but PNI/GNRI don't

### H4: Food security interaction — INCONCLUSIVE
- Food secure: β=0.040, p=0.039 (significant)
- Food insecure: β=0.038, p=0.406 (NS, N=284)
- Similar point estimate but wide CI in insecure group
- Power limitation: only 284 food insecure with complete data

### Key Subgroup Patterns
- E-DII effect strongest in: age 60-69, higher income (>130% PIR), non-diabetic
- E-DII effect null in: lower income, diabetic patients
- This suggests **effect modification by SES and metabolic health**

### Sensitivity
- Primary analysis robust to: excluding stroke, cancer, depression
- Attenuated at age ≥65 (p=0.068) — effect driven by 60-69 age group
- Consistent across sex, education, CVD strata

## Completed Arms
### GI Cancer + PNI/CONUT/GNRI (manuscript complete)
- Dual-cohort NHANES III + Continuous NHANES (1988-2019), 313 GI cancer patients
- PNI: HR 0.78 (0.63-0.96, p=0.020), GNRI: HR 0.67 (0.57-0.80, p<0.001)
- Time-dependent effect: PNI strongest in first 2 years (HR 0.34)
- Competing risks: PNI protects against non-cancer death, not cancer-specific death
- Manuscript targeting *Clinical Nutrition*

### Cognition + Nutrition (baseline complete)
- NHANES 2011-2014, N=2,313 complete cases ≥60 years
- CONUT strongest predictor of cognitive function (β=0.052, p=0.001)
- Depression mediation NOT supported for DII→cognition
- Food security moderation unclear due to power limitations

## Open Questions
1. Does depression mediate PNI/CONUT/GNRI→cognition? (Only tested for DII so far)
2. Would a multi-index composite score outperform individual indices?
3. Is the SES gradient in DII→cognition real or an artifact of differential recall bias?
4. Should we pursue the null mediation finding as a negative results paper?
5. Do kidney function (CKD) affect the associations? (eGFR data not yet integrated)

## Lessons and Constraints
- `nhanesA` package has timeout issues with CDC server — use direct XPT download instead
- NHANES 2011-2014 CFQ variable names differ from documentation (use CFDCST1/2/3 for CERAD, CFDCSR for delayed recall)
- DPQ variables are DPQ010-DPQ090 (not DPQ010-DPQ019 as in older codebooks)
- `svyglm` predict returns `svystat` object (named vector with `var` attribute), not a list
- `data_frame()` is deprecated in dplyr — use `tibble()` or `data.frame()`
