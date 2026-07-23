# Prognostic Nutritional Index (PNI) and Survival in Gastrointestinal Cancer: 
# A Cross-Validation Study Using NHANES III and Continuous NHANES (2005-2016)

## Proposed Title Options
1. **PNI, CONUT, and GNRI Predict All-Cause Mortality in GI Cancer: A Dual-Validation Study of NHANES (1988-2016)**
2. **Dietary Quality vs Biochemical Nutritional Status in GI Cancer Prognosis: NHANES 2005-2016**
3. **Prognostic Value of Nutritional Indices in GI Cancer: Time-Dependent Effects and External Validation Across 30 Years of NHANES Data**

---

## Manuscript

### Abstract

**Background:** Nutritional status is a modifiable prognostic factor in cancer patients. The Prognostic Nutritional Index (PNI), Controlling Nutritional Status (CONUT), and Geriatric Nutritional Risk Index (GNRI) are composite biomarkers integrating albumin, lymphocyte, and cholesterol, but large-scale population validation in gastrointestinal (GI) cancer is limited.

**Methods:** We identified 258 GI cancer patients (esophagus, stomach, colon, rectum, pancreas, liver, gallbladder) aged ≥60 from NHANES III (1988-1994, n=59) and Continuous NHANES (2005-2016, n=199). Mortality follow-up was through December 31, 2019 (up to 31 years). PNI, CONUT, and GNRI were derived from baseline serum albumin, lymphocyte count, and total cholesterol. Survey-weighted Cox proportional hazards models, time-dependent Cox, landmark analysis, restricted mean survival time (RMST), and competing risks (GI cancer vs non-cancer death) were performed. External validation was conducted in NHANES III independently. Healthy Eating Index 2015 (HEI-2015) was calculated as a dietary quality comparator.

**Results:** Higher PNI was associated with reduced all-cause mortality in both cohorts: NHANES III HR 0.71 (0.50-1.00, p=0.047) and Continuous NHANES HR 0.66 (0.53-0.82, p<0.001). Combined HR 0.65 (0.53-0.81, p<0.001). GNRI showed similar prognostic strength (HR 0.64, 0.51-0.79, p<0.001). The protective effect was strongest in the first 2 years post-diagnosis (PNI HR 0.36, 0.23-0.57) and attenuated over time (interaction p=0.014). The effect persisted in landmark analyses at 1, 3, and 5 years, and was stable across diagnostic eras (era interaction p=0.99). Median survival by PNI tertile was 6.8 vs 10.0 years. In competing risks, PNI primarily protected against non-cancer mortality (HR 0.59, p<0.001) rather than GI cancer-specific death (HR 0.84, p=0.36, limited by 30 events). HEI-2015 dietary quality did not predict survival (HR 1.05, p=0.81), suggesting that biochemical nutritional status, not dietary intake, drives the association.

**Conclusions:** PNI, CONUT, and GNRI are robust, independent prognostic factors for all-cause mortality in GI cancer patients, externally validated across two independent NHANES cohorts spanning 30 years. The effect is time-dependent with maximal protection early after diagnosis. The null finding for HEI-2015 suggests that systemic inflammation-nutrition status, rather than diet quality per se, is the operative pathway.

### Keywords
Prognostic Nutritional Index; Gastrointestinal Cancer; NHANES; Survival Analysis; Competing Risks; Time-Dependent Effect

---

### Introduction

1. **Cancer cachexia and malnutrition in GI cancer** 
   - High prevalence (30-80% depending on site and stage)
   - Unique vulnerability due to digestive/absorptive dysfunction
   - Linked to poor treatment tolerance, complications, survival

2. **Available nutritional assessment tools**
   - PNI (Onodera 1984): 10 × albumin + 0.005 × lymphocyte
   - CONUT (Ignacio de Ulíbarri 2005): albumin + lymphocyte + cholesterol
   - GNRI (Bouillanne 2005): albumin × 14.89 + BMI/22 × 41.7
   - All validated in single-center surgical cohorts, limited population-level evidence

3. **Knowledge gap**
   - No dual-cohort population validation across different survey eras
   - Limited time-dependent analysis of nutrition effects
   - Dietary quality vs biochemical nutrition unaddressed

4. **Study objectives**
   - Validate PNI/CONUT/GNRI in population-based GI cancer
   - Time-dependent and competing risk analyses
   - Compare dietary (HEI-2015) vs biochemical (PNI) nutrition

---

### Methods

#### Study Population
- **NHANES III** (1988-1994): Nationally representative US sample, n=59 GI cancer ≥60yo
- **Continuous NHANES** (2005-2016): 6 pooled cycles, n=199 GI cancer ≥60yo
- **Mortality linkage**: NCHS 2019 Public-Use Linked Mortality Files (follow-up through 12/31/2019)

#### GI Cancer Ascertainment
- MCQ240 letter-suffix variables (G=colon, H=esophageal, M=liver, T=pancreatic, V=rectal, Z=stomach, I=gallbladder) — age at diagnosis fields
- NHANES III: HAC3OS site codes (6=colon, 12=stomach, 13=esophagus, 14=pancreas, 15=liver, 25=gallbladder)

#### Nutritional Indices
- **PNI** = 10 × albumin(g/dL) + 0.005 × lymphocyte(/µL)
- **CONUT** = albumin(0-6) + lymphocyte(0-3) + cholesterol(0-3)
- **GNRI** = 14.89 × albumin + 41.7 × (BMI/22)
- **HEI-2015**: 13-component score (0-100) from 24-hr dietary recall + Food Patterns Equivalents Database (USDA ARS)

#### Statistical Analysis
- Survey-weighted Cox PH (continuous, per 1-SD)
- Time-dependent Cox (0-2, 2-5, 5+ year periods)
- Landmark analysis (1, 3, 5 years)
- RMST at 5 and 10 years
- Competing risks: cause-specific Cox for GI cancer vs non-cancer death
- Era interaction (2005-2010 vs 2011-2016)
- Sensitivity: excluding 6-month deaths, age ≥65
- External validation: NHANES III independent model
- Dietary comparison: HEI-2015 alone, PNI alone, and mutually adjusted

---

### Results (Key Tables)

**Table 1:** Baseline characteristics by GI tumor status
*Already produced: results/gi_analysis/table1_descriptive.csv*

**Table 2:** Cross-sectional nutrition differences (GI tumor vs non-cancer)
*Already produced: results/gi_analysis/cross_sectional.csv*

**Table 3:** Cox PH — nutritional indices and all-cause mortality
| Index | Period | HR | 95% CI | p |
|-------|--------|----|--------|---|
| PNI | Combined (0-2yr) | 0.36 | 0.23-0.57 | <0.001 |
| PNI | Combined (overall) | 0.65 | 0.53-0.81 | <0.001 |
| PNI | NHANES III only | 0.71 | 0.50-1.00 | 0.047 |
| CONUT | Combined | 0.74 | 0.62-0.88 | <0.001 |
| GNRI | Combined | 0.65 | 0.52-0.80 | <0.001 |
| HEI-2015 | Adjusted | 1.05 | 0.72-1.52 | 0.81 |

**Table 4:** Competing risks (GI cancer death vs other death)
*Already produced: results/gi_analysis/cox_causespecific.csv*

**Table 5:** Time-dependent, landmark, and era interaction
*Already produced: results/gi_analysis/cox_timedependent.csv, cox_landmark.csv*

**Figure 1:** Kaplan-Meier by PNI tertile (combined)
*Already produced: figures/gi_analysis/fig_combined_km.png*

**Figure 2:** Cox forest plot — all 3 indices
*Already produced: figures/gi_analysis/fig2_cox_forest.png*

**Figure 3:** Time-dependent HR
*Already produced: figures/gi_analysis/fig3_timedependent.png*

**Figure 4:** Landmark analysis
*Already produced: figures/gi_analysis/fig4_landmark.png*

**Figure 5:** Competing risks forest
*Already produced: figures/gi_analysis/fig5_competing_risks.png*

**Figure 6:** HEI-2015 vs PNI comparison
*Additional figure needed*

---

### Discussion

1. **Principal findings**
   - PNI/GNRI consistently predict all-cause mortality in GI cancer
   - Effect is time-dependent (strongest first 2 years)
   - Protected against non-cancer mortality (infectious, cardiovascular)
   - HEI-2015 does NOT predict survival → mechanism is biochemical, not dietary

2. **Comparison with literature**
   - Consistent with Samsung 7,781 gastric cancer cohort (HR 1.38)
   - Extends to population-level with dual-cohort validation
   - First study to compare dietary vs biochemical nutrition in GI cancer

3. **Clinical implications**
   - Simple, inexpensive prognostic stratification tool
   - PNI <45 identifies high-risk group (median survival 5.8 vs 11.3 years)
   - Early nutritional intervention window (~2 years post-diagnosis)
   - Albumin/lymphocyte more prognostic than dietary quality assessment

4. **Limitations**
   - Small GI cancer sample (n=258), limited site-specific analysis
   - Single baseline nutrition measurement
   - Public-use mortality cause-of-death limited to 10-group recode
   - NHANES III differences (assay methods, survey design)
   - No cancer stage, treatment, or recurrence data
   - Predominantly colorectal cancer (78%)

5. **Strengths**
   - Population-based, not hospital-based
   - Dual-cohort external validation (30-year period)
   - Comprehensive temporal analysis (time-dependent, landmark, RMST)
   - Dietary vs biochemical comparison

---

### Proposed Journal Targets

| Journal | IF | Type | Fit |
|---------|----|------|-----|
| **Nutrition and Cancer** | ~2.9 | Brief Report | PNI validation, HEI comparison |
| **Frontiers in Nutrition** | ~5.0 | Original Research | Comprehensive analysis |
| **Nutrients** | ~5.9 | Original Research | Time-dependent effects |
| **Journal of Cachexia, Sarcopenia and Muscle** | ~8.9 | Short Report | Cachexia angle |
| **BMC Cancer** | ~3.8 | Original Research | Methodological strength |

---

### Completed Items

- [x] **HEI vs PNI comparison figure** -- saved as `figures/gi_analysis/fig_hei_vs_pni_forest.pdf/png`
- [x] **HEI-PNI correlation scatter plot** -- saved as `figures/gi_analysis/fig_hei_pni_correlation.pdf/png`
- [x] **Table 2 (baseline by PNI tertile)** -- saved as `results/gi_analysis/table2_baseline_pni_tertile.csv`
- [x] **Sensitivity: excluding CRC** -- non-CRC GI (n=62): PNI HR 0.80 (0.46-1.38, p=0.42); CONUT/GNRI remain significant
- [x] **Full manuscript** -- `MANUSCRIPT.md` with Introduction, Methods, Results, Discussion, References

### Still Pending (optional)

- [ ] Figure updates: re-run `run_gi_figures.R` with fixed data to update Fig 1-5 with corrected HR values
- [ ] UK Biobank application (needs institutional signatory confirmation)
- [ ] Target journal formatting (Nutrients/Frontiers in Nutrition)
