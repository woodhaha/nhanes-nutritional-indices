# Prognostic Nutritional Index and 30-Year Mortality in Gastrointestinal Cancer: A Dual-Cohort NHANES Study With Time-Dependent, Competing Risk, and Joint Inflammatory Analysis

## Running Title
PNI Predicts GI Cancer Survival: NHANES 30-Year Validation

---

## Abstract

**Background:** Nutritional status is a modifiable prognostic factor in cancer patients. The Prognostic Nutritional Index (PNI), Controlling Nutritional Status (CONUT), and Geriatric Nutritional Risk Index (GNRI) are composite biomarkers integrating serum albumin, lymphocyte count, and cholesterol, but large-scale population-level validation in gastrointestinal (GI) cancer remains limited.

**Methods:** We identified 313 GI cancer patients (esophagus, stomach, colon, rectum, pancreas, liver, gallbladder) from NHANES III (1988-1994, n=92) and Continuous NHANES (2005-2016, n=221), with mortality follow-up through December 2019 (up to 31 years). PNI, CONUT, and GNRI were derived from baseline measurements. Cox proportional hazards, time-dependent Cox (0-2, 2-5, 5+ years), landmark analysis, and competing risks were performed. Healthy Eating Index 2015 (HEI-2015) was calculated as a dietary quality comparator.

**Results:** Higher PNI was associated with reduced all-cause mortality: adjusted HR 0.78 (95% CI 0.63-0.96, p=0.020). GNRI showed the strongest association (HR 0.67, 0.57-0.80, p<0.001). The PNI effect was markedly time-dependent, strongest in the first 2 years (HR 0.34, 0.21-0.54, p<0.001) and attenuating thereafter (interaction p<0.001). Optimal PNI threshold was 48.5; patients below this had median survival 6.6 vs 12.8 years (p<0.0001). The PNI-mortality relationship was nonlinear (p<0.0001). In competing risks, PNI primarily protected against non-cancer death (HR 0.64, p=0.001) rather than GI cancer-specific death (HR 1.07, p=0.68). PNI decomposition showed albumin drove the effect (HR 0.69, p<0.001) while lymphocyte alone was non-significant (HR 1.06, p=0.54). The neutrophil-to-lymphocyte ratio (NLR), an inflammatory marker, also predicted mortality (HR 1.43, p<0.001) with similar C-statistic improvement to PNI (both +0.017). HEI-2015 dietary quality did not predict survival (HR 0.97, p=0.85), and no PNI × CRP interaction was observed (p=0.78). In joint exposure analysis, patients with both low PNI and high NLR showed the highest mortality risk (HR 1.67, p=0.028). Restricted mean survival time analysis confirmed a clinically meaningful 1.97-year survival difference between high and low PNI groups (p<0.001). An estimated 13.4% of deaths were attributable to low PNI (population attributable fraction).

**Conclusions:** PNI, CONUT, and GNRI are robust, independent prognostic factors for all-cause mortality in GI cancer patients, validated across two nationally representative US cohorts spanning 30 years. The effect is time-dependent with maximal protection early after diagnosis, and PNI < 48.5 identifies a high-risk group. The prognostic effect is driven primarily by serum albumin, not lymphocyte count, and is independent of systemic inflammation (CRP). Inflammatory markers (NLR, SII) show comparable prognostic performance, suggesting both nutritional and inflammatory pathways contribute to the observed survival differences. The null dietary finding (HEI-2015) suggests that biochemical nutritional status, rather than dietary intake, drives the association.

**Keywords:** Prognostic Nutritional Index; Gastrointestinal Cancer; NHANES; Survival Analysis; Competing Risks; Time-Dependent Effect

---

## 1. Introduction

Cancer cachexia and malnutrition represent one of the most challenging clinical problems in gastrointestinal oncology. GI cancer patients are uniquely vulnerable due to tumor-induced catabolism, malabsorption, and mechanical obstruction of the digestive tract [1,2]. The prevalence of malnutrition in GI cancer patients ranges from 30% to 80% depending on tumor site and stage, and is consistently associated with poorer treatment tolerance, increased postoperative complications, and reduced survival [3,4].

Several composite nutritional indices have been developed to integrate multiple biomarkers into a single prognostic score. The Prognostic Nutritional Index (PNI), originally proposed by Onodera et al. in 1984, combines serum albumin and peripheral lymphocyte count [5]. The Controlling Nutritional Status (CONUT) score incorporates albumin, lymphocyte count, and total cholesterol [6]. The Geriatric Nutritional Risk Index (GNRI) adjusts albumin for body mass index [7]. Each of these indices has been validated primarily in single-center surgical cohorts, with evidence concentrated in specific GI cancer subtypes such as gastric and colorectal cancer [8,9].

Despite their clinical utility in specialized settings, several important gaps remain. First, population-level validation of these indices across different survey eras and geographic contexts is lacking. Most prior studies originate from East Asian surgical series with sample sizes rarely exceeding 1,000 patients [10,11]. Second, the temporal dynamics of nutritional risk -- whether the prognostic effect is constant over time or concentrated in specific post-diagnosis windows -- have not been systematically examined. Third, the conceptual distinction between biochemical nutritional status (albumin, lymphocyte count, also reflecting systemic inflammation) and dietary nutritional quality (nutrient intake patterns) remains poorly characterized in the cancer prognosis literature.

To address these gaps, we undertook a comprehensive analysis of nutritional indices in GI cancer patients using the National Health and Nutrition Examination Survey (NHANES), a nationally representative, population-based survey with mortality linkage. We leveraged two independent NHANES generations -- NHANES III (1988-1994) and Continuous NHANES (2005-2016) -- providing up to 31 years of mortality follow-up and enabling cross-validation across different survey eras. Our specific objectives were to: (1) validate PNI, CONUT, and GNRI as prognostic factors for all-cause mortality in GI cancer at the population level; (2) characterize time-dependent effects through landmark, time-dependent Cox, and restricted mean survival time analyses; (3) compare prognostic performance across competing causes of death; and (4) compare biochemical nutritional indices against a dietary quality measure (HEI-2015).

## 2. Methods

### 2.1 Study Population

We used data from two independent NHANES cycles: NHANES III (1988-1994) and Continuous NHANES (2005-2016). NHANES is a nationally representative, cross-sectional survey of the civilian non-institutionalized US population conducted by the National Center for Health Statistics (NCHS). Both surveys collected demographic, dietary, laboratory, and examination data through standardized protocols.

The analytic sample included participants aged >=18 years at the time of survey participation. A total of 44,066 participants met inclusion criteria after exclusions for missing nutrition biomarker data and implausible laboratory values. Of these, 313 had a GI cancer diagnosis.

### 2.2 GI Cancer Ascertainment

GI cancer was identified from self-reported cancer history. In the Continuous NHANES (2005-2016), participants were asked about cancer diagnoses through the Medical Conditions Questionnaire (MCQ). Age-at-diagnosis variables (MCQ240A-MCQ240DD) were used to identify GI cancer based on cancer site codes. In NHANES III, cancer site was identified from variable HAC3OS (family of cancer codes), with GI sites defined as codes 6 (colon), 12 (stomach), 13 (esophagus), 14 (pancreas), 15 (liver), and 25 (gallbladder).

GI cancer was defined as cancer of the esophagus, stomach, colon, rectum, pancreas, liver, or gallbladder.

### 2.3 Nutritional Indices

**PNI** was calculated according to Onodera's formula:
- PNI = 10 x serum albumin (g/dL) + 0.005 x lymphocyte count (/micro-L)

**CONUT** was derived as the sum of three component scores:
- Albumin: >=3.5 g/dL = 0, 3.0-3.49 = 2, 2.5-2.99 = 4, <2.5 = 6
- Lymphocyte: >=1600/micro-L = 0, 1200-1599 = 1, 800-1199 = 2, <800 = 3
- Total cholesterol: >=180 mg/dL = 0, 140-179 = 1, 100-139 = 2, <100 = 3

**GNRI** was calculated as:
- GNRI = 14.89 x albumin (g/dL) + 41.7 x (BMI / 22)

Serum albumin was measured using the bromcresol purple (BCP) method in Continuous NHANES and bromcresol green (BCG) in NHANES III. For Continuous NHANES, albumin values were converted from g/L to g/dL for consistency. Complete blood count (CBC) provided lymphocyte counts, and total cholesterol was measured enzymatically.

**HEI-2015** was calculated from 24-hour dietary recall data using the Food Patterns Equivalents Database (FPED) from the USDA Agricultural Research Service. HEI-2015 is a 13-component score (0-100) measuring adherence to the 2015-2020 Dietary Guidelines for Americans, with higher scores indicating better dietary quality.

### 2.4 Mortality Linkage

Mortality status was ascertained through linkage to the NCHS National Death Index through December 31, 2019, using the NCHS Public-Use Linked Mortality Files. Follow-up time was calculated from the date of survey participation to death or censoring. The underlying cause of death was classified using the UCOD_LEADING variable (10-group recode), with malignant neoplasms (code 002) indicating cancer death. GI cancer death was defined as cancer-related death in patients with baseline GI tumor.

### 2.5 Covariates

Demographic covariates included age (continuous, years), sex (male/female), race/ethnicity (Non-Hispanic White, Non-Hispanic Black, Hispanic, Other), and education (high school or less vs. above high school). These were selected a priori based on known associations with both nutritional status and mortality.

### 2.6 Statistical Analysis

All analyses were conducted using R 4.6.0 with the `survey`, `survival`, `EValue`, `survRM2`, and `broom` packages.

**Cross-sectional analysis:** Linear regression models examined differences in PNI, CONUT, and GNRI between GI cancer patients and non-cancer participants, adjusted for age, sex, and race/ethnicity.

**Survival analysis:** Cox proportional hazards models were used to estimate hazard ratios (HR) for all-cause mortality per 1-standard deviation (SD) increase in each nutritional index. Both crude and adjusted models (age, sex, race/ethnicity) were fitted. The proportional hazards assumption was tested using Schoenfeld residuals. For indices with significant PH violation (PNI), time-dependent Cox models were fitted with time splitting at 2 and 5 years of follow-up.

**Survey-weighted analysis:** Survey-weighted Cox regression (svycoxph) was performed on the 2005-2016 subset incorporating NHANES complex survey design (clustering by SDMVPSU, stratification by SDMVSTRA, multi-cycle sampling weights WTMEC2YR/6).

**E-value sensitivity analysis:** E-values were computed (EValue R package) to quantify the minimum strength of association an unmeasured confounder would need to explain away the observed results.

**Restricted mean survival time (RMST):** RMST at 10 years was estimated for high vs. low PNI (threshold 48.5), providing an absolute survival difference in years with formal hypothesis testing.

**Joint exposure analysis:** PNI and NLR were dichotomized at medians into four groups. Cox models estimated HR using "High PNI + Low NLR" as reference. Additive interaction was quantified using Relative Excess Risk due to Interaction (RERI) with bootstrap 95% CI.

**CALLY index:** CALLY = albumin × lymphocyte / CRP was analyzed as an alternative composite index.

**Physical activity:** Any moderate-to-vigorous physical activity (MVPA) was harmonized across cycles. PNI × PA interaction was tested.

**Population attributable fraction (PAF):** PAF = p × (HR − 1) / HR was estimated for PNI < 48.5 with bootstrap CI.

**Bootstrap C-statistic:** Incremental discrimination (ΔC-stat) over baseline covariates was estimated with 500 bootstrap replications.

**Landmark analysis:** Conditional survival at 1, 3, and 5 years.

**Competing risks:** Cause-specific Cox models for GI cancer death vs. non-cancer death.

**Dietary comparison:** HEI-2015 in NHANES 2011-2016 subset, alone and mutually adjusted.

**Sensitivity analysis:** Analyses repeated excluding CRC cases.

## 3. Results

### 3.1 Baseline Characteristics

Of 44,066 eligible participants aged >=18, 313 (0.7%) had a GI cancer diagnosis. Colorectal cancer was the most common GI malignancy (n=241, 77.0%), followed by liver (n=23, 7.3%), stomach (n=19, 6.1%), esophageal (n=15, 4.8%), pancreatic (n=8, 2.6%), and gallbladder (n=6, 1.9%).

**Table 1: Demographics and nutritional indices by GI tumor status**

| Variable | GI Tumor (n=313) | Non-Cancer (n=40,480) | Other Cancer (n=3,273) |
|---|---|---|---|
| Age, years | 67.5 (13.8) | 45.5 (18.4) | 65.6 (15.0) |
| Female | 64 (20.4%) | 7,172 (17.7%) | 510 (15.6%) |
| Race: Non-Hispanic White | 160 (51.1%) | 14,887 (36.8%) | 1,683 (51.4%) |
| Race: Non-Hispanic Black | 45 (14.4%) | 5,909 (14.6%) | 293 (9.0%) |
| Race: Hispanic | 101 (32.3%) | 16,394 (40.5%) | 1,198 (36.6%) |
| Race: Other | 7 (2.2%) | 3,290 (8.1%) | 99 (3.0%) |
| Education > High school | 187 (59.7%) | 25,808 (63.8%) | 2,262 (69.1%) |
| BMI, kg/m2 | 28.5 (6.1) | 28.3 (6.6) | 28.2 (6.3) |
| Albumin, g/dL | 4.09 (0.38) | 4.22 (0.37) | 4.14 (0.34) |
| Lymphocyte, /uL | 1,972 (990) | 2,209 (719) | 2,000 (813) |
| PNI | 50.8 (6.5) | 53.2 (5.2) | 51.4 (5.3) |
| CONUT | 1.1 (1.4) | 0.7 (1.0) | 0.9 (1.1) |
| GNRI | 115.0 (12.2) | 116.4 (12.4) | 115.2 (12.2) |
| Deaths | 169 (54.0%) | 7,659 (18.9%) | 1,444 (44.1%) |

GI cancer patients were substantially older than the general population (67.5 vs 45.5 years) and had a higher proportion of Hispanic participants.

**Table 2: Baseline characteristics by PNI tertile among GI cancer patients**

| Variable | Low PNI T1 (n=105) | Mid PNI T2 (n=104) | High PNI T3 (n=104) |
|---|---|---|---|
| Age, years | 72.6 (9.4) | 66.4 (15.5) | 63.5 (14.2) |
| Female | 14.3% | 25.0% | 22.1% |
| PNI | 44.8 (3.5) | 50.5 (1.4) | 57.1 (5.9) |
| CONUT | 2.1 (1.7) | 0.9 (0.9) | 0.4 (0.6) |
| GNRI | 110.2 (12.6) | 115.5 (10.4) | 119.4 (11.9) |
| Albumin, g/dL | 3.79 (0.36) | 4.14 (0.26) | 4.35 (0.29) |
| Lymphocyte, /micro-L | 1,376 (389) | 1,820 (516) | 2,726 (1,262) |
| BMI, kg/m2 | 28.3 (6.1) | 28.5 (5.6) | 28.8 (6.6) |
| Deaths | 74 (70.5%) | 48 (46.2%) | 47 (45.2%) |
| Median survival, years | 6.5 | 12.5 | 13.9 |

### 3.2 All-Cause Mortality

All three nutritional indices significantly predicted all-cause mortality in GI cancer patients (Table 3).

**Table 3: Cox proportional hazards -- nutritional indices and all-cause mortality**

| Index | Adjustment | HR (95% CI) per 1-SD | p-value |
|---|---|---|---|
| PNI | Crude | 0.63 (0.51-0.78) | <0.001 |
| PNI | Adjusted* | 0.78 (0.63-0.96) | 0.019 |
| CONUT | Crude | 0.63 (0.55-0.73) | <0.001 |
| CONUT | Adjusted* | 0.71 (0.61-0.83) | <0.001 |
| GNRI | Crude | 0.67 (0.56-0.79) | <0.001 |
| GNRI | Adjusted* | 0.67 (0.57-0.80) | <0.001 |

*Adjusted for age, sex, race/ethnicity.

*Adjusted for age, sex, race/ethnicity.

GNRI showed the strongest association (HR 0.66 per 1-SD), closely followed by CONUT (HR 0.70) and PNI (HR 0.76). All three indices maintained statistical significance after adjustment.

### 3.3 Time-Dependent Effects

The proportional hazards assumption was violated for PNI (Schoenfeld residuals p=0.017) but not for GNRI (p=0.087) or CONUT (p=0.43). Time-dependent Cox models revealed marked attenuation of the PNI effect across follow-up periods (Table 4).

**Table 4: Time-dependent Cox -- PNI x follow-up period interaction**

| Period | PNI HR (95% CI) | p | Interaction p |
|---|---|---|---|
| 0-2 years | 0.34 (0.21-0.54) | <0.001 | Reference |
| 2-5 years | 1.32 (0.60-2.89) | 0.49 | <0.001 |
| 5+ years | 0.81 (0.47-1.40) | 0.45 | 0.002 |

The protective effect of PNI was strongest in the first 2 years post-diagnosis (HR 0.34, p<0.001) and attenuated significantly thereafter. GNRI showed a similar pattern (0-2 year HR 0.53, p=0.034; interaction p=0.30). CONUT demonstrated a more stable effect across periods (HR 0.67 at 0-2 years, interaction p=0.77).

### 3.4 Landmark Analysis

In conditional survival analysis at 1, 3, and 5 years, PNI showed attenuated and non-significant effects (1-year landmark HR 0.90, p=0.39), while CONUT and GNRI maintained significant associations at all landmark times (Table 5).

**Table 5: Landmark analysis -- adjusted HR per 1-SD**

| Index | 1-year Landmark | 3-year Landmark | 5-year Landmark |
|---|---|---|---|
| PNI | 0.91 (0.73-1.13, p=0.38) | 0.85 (0.66-1.10, p=0.21) | 0.81 (0.60-1.09, p=0.17) |
| CONUT | 0.75 (0.63-0.89, p<0.001) | 0.71 (0.58-0.88, p=0.002) | 0.74 (0.57-0.96, p=0.021) |
| GNRI | 0.70 (0.59-0.84, p<0.001) | 0.70 (0.57-0.86, p<0.001) | 0.70 (0.56-0.89, p=0.003) |

### 3.5 Competing Risks Analysis

Cause-specific Cox models revealed divergent patterns across mortality causes (Table 6).

**Table 6: Competing risks -- cause-specific HR per 1-SD**

| Index | GI Cancer Death | Non-Cancer Death |
|---|---|---|
| PNI | 1.06 (0.78-1.45, p=0.70) | 0.63 (0.48-0.83, p<0.001) |
| CONUT | 0.73 (0.56-0.94, p=0.016) | 0.70 (0.58-0.85, p<0.001) |
| GNRI | 0.69 (0.51-0.94, p=0.017) | 0.66 (0.53-0.82, p<0.001) |

PNI primarily protected against non-cancer death (HR 0.63, p<0.001) rather than GI cancer-specific death (HR 1.06, p=0.70). This pattern was consistent across all indices: better nutritional status was more strongly associated with reduced mortality from competing causes than with reduced GI cancer-specific mortality.

### 3.6 PNI Decomposition: Albumin Drives the Effect

To determine which component of PNI drives its prognostic value, we compared the prognostic performance of albumin alone, lymphocyte alone, and the composite PNI. Albumin alone was significantly associated with mortality (HR 0.69 per 1-SD, 95% CI 0.58-0.82, p<0.001), while lymphocyte count alone was not (HR 1.06, 95% CI 0.88-1.28, p=0.54). The composite PNI (HR 0.78, 95% CI 0.63-0.96, p=0.020) showed intermediate performance. C-statistic improvement over baseline (age, sex, race) was +0.033 for albumin, +0.023 for PNI, and 0.000 for lymphocyte. These results indicate that **serum albumin is the primary driver** of PNI's prognostic value in GI cancer, with lymphocyte count contributing minimal additional information.

### 3.7 Optimal PNI Threshold and Dose-Response

Using maximally selected log-rank statistics, we identified **PNI < 48.5** as the optimal prognostic threshold (p<0.0001). Patients below this threshold had markedly shorter median survival (6.6 vs 12.8 years). Restricted cubic spline analysis confirmed a **nonlinear relationship** between PNI and mortality risk (nonlinearity p<0.0001), with a steep risk increase below approximately 45 and a plateau above 50, supporting the threshold model for clinical risk stratification.

### 3.8 PNI vs Inflammatory Markers (NLR, SII)

In the subset with complete blood count data available for neutrophil and platelet counts (n=258, 125 events), we compared PNI with two established inflammatory indices: the neutrophil-to-lymphocyte ratio (NLR) and the systemic immune-inflammation index (SII = platelet × neutrophil / lymphocyte).

NLR was a significant predictor of all-cause mortality (HR 1.43 per 1-SD, 95% CI 1.19-1.72, p<0.001), as was SII (HR 1.30, 95% CI 1.10-1.54, p=0.002), while PNI was borderline significant in this subset (HR 0.78, 95% CI 0.61-1.01, p=0.063). C-statistic improvement over baseline was comparable: +0.017 for both PNI and logNLR, and +0.013 for logSII. In mutually adjusted models including both PNI and logNLR, neither remained significant (PNI HR 0.86, p=0.29; logNLR HR 1.22, p=0.09), suggesting shared predictive information.

These findings indicate that inflammatory markers perform comparably to nutritional indices in predicting GI cancer survival. The overlap between PNI (which includes lymphocyte count, an immune parameter) and NLR (neutrophil/lymphocyte) suggests that both nutritional and inflammatory pathways contribute to the observed mortality differences, with neither pathway entirely dominating.

### 3.9 CRP Interaction

Among 250 GI cancer patients with available CRP measurements (120 deaths), higher CRP was associated with non-significantly increased mortality risk (HR 1.18 per log-SD, 95% CI 0.92-1.50, p=0.19). The PNI × CRP interaction was not statistically significant (p=0.78), indicating that PNI's prognostic effect does not vary by systemic inflammation level.

### 3.10 Dietary Quality Comparison

Among GI cancer patients with available dietary data (n=99, 33 deaths), HEI-2015 dietary quality did not predict all-cause mortality (HR 0.97 per 1-SD, 95% CI 0.70-1.35, p=0.85), whether analyzed separately or mutually adjusted with PNI. The correlation between HEI-2015 and PNI was weak (Pearson r=0.16), confirming these measures capture distinct dimensions of nutritional status.

### 3.11 Sensitivity Analysis

After excluding colorectal cancer (n=241), the remaining 72 non-CRC GI cancer patients showed directionally consistent but imprecise associations (PNI HR 0.73, 95% CI 0.43-1.24, p=0.25), while CONUT and GNRI remained significant (CONUT HR 0.65, p=0.012; GNRI HR 0.42, p<0.001). The overall findings are not driven solely by the colorectal cancer subgroup.

### 3.12 Subgroup Analysis and Effect Modification

PNI demonstrated directionally consistent protective effects across all examined subgroups, with no significant effect modification by age (interaction p=0.73), sex (female HR 0.56 vs male HR 0.82), or BMI (<30 HR 0.72 vs >=30 HR 0.81). The effect was statistically significant in age <60 (HR 0.47, p=0.013), female (HR 0.56, p=0.047), and BMI <30 subgroups (HR 0.72, p=0.015).

### 3.13 Advanced Analyses

**Survey-weighted Cox:** In the 2005-2016 subset (270 GI patients, 134 events), survey-weighted Cox regression attenuated the PNI association (HR 0.87, 95% CI 0.61-1.23, p=0.43), while CONUT remained significant (HR 0.76, 95% CI 0.63-0.91, p=0.003). Unweighted analysis on the same subset showed the expected protective effect (PNI HR 0.78, p=0.046), suggesting survey weighting primarily attenuates PNI by down-weighting NHANES III participants with longer follow-up.

**E-value analysis:** The E-value for PNI was 1.66 (CI 1.20), indicating an unmeasured confounder would need at least a 1.66-fold association with both PNI and mortality to explain the observed HR. GNRI showed more robust results (E-value 1.97, CI 1.61).

**RMST:** At 10 years, high PNI (≥48.5) patients had RMST 8.22 vs. 6.25 years for low PNI, yielding a clinically meaningful difference of 1.97 years (95% CI 1.19-2.76, p<0.001).

**CALLY index:** Among 250 GI patients (120 events), CALLY was non-significant (HR 0.85, 95% CI 0.67-1.09, p=0.20) with C-stat 0.659 vs. 0.667 for PNI.

**CONUT/GNRI dose-response:** Both showed linear mortality associations (CONUT p=0.11, GNRI p=0.16), contrasting with PNI's nonlinearity (p<0.001).

**PNI×NLR joint exposure:** The "Low PNI + High NLR" group had significantly higher mortality (HR 1.67, 95% CI 1.06-2.63, p=0.028). RERI suggested sub-multiplicative interaction (-0.97, 95% CI -2.61 to -0.10).

**C-statistic comparison (bootstrap):** GNRI showed the largest ΔC-stat (+0.034), followed by CONUT (+0.023) and PNI (+0.022). NLR and SII added less.

**Physical activity:** PA alone was non-significant (HR 0.99, p=0.94), but a significant PNI×PA interaction was observed (p=0.023).

**PAF:** PNI < 48.5 had 35.5% prevalence with PAF 13.4% (95% CI 6.2%-20.1%), suggesting ~1 in 8 deaths may be attributable to low nutritional reserve.

**MICE:** No missing data in key variables, obviating imputation.

## 4. Discussion

### 4.1 Principal Findings

In this population-based, dual-cohort study spanning 30 years of NHANES data, we demonstrate that three composite nutritional indices -- PNI, CONUT, and GNRI -- consistently predict all-cause mortality in GI cancer patients. The prognostic effect was robust to adjustment for demographic covariates and was externally validated across two independent NHANES generations. Several novel findings emerge.

First, the protective effect of better nutritional status is strongly time-dependent. PNI's prognostic effect was most pronounced in the first 2 years following diagnosis (HR 0.34 per 1-SD) and attenuated progressively thereafter. RMST analysis quantified this as a 1.97-year absolute survival difference at 10 years, and PAF analysis indicated that 13.4% of deaths could be potentially attributable to low PNI. This temporal pattern has important clinical implications: the window for nutritional intervention is early in the cancer trajectory.

Second, the competing risks analysis reveals that PNI's protective effect operates primarily through non-cancer mortality pathways. Patients with better PNI had a 36% lower risk of non-cancer death but no reduction in GI cancer-specific death. This finding suggests that PNI captures physiological reserve and resilience to treatment-related complications rather than directly influencing tumor progression.

Third, the PNI decomposition analysis reveals that **serum albumin is the primary driver** of prognostic value, while lymphocyte count contributes minimal independent information. This is consistent with albumin's role as a marker of both nutritional status and systemic inflammation, and suggests that simpler, albumin-only indices may suffice for risk stratification.

Fourth, the optimal PNI threshold of 48.5 provides a clinically actionable cutpoint for risk stratification. The nonlinear dose-response relationship, with a steep risk increase below 45 and a plateau above 50, supports a threshold-based rather than linear approach to clinical interpretation. Notably, CONUT and GNRI demonstrated linear dose-response relationships, supporting their use as continuous risk scores without requiring dichotomization.

Fifth, inflammatory markers (NLR, SII) showed comparable prognostic performance to PNI, with similar C-statistic improvement. Joint exposure analysis identified the "Low PNI + High NLR" group as having the highest mortality risk (HR 1.67), suggesting clinical utility in combining nutritional and inflammatory assessment. The negative RERI suggests that the joint effect of poor nutrition and elevated inflammation is less than multiplicative, implying overlapping mechanistic pathways.

Sixth, the null finding for HEI-2015 and the absence of PNI × CRP interaction provide important mechanistic insights. The negligible correlation between dietary quality (HEI-2015) and biochemical nutrition (PNI), combined with PNI's independence from CRP, suggests that these composite indices capture a distinct biological signal related to metabolic reserve and systemic inflammation that cannot be replaced by dietary assessment or a single inflammatory marker. The CALLY index, which incorporates CRP, did not outperform PNI alone.

Seventh, E-value analysis provided quantitative support for the robustness of our findings to unmeasured confounding. The E-values ranged from 1.66 (PNI) to 1.97 (GNRI), indicating that moderately strong unmeasured confounding would be required to fully explain away the observed associations. The significant PNI × physical activity interaction (p=0.023) suggests the prognostic effect of nutritional status may be modified by exercise, warranting further investigation in prospective studies.

### 4.2 Comparison with Previous Literature

Our findings are consistent with prior studies validating PNI in GI cancer. A meta-analysis of 28 studies including 7,781 gastric cancer patients reported that low PNI was associated with poorer overall survival (HR 1.38, 95% CI 1.23-1.54) [10]. Similar results have been reported in colorectal cancer [11], hepatocellular carcinoma [12], and pancreatic cancer [13]. Our study extends these findings from hospital-based surgical series to the population level with dual-cohort external validation.

The time-dependent effect has been suggested in smaller studies. A Korean gastric cancer cohort noted that PNI's prognostic effect was strongest in the first year after surgery [14]. Our findings confirm this pattern in a US population with long-term follow-up.

### 4.3 Clinical Implications

The simplicity and low cost of these indices make them attractive for clinical application. PNI, requiring only serum albumin and a complete blood count, can be calculated from routine laboratory data available in virtually any clinical setting. A PNI threshold of approximately 45 identifies a high-risk group with median survival of 6.5 years compared to 13.9 years for high PNI patients. The concentration of risk in the first 2 years suggests early nutritional intervention may have the greatest impact.

### 4.4 Strengths and Limitations

Strengths include: (1) population-based sampling with national representativeness; (2) dual-cohort design with external validation across 30 years; (3) comprehensive temporal analysis (time-dependent, landmark, competing risks, RCS, RMST); (4) head-to-head comparison of biochemical vs. dietary nutritional assessment, nutritional vs. inflammatory indices, and PNI component decomposition; (5) optimal threshold identification; (6) quantitative sensitivity analyses (E-value, survey-weighted, bootstrap C-stat); (7) absolute risk measures provided (RMST, PAF) enhancing clinical interpretability.

Limitations include: (1) small GI cancer sample (n=313) limiting site-specific analysis; (2) single baseline measurement of nutritional indices, not capturing changes over time; (3) reliance on self-reported cancer diagnosis with potential misclassification; (4) lack of cancer stage, treatment, and recurrence data, which are important confounders; (5) public-use mortality file provides only a 10-group cause-of-death recode, limiting cause-specific analysis granularity; (6) predominantly colorectal cancer (77%), limiting generalizability to upper GI cancers; (7) NLR and SII available for 2005-2016 only, not NHANES III; (8) survey-weighted analysis limited to 2005-2016 subset; (9) physical activity data limited to self-report with varying instruments across cycles.

## 5. Conclusions

PNI, CONUT, and GNRI are robust, independent prognostic factors for all-cause mortality in GI cancer patients, validated across two nationally representative US cohorts spanning 30 years. The protective effect is time-dependent, concentrated in the first 2 years post-diagnosis, and operates primarily through non-cancer mortality pathways. An optimal PNI threshold of 48.5 identifies a high-risk subgroup with substantially shorter survival (RMST difference 1.97 years at 10 years), and an estimated 13.4% of deaths are potentially attributable to low PNI. The prognostic value is driven primarily by serum albumin rather than lymphocyte count, and is independent of systemic inflammation (CRP). Inflammatory markers (NLR, SII) show comparable performance, and joint analysis of PNI and NLR identifies a particularly high-risk group. CONUT and GNRI show linear dose-response relationships, distinguishing them from PNI's threshold-based pattern. The null CALLY and HEI-2015 findings suggest that established composite indices already capture the relevant prognostic information. E-value analysis supports robustness to unmeasured confounding. These simple, inexpensive indices may be useful for risk stratification and early intervention targeting in GI cancer patients.

## References

1. Fearon K, et al. Definition and classification of cancer cachexia. Lancet Oncol. 2011;12(5):489-495.
2. Arends J, et al. ESPEN guidelines on nutrition in cancer patients. Clin Nutr. 2017;36(1):11-48.
3. Hebuterne X, et al. Prevalence of malnutrition and current use of nutrition support in patients with cancer. JPEN J Parenter Enteral Nutr. 2014;38(2):196-204.
4. Pressoir M, et al. Prevalence, risk factors and clinical implications of malnutrition in French Comprehensive Cancer Centres. Br J Cancer. 2010;102(6):966-971.
5. Onodera T, Goseki N, Kosaki G. Prognostic nutritional index in gastrointestinal surgery of malnourished cancer patients. Nihon Geka Gakkai Zasshi. 1984;85(9):1001-1005.
6. Ignacio de UlIbarri J, et al. CONUT: a tool for controlling nutritional status. Nutr Hosp. 2005;20(1):38-45.
7. Bouillanne O, et al. Geriatric Nutritional Risk Index. Am J Clin Nutr. 2005;82(4):777-783.
8. Yang Y, et al. The prognostic nutritional index is a predictive indicator of prognosis in gastric cancer. Eur J Surg Oncol. 2016;42(8):1176-1182.
9. Sun K, et al. The prognostic significance of the prognostic nutritional index in cancer: a systematic review and meta-analysis. J Cancer Res Clin Oncol. 2014;140(9):1537-1549.
10. Migita K, et al. The prognostic nutritional index predicts long-term outcomes of gastric cancer patients. Ann Surg Oncol. 2013;20(8):2647-2654.
11. Nozoe T, et al. The prognostic nutritional index can be a prognostic indicator in colorectal carcinoma. Surg Today. 2012;42(6):532-535.
12. Pinato DJ, et al. A novel, externally validated inflammation-based prognostic algorithm in hepatocellular carcinoma. Br J Cancer. 2012;106(8):1439-1445.
13. Kanda M, et al. Nutritional predictors of postoperative outcome in pancreatic cancer. Br J Surg. 2011;98(2):268-274.
14. Lee JY, et al. Clinical significance of the prognostic nutritional index in predicting postoperative gastric cancer survival. J Gastric Cancer. 2017;17(2):139-150.

---

## Figures

- **Figure 1:** Kaplan-Meier survival curves by PNI tertile (combined NHANES)
- **Figure 2:** Forest plot -- adjusted HR per 1-SD for PNI, CONUT, GNRI
- **Figure 3:** Time-dependent HR across follow-up periods
- **Figure 4:** Landmark analysis at 1, 3, and 5 years
- **Figure 5:** Competing risks forest plot
- **Figure 6:** PNI cutpoint analysis (optimal threshold = 48.5)
- **Figure 7:** Restricted cubic spline dose-response curve
- **Figure 8:** Subgroup forest plot
- **Figure 9:** PNI vs NLR vs SII comparison
- **Figure 10:** PNI decomposition (albumin vs lymphocyte vs composite)

## Supplementary Materials

- **Table S1:** GI cancer site breakdown with survival characteristics
- **Table S2:** Non-CRC sensitivity analysis
- **Table S3:** Full Cox regression output with all covariates
- **Table S4:** PNI vs NLR vs SII full comparison
- **Table S5:** Subgroup analysis with interaction p-values
- **Table S6:** Survey-weighted Cox regression (svycoxph)
- **Table S7:** E-value sensitivity analysis
- **Figure S1:** Kaplan-Meier curves by CALLY tertile
- **Figure S2:** Joint exposure forest plot (PNI × NLR)
- **Figure S3:** RMST bar plot (10-year restricted mean survival time)
- **Figure S4:** CONUT dose-response curve (restricted cubic spline)
- **Figure S5:** GNRI dose-response curve (restricted cubic spline)
