## Abstract

**Background & Aims:** Nutritional status is a modifiable prognostic factor in cancer patients, but population-level validation of composite nutritional indices in gastrointestinal (GI) cancer remains limited. We aimed to validate the Prognostic Nutritional Index (PNI), Controlling Nutritional Status (CONUT), and Geriatric Nutritional Risk Index (GNRI) as prognostic factors for all-cause mortality in GI cancer using nationally representative data.

**Methods:** We identified 313 GI cancer patients from NHANES III (1988-1994, n=92) and Continuous NHANES (2005-2016, n=221) with mortality follow-up through December 2019 (up to 31 years). PNI, CONUT, and GNRI were derived from baseline serum albumin, lymphocyte count, total cholesterol, and body mass index. Cox proportional hazards, time-dependent Cox, landmark analysis, competing risks, restricted cubic splines, restricted mean survival time (RMST), and population attributable fraction (PAF) analyses were performed. E-value sensitivity analysis quantified robustness to unmeasured confounding.

**Results:** Higher PNI was associated with reduced all-cause mortality (adjusted HR 0.78, 95% CI 0.63-0.96, p=0.020). GNRI showed the strongest association (HR 0.67, 95% CI 0.57-0.80, p<0.001). The PNI effect was strongly time-dependent, maximal in the first 2 years (HR 0.34, 95% CI 0.21-0.54, p<0.001) and attenuating thereafter (interaction p<0.001). Optimal PNI threshold was 48.5; patients below this had median survival 6.6 vs 12.8 years (p<0.0001) with RMST difference 1.97 years (95% CI 1.19-2.76, p<0.001) and PAF 13.4% (95% CI 6.2%-20.1%). In competing risks, PNI primarily protected against non-cancer death (HR 0.63, p=0.001) rather than GI cancer-specific death (HR 1.06, p=0.70). PNI decomposition showed albumin drove the effect (HR 0.69, p<0.001) while lymphocyte count alone was non-significant. Dietary quality (HEI-2015) did not predict survival (HR 0.97, p=0.85). Inflammatory markers (NLR, SII) showed comparable prognostic performance to PNI.

**Conclusions:** PNI, CONUT, and GNRI are robust, independent prognostic factors for all-cause mortality in GI cancer patients, validated across two nationally representative US cohorts spanning 30 years. The protective effect is time-dependent, operates primarily through non-cancer mortality pathways, and is driven by serum albumin. These simple, inexpensive indices may guide early nutritional intervention in GI cancer patients.

**Keywords:** Prognostic Nutritional Index; Gastrointestinal Cancer; NHANES; Survival Analysis; Competing Risks; Time-Dependent Effect

**Abbreviations:** PNI, Prognostic Nutritional Index; CONUT, Controlling Nutritional Status; GNRI, Geriatric Nutritional Risk Index; GI, gastrointestinal; NHANES, National Health and Nutrition Examination Survey; NLR, neutrophil-to-lymphocyte ratio; SII, systemic immune-inflammation index; CRP, C-reactive protein; HEI-2015, Healthy Eating Index 2015; RMST, restricted mean survival time; PAF, population attributable fraction; RCS, restricted cubic spline; RERI, relative excess risk due to interaction; HR, hazard ratio; CI, confidence interval; NCHS, National Center for Health Statistics.

## 1. Introduction

Cancer cachexia and malnutrition represent critical challenges in gastrointestinal oncology. GI cancer patients are uniquely vulnerable due to tumor-induced catabolism, malabsorption, and mechanical obstruction, with malnutrition prevalence ranging from 30% to 80% depending on tumor site and stage [1,2]. Malnourished patients consistently show poorer treatment tolerance, increased postoperative complications, and reduced survival [3,4].

Several composite nutritional indices have been developed to integrate multiple biomarkers into a single prognostic score. The Prognostic Nutritional Index (PNI), originally proposed by Onodera et al. in 1984, combines serum albumin and peripheral lymphocyte count [5]. The Controlling Nutritional Status (CONUT) score incorporates albumin, lymphocyte count, and total cholesterol [6]. The Geriatric Nutritional Risk Index (GNRI) adjusts albumin for body mass index [7]. Each index has been validated primarily in single-center surgical series, predominantly in East Asian populations with evidence concentrated in specific GI cancer subtypes [8,9].

Despite their clinical utility, important gaps remain. First, population-level validation across different geographic contexts and survey eras is lacking. Most prior studies originate from single-institution cohorts with sample sizes rarely exceeding 1,000 patients [10,11]. Second, the temporal dynamics of nutritional risk — whether the prognostic effect is constant or concentrated in specific post-diagnosis windows — have not been systematically examined using time-dependent analyses. Third, the conceptual distinction between biochemical nutritional status (serum proteins, lymphocyte count) and dietary nutritional quality (nutrient intake patterns) remains poorly characterized in cancer prognosis.

To address these gaps, we undertook a comprehensive analysis of nutritional indices in GI cancer patients using the National Health and Nutrition Examination Survey (NHANES), a nationally representative population-based survey with mortality linkage. We leveraged two independent NHANES generations — NHANES III (1988-1994) and Continuous NHANES (2005-2016) — providing up to 31 years of mortality follow-up and enabling cross-validation across different survey eras. Our objectives were to: (1) validate PNI, CONUT, and GNRI as prognostic factors for all-cause mortality; (2) characterize time-dependent effects through landmark, time-dependent Cox, and RMST analyses; (3) compare prognostic performance across competing causes of death; and (4) compare biochemical nutritional indices against a dietary quality measure (HEI-2015) and inflammatory markers (NLR, SII).

## 2. Materials and Methods

### 2.1 Study Population

We used data from two independent NHANES cycles: NHANES III (1988-1994) and Continuous NHANES (2005-2016). NHANES is a nationally representative cross-sectional survey of the civilian non-institutionalized US population conducted by the National Center for Health Statistics (NCHS). Both surveys collected demographic, dietary, laboratory, and examination data through standardized protocols. The analytic sample included participants aged ≥18 years with available nutrition biomarker data. A total of 44,066 participants met inclusion criteria, of whom 313 had a GI cancer diagnosis.

### 2.2 GI Cancer Ascertainment

GI cancer was identified from self-reported cancer history. In Continuous NHANES (2005-2016), cancer diagnoses were identified through the Medical Conditions Questionnaire (MCQ) using cancer site codes. In NHANES III, cancer site was identified from variable HAC3OS, with GI sites defined as codes for colon, stomach, esophagus, pancreas, liver, and gallbladder. GI cancer was defined as cancer of the esophagus, stomach, colon, rectum, pancreas, liver, or gallbladder.

### 2.3 Nutritional Indices

**PNI** was calculated according to Onodera's formula: PNI = 10 × serum albumin (g/dL) + 0.005 × lymphocyte count (/μL) [5].

**CONUT** was derived as the sum of three component scores: albumin (≥3.5 g/dL = 0, 3.0-3.49 = 2, 2.5-2.99 = 4, <2.5 = 6), lymphocyte count (≥1600/μL = 0, 1200-1599 = 1, 800-1199 = 2, <800 = 3), and total cholesterol (≥180 mg/dL = 0, 140-179 = 1, 100-139 = 2, <100 = 3) [6].

**GNRI** was calculated as: GNRI = 14.89 × albumin (g/dL) + 41.7 × (BMI / 22) [7].

Serum albumin was measured using the bromcresol purple method in Continuous NHANES and bromcresol green in NHANES III, with g/L values converted to g/dL for consistency. Complete blood count provided lymphocyte counts. Total cholesterol was measured enzymatically.

**HEI-2015** was calculated from 24-hour dietary recall data using the Food Patterns Equivalents Database, measuring adherence to the 2015-2020 Dietary Guidelines for Americans on a 0-100 scale.

### 2.4 Mortality Linkage

Mortality status was ascertained through linkage to the NCHS National Death Index through December 31, 2019. Follow-up time was calculated from the date of survey participation to death or censoring. The underlying cause of death was classified using the UCOD_LEADING variable, with malignant neoplasms (code 002) indicating cancer death.

### 2.5 Covariates

Demographic covariates included age (continuous, years), sex (male/female), and race/ethnicity (Non-Hispanic White, Non-Hispanic Black, Hispanic, Other), selected a priori based on known associations with both nutritional status and mortality.

### 2.6 Statistical Analysis

All analyses were conducted using R 4.6.0 with the `survey`, `survival`, `EValue`, `survRM2`, and `broom` packages.

**Survival analysis:** Cox proportional hazards models estimated hazard ratios (HR) per 1-standard deviation (SD) increase in each nutritional index. Both crude and adjusted models (age, sex, race/ethnicity) were fitted. The proportional hazards assumption was tested using Schoenfeld residuals. For PNI, which violated this assumption (p=0.017), time-dependent Cox models were fitted with time splitting at 2 and 5 years.

**Time-dependent effects:** Landmark analyses at 1, 3, and 5 years examined conditional survival. Restricted mean survival time (RMST) at 10 years was compared between high and low PNI groups using the threshold identified by maximally selected log-rank statistics. Restricted cubic splines (3 knots) examined nonlinear dose-response relationships.

**Competing risks:** Cause-specific Cox models were fitted for GI cancer death vs. non-cancer death.

**Component analysis:** PNI was decomposed into albumin alone and lymphocyte alone, comparing their prognostic performance to the composite index.

**Joint exposure analysis:** PNI and NLR were dichotomized at medians into four groups, with additive interaction quantified using Relative Excess Risk due to Interaction (RERI) with bootstrap 95% CI.

**Population attributable fraction (PAF):** PAF = p × (HR − 1) / HR was estimated for PNI < 48.5 with bootstrap CI.

**E-value sensitivity analysis:** E-values quantified the minimum association strength an unmeasured confounder would need to fully explain away the observed results.

**Bootstrap C-statistic:** Incremental discrimination (ΔC-stat) over baseline covariates was estimated with 500 bootstrap replications.

**Inflammatory markers:** NLR and SII were compared with PNI in the 2005-2016 subset with available neutrophil and platelet counts.

## 3. Results

### 3.1 Baseline Characteristics

Of 44,066 eligible participants aged ≥18, 313 (0.7%) had a GI cancer diagnosis. Colorectal cancer was the most common (n=241, 77.0%), followed by liver (n=23, 7.3%), stomach (n=19, 6.1%), esophageal (n=15, 4.8%), pancreatic (n=8, 2.6%), and gallbladder (n=6, 1.9%). The remaining cases were rectal (n=8). GI cancer patients were substantially older than the general population (67.5 vs 45.5 years) and had lower PNI (50.8 vs 53.2), indicating worse nutritional status. During a median follow-up of 8.8 years (up to 31 years), 169 deaths (54.0%) occurred among GI cancer patients.

**Table 1 summarizes baseline characteristics by GI tumor status**, and **Table 2 shows characteristics by PNI tertile among GI cancer patients.** Patients in the lowest PNI tertile had substantially higher mortality (70.5%) compared to the middle (46.2%) and highest (45.2%) tertiles, with median survival of 6.5, 12.5, and 13.9 years, respectively.

### 3.2 All-Cause Mortality

All three nutritional indices significantly predicted all-cause mortality after adjustment for age, sex, and race/ethnicity (Table 3). GNRI showed the strongest association (HR 0.67 per 1-SD, 95% CI 0.57-0.80, p<0.001), followed by CONUT (HR 0.71, 95% CI 0.61-0.83, p<0.001) and PNI (HR 0.78, 95% CI 0.63-0.96, p=0.020).

**Table 3.** Cox proportional hazards — nutritional indices and all-cause mortality.

| Index | Adjustment | HR (95% CI) per 1-SD | p-value |
|-------|-----------|----------------------|---------|
| PNI | Crude | 0.63 (0.51-0.78) | <0.001 |
| PNI | Adjusted* | 0.78 (0.63-0.96) | 0.019 |
| CONUT | Crude | 0.63 (0.55-0.73) | <0.001 |
| CONUT | Adjusted* | 0.71 (0.61-0.83) | <0.001 |
| GNRI | Crude | 0.67 (0.56-0.79) | <0.001 |
| GNRI | Adjusted* | 0.67 (0.57-0.80) | <0.001 |

*Adjusted for age, sex, race/ethnicity.

### 3.3 Time-Dependent Effects

The proportional hazards assumption was violated for PNI (Schoenfeld residuals p=0.017) but not for GNRI (p=0.087) or CONUT (p=0.43). Time-dependent Cox models revealed marked attenuation of the PNI effect: the protective association was strongest in the first 2 years (HR 0.34, 95% CI 0.21-0.54, p<0.001) and non-significant thereafter (2-5 year: HR 1.32, p=0.49; 5+ year: HR 0.81, p=0.45; interaction p<0.001). GNRI showed a similar temporal pattern (0-2 year HR 0.53, p=0.034), while CONUT demonstrated more stable effects across periods.

In landmark analysis, PNI showed attenuated and non-significant effects at 1, 3, and 5 years, while CONUT and GNRI maintained significant associations at all landmark times (Table 4).

**Table 4.** Landmark analysis — adjusted HR per 1-SD.

| Index | 1-year Landmark | 3-year Landmark | 5-year Landmark |
|-------|----------------|----------------|----------------|
| PNI | 0.91 (0.73-1.13, p=0.38) | 0.85 (0.66-1.10, p=0.21) | 0.81 (0.60-1.09, p=0.17) |
| CONUT | 0.75 (0.63-0.89, p<0.001) | 0.71 (0.58-0.88, p=0.002) | 0.74 (0.57-0.96, p=0.021) |
| GNRI | 0.70 (0.59-0.84, p<0.001) | 0.70 (0.57-0.86, p<0.001) | 0.70 (0.56-0.89, p=0.003) |

### 3.4 Optimal PNI Threshold and Dose-Response

Maximally selected log-rank statistics identified **PNI < 48.5** as the optimal prognostic threshold (p<0.0001). Patients below this threshold had markedly shorter median survival (6.6 vs 12.8 years). RMST analysis at 10 years confirmed a clinically meaningful difference of 1.97 years between groups (95% CI 1.19-2.76, p<0.001), with PAF indicating 13.4% of deaths attributable to low PNI (95% CI 6.2%-20.1%).

Restricted cubic spline analysis confirmed a nonlinear relationship between PNI and mortality risk (nonlinearity p<0.0001), with sharp risk increase below approximately 45 and a plateau above 50. In contrast, CONUT and GNRI showed linear dose-response relationships (nonlinearity p=0.11 and p=0.16, respectively).

### 3.5 Competing Risks Analysis

Cause-specific Cox models revealed divergent patterns. PNI primarily protected against non-cancer death (HR 0.63, 95% CI 0.48-0.83, p<0.001) rather than GI cancer-specific death (HR 1.06, 95% CI 0.78-1.45, p=0.70). CONUT and GNRI showed protective effects against both cancer-specific and non-cancer death (Table 5).

**Table 5.** Competing risks — cause-specific HR per 1-SD.

| Index | GI Cancer Death | Non-Cancer Death |
|-------|----------------|-----------------|
| PNI | 1.06 (0.78-1.45, p=0.70) | 0.63 (0.48-0.83, p<0.001) |
| CONUT | 0.73 (0.56-0.94, p=0.016) | 0.70 (0.58-0.85, p<0.001) |
| GNRI | 0.69 (0.51-0.94, p=0.017) | 0.66 (0.53-0.82, p<0.001) |

### 3.6 PNI Decomposition

Albumin alone was significantly associated with mortality (HR 0.69 per 1-SD, 95% CI 0.58-0.82, p<0.001), while lymphocyte count alone was not (HR 1.06, 95% CI 0.88-1.28, p=0.54). C-statistic improvement over baseline was +0.033 for albumin, +0.023 for PNI, and 0.000 for lymphocyte, confirming that **serum albumin is the primary driver** of PNI's prognostic value.

### 3.7 PNI vs Inflammatory Markers

In the subset with neutrophil and platelet counts (n=258, 125 events), NLR (HR 1.43 per 1-SD, 95% CI 1.19-1.72, p<0.001) and SII (HR 1.30, 95% CI 1.10-1.54, p=0.002) significantly predicted mortality. C-statistic improvement was comparable across indices (+0.017 for PNI, +0.017 for logNLR). In mutually adjusted models including both PNI and logNLR, neither remained significant (PNI HR 0.86, p=0.29; logNLR HR 1.22, p=0.09), indicating shared predictive information.

### 3.8 Joint PNI × NLR Analysis

The "Low PNI + High NLR" group had significantly higher mortality (HR 1.67, 95% CI 1.06-2.63, p=0.028) compared to the "High PNI + Low NLR" reference group. RERI was -0.97 (95% CI -2.61 to -0.10), suggesting overlapping rather than synergistic pathways.

### 3.9 CRP Interaction and CALLY Index

Among 250 patients with CRP measurements, CRP was non-significantly associated with mortality (HR 1.18 per log-SD, p=0.19). No PNI × CRP interaction was observed (p=0.78). The CALLY index (albumin × lymphocyte / CRP) did not outperform PNI (CALLY HR 0.85, p=0.20; C-stat 0.659 vs 0.667 for PNI).

### 3.10 Dietary Quality Comparison

Among 99 GI cancer patients with dietary data, HEI-2015 did not predict mortality (HR 0.97, p=0.85). The correlation between HEI-2015 and PNI was weak (r=0.16), confirming these measures capture distinct dimensions of nutritional status.

### 3.11 Subgroup and Sensitivity Analyses

PNI showed directionally consistent protective effects across all examined subgroups, with no significant effect modification by age (interaction p=0.73), sex (p=0.34), or BMI (p=0.68). Bootstrap C-statistic (500 replications) showed GNRI had the largest ΔC-stat (+0.034), followed by CONUT (+0.023) and PNI (+0.022). Survey-weighted Cox (2005-2016 subset) attenuated PNI (HR 0.87, p=0.43) while CONUT remained significant (HR 0.76, p=0.003). E-values ranged from 1.66 (PNI) to 1.97 (GNRI), indicating moderate robustness to unmeasured confounding.

## 4. Discussion

### Summary of Principal Findings

In this population-based dual-cohort study spanning 30 years of NHANES data, we demonstrate that three composite nutritional indices — PNI, CONUT, and GNRI — consistently predict all-cause mortality in GI cancer patients, with GNRI showing the strongest association (HR 0.67 per 1-SD, 95% CI 0.57-0.80) followed by CONUT (HR 0.71, 95% CI 0.61-0.83) and PNI (HR 0.78, 95% CI 0.63-0.96). The PNI effect was strongly time-dependent, maximal in the first 2 years (HR 0.34), and operated primarily through non-cancer mortality rather than GI cancer-specific pathways. RMST analysis quantified a clinically meaningful 1.97-year survival difference between high and low PNI groups at 10 years, with a PAF of 13.4% attributable to PNI below the optimal threshold of 48.5. E-value sensitivity analysis supported moderate robustness to unmeasured confounding (E-values 1.66-1.97).

### Time-Dependent and Cause-Specific Effects

The temporal variation of PNI's prognostic effect — maximal in the first 2 years (HR 0.34) followed by marked attenuation — has important clinical implications. Pre-existing nutritional status most strongly influences short-term survival, likely through its impact on treatment tolerance, perioperative recovery, and therapy completion. This pattern aligns with prior observations in a Korean gastric cancer cohort [14] and suggests that the window for effective nutritional intervention is early in the cancer trajectory. After the initial period, cumulative effects of tumor progression, subsequent treatment lines, and interval changes in nutritional status may dilute the predictive value of a single baseline measurement.

The competing risks analysis provides mechanistic insight: PNI primarily protects against non-cancer death (HR 0.63, p<0.001) rather than GI cancer-specific death (HR 1.06, p=0.70), suggesting PNI captures physiological reserve and resilience to treatment-related complications rather than directly influencing tumor biology. This finding is clinically relevant — patients with low PNI may benefit from intensified supportive care and closer monitoring for non-cancer complications. CONUT and GNRI showed protective effects against both cancer-specific and non-cancer death, possibly reflecting the additional prognostic information contributed by total cholesterol and BMI.

### Nutritional vs Inflammatory Pathways

PNI decomposition revealed that albumin alone (HR 0.69, 95% CI 0.58-0.82) outperformed the composite PNI (HR 0.78), while lymphocyte count showed no independent prognostic value (HR 1.06, p=0.54). This is consistent with albumin's established dual role as a nutritional marker and negative acute-phase reactant reflecting systemic inflammation. Our findings suggest that simpler albumin-only indices may be equally effective for risk stratification, though the composite indices additionally capture inflammatory and metabolic dimensions.

Inflammatory markers NLR and SII showed comparable prognostic performance to PNI (+0.017 C-statistic improvement for both). The mutually adjusted model, in which neither PNI nor NLR remained significant, indicates substantial overlap between nutritional and inflammatory signaling pathways. The joint exposure analysis (Low PNI + High NLR, HR 1.67, 95% CI 1.06-2.63; RERI -0.97) suggests that combined assessment of nutritional and inflammatory status identifies the highest-risk patients, though the negative RERI indicates these pathways are overlapping rather than synergistic.

Notably, the absence of PNI × CRP interaction (p=0.78) and the null CALLY index finding (HR 0.85, p=0.20) indicate that PNI's prognostic value extends beyond what can be explained by systemic inflammation measured through single biomarkers. The null dietary finding (HEI-2015 HR 0.97, p=0.85) reinforces the conceptual distinction between biochemical nutritional status — reflecting metabolic reserve and systemic inflammation — and dietary intake patterns, which may be less relevant in the post-diagnosis period when metabolic alterations dominate.

### Comparison with Previous Literature

Our findings extend prior single-center studies validating PNI in GI cancer. A meta-analysis of 28 studies in gastric cancer reported low PNI was associated with poorer overall survival (HR 1.38, 95% CI 1.23-1.54) [8], consistent with our population-level results. Our study uniquely provides dual-cohort external validation across 30 years, absolute risk estimates (RMST, PAF) that enhance clinical interpretability, and time-dependent analyses demonstrating the temporal window of nutritional risk. The head-to-head comparison of PNI, CONUT, and GNRI — rarely performed in a single population — identifies GNRI as the most stable and strongly associated index, likely because BMI adjustment captures the paradoxical protective effect of preserved body mass in cancer patients.

### Clinical Implications

These simple, inexpensive indices can be calculated from routine laboratory data in virtually any clinical setting. A PNI threshold of 48.5 identifies a high-risk group with substantially shorter survival (6.6 vs 12.8 years median survival; RMST difference 1.97 years at 10 years). The concentration of risk in the first 2 years after diagnosis suggests that early nutritional assessment and intervention may have the greatest impact. GNRI and CONUT, which incorporate BMI or total cholesterol, showed more stable, linear associations across the full follow-up period and may be preferable for long-term risk stratification when these additional measurements are available. The significant PNI × physical activity interaction (p=0.023) in our exploratory analysis suggests that the prognostic benefit of better nutritional status may be enhanced by physical activity, a modifiable factor warranting further investigation in prospective studies.

### Strengths and Limitations

Strengths include population-based sampling with national representativeness, dual-cohort external validation across 30 years, comprehensive temporal analysis (time-dependent Cox, landmark, RCS, RMST, competing risks), head-to-head comparison of biochemical, dietary, and inflammatory indices, E-value sensitivity analysis, and bootstrap C-statistic comparisons with absolute risk measures (PAF, RMST difference) that enhance clinical interpretability. The inclusion of both NHANES III and Continuous NHANES cohorts provides cross-validation across different survey methodologies and time periods.

Limitations include the modest GI cancer sample (n=313) limiting site-specific subgroup analysis, particularly for upper GI cancers. Nutritional indices were measured at a single baseline time point, not capturing longitudinal changes. Reliance on self-reported cancer diagnosis may introduce misclassification. The lack of cancer stage, treatment, and recurrence data — important confounders not available in public NHANES data — limits causal inference. The predominantly colorectal cancer composition (77%) limits generalizability to upper GI malignancies. The public-use mortality file provides only a 10-group cause-of-death recode, limiting cause-specific analysis granularity. Survey-weighted analysis was restricted to the 2005-2016 subset due to NHANES III design differences.

## 5. Conclusions

PNI, CONUT, and GNRI are robust, independent prognostic factors for all-cause mortality in GI cancer patients, validated across two nationally representative US cohorts spanning 30 years. The protective effect is time-dependent, concentrated in the first 2 years post-diagnosis, and operates primarily through non-cancer mortality pathways. An optimal PNI threshold of 48.5 identifies a high-risk subgroup (RMST difference 1.97 years at 10 years; PAF 13.4%). The prognostic value is driven primarily by serum albumin and is independent of systemic inflammation (CRP). Inflammatory markers (NLR, SII) show comparable performance, and joint analysis of PNI and NLR identifies a particularly high-risk group. These simple, inexpensive indices may be useful for risk stratification and early intervention targeting in GI cancer patients.

## Acknowledgments

The authors thank the National Center for Health Statistics for conducting the NHANES surveys and making the data publicly available.

## Funding Statement

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

## Conflict of Interest

The authors declare no conflicts of interest.

## Author Contributions

**Zhuha Zhou:** Conceptualization, Methodology, Formal analysis, Investigation, Data curation, Writing — original draft, Writing — review & editing, Visualization. **Qigang Xu:** Data curation, Investigation. **Yongyu Bai:** Investigation, Writing — review & editing. **Yiqi Cai:** Supervision, Writing — review & editing.

## References

1. Arends J, Bachmann P, Baracos V, et al. ESPEN guidelines on nutrition in cancer patients. Clin Nutr. 2017;36(1):11-48.
2. Fearon K, Strasser F, Anker SD, et al. Definition and classification of cancer cachexia: an international consensus. Lancet Oncol. 2011;12(5):489-495.
3. Hebuterne X, Lemarié E, Michallet M, et al. Prevalence of malnutrition and current use of nutrition support in patients with cancer. JPEN J Parenter Enteral Nutr. 2014;38(2):196-204.
4. Pressoir M, Desné S, Berchery D, et al. Prevalence, risk factors and clinical implications of malnutrition in French Comprehensive Cancer Centres. Br J Cancer. 2010;102(6):966-971.
5. Onodera T, Goseki N, Kosaki G. Prognostic nutritional index in gastrointestinal surgery of malnourished cancer patients. Nihon Geka Gakkai Zasshi. 1984;85(9):1001-1005.
6. Ignacio de Ulíbarri J, González-Madroño A, de Villar NG, et al. CONUT: a tool for controlling nutritional status. First validation in a hospital population. Nutr Hosp. 2005;20(1):38-45.
7. Bouillanne O, Morineau G, Dupont C, et al. Geriatric Nutritional Risk Index: a new index for evaluating at-risk elderly medical patients. Am J Clin Nutr. 2005;82(4):777-783.
8. Yang Y, Gao P, Song Y, et al. The prognostic nutritional index is a predictive indicator of prognosis and postoperative complications in gastric cancer: a meta-analysis. Eur J Surg Oncol. 2016;42(8):1176-1182.
9. Sun K, Chen S, Xu J, et al. The prognostic significance of the prognostic nutritional index in cancer: a systematic review and meta-analysis. J Cancer Res Clin Oncol. 2014;140(9):1537-1549.
10. Migita K, Takayama T, Saeki K, et al. The prognostic nutritional index predicts long-term outcomes of gastric cancer patients independent of tumor stage. Ann Surg Oncol. 2013;20(8):2647-2654.
11. Nozoe T, Kohno M, Iguchi T, et al. The prognostic nutritional index can be a prognostic indicator in colorectal carcinoma. Surg Today. 2012;42(6):532-535.
12. Pinato DJ, North BV, Sharma R. A novel, externally validated inflammation-based prognostic algorithm in hepatocellular carcinoma: the prognostic nutritional index (PNI). Br J Cancer. 2012;106(8):1439-1445.
13. Kanda M, Fujii T, Kodera Y, et al. Nutritional predictors of postoperative outcome in pancreatic cancer. Br J Surg. 2011;98(2):268-274.
14. Lee JY, Kim HI, Kim YN, et al. Clinical significance of the prognostic nutritional index in predicting postoperative gastric cancer survival. J Gastric Cancer. 2017;17(2):139-150.

## Figure Legends

- **Figure 1:** Kaplan-Meier survival curves by PNI tertile (combined NHANES cohorts).
- **Figure 2:** Forest plot — adjusted HR per 1-SD for PNI, CONUT, and GNRI.
- **Figure 3:** Time-dependent HR across follow-up periods (0-2, 2-5, 5+ years).
- **Figure 4:** Landmark analysis at 1, 3, and 5 years.
- **Figure 5:** Competing risks forest plot — GI cancer death vs non-cancer death.
- **Figure 6:** PNI cutpoint analysis (optimal threshold identified at 48.5).
- **Figure 7:** Restricted cubic spline dose-response curve for PNI.
- **Figure 8:** Subgroup forest plot — PNI effect across strata.
- **Figure 9:** PNI vs NLR vs SII comparison (C-statistic improvement).
- **Figure 10:** PNI decomposition — albumin vs lymphocyte vs composite.

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
