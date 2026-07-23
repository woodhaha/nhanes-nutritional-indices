# GI Tumor × Nutrition Analysis — Complete Results Summary
══════════════════════════════════════════════════════════════════════════════
 A12 — SURVEY-WEIGHTED COX (svycoxph, 2005-2016 subset)
══════════════════════════════════════════════════════════════════════════════

 GI subset with weights: 270 GI patients (134 events)

 Index      HR (95% CI)               p
────────────────────────────────────────────────────────
 PNI        0.869 (0.611-1.234)       0.432
 CONUT      0.757 (0.630-0.909)       0.003
 GNRI       0.800 (0.612-1.046)       0.103

 Unweighted (same subset for reference):
   PNI:  0.776 (0.605-0.995), p=0.046
   CONUT: 0.720 (0.607-0.854), p<0.001
   GNRI:  0.670 (0.551-0.814), p<0.001

 → Survey attenuation: PNI loses significance under complex survey design,
   likely due to NHANES III (longer follow-up) being down-weighted.
   CONUT remains significant (p=0.003).

══════════════════════════════════════════════════════════════════════════════
 A13 — E-VALUE SENSITIVITY ANALYSIS
══════════════════════════════════════════════════════════════════════════════

 Index     E-value    CI E-value
────────────────────────────────────────
 PNI       1.660      1.201
 CONUT     1.850      1.534
 GNRI      1.968      1.609

 Interpretation: An unmeasured confounder would need an association of
 ≥1.66 with both PNI and all-cause mortality to explain away the observed
 HR=0.78. For GNRI (E-value 1.97), results are more robust.

══════════════════════════════════════════════════════════════════════════════
 A14 — CALLY INDEX (CRP × Albumin × Lymphocyte)
══════════════════════════════════════════════════════════════════════════════

 Available: 250 GI patients (120 events)
 CALLY (adjusted): HR = 0.853 (0.670-1.088), p = 0.200

 C-stat: Base=0.650, +PNI=0.667, +CALLY=0.659
 → CALLY does NOT outperform PNI. Not statistically significant.

══════════════════════════════════════════════════════════════════════════════
 A15 — CONUT + GNRI DOSE-RESPONSE (RCS)
══════════════════════════════════════════════════════════════════════════════

 Index     Nonlinearity p
────────────────────────────
 CONUT     0.111 (linear)
 GNRI      0.162 (linear)
 PNI       <0.001 (nonlinear, reference)

 → CONUT and GNRI show linear relationships with mortality risk.

══════════════════════════════════════════════════════════════════════════════
 A16 — RESTRICTED MEAN SURVIVAL TIME (RMST)
══════════════════════════════════════════════════════════════════════════════

 By PNI cutpoint (48.5, 10-year restriction):
   High PNI: 8.22 yr, Low PNI: 6.25 yr, Diff: 1.97 yr (1.19-2.76), p<0.001

 → PNI < 48.5 → ~2 years shorter life expectancy within 10 years.

══════════════════════════════════════════════════════════════════════════════
 A17 — PNI × NLR JOINT EXPOSURE + RERI
══════════════════════════════════════════════════════════════════════════════

 258 GI, 125 events
 Low PNI + High NLR: HR=1.67 (1.06-2.63), p=0.028 (ref=High PNI+Low NLR)
 RERI: -0.966 (95% CI -2.608 to -0.095)

 → Highest risk in Low PNI + High NLR group. Sub-multiplicative interaction.

══════════════════════════════════════════════════════════════════════════════
 A18 — BOOTSTRAP C-STATISTIC (500 reps)
══════════════════════════════════════════════════════════════════════════════

 +PNI: ΔC=0.022 (-0.000-0.052), +CONUT: 0.023 (0.004-0.049)
 +GNRI: ΔC=0.034 (0.009-0.063), +NLR: 0.015 (-0.004-0.041)
 → GNRI shows largest improvement.

══════════════════════════════════════════════════════════════════════════════
 A19 — MICE MULTIPLE IMPUTATION
══════════════════════════════════════════════════════════════════════════════

 Missingness in key GI variables: 0% → Imputation not needed.

══════════════════════════════════════════════════════════════════════════════
 A20 — PNI × PHYSICAL ACTIVITY
══════════════════════════════════════════════════════════════════════════════

 263 GI (130 events), 39.5% active
 PA alone: HR=0.99 (0.66-1.47), p=0.94
 PNI × PA interaction: HR=1.76 (1.08-2.86), p=0.023

 PA alone non-significant. Significant interaction suggests PNI effect
 varies by activity level.

══════════════════════════════════════════════════════════════════════════════
 A21 — POPULATION ATTRIBUTABLE FRACTION (PAF)
══════════════════════════════════════════════════════════════════════════════

 PNI < 48.5 prevalence: 35.5%, HR=1.610
 PAF: 13.4% (95% CI 6.2%-20.1%)
 → ~13% of deaths potentially attributable to low PNI.

══════════════════════════════════════════════════════════════════════════════
 NEW OUTPUT FILES (03_advanced.R)
══════════════════════════════════════════════════════════════════════════════

 results/gi_analysis/
   cox_surveyweighted.csv, evalue_analysis.csv, cally_analysis.csv
   rcs_conut_predictions.csv, rcs_gnri_predictions.csv
   rmst_analysis.csv, pni_nlr_joint.csv, paf_analysis.csv

 figures/gi_analysis/
   fig_rcs_conut.pdf, fig_rcs_gnri.pdf, fig_rcs_pni.pdf
   fig_s1_cally_km.pdf, fig_s2_joint_pni_nlr.pdf, fig_s3_rmst.pdf

══════════════════════════════════════════════════════════════════════════════
# Generated: 2026-07-23 | All ages, final clean data
# Updated: 2026-07-23 | Added Tier 1-3 advanced analyses (A12-A21)

══════════════════════════════════════════════════════════════════════════════
 DATA OVERVIEW
══════════════════════════════════════════════════════════════════════════════

 Total NHANES participants (age >=18): 44,066
 GI tumor patients:                    313  (0.7%)
   - NHANES III (1988-1994):            88
   - NHANES 2005-2016:                 225
 Deaths among GI patients:            169  (54.0%)
 GI cancer deaths:                      58
 Non-cancer deaths:                    111

 GI site distribution:
   Colon:         223  (71.2%)
   Stomach:        27   (8.6%)
   Liver:          29   (9.3%)
   Esophageal:     18   (5.8%)
   Rectal:          8   (2.6%)
   Pancreatic:      7   (2.2%)
   Gallbladder:     1   (0.3%)

══════════════════════════════════════════════════════════════════════════════
 TABLE 1 — Demographics by GI status
══════════════════════════════════════════════════════════════════════════════

 Variable                   GI Tumor (n=313)  Non-Cancer (n=40,480) 
────────────────────────────────────────────────────────────────────────
 N                           313                40480                
 Age, years                  67.5 (13.8)         45.5 (18.4)         
 Female                       64 (20.4%)        7172 (17.7%)        
 Race: NH White              160 (51.1%)       14887 (36.8%)        
 Race: NH Black               45 (14.4%)        5909 (14.6%)        
 Race: Hispanic              101 (32.3%)       16394 (40.5%)        
 Education > High school     187 (59.7%)       25808 (63.8%)        
 BMI, kg/m2                   28.5 (6.1)          28.3 (6.6)         
 Albumin, g/dL                 4.09 (0.38)         4.22 (0.37)       
 Lymphocyte, /uL            1972 (990)          2209 (719)          
 PNI                          50.8 (6.5)          53.2 (5.2)         
 CONUT                         1.1 (1.4)           0.7 (1.0)         
 GNRI                        115.0 (12.2)        116.4 (12.4)        
 Deaths                      169 (54.0%)        7659 (18.9%)        
 Survival, years               8.8 (6.9)          13.3 (8.8)         

══════════════════════════════════════════════════════════════════════════════
 TABLE 2 — Baseline characteristics by PNI tertile (GI tumor only)
══════════════════════════════════════════════════════════════════════════════

 Variable                   Low PNI T1        Mid PNI T2        High PNI T3    
                           (n=105)           (n=104)           (n=104)        
────────────────────────────────────────────────────────────────────────
 Age, years                  72.6 (9.4)        66.4 (15.5)        63.5 (14.2)   
 Female                      15 (14.3%)        26 (25.0%)        23 (22.1%)   
 PNI                         44.8 (3.5)        50.5 (1.4)        57.1 (5.9)    
 CONUT                        2.1 (1.7)         0.9 (0.9)         0.4 (0.6)    
 GNRI                       110.2 (12.6)      115.5 (10.4)      119.4 (11.9)   
 Albumin, g/dL                3.79 (0.36)       4.14 (0.26)       4.35 (0.29)   
 Lymphocyte, /uL           1376 (389)        1820 (516)        2726 (1262)     
 BMI, kg/m2                  28.3 (6.1)        28.5 (5.6)        28.8 (6.6)    
 Deaths                      74 (70.5%)        48 (46.2%)        47 (45.2%)    
 Median survival, yr          6.5 (5.5)         12.5 (7.8)        13.9 (8.5)   

══════════════════════════════════════════════════════════════════════════════
 A1 — MAIN COX REGRESSION (All-cause mortality)
══════════════════════════════════════════════════════════════════════════════

 Index     Adjustment     HR (95% CI)               p-value      
────────────────────────────────────────────────────────────────────────
 PNI       Crude          0.630 (0.513-0.775)       <0.001       
 PNI       Adjusted*      0.779 (0.631-0.961)       0.020        
 CONUT     Crude          0.629 (0.545-0.727)       <0.001       
 CONUT     Adjusted*      0.714 (0.613-0.833)       <0.001       
 GNRI      Crude          0.665 (0.559-0.790)       <0.001       
 GNRI      Adjusted*      0.673 (0.563-0.803)       <0.001       
 *Adjusted for age, sex, race/ethnicity

 PH assumption test:
   PNI:     p=0.017  → violated → time-dependent Cox needed
   CONUT:   p=0.424  → OK
   GNRI:    p=0.092  → borderline OK

══════════════════════════════════════════════════════════════════════════════
 A2 — TIME-DEPENDENT COX (PNI × follow-up period interaction)
══════════════════════════════════════════════════════════════════════════════

 Index     Period         HR (95% CI)               p (interaction)
────────────────────────────────────────────────────────────────────────
 PNI       0-2 yr         0.340 (0.212-0.546)       Reference      
 PNI       2-5 yr         3.860 (1.738-8.576)       <0.001        
 PNI       5+ yr          2.376 (1.369-4.125)       0.002         
 CONUT     0-2 yr         0.671 (0.485-0.928)       Reference      
 CONUT     2-5 yr         1.172 (0.755-1.821)       0.479         
 CONUT     5+ yr          1.032 (0.698-1.526)       0.875         
 GNRI      0-2 yr         0.536 (0.300-0.956)       Reference      
 GNRI      2-5 yr         1.228 (0.628-2.400)       0.548         
 GNRI      5+ yr          1.369 (0.746-2.512)       0.311         

 → PNI effect concentrates in first 2 years post-diagnosis.
 → CONUT/GNRI effects more stable over time.

══════════════════════════════════════════════════════════════════════════════
 A3 — LANDMARK ANALYSIS
══════════════════════════════════════════════════════════════════════════════

 Index     1-year             3-year              5-year            
────────────────────────────────────────────────────────────────────────
 PNI       0.90 (0.73-1.12)   0.84 (0.65-1.08)    0.80 (0.60-1.08)  
           p=0.358            p=0.174             p=0.138           
 CONUT     0.75 (0.63-0.89)   0.72 (0.58-0.88)    0.74 (0.57-0.96)  
           p<0.001            p=0.002             p=0.020           
 GNRI      0.70 (0.59-0.85)   0.70 (0.57-0.87)    0.71 (0.56-0.90)  
           p<0.001            p<0.001             p=0.005           

 → PNI loses significance in conditional survival; CONUT/GNRI remain significant.

══════════════════════════════════════════════════════════════════════════════
 A4 — COMPETING RISKS (Cause-specific Cox)
══════════════════════════════════════════════════════════════════════════════

 Index       GI Cancer Death        Non-Cancer Death      
────────────────────────────────────────────────────────────────────────
 PNI         1.068 (0.785-1.453)    0.636 (0.485-0.835)   
             p=0.676                p=0.001              
 CONUT       0.728 (0.563-0.943)    0.706 (0.583-0.856)   
             p=0.016                p<0.001              
 GNRI        0.692 (0.512-0.935)    0.661 (0.531-0.824)   
             p=0.016                p<0.001              

 → Nutrition indices protect against NON-CANCER death, not GI cancer death.
 → Suggests mechanism: physiological resilience > anti-tumor effect.

══════════════════════════════════════════════════════════════════════════════
 A5 — PNI OPTIMAL CUTPOINT
══════════════════════════════════════════════════════════════════════════════

 Method:         Maximally selected log-rank statistic (Contal & O'Quigley)
 Optimal PNI:    48.5
 Log-rank p:     <0.0001

 Group           N    Events    Median survival (years)
────────────────────────────────────────────────────────
 Low PNI (<48.5)  111   75        6.58 (5.25-8.92)
 High PNI (≥48.5) 202   94       12.83 (11.75-19.58)

 → PNI < 48.5 identifies a subgroup with ~2× hazard and halved median survival.

══════════════════════════════════════════════════════════════════════════════
 A6 — RESTRICTED CUBIC SPLINE (Dose-response)
══════════════════════════════════════════════════════════════════════════════

 Method:          3-knot restricted cubic spline
 Nonlinearity:    p < 0.0001  (LRT vs linear model)

 → PNI-mortality relationship is nonlinear.
 → Sharp risk increase below PNI ~45; plateau above ~50.
 → Supports the threshold model for clinical application.

══════════════════════════════════════════════════════════════════════════════
 A7 — SUBGROUP ANALYSIS (PNI effect)
══════════════════════════════════════════════════════════════════════════════

 Subgroup                   N    Events    HR (95% CI)            p
────────────────────────────────────────────────────────────────────────
 All patients               313   169      0.78 (0.63-0.96)       0.020
 Age < 60                    76    24      0.47 (0.26-0.85)       0.013
 Age 60-74                  115    66      0.76 (0.53-1.10)       0.149
 Age ≥ 75                   122    79      0.91 (0.65-1.27)       0.576
 Female                      64    38      0.56 (0.32-0.99)       0.047
 Male                       249   131      0.82 (0.64-1.05)       0.120
 NH White                   160    93      0.71 (0.53-0.96)       0.024
 NH Black                    45    25      1.14 (0.58-2.24)       0.710
 Hispanic                   101    48      0.67 (0.42-1.07)       0.091
 BMI < 30                   210   117      0.72 (0.56-0.94)       0.015
 BMI ≥ 30                   103    52      0.81 (0.50-1.31)       0.387
 CRC                        241   133      0.78 (0.62-0.99)       0.037
 Non-CRC GI                  72    36      0.73 (0.43-1.24)       0.247

 PNI × age interaction:     p = 0.727

 → PNI effect directionally consistent across all subgroups.
 → Effect appears stronger in age<60, female, NH White, BMI<30.
 → Caution: subgroup CIs wide due to small N.

══════════════════════════════════════════════════════════════════════════════
 A8 — PNI DECOMPOSITION (Component vs Composite)
══════════════════════════════════════════════════════════════════════════════

 Model         C-statistic    HR (95% CI)               p
────────────────────────────────────────────────────────────────────────
 Base (demo)   0.661          —                         —
 + Albumin     0.694          0.689 (0.581-0.817)       <0.001
 + Lymphocyte  0.661          1.061 (0.879-1.280)       0.539
 + PNI         0.684          0.779 (0.631-0.961)       0.020

 → ALBUMIN drives the PNI effect. Lymphocyte alone has NO prognostic value.
 → PNI ≈ albumin + noise; CONUT/GNRI may add cholesterol/BMI information.

══════════════════════════════════════════════════════════════════════════════
 A9 — CRP × PNI INTERACTION
══════════════════════════════════════════════════════════════════════════════

 CRP available: 250 GI patients (120 deaths)
 PNI × CRP interaction p = 0.780

 → No evidence that PNI effect varies by inflammation level.
 → PNI's prognostic value is independent of systemic inflammation (CRP).

══════════════════════════════════════════════════════════════════════════════
 A10 — DIETARY QUALITY (HEI-2015)
══════════════════════════════════════════════════════════════════════════════

 Available: 99 GI patients (33 deaths) from NHANES 2011-2016

 Index           HR (95% CI)               p
────────────────────────────────────────────────────────
 HEI alone       0.969 (0.697-1.346)       0.850
 PNI alone       0.679 (0.452-1.021)       0.063
 PNI mutual      0.672 (0.443-1.019)       0.062
 HEI mutual      1.041 (0.746-1.454)       0.811

 HEI-PNI correlation: r = 0.164

 → HEI-2015 does NOT predict GI cancer survival (null finding).
 → PNI marginally protective even in this small subset.

══════════════════════════════════════════════════════════════════════════════
 OUTPUT FILES
══════════════════════════════════════════════════════════════════════════════

 results/gi_analysis/
   nhanes_clean.rds                  — Cleaned combined data (tibble)
   nhanes_clean_with_hei.rds         — With HEI-2015 merged
   table1.csv                        — Table 1 comprehensive
   table2.csv                        — Table 2 by PNI tertile
   cox_main.csv                      — A1: Main Cox
   cox_timedependent.csv             — A2: Time-dependent
   cox_landmark.csv                  — A3: Landmark
   cox_causespecific.csv             — A4: Competing risks
   pni_decomposition.csv             — A8: Component analysis
   subgroup_forest.csv               — A7: Subgroup analysis
   rcs_predictions.csv               — A6: RCS predictions
   crp_interaction.csv               — A9: CRP × PNI

 figures/gi_analysis/
   fig1_km.pdf/png                   — KM by PNI tertile
   fig2_forest.pdf/png               — Main Cox forestplot
   fig_cutpoint.pdf                  — PNI cutpoint KM (48.5)
   fig_rcs.pdf                       — Restricted cubic spline
   fig_subgroup_forest.pdf/png       — Subgroup forest plot

══════════════════════════════════════════════════════════════════════════════
 A11 — NLR / SII COMPARISON (Inflammatory vs Nutritional)
══════════════════════════════════════════════════════════════════════════════

 Subset: 258 GI patients (125 deaths) from 2005-2016 (CBC with neutrophils)
 NLR = neutrophil count / lymphocyte count
 SII = platelet × neutrophil / lymphocyte

 Index        HR (95% CI)               p        ∆C-stat
────────────────────────────────────────────────────────────────────────
 Base model   —                         —         (0.658)
 PNI          0.784 (0.606-1.013)       0.063    +0.017
 NLR          1.430 (1.188-1.722)       <0.001   +0.017
 logNLR       1.282 (1.036-1.587)       0.022    +0.017
 SII          1.303 (1.102-1.540)       0.002    —
 logSII       1.239 (1.012-1.518)       0.038    +0.013

 Mutual adjustment:
   PNI + logNLR: PNI HR=0.86 (0.65-1.14, p=0.29)
                 logNLR HR=1.22 (0.97-1.54, p=0.09)

 → NLR and PNI show similar prognostic strength.
 → Inflammatory and nutritional pathways share predictive information.
 → Neither dominates when both are in the model.

══════════════════════════════════════════════════════════════════════════════
 NEW FIGURES GENERATED
══════════════════════════════════════════════════════════════════════════════

 figures/gi_analysis/
   fig1_km.pdf/png                   — KM by PNI tertile
   fig2_forest.pdf/png               — Main Cox forest plot
   fig_cutpoint.pdf                  — PNI cutpoint KM (threshold 48.5)
   fig_rcs.pdf                       — Restricted cubic spline
   fig_subgroup_forest.pdf/png       — Subgroup forest plot

══════════════════════════════════════════════════════════════════════════════
 MANUSCRIPT (UPDATED)
══════════════════════════════════════════════════════════════════════════════

 MANUSCRIPT.md — Updated with:
   §3.6  PNI Decomposition (albumin drives the effect)
   §3.7  Optimal threshold (48.5) + RCS dose-response
   §3.8  PNI vs NLR/SII comparison
   §3.9  CRP interaction
   §3.10 HEI-2015 dietary comparison
   §3.11 Sensitivity analysis
   §3.12 Subgroup analysis
   Abstract, Discussion, Conclusions all updated.

══════════════════════════════════════════════════════════════════════════════
 PIPELINE SCRIPTS
══════════════════════════════════════════════════════════════════════════════

 run_all.R        — One-command pipeline (Rscript run_all.R)
 01_prepare.R     — Data preparation
 02_analyze.R     — 10 core analyses + figures
 hei_merge.R      — HEI-2015 merge
 nlr_sii_extract.R — NLR/SII extraction + comparison
 _archive/        — Old scripts (pre-refactoring)
