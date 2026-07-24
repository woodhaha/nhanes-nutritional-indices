# CKD stratification — final robustness check
source(here::here("R_scripts", "00_config.R"))
library(haven); library(survey); library(dplyr)

load(file.path(RESULTS_DIR, "analysis_workspace.RData"))

# ── Compute eGFR (CKD-EPI 2021) ──────────────────────────────────────────
# Load creatinine from BIOPRO
scr_list <- list()
for (cyc in c("G", "H")) {
  xpt <- read_xpt(file.path(DATA_DIR, "nhanes_xpt", sprintf("BIOPRO_%s.XPT", cyc)))
  scr_list[[cyc]] <- xpt[, c("SEQN", "LBXSCR")]
}
scr <- do.call(rbind, scr_list)

# Merge creatinine, compute eGFR
df_cc <- df_cc %>%
  left_join(scr, by = "SEQN") %>%
  mutate(
    scr = as.numeric(LBXSCR),
    # CKD-EPI 2021: κ=0.7(F)/0.9(M), α=-0.241(F)/-0.302(M)
    kappa = ifelse(sex_binary == 1, 0.7, 0.9),
    alpha = ifelse(sex_binary == 1, -0.241, -0.302),
    egfr = 142 * pmin(scr / kappa, 1)^alpha * pmax(scr / kappa, 1)^(-1.200) *
           0.9938^age * ifelse(sex_binary == 1, 1.012, 1),
    ckd_stage = case_when(
      egfr >= 90 ~ "Stage 1 (≥90)",
      egfr >= 60 ~ "Stage 2 (60-89)",
      egfr >= 45 ~ "Stage 3a (45-59)",
      egfr >= 30 ~ "Stage 3b (30-44)",
      egfr >= 15 ~ "Stage 4 (15-29)",
      egfr < 15  ~ "Stage 5 (<15)",
      TRUE ~ NA_character_
    ),
    ckd_3bplus = ifelse(egfr < 45, 1, 0)
  )

cat(sprintf("eGFR available for %d participants\n", sum(!is.na(df_cc$egfr))))
cat(sprintf("CKD Stage 3b+: %d (%.1f%%)\n", sum(df_cc$ckd_3bplus == 1, na.rm=TRUE),
            mean(df_cc$ckd_3bplus == 1, na.rm=TRUE)*100))
print(table(df_cc$ckd_stage, useNA="ifany"))

# Redefine design
design_ckd <- svydesign(id=~SDMVPSU, strata=~strata_pool, weights=~wt_mec_4yr,
                         nest=TRUE, data=df_cc)

# ── Main models with eGFR adjustment ─────────────────────────────────────
cat("\n── Model with eGFR adjustment ──\n")
exposures <- c("z_CONUT", "z_EDII", "z_PNI", "z_GNRI")
# Create standardized versions (higher = better for all)
df_cc$z_CONUT <- -as.numeric(scale(df_cc$CONUT)[,1])
df_cc$z_EDII  <- -as.numeric(scale(df_cc$E_DII)[,1])
df_cc$z_PNI   <- as.numeric(scale(df_cc$PNI)[,1])
df_cc$z_GNRI  <- as.numeric(scale(df_cc$GNRI)[,1])

design_ckd <- svydesign(id=~SDMVPSU, strata=~strata_pool, weights=~wt_mec_4yr,
                         nest=TRUE, data=df_cc)

# Compare models with vs without eGFR
for (exp_var in exposures) {
  base_form <- as.formula(paste("cog_z_composite ~", exp_var, "+",
                                paste(covars_m3, collapse=" + ")))
  egfr_form <- as.formula(paste("cog_z_composite ~", exp_var, "+ egfr +",
                                paste(covars_m3, collapse=" + ")))
  b <- svyglm(base_form, design=design_ckd)
  e <- svyglm(egfr_form, design=design_ckd)
  bs <- summary(b)$coefficients[2,]
  es <- summary(e)$coefficients[2,]
  cat(sprintf("%-10s base: β=%.4f p=%.4f | +eGFR: β=%.4f p=%.4f\n",
              exp_var, bs[1], bs[4], es[1], es[4]))
}

# ── Stratified: exclude CKD 3b+ ──────────────────────────────────────────
cat("\n── Excluding CKD Stage 3b+ (eGFR <45) ──\n")
df_no_ckd <- filter(df_cc, is.na(ckd_3bplus) | ckd_3bplus == 0)
design_no_ckd <- svydesign(id=~SDMVPSU, strata=~strata_pool, weights=~wt_mec_4yr,
                            nest=TRUE, data=df_no_ckd)
for (exp_var in exposures) {
  form <- as.formula(paste("cog_z_composite ~", exp_var, "+",
                           paste(covars_m3, collapse=" + ")))
  fit <- svyglm(form, design=design_no_ckd)
  s <- summary(fit)$coefficients[2,]
  cat(sprintf("%-10s β=%.4f p=%.4f (N=%d)\n", exp_var, s[1], s[4], nrow(df_no_ckd)))
}

cat("\n── CKD analysis complete ──\n")
