# ==============================================================================
# 06_mediation_and_composite.R — Extended mediation + multi-index composite
# Inner Loop 2: (1) SEM mediation for all 4 indices
#              (2) Multi-index composite score
# ==============================================================================

source(here::here("R_scripts", "00_config.R"))
library(lavaan)

# Load workspace
load(file.path(RESULTS_DIR, "analysis_workspace.RData"))

# ── SECTION A: Mediation for all 4 indices ──────────────────────────────────
cat("── Extended Mediation: All 4 Indices ──\n")

exposures_med <- c("E_DII_raw", "PNI", "CONUT", "GNRI")
exposure_labels_med <- c("E-DII (pro-inflam → anti)", "PNI", "CONUT (high = good)", "GNRI")

# For mediation: E-DII higher = pro-inflammatory (direction matters)
# For PNI/CONUT/GNRI: higher = better nutrition
# Reverse PNI/CONUT/GNRI so consistent direction (higher = better function)
# CONUT: reverse so higher = better (original CONUT higher = worse)
sem_df <- df_cc %>%
  mutate(
    E_DII_med = E_DII_raw,             # original direction (higher = more inflammatory)
    PNI_med = PNI,                     # higher = better
    CONUT_med = -CONUT + max(CONUT, na.rm = TRUE),  # reversed so higher = better
    GNRI_med = GNRI                    # higher = better
  ) %>%
  dplyr::select(SEQN, cog_z_composite, phq9_total,
                E_DII_med, PNI_med, CONUT_med, GNRI_med,
                age, sex_binary, edu_binary, income_pir, bmi, comorb_count) %>%
  filter(complete.cases(.))

cat(sprintf("SEM complete-case N = %d\n", nrow(sem_df)))

mediation_results <- list()
med_exposures <- c("E_DII_med", "PNI_med", "CONUT_med", "GNRI_med")
med_labels <- c("E-DII (pro-inflammatory)", "PNI", "CONUT (reversed)", "GNRI")

for (i in seq_along(med_exposures)) {
  exp_var <- med_exposures[i]
  exp_lab <- med_labels[i]

  cat(sprintf("\n  Mediation: %s\n", exp_lab))

  sem_model <- sprintf('
    cog_z_composite ~ c*%s + age + sex_binary + edu_binary + income_pir + bmi + comorb_count
    phq9_total ~ a*%s + age + sex_binary + edu_binary + income_pir + bmi + comorb_count
    cog_z_composite ~ b*phq9_total
    ab := a * b
    total := c + (a * b)
  ', exp_var, exp_var)

  sem_fit <- tryCatch({
    lavaan::sem(sem_model, data = sem_df, estimator = "MLR",
                missing = "fiml", bootstrap = 1000)
  }, error = function(e) {
    warning(sprintf("  SEM failed for %s: %s", exp_lab, e$message))
    return(NULL)
  })

  if (!is.null(sem_fit)) {
    pe <- parameterEstimates(sem_fit, boot.ci.type = "bca.simple", standardized = TRUE)
    pe$exposure <- exp_lab
    pe$N <- lavInspect(sem_fit, "ntotal")
    pe$converged <- lavInspect(sem_fit, "converged")
    mediation_results[[exp_lab]] <- pe

    # Extract key paths
    a_path <- pe$est[pe$label == "a"]
    b_path <- pe$est[pe$label == "b"]
    c_path <- pe$est[pe$label == "c"]
    ab_est <- pe$est[pe$label == "ab"]
    total_est <- pe$est[pe$label == "total"]
    a_p <- pe$pvalue[pe$label == "a"]

    cat(sprintf("    a (exposure→PHQ9): %.4f (p=%.4f)\n", a_path, a_p))
    cat(sprintf("    b (PHQ9→cognition): %.4f (p=%s)\n", b_path, format(pe$pvalue[pe$label == "b"], digits=4)))
    cat(sprintf("    c (direct): %.4f (p=%s)\n", c_path, format(pe$pvalue[pe$label == "c"], digits=4)))
    cat(sprintf("    Indirect (a*b): %.4f\n", ab_est))
    cat(sprintf("    Total: %.4f\n", total_est))
    cat(sprintf("    Prop mediated: %.1f%%\n", ab_est / total_est * 100))
  } else {
    cat("    SEM failed to converge\n")
  }
}
all_mediation <- do.call(rbind, mediation_results)
write.csv(all_mediation, file.path(RESULTS_DIR, "mediation_all_indices.csv"), row.names = FALSE)
cat(sprintf("\nExtended mediation saved: %d rows\n", nrow(all_mediation)))

# ── SECTION B: Multi-index composite score ──────────────────────────────────
cat("\n── Multi-Index Composite Score ──\n")

# Create standardized composite from the two strongest indices
# We standardize each and average them
composite_df <- df_cc %>%
  mutate(
    z_CONUT = -scale(CONUT)[, 1],
    z_EDII  = -scale(E_DII)[, 1],
    z_PNI   = scale(PNI)[, 1],
    z_GNRI  = scale(GNRI)[, 1]
  )
composite_df$nutri_composite_2 <- rowMeans(cbind(composite_df$z_CONUT, composite_df$z_EDII), na.rm = FALSE)
composite_df$nutri_composite_4 <- rowMeans(cbind(composite_df$z_CONUT, composite_df$z_EDII, composite_df$z_PNI, composite_df$z_GNRI), na.rm = FALSE)

# Update design
design_comp <- svydesign(
  id = ~SDMVPSU, strata = ~strata_pool, weights = ~wt_mec_4yr,
  nest = TRUE, data = composite_df
)

# Compare composites vs individual indices in Model 3
compare_vars <- c("z_CONUT", "z_EDII", "z_PNI", "z_GNRI", "nutri_composite_2", "nutri_composite_4")
compare_labels <- c("CONUT", "E-DII", "PNI", "GNRI", "Composite (CONUT+E-DII)", "Composite (All 4)")

compare_results <- list()
for (i in seq_along(compare_vars)) {
  form <- as.formula(paste("cog_z_composite ~", compare_vars[i], "+",
                           paste(covars_m3, collapse = " + ")))
  fit <- tryCatch(svyglm(form, design = design_comp), error = function(e) NULL)
  if (!is.null(fit)) {
    s <- summary(fit)$coefficients[2, ]
    compare_results[[compare_labels[i]]] <- data.frame(
      exposure = compare_labels[i],
      beta = s[1], se = s[2], t = s[3], p = s[4],
      stringsAsFactors = FALSE
    )
  }
}
comp_table <- do.call(rbind, compare_results)
comp_table$beta_abs <- abs(comp_table$beta)
comp_table <- comp_table[order(comp_table$beta_abs, decreasing = TRUE), ]
write.csv(comp_table, file.path(RESULTS_DIR, "composite_score_comparison.csv"), row.names = FALSE)
cat("Composite score comparison:\n")
print(comp_table[, c("exposure", "beta", "p")], row.names = FALSE)

cat("\n── 06_mediation_and_composite.R complete ──\n")
