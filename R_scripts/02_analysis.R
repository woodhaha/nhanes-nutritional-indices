# ==============================================================================
# 02_analysis.R — Survey-weighted analysis
# Main models, RCS dose-response, subgroup, mediation, sensitivity
# ==============================================================================

source(here::here("R_scripts", "00_config.R"))

# ── Load derived data ─────────────────────────────────────────────────────────
df <- readRDS(file.path(DATA_DIR, "nhanes_2011_2014_derived.rds"))
cat(sprintf("Loaded: N = %d participants\n", nrow(df)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION A: Survey Design Setup
# ═══════════════════════════════════════════════════════════════════════════════

# Pooled 4-year MEC weight: WTMEC2YR / 2 (2 cycles)
df <- df %>%
  mutate(
    wt_mec_4yr = WTMEC2YR / 2,
    # For pooled strata: combine cycle + original strata to ensure uniqueness
    strata_pool = paste0(SDDSRVYR, "_", SDMVSTRA)
  )

design <- svydesign(
  id      = ~SDMVPSU,
  strata  = ~strata_pool,
  weights = ~wt_mec_4yr,
  nest    = TRUE,
  data    = df
)

cat(sprintf("Survey design: %d strata, %d PSUs\n",
            length(unique(df$strata_pool)),
            length(unique(df$SDMVPSU))))

# ── Analytic subset: complete cases for regression models ─────────────────────
# Define model variables
model_vars <- c("cog_z_composite", "phq9_total", "E_DII", "PNI", "CONUT",
                "age", "sex_binary", "race_eth", "edu_binary", "income_pir",
                "income_low", "food_insecure",
                "bmi", "smoking", "alcohol",
                "dm_diagnosed", "htn", "cvd_any", "obesity", "comorb_count",
                "crp_mgdl", "pa_total_mets", "pa_active",
                "wt_mec_4yr", "strata_pool", "SDMVPSU")

df_cc <- df %>%
  filter(
    !is.na(cog_z_composite),
    !is.na(E_DII),
    !is.na(PNI),
    !is.na(phq9_total),
    !is.na(age),
    !is.na(edu_binary),
    !is.na(income_pir),
    !is.na(bmi),
    !is.na(smoking)
  )

design_cc <- svydesign(
  id      = ~SDMVPSU,
  strata  = ~strata_pool,
  weights = ~wt_mec_4yr,
  nest    = TRUE,
  data    = df_cc
)

cat(sprintf("Complete-case analytic sample: N = %d\n", nrow(df_cc)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION B: Table 1 — Baseline by PNI tertile
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Table 1: Baseline characteristics ──\n")

# Create PNI tertiles in the data
df_cc$pni_tertile <- factor(ntile(df_cc$PNI, 3), levels = 1:3,
                            labels = c("T1 (Low)", "T2 (Mid)", "T3 (High)"))

# Update design with new variable
design_cc <- svydesign(
  id = ~SDMVPSU, strata = ~strata_pool, weights = ~wt_mec_4yr,
  nest = TRUE, data = df_cc
)

# Weighted means/percentages by PNI tertile
table1_vars <- c("age", "sex_binary", "race_eth", "edu_binary",
                 "income_pir", "income_low", "food_insecure",
                 "bmi", "smoking", "alcohol",
                 "dm_diagnosed", "htn", "cvd_any", "obesity",
                 "comorb_count", "crp_mgdl", "pa_active",
                 "cerad_imm", "cerad_del", "animal_flu", "dsst",
                 "cog_z_composite", "phq9_total", "phq9_depression",
                 "E_DII", "PNI", "CONUT", "albumin_gdl")

# Continuous variables: weighted mean (SD)
continuous_vars <- c("age", "bmi", "income_pir", "comorb_count",
                     "crp_mgdl", "cerad_imm", "cerad_del", "animal_flu",
                     "dsst", "cog_z_composite", "phq9_total",
                     "E_DII", "PNI", "CONUT", "albumin_gdl")

cat_table1 <- function(var, design_obj, group_var = "pni_tertile") {
  form <- as.formula(paste0("~", var, " + ", group_var))

  if (var %in% continuous_vars) {
    # Weighted mean and SE
    means <- svyby(as.formula(paste0("~", var)), as.formula(paste0("~", group_var)),
                   design_obj, svymean, na.rm = TRUE)
    vars  <- svyby(as.formula(paste0("~", var)), as.formula(paste0("~", group_var)),
                   design_obj, svyvar, na.rm = TRUE)

    result <- data.frame(
      variable = var,
      level    = "",
      stringsAsFactors = FALSE
    )
    for (i in 1:nrow(means)) {
      m  <- means[i, 2]
      se <- sqrt(vars[i, 2])
      result[[paste0("T", i, "_mean")]] <- sprintf("%.2f (%.2f)", m, se)
    }
    result$p_value <- NA_real_
    return(result)

  } else {
    # Categorical: weighted proportions
    props <- svyby(as.formula(paste0("~", var)), as.formula(paste0("~", group_var)),
                   design_obj, svymean, na.rm = TRUE)

    # Chi-square test
    tryCatch({
      chi_test <- svychisq(as.formula(paste0("~", var, "+", group_var)),
                          design_obj, statistic = "F")
      p_val <- chi_test$p.value
    }, error = function(e) p_val <<- NA_real_)

    # Build result rows
    prop_cols <- grep("^se\\.", names(props), invert = TRUE, value = TRUE)
    prop_cols <- prop_cols[prop_cols != group_var]

    result_list <- list()
    for (col in prop_cols) {
      result_list[[length(result_list) + 1]] <- data.frame(
        variable = var,
        level = gsub(paste0("^", var), "", col),
        stringsAsFactors = FALSE
      )
    }

    # Fill in values per tertile
    for (i in 1:nrow(props)) {
      for (j in seq_along(prop_cols)) {
        col <- prop_cols[j]
        val <- props[i, col]
        se_col <- grep(paste0("^se\\.", col), names(props), value = TRUE)
        se <- if (length(se_col) > 0) props[i, se_col] else NA_real_
        result_list[[j]][[paste0("T", i, "_prop")]] <-
          sprintf("%.1f%%", val * 100)
      }
    }

    result <- do.call(rbind, result_list)
    result$p_value <- c(p_val, rep(NA_real_, nrow(result) - 1))
    return(result)
  }
}

table1_list <- lapply(table1_vars, function(v) {
  tryCatch(cat_table1(v, design_cc), error = function(e) NULL)
})

table1 <- do.call(rbind, Filter(Negate(is.null), table1_list))
write.csv(table1, file.path(RESULTS_DIR, "table1_baseline.csv"), row.names = FALSE)

cat(sprintf("Table 1 saved: %d variables, %d rows\n",
            length(unique(table1$variable)), nrow(table1)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION C: Exposures — standardize for comparability
# ═══════════════════════════════════════════════════════════════════════════════

# Standardize all four nutritional exposures so betas are comparable
# E_DII: reverse-coded so higher = better nutrition (anti-inflammatory)
df_cc <- df_cc %>%
  mutate(
    E_DII_scaled = -scale(E_DII)[, 1],    # reverse: higher = anti-inflammatory
    PNI_scaled    = scale(PNI)[, 1],
    CONUT_scaled  = -scale(CONUT)[, 1],    # reverse: higher CONUT = worse nutrition
    GNRI_scaled   = scale(GNRI)[, 1]
  )

design_cc <- svydesign(
  id = ~SDMVPSU, strata = ~strata_pool, weights = ~wt_mec_4yr,
  nest = TRUE, data = df_cc
)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION D: Main Analysis — Sequential Weighted Linear Regression
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Main Analysis: Sequential models ──\n")

# ── Base covariates ───────────────────────────────────────────────────────────
covars_m1 <- c("age", "sex_binary")
covars_m2 <- c(covars_m1, "race_eth", "edu_binary", "income_pir")
covars_m3 <- c(covars_m2, "bmi", "smoking", "alcohol", "pa_active",
               "comorb_count", "crp_mgdl")

#' Fit and extract one model
run_model <- function(exposure, covariates, outcome = "cog_z_composite",
                     design_obj = design_cc, label = "") {
  form <- as.formula(
    paste(outcome, "~", exposure, "+", paste(covariates, collapse = " + "))
  )

  fit <- svyglm(form, design = design_obj)

  coef_row <- tidy(fit, conf.int = TRUE) %>%
    filter(term == exposure)

  data.frame(
    exposure  = exposure,
    model     = label,
    beta      = coef_row$estimate,
    se        = coef_row$std.error,
    lower     = coef_row$conf.low,
    upper     = coef_row$conf.high,
    p_value   = coef_row$p.value,
    N         = nobs(fit),
    R2        = 1 - fit$deviance / fit$null.deviance,
    stringsAsFactors = FALSE
  )
}

# ── Run all models for all four exposures ─────────────────────────────────────
exposures <- c("E_DII_scaled", "PNI_scaled", "CONUT_scaled", "GNRI_scaled")
exposure_labels <- c("E-DII (anti-inflam)", "PNI", "CONUT (reversed)",
                     "GNRI")

all_results <- list()

for (i in seq_along(exposures)) {
  exp_var <- exposures[i]
  exp_lab <- exposure_labels[i]

  all_results[[length(all_results) + 1]] <-
    run_model(exp_var, covars_m1, label = "Model 1: Age+Sex")
  all_results[[length(all_results) + 1]] <-
    run_model(exp_var, covars_m2, label = "Model 2: +SES+Race")
  all_results[[length(all_results) + 1]] <-
    run_model(exp_var, covars_m3, label = "Model 3: +Lifestyle+Comorb")

  # Model 4: additionally adjust for depression
  m4_covars <- c(covars_m3, "phq9_total")
  all_results[[length(all_results) + 1]] <-
    run_model(exp_var, m4_covars,
              label = "Model 4: +Depression")
}

main_results <- do.call(rbind, all_results)
write.csv(main_results, file.path(RESULTS_DIR, "main_results.csv"), row.names = FALSE)

cat(sprintf("Main analysis complete: %d models run\n", nrow(main_results)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION E: Head-to-head comparison — standardized betas
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Head-to-head nutritional index comparison ──\n")

# Model 3 (fully adjusted) comparison for all 4 indices
head2head <- main_results %>%
  filter(model == "Model 3: +Lifestyle+Comorb") %>%
  mutate(
    exposure_label = factor(exposure,
      levels = exposures,
      labels = exposure_labels),
    beta_abs = abs(beta)
  ) %>%
  arrange(desc(beta_abs))

write.csv(head2head, file.path(RESULTS_DIR, "head2head_comparison.csv"),
          row.names = FALSE)

# ── Incremental R²: test if nutrition adds to base model ──────────────────────
base_form <- as.formula(paste("cog_z_composite ~", paste(covars_m3, collapse = " + ")))
base_model <- svyglm(base_form, design = design_cc)

r2_increment <- data.frame(
  exposure = exposures,
  label    = exposure_labels,
  base_R2  = 1 - base_model$deviance / base_model$null.deviance,
  full_R2  = NA_real_,
  delta_R2 = NA_real_,
  F_stat   = NA_real_,
  p_Ftest  = NA_real_
)

for (i in seq_along(exposures)) {
  full_form <- as.formula(paste("cog_z_composite ~", exposures[i], "+",
                                paste(covars_m3, collapse = " + ")))
  full_model <- svyglm(full_form, design = design_cc)
  r2_increment$full_R2[i] <- 1 - full_model$deviance / full_model$null.deviance
  r2_increment$delta_R2[i] <- r2_increment$full_R2[i] - r2_increment$base_R2[i]

  # Wald test for added term
  wald_test <- regTermTest(full_model, exposures[i])
  r2_increment$F_stat[i]  <- wald_test$Ftest
  r2_increment$p_Ftest[i] <- wald_test$p
}

write.csv(r2_increment, file.path(RESULTS_DIR, "r2_increment.csv"), row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION F: RCS — Dose-Response (E-DII ~ cognitive function)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── RCS dose-response analysis ──\n")

# Use raw (non-standardized) E-DII for RCS
# Negative E_DII = anti-inflammatory, positive = pro-inflammatory
df_cc$E_DII_raw <- df_cc$E_DII

# Re-fit design
design_cc <- svydesign(
  id = ~SDMVPSU, strata = ~strata_pool, weights = ~wt_mec_4yr,
  nest = TRUE, data = df_cc
)

# Hmisc rcs with 3 knots
knots <- quantile(df_cc$E_DII_raw, probs = c(0.10, 0.50, 0.90), na.rm = TRUE)

rcs_form <- as.formula(
  paste("cog_z_composite ~ rcs(E_DII_raw, parms = knots) +",
        paste(covars_m3, collapse = " + "))
)

rcs_fit <- svyglm(rcs_form, design = design_cc)

# ── Generate predictions for smooth curve ─────────────────────────────────────
E_DII_seq <- seq(
  quantile(df_cc$E_DII_raw, 0.01, na.rm = TRUE),
  quantile(df_cc$E_DII_raw, 0.99, na.rm = TRUE),
  length.out = 100
)

# Prediction data frame
pred_df <- data.frame(
  E_DII_raw = E_DII_seq,
  age       = median(df_cc$age, na.rm = TRUE),
  sex_binary = 0.5,
  race_eth  = factor("Non-Hispanic White",
                     levels = levels(as.factor(df_cc$race_eth))),
  edu_binary = median(df_cc$edu_binary, na.rm = TRUE),
  income_pir = median(df_cc$income_pir, na.rm = TRUE),
  bmi       = median(df_cc$bmi, na.rm = TRUE),
  smoking   = factor("Never", levels = levels(as.factor(df_cc$smoking))),
  alcohol   = factor("Never", levels = levels(as.factor(df_cc$alcohol))),
  pa_active = median(df_cc$pa_active, na.rm = TRUE),
  comorb_count = median(df_cc$comorb_count, na.rm = TRUE),
  crp_mgdl  = median(df_cc$crp_mgdl, na.rm = TRUE)
)

preds <- predict(rcs_fit, newdata = pred_df, se.fit = TRUE)
pred_df$y_hat   <- preds$link
pred_df$y_se    <- preds$SE
pred_df$y_lower <- preds$link - 1.96 * preds$SE
pred_df$y_upper <- preds$link + 1.96 * preds$SE

# ── Non-linearity test ────────────────────────────────────────────────────────
# Compare to linear model
lin_form <- as.formula(
  paste("cog_z_composite ~ E_DII_raw +", paste(covars_m3, collapse = " + "))
)
lin_fit <- svyglm(lin_form, design = design_cc)

# Likelihood-ratio-type test using anova
nonlin_p <- tryCatch({
  anova(lin_fit, rcs_fit, method = "Wald")$p[2]
}, error = function(e) NA_real_)

cat(sprintf("Non-linearity test p-value: %.4f\n", nonlin_p))

# Save RCS predictions
write.csv(pred_df, file.path(RESULTS_DIR, "rcs_predictions.csv"), row.names = FALSE)
cat(sprintf("RCS predictions saved (%d points, nonlin p = %.4f)\n",
            nrow(pred_df), nonlin_p))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION G: Subgroup Analysis
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Subgroup analysis ──\n")

subgroup_vars <- list(
  "Sex"        = list(var = "sex_binary", strata = c(0, 1),
                      labels = c("Male", "Female")),
  "Age Group"  = list(var = "age_group",
                      strata = c("60-69", "70-79", "80+"),
                      labels = c("60-69", "70-79", "80+")),
  "Education"  = list(var = "edu_binary", strata = c(0, 1),
                      labels = c("≤HS", ">HS")),
  "Income"     = list(var = "income_low", strata = c(0, 1),
                      labels = c("PIR>130%", "PIR≤130%")),
  "Food Security" = list(var = "food_insecure", strata = c(0, 1),
                        labels = c("Secure", "Insecure")),
  "Depression" = list(var = "phq9_depression", strata = c(0, 1),
                      labels = c("PHQ-9<10", "PHQ-9≥10")),
  "CVD"        = list(var = "cvd_any", strata = c(0, 1),
                      labels = c("No CVD", "CVD")),
  "Diabetes"   = list(var = "dm_diagnosed", strata = c(0, 1),
                      labels = c("No DM", "DM"))
)

# Using the best-performing exposure: E_DII_scaled (from head-to-head)
# Adapt based on results — default to E_DII

run_subgroup <- function(sg_var, sg_level, sg_label, exposure_var,
                         covariates, design_obj) {
  sub_data <- subset(design_obj, get(sg_var) == sg_level)
  n_sub <- nrow(sub_data$variables)

  if (n_sub < 30) return(NULL)

  tryCatch({
    form <- as.formula(
      paste("cog_z_composite ~", exposure_var, "+",
            paste(setdiff(covariates, sg_var), collapse = " + "))
    )
    fit <- svyglm(form, design = sub_data)
    coef_row <- tidy(fit, conf.int = TRUE) %>% filter(term == exposure_var)

    data.frame(
      subgroup_var   = sg_var,
      subgroup_level = sg_label,
      N              = n_sub,
      beta           = coef_row$estimate,
      se             = coef_row$std.error,
      lower          = coef_row$conf.low,
      upper          = coef_row$conf.high,
      p_value        = coef_row$p.value,
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
}

# Run subgroups for E-DII (primary exposure)
subgroup_list <- list()
for (sg_name in names(subgroup_vars)) {
  sg <- subgroup_vars[[sg_name]]
  for (j in seq_along(sg$strata)) {
    res <- run_subgroup(sg$var, sg$strata[j], sg$labels[j],
                       "E_DII_scaled", covars_m3, design_cc)
    if (!is.null(res)) subgroup_list[[length(subgroup_list) + 1]] <- res
  }
}

subgroup_results <- do.call(rbind, subgroup_list)

# ── Interaction tests ─────────────────────────────────────────────────────────
interaction_results <- data.frame(
  subgroup_var = character(),
  p_interaction = numeric(),
  stringsAsFactors = FALSE
)

for (sg_name in names(subgroup_vars)) {
  sg <- subgroup_vars[[sg_name]]
  int_form <- as.formula(
    paste("cog_z_composite ~ E_DII_scaled *", sg$var, "+",
          paste(setdiff(covars_m3, sg$var), collapse = " + "))
  )
  int_fit <- svyglm(int_form, design = design_cc)

  # Get interaction term p-value
  int_term <- grep("E_DII_scaled:", names(coef(int_fit)),
                   value = TRUE, fixed = TRUE)
  if (length(int_term) > 0) {
    int_p <- tidy(int_fit) %>% filter(term %in% int_term) %>% pull(p.value)
    interaction_results <- rbind(interaction_results, data.frame(
      subgroup_var = sg_name,
      p_interaction = int_p[1],
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(subgroup_results, file.path(RESULTS_DIR, "subgroup_results.csv"),
          row.names = FALSE)
write.csv(interaction_results, file.path(RESULTS_DIR, "interaction_results.csv"),
          row.names = FALSE)

cat(sprintf("Subgroup analysis: %d strata, %d interactions tested\n",
            nrow(subgroup_results), nrow(interaction_results)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION H: Mediation Analysis
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Mediation analysis (SEM) ──\n")

# Path: E_DII → depression (phq9_total) → cognitive function
# Using lavaan with survey-corrected SEs (bootstrapped)
# Note: lavaan doesn't natively support survey weights,
# so we use bootstrap or report unweighted SEM as sensitivity

# Prepare data frame for SEM
sem_df <- df_cc %>%
  dplyr::select(cog_z_composite, phq9_total, E_DII_raw, age, sex_binary,
                edu_binary, income_pir, bmi, smoking, comorb_count) %>%
  filter(complete.cases(.))

# SEM model specification
# Model: E_DII → cognition (direct)
#        E_DII → depression → cognition (indirect via depression)
#        E_DII → cognition → depression (reverse path for comparison)
sem_model <- '
  # Direct effects
  cog_z_composite ~ c*E_DII_raw + age + sex_binary + edu_binary + income_pir + bmi + comorb_count
  phq9_total ~ a*E_DII_raw + age + sex_binary + edu_binary + income_pir + bmi + comorb_count

  # Path from depression to cognition
  cog_z_composite ~ b*phq9_total

  # Indirect effect
  ab := a * b

  # Total effect
  total := c + (a * b)
'

sem_fit <- tryCatch({
  lavaan::sem(sem_model, data = sem_df, estimator = "MLR",
              missing = "fiml", bootstrap = 1000)
}, error = function(e) {
  warning("SEM failed: ", e$message)
  return(NULL)
})

if (!is.null(sem_fit)) {
  sem_summary <- parameterEstimates(sem_fit, boot.ci.type = "bca.simple",
                                    standardized = TRUE)
  write.csv(sem_summary, file.path(RESULTS_DIR, "mediation_sem.csv"),
            row.names = FALSE)

  cat(sprintf("SEM convergence: %s\n", lavInspect(sem_fit, "converged")))
  cat(sprintf("  N used: %d\n", lavInspect(sem_fit, "ntotal")))
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION I: Sensitivity Analyses
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Sensitivity analyses ──\n")

sensitivity_list <- list(
  "Primary (Model 3)"   = list(subset = NULL, label = "Primary analysis"),
  "Exclude stroke"     = list(subset = expression(stroke == 0),
                              label = "Excluding stroke history"),
  "Exclude cancer"     = list(subset = expression(cancer_hx == 0),
                              label = "Excluding cancer history"),
  "Exclude PHQ9 ≥ 10"  = list(subset = expression(phq9_depression == 0),
                              label = "Excluding depression"),
  "Exclude albumin<3.5" = list(subset = expression(albumin_gdl >= 3.5),
                               label = "Excluding hypoalbuminemia"),
  "Age 65+ only"       = list(subset = expression(age >= 65),
                              label = "Restricted to age ≥ 65"),
  "Age 60-75"          = list(subset = expression(age >= 60 & age <= 75),
                              label = "Restricted to age 60-75")
)

sensitivity_results <- data.frame()

for (sens_name in names(sensitivity_list)) {
  sens <- sensitivity_list[[sens_name]]

  if (is.null(sens$subset)) {
    sub_design <- design_cc
  } else {
    sub_design <- subset(design_cc, eval(sens$subset))
  }

  n_sens <- nrow(sub_design$variables)

  if (n_sens < 100) {
    warning("Sensitivity [", sens$label, "]: N=", n_sens, " — too small, skipping")
    next
  }

  # Fit Model 3 on subset
  sens_form <- as.formula(
    paste("cog_z_composite ~ E_DII_scaled +", paste(covars_m3, collapse = " + "))
  )
  sens_fit <- svyglm(sens_form, design = sub_design)
  coef_row <- tidy(sens_fit, conf.int = TRUE) %>% filter(term == "E_DII_scaled")

  sensitivity_results <- rbind(sensitivity_results, data.frame(
    analysis = sens$label,
    N        = n_sens,
    beta     = coef_row$estimate,
    se       = coef_row$std.error,
    lower    = coef_row$conf.low,
    upper    = coef_row$conf.high,
    p_value  = coef_row$p.value,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  [%s] N=%d, beta=%.3f (%.3f, %.3f), p=%.4f\n",
              sens$label, n_sens,
              coef_row$estimate, coef_row$conf.low, coef_row$conf.high,
              coef_row$p.value))
}

write.csv(sensitivity_results, file.path(RESULTS_DIR, "sensitivity_analysis.csv"),
          row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION J: Depression as outcome — secondary analysis
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Secondary: Nutrition → Depression ──\n")

dep_covars <- setdiff(covars_m3, "comorb_count")
dep_covars <- c(dep_covars, "comorb_count")

dep_results <- data.frame()
for (i in seq_along(exposures)) {
  exp_var <- exposures[i]
  exp_lab <- exposure_labels[i]

  # Linear: PHQ-9 continuous
  dep_lin <- run_model(exp_var, dep_covars, outcome = "phq9_total",
                       label = "PHQ-9 (continuous)")

  # Logistic: PHQ-9 ≥ 10
  log_form <- as.formula(
    paste("phq9_depression ~", exp_var, "+", paste(dep_covars, collapse = " + "))
  )
  log_fit <- svyglm(log_form, design = design_cc, family = quasibinomial())
  log_coef <- tidy(log_fit, conf.int = TRUE) %>% filter(term == exp_var)

  dep_log <- data.frame(
    exposure = exp_var,
    model    = "PHQ-9 ≥ 10 (logistic)",
    beta     = log_coef$estimate,
    se       = log_coef$std.error,
    lower    = log_coef$conf.low,
    upper    = log_coef$conf.high,
    p_value  = log_coef$p.value,
    N        = nobs(log_fit),
    R2       = NA_real_,
    stringsAsFactors = FALSE
  )

  dep_results <- rbind(dep_results, dep_lin, dep_log)
}

write.csv(dep_results, file.path(RESULTS_DIR, "depression_outcome.csv"),
          row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION K: Save workspace
# ═══════════════════════════════════════════════════════════════════════════════

save.image(file.path(RESULTS_DIR, "analysis_workspace.RData"))
cat("\n── 02_analysis.R complete ──\n")
cat(sprintf("All results saved to %s\n", RESULTS_DIR))
