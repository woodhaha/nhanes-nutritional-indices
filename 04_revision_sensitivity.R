# ─── 04_revision_sensitivity.R — Reviewer-requested sensitivity analyses ───
# Addresses 4-panel review consensus:
#   R1+R3+R4: cohort pooling, confounding, PNI vs albumin, threshold bias
#   R2: BCG/BCP, GLIM discussion (text only)
#   R4: cross-validation, cohort-stratified, extended covariates
#
# Rscript 04_revision_sensitivity.R
# Requires: 01_prepare.R + 02_analyze.R + 03_advanced.R already run

library(dplyr); library(survival); library(broom); library(ggplot2)

PROJ <- normalizePath(".")
RES_DIR <- file.path(PROJ, "results/gi_analysis")
FIG_DIR <- file.path(PROJ, "figures/gi_analysis")

# -- Load -----------------------------------------------------------------
df <- readRDS(file.path(RES_DIR, "nhanes_clean.rds"))
for (nm in names(df)) {
  if (is.data.frame(df[[nm]])) df[[nm]] <- as.numeric(df[[nm]][[1]])
  if (is.list(df[[nm]])) df[[nm]] <- as.numeric(unlist(df[[nm]]))
}

gi <- df %>% filter(gi_tumor==1, surv_years>0) %>%
  mutate(PNI_s = as.numeric(scale(PNI)),
         CONUT_s = as.numeric(scale(-CONUT)),
         GNRI_s = as.numeric(scale(GNRI)),
         age_s = as.numeric(scale(age)),
         sex_b = ifelse(sex=="Female", 1L, 0L),
         race_eth_b = ifelse(race_eth=="Non-Hispanic White", 1L, 0L),
         cohort = ifelse(grepl("III", cycle), "NHANES_III", "Continuous_NHANES"))
gi <- gi %>% filter(!is.na(death), !is.na(PNI))
cat(sprintf("GI N=%d, events=%d\n", nrow(gi), sum(gi$death)))

# ═══════════════════════════════════════════════════════════════════════════
# 1. COHORT-STRATIFIED ANALYSIS (R1, R3, R4)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== 1. Cohort-stratified Cox ===\n")

adj_cov <- "age + sex_b + race_eth_b"

strat_res <- bind_rows(lapply(c("PNI_s","CONUT_s","GNRI_s"), function(idx) {
  bind_rows(lapply(c("NHANES_III","Continuous_NHANES"), function(co) {
    s <- gi %>% filter(cohort == co)
    if (sum(s$death) < 10) return(data.frame())
    f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", adj_cov))
    hr <- tidy(coxph(f, data=s), conf.int=TRUE) %>% filter(term==idx)
    data.frame(Index=gsub("_s","",idx), Cohort=co,
               N=nrow(s), Events=sum(s$death),
               HR=exp(hr$estimate), Lower=exp(hr$conf.low),
               Upper=exp(hr$conf.high), P=hr$p.value)
  }))
}))
strat_res <- strat_res %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
cat("Cohort-stratified:\n"); print(strat_res[, c("Index","Cohort","N","Events","HR_CI","P")])
write.csv(strat_res, file.path(RES_DIR, "sensitivity_cohort_stratified.csv"), row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════
# 2. COHORT × PNI INTERACTION (R1, R4)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== 2. Cohort × PNI interaction ===\n")
gi <- gi %>% mutate(cohort_b = ifelse(cohort == "Continuous_NHANES", 1L, 0L))
for (idx in c("PNI_s","CONUT_s","GNRI_s")) {
  f <- as.formula(paste("Surv(surv_years, death) ~", idx, "* cohort_b +", adj_cov))
  h <- tidy(coxph(f, data=gi), conf.int=TRUE)
  int <- h %>% filter(grepl(":", term))
  cat(sprintf("  %s × cohort: HR=%.3f (%.3f-%.3f), p=%.4f\n",
      gsub("_s","",idx), exp(int$estimate), exp(int$conf.low), exp(int$conf.high), int$p.value))
}
write.csv(data.frame(Sensitivity="cohort_interaction", p_PNI=NA),  # placeholder
          file.path(RES_DIR, "sensitivity_cohort_interaction.csv"), row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════
# 3. EXTENDED COVARIATE MODEL (R4: SES + physical activity)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== 3. Extended covariate model ===\n")
gi_ext <- gi %>% filter(!is.na(edu_binary))
cat(sprintf("  Extended model N=%d, events=%d\n", nrow(gi_ext), sum(gi_ext$death)))

adj_ext <- "age + sex_b + race_eth_b + edu_binary + income_pir"
adj_ext_pa <- "age + sex_b + race_eth_b + edu_binary + income_pir + any_mvpa"

ext_res <- bind_rows(lapply(c("PNI_s","CONUT_s","GNRI_s"), function(idx) {
  bind_rows(list(
    run1 = {f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+ age + sex_b + race_eth_b"));
            hr <- tidy(coxph(f, data=gi_ext), conf.int=TRUE) %>% filter(term==idx);
            data.frame(Index=gsub("_s","",idx), Model="Basic (age+sex+race)",
                       N=nrow(gi_ext), Events=sum(gi_ext$death),
                       HR=exp(hr$estimate), Lower=exp(hr$conf.low),
                       Upper=exp(hr$conf.high), P=hr$p.value)},
    run2 = {f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", adj_ext));
            hr <- tidy(coxph(f, data=gi_ext), conf.int=TRUE) %>% filter(term==idx);
            data.frame(Index=gsub("_s","",idx), Model="+ Education + PIR",
                       N=nrow(gi_ext), Events=sum(gi_ext$death),
                       HR=exp(hr$estimate), Lower=exp(hr$conf.low),
                       Upper=exp(hr$conf.high), P=hr$p.value)},
    run3 = {s <- gi_ext %>% filter(!is.na(any_mvpa));
            f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", adj_ext_pa));
            hr <- tidy(coxph(f, data=s), conf.int=TRUE) %>% filter(term==idx);
            data.frame(Index=gsub("_s","",idx), Model="+ Education + PIR + PA",
                       N=nrow(s), Events=sum(s$death),
                       HR=exp(hr$estimate), Lower=exp(hr$conf.low),
                       Upper=exp(hr$conf.high), P=hr$p.value)}
  ))
}))
ext_res <- ext_res %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
cat("Extended covariate models:\n"); print(ext_res[, c("Index","Model","N","Events","HR_CI","P")])
write.csv(ext_res, file.path(RES_DIR, "sensitivity_extended_covariates.csv"), row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════
# 4. PNI vs ALBUMIN vs LYMPHOCYTE — formal comparison (R1, R2)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== 4. PNI vs albumin vs lymphocyte (decomposition) ===\n")
gi <- gi %>% mutate(alb_s = as.numeric(scale(albumin_gdl)),
                     lym_s = as.numeric(scale(lymph_abs)))

# C-stat with bootstrap CI
boot_cstat_compare <- function(B=2000) {
  set.seed(42)
  res <- replicate(B, {
    idx <- sample(nrow(gi), nrow(gi), replace=TRUE)
    b <- gi[idx, ]
    f0 <- coxph(Surv(surv_years, death) ~ age + sex_b + race_eth_b, data=b)
    f_a <- coxph(Surv(surv_years, death) ~ alb_s + age + sex_b + race_eth_b, data=b)
    f_l <- coxph(Surv(surv_years, death) ~ lym_s + age + sex_b + race_eth_b, data=b)
    f_p <- coxph(Surv(surv_years, death) ~ PNI_s + age + sex_b + race_eth_b, data=b)
    c0 <- as.numeric(f0$concordance["concordance"])
    c(Albumin=as.numeric(f_a$concordance["concordance"]) - c0,
      Lymphocyte=as.numeric(f_l$concordance["concordance"]) - c0,
      PNI=as.numeric(f_p$concordance["concordance"]) - c0)
  })
  t(apply(res, 1, function(x) c(mean=mean(x, na.rm=TRUE),
                                 lower=quantile(x, 0.025, na.rm=TRUE),
                                 upper=quantile(x, 0.975, na.rm=TRUE))))
}
bc <- boot_cstat_compare(2000)
cat("Bootstrap ΔC-stat (2000 reps):\n")
for (nm in rownames(bc)) cat(sprintf("  %s: %.4f (%.4f-%.4f)\n", nm, bc[nm,1], bc[nm,2], bc[nm,3]))

# Likelihood ratio test: PNI vs albumin nested models
f_alb <- coxph(Surv(surv_years, death) ~ alb_s + age + sex_b + race_eth_b, data=gi)
f_pni <- coxph(Surv(surv_years, death) ~ PNI_s + age + sex_b + race_eth_b, data=gi)
f_both <- coxph(Surv(surv_years, death) ~ alb_s + lym_s + age + sex_b + race_eth_b, data=gi)
lrt_alb_vs_pni <- anova(f_alb, f_pni)
lrt_alb_vs_both <- anova(f_alb, f_both)
cat(sprintf("  LRT: albumin (AIC=%.1f) vs PNI (AIC=%.1f): p=%.4f\n",
    AIC(f_alb), AIC(f_pni), lrt_alb_vs_pni$`Pr(>|Chi|)`[2]))
cat(sprintf("  LRT: albumin+lymph (AIC=%.1f) vs albumin (AIC=%.1f): p=%.4f\n",
    AIC(f_both), AIC(f_alb), lrt_alb_vs_both$`Pr(>|Chi|)`[2]))
write.csv(data.frame(
  Model=c("Albumin","PNI","Albumin+Lymph"),
  AIC=round(c(AIC(f_alb), AIC(f_pni), AIC(f_both)), 1),
  DeltaCstat=paste0(sprintf("%.4f", bc[c("Albumin","PNI","Lymphocyte"),1]),
                     " (", sprintf("%.4f", bc[c("Albumin","PNI","Lymphocyte"),2]),
                     "-", sprintf("%.4f", bc[c("Albumin","PNI","Lymphocyte"),3]), ")")),
  file.path(RES_DIR, "sensitivity_pni_decomposition.csv"), row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════
# 5. CROSS-VALIDATED PNI THRESHOLD (R1, R3)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== 5. Cross-validated PNI threshold ===\n")
suppressPackageStartupMessages(library(maxstat))

# 5-fold cross-validated threshold
set.seed(42)
folds <- sample(rep(1:5, length.out=nrow(gi)))
cv_thresholds <- sapply(1:5, function(k) {
  train <- gi[folds != k, ]
  mt <- maxstat.test(Surv(surv_years, death) ~ PNI, data=train,
                     smethod="LogRank", pmethod="Lau92", minprop=0.2, maxprop=0.8)
  as.numeric(mt$estimate)
})
cat(sprintf("  5-fold CV thresholds: %s\n", paste(sprintf("%.1f", cv_thresholds), collapse=", ")))
cat(sprintf("  Mean CV threshold: %.1f (SD=%.1f)\n", mean(cv_thresholds), sd(cv_thresholds)))

# Bootstrap threshold stability
boot_thresholds <- replicate(500, {
  idx <- sample(nrow(gi), nrow(gi), replace=TRUE)
  b <- gi[idx, ]
  mt <- tryCatch(maxstat.test(Surv(surv_years, death) ~ PNI, data=b,
                              smethod="LogRank", pmethod="Lau92", minprop=0.2, maxprop=0.8),
                 error=function(e) NULL)
  if (is.null(mt)) return(NA)
  as.numeric(mt$estimate)
})
boot_thresholds <- boot_thresholds[!is.na(boot_thresholds)]
cat(sprintf("  Bootstrap threshold (500 reps): median=%.1f, 95%% CI (%.1f-%.1f)\n",
    median(boot_thresholds), quantile(boot_thresholds, 0.025), quantile(boot_thresholds, 0.975)))
write.csv(data.frame(
  Optimal_threshold=48.5,
  CV_mean=mean(cv_thresholds), CV_sd=sd(cv_thresholds),
  CV_values=paste(sprintf("%.1f", cv_thresholds), collapse=","),
  Boot_median=median(boot_thresholds),
  Boot_CI_lower=quantile(boot_thresholds, 0.025),
  Boot_CI_upper=quantile(boot_thresholds, 0.975)),
  file.path(RES_DIR, "sensitivity_threshold_validation.csv"), row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════
# 6. COHORT AS FIXED EFFECT (BCG/BCP proxy) (R2)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== 6. Cohort as fixed effect (BCG/BCP sensitivity) ===\n")
for (idx in c("PNI_s","CONUT_s","GNRI_s")) {
  f_base <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", adj_cov))
  f_cohort <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", adj_cov, "+ cohort_b"))
  h_base <- tidy(coxph(f_base, data=gi), conf.int=TRUE) %>% filter(term==idx)
  h_cohort <- tidy(coxph(f_cohort, data=gi), conf.int=TRUE) %>% filter(term==idx)
  cat(sprintf("  %s: w/o cohort HR=%.3f (%.3f-%.3f), w/ cohort HR=%.3f (%.3f-%.3f)\n",
      gsub("_s","",idx),
      exp(h_base$estimate), exp(h_base$conf.low), exp(h_base$conf.high),
      exp(h_cohort$estimate), exp(h_cohort$conf.low), exp(h_cohort$conf.high)))
}

# ═══════════════════════════════════════════════════════════════════════════
# 7. REVERSE CAUSATION: EXCLUDE DEATHS < 2 YEARS (R4)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== 7. Excluding early deaths (<2 yr) ===\n")
gi_2yr <- gi %>% filter(surv_years >= 2)
cat(sprintf("  Excluding <2yr: N=%d, events=%d\n", nrow(gi_2yr), sum(gi_2yr$death)))

for (idx in c("PNI_s","CONUT_s","GNRI_s")) {
  f <- as.formula(paste("Surv(surv_years-2, death) ~", idx, "+", adj_cov))
  hr <- tidy(coxph(f, data=gi_2yr), conf.int=TRUE) %>% filter(term==idx)
  cat(sprintf("  %s: HR=%.3f (%.3f-%.3f), p=%.4f\n",
      gsub("_s","",idx), exp(hr$estimate), exp(hr$conf.low), exp(hr$conf.high), hr$p.value))
}

# ═══════════════════════════════════════════════════════════════════════════
# 8. ADDITIONAL: COHORT STRATIFIED + EXTENDED COVARIATES (R4)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== 8. Continuous NHANES only: extended covariates ===\n")
gi_cont <- gi %>% filter(cohort == "Continuous_NHANES") %>% filter(!is.na(edu_binary))
cat(sprintf("  Cont NHANES N=%d, events=%d\n", nrow(gi_cont), sum(gi_cont$death)))

adj_full <- "age + sex_b + race_eth_b + edu_binary + income_pir"
if (sum(!is.na(gi_cont$any_mvpa)) > 50) adj_full <- paste0(adj_full, " + any_mvpa")

for (idx in c("PNI_s","CONUT_s","GNRI_s")) {
  f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", adj_full))
  hr <- tryCatch({h <- tidy(coxph(f, data=gi_cont), conf.int=TRUE) %>% filter(term==idx);
                  h}, error=function(e) NULL)
  if (!is.null(hr)) cat(sprintf("  %s: HR=%.3f (%.3f-%.3f), p=%.4f [+edu+PIR+PA]\n",
      gsub("_s","",idx), exp(hr$estimate), exp(hr$conf.low), exp(hr$conf.high), hr$p.value))
}

cat("\n=== 04_revision_sensitivity.R COMPLETE ===\n")
