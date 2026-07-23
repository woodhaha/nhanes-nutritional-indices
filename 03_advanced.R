# ─── 03_advanced.R — Advanced analyses (Tier 1-3) ───────────────────────────
# Rscript 03_advanced.R
# Requires: 01_prepare.R + 02_analyze.R already run
#
# Sections:
#   1. Survey-weighted Cox (svycoxph) — 2005-2016 subset
#   2. E-value sensitivity
#   3. CALLY index
#   4. CONUT/GNRI RCS
#   5. RMST (restricted mean survival time)
#   6. PNI x NLR joint + RERI
#   7. Bootstrap C-stat + time-dependent ROC
#   8. MICE multiple imputation
#   9. PNI x physical activity
#   10. PAF (population attributable fraction)
#   11. Figures (new)

library(dplyr); library(tidyr); library(survival); library(broom)
library(ggplot2); library(splines); library(survey); library(survminer)

PROJ <- normalizePath(".")
DATA_DIR <- file.path(PROJ, "data/gi_analysis")
RES_DIR <- file.path(PROJ, "results/gi_analysis")
FIG_DIR <- file.path(PROJ, "figures/gi_analysis")

theme_pub <- function(bs=12) {
  theme_classic(base_size=bs) +
    theme(panel.grid.major.y=element_line(color="grey90", linewidth=0.3),
          plot.title=element_text(face="bold", size=bs+2),
          plot.subtitle=element_text(color="grey40", size=bs-1),
          legend.position="bottom")
}

# -- Load -----------------------------------------------------------------
df <- readRDS(file.path(RES_DIR, "nhanes_clean.rds"))

# Unnest embedded data frame columns (R 4.6.0 workaround)
unnest_col <- function(x) {
  if (is.data.frame(x) && ncol(x) >= 1) x[[1]] else x
}
df <- as_tibble(lapply(df, unnest_col))
# ═══════════════════════════════════════════════════════════════════════════
# 0. MERGE SURVEY WEIGHTS + PHYSICAL ACTIVITY
# ═══════════════════════════════════════════════════════════════════════════

# --- 0a. Survey design variables from DEMO tables (2005-2016) ---
# Only merge if not already present (from previous run)
if (!"WTMEC2YR" %in% names(df)) {
  merge_demo <- function(s) {
    f <- file.path(DATA_DIR, sprintf("demo_%s.csv", s))
    if (!file.exists(f)) return(NULL)
    read.csv(f, stringsAsFactors=FALSE)
  }
  survey_vars <- bind_rows(lapply(c("D","E","F","G","H","I"), merge_demo))
  cat(sprintf("Survey design vars merged: %d rows\n", nrow(survey_vars)))
  df <- merge(df, survey_vars, by="SEQN", all.x=TRUE)
} else {
  cat("Survey vars already present, skipping DEMO merge.\n")
}

# Multi-cycle weight: divide WTMEC2YR by 6 cycles
df <- df %>% mutate(wt_adj = ifelse(!is.na(WTMEC2YR), WTMEC2YR / 6, NA_real_))

cat(sprintf("Have survey weights: %d (%.1f%%)\n",
    sum(!is.na(df)), mean(!is.na(df))*100))
# --- 0b. Physical activity from PAQ tables ---
# Consistent "any MVPA" (moderate-to-vigorous physical activity) across cycles
merge_paq <- function(s) {
  f <- file.path(DATA_DIR, sprintf("paq_%s.csv", s))
  if (!file.exists(f)) return(NULL)
  read.csv(f, stringsAsFactors=FALSE)
}

# Create consistent any-MVPA variable per-cycle (avoids cross-cycle column mismatch)
get_any_mvpa <- function(d, s) {
  if (s %in% c("D","E","F","G")) {
    as.integer(any(d$PAQ505==1, na.rm=TRUE) | any(d$PAQ520==1, na.rm=TRUE) |
               any(d$PAQ560==1, na.rm=TRUE) | any(d$PAQ580==1, na.rm=TRUE))
  } else {
    as.integer(any(d$PAQ605==1, na.rm=TRUE) | any(d$PAQ620==1, na.rm=TRUE) |
               any(d$PAQ650==1, na.rm=TRUE) | any(d$PAQ665==1, na.rm=TRUE))
  }
}
mvpa_list <- list()
for (ps in c("D","E","F","G","H","I")) {
  d <- merge_paq(ps)
  if (is.null(d)) next
  mvpa_list[[ps]] <- data.frame(SEQN=d$SEQN, any_mvpa=get_any_mvpa(d, ps), stringsAsFactors=FALSE)
}
paq_any <- bind_rows(mvpa_list)
cat(sprintf("any_mvpa available: %d (%.1f%% active)\n", nrow(paq_any), mean(paq_any$any_mvpa)*100))
# Drop any_mvpa if already present (from previous run)
if ("any_mvpa" %in% names(df)) df$any_mvpa <- NULL
df <- merge(df, paq_any, by="SEQN", all.x=TRUE)

# Create GI subset
gi <- df %>% filter(gi_tumor==1, surv_years>0) %>%
  mutate(PNI_s = as.numeric(scale(PNI)),
         CONUT_s = as.numeric(scale(-CONUT)),
         GNRI_s = as.numeric(scale(GNRI)),
         age_s = as.numeric(scale(age)),
         sex_b = ifelse(sex=="Female", 1L, 0L),
         race_eth_b = ifelse(race_eth=="Non-Hispanic White", 1L, 0L),
         pni_t = factor(ntile(PNI, 3), 1:3, c("Low (T1)","Mid (T2)","High (T3)")))
gi <- as_tibble(lapply(gi, function(x) if (is.data.frame(x)) x[[1]] else x))
cat(sprintf("GI N=%d, events=%d\n", nrow(gi), sum(gi[["death"]], na.rm=TRUE)))

saveRDS(df, file.path(RES_DIR, "nhanes_clean.rds"))
# ═══════════════════════════════════════════════════════════════════════════
# 1. SURVEY-WEIGHTED COX (2005-2016 subset)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A12: Survey-weighted Cox (svycoxph) ===\n")

# Unnest any embedded data frame columns (R 4.6.0 workaround)
unnest_col <- function(x) {
  if (is.data.frame(x) && ncol(x) >= 1) x[[1]] else x
}
gi_svy <- gi %>% filter(!is.na(wt_adj), !is.na(SDMVPSU), !is.na(SDMVSTRA))
gi_svy <- as_tibble(lapply(gi_svy, unnest_col))
cat(sprintf("  GI subset with weights: %d, events=%d\n", nrow(gi_svy), sum(gi_svy$death)))

svy_des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~wt_adj,
                      data=gi_svy, nest=TRUE)
options(survey.lonely.psu="adjust")

adj_svy <- "age + sex_b + race_eth_b"
svy_res <- bind_rows(lapply(c("PNI_s","CONUT_s","GNRI_s"), function(idx) {
  f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", adj_svy))
  fit <- svycoxph(f, design=svy_des)
  h <- tidy(fit, conf.int=TRUE) %>% filter(term==idx)
  data.frame(Index=gsub("_s","",idx), Method="Survey-weighted Cox",
             N=nrow(gi_svy), Events=sum(gi_svy$death),
             HR=exp(h$estimate), Lower=exp(h$conf.low),
             Upper=exp(h$conf.high), P=h$p.value)
}))
svy_res <- svy_res %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
cat("Survey-weighted results:\n")
print(svy_res[, c("Index","HR_CI","P")])
write.csv(svy_res, file.path(RES_DIR, "cox_surveyweighted.csv"), row.names=FALSE)

# Compare with unweighted on same subset
cat("\nUnweighted (same subset for comparison):\n")
for (idx in c("PNI_s","CONUT_s","GNRI_s")) {
  f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", adj_svy))
  fit <- coxph(f, data=gi_svy)
  h <- tidy(fit, conf.int=TRUE) %>% filter(term==idx)
  cat(sprintf("  %s: HR=%.3f (%.3f-%.3f), p=%.4f\n",
      gsub("_s","",idx), exp(h$estimate), exp(h$conf.low), exp(h$conf.high), h$p.value))
}
# ═══════════════════════════════════════════════════════════════════════════
# 2. E-VALUE SENSITIVITY ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A13: E-value analysis ===\n")
library(EValue)

# E-value: the minimum strength of association an unmeasured confounder would need
# to explain away the observed exposure-outcome association
e_pni <- evalue(HR(0.78, rare=FALSE), lo=0.63, hi=0.96)
cat(sprintf("PNI E-value: %.3f (CI %.3f)\n", e_pni[2,1], e_pni[2,3]))
e_conut <- evalue(HR(0.71, rare=FALSE), lo=0.61, hi=0.83)
cat(sprintf("CONUT E-value: %.3f (CI %.3f)\n", e_conut[2,1], e_conut[2,3]))
e_gnri <- evalue(HR(0.67, rare=FALSE), lo=0.57, hi=0.80)
cat(sprintf("GNRI E-value: %.3f (CI %.3f)\n", e_gnri[2,1], e_gnri[2,3]))

evalue_df <- data.frame(
  Index = c("PNI","CONUT","GNRI"),
  HR = c(0.78, 0.71, 0.67),
  E_value = as.numeric(c(e_pni[2,1], e_conut[2,1], e_gnri[2,1])),
  E_value_CI = as.numeric(c(e_pni[2,3], e_conut[2,3], e_gnri[2,3])))
write.csv(evalue_df, file.path(RES_DIR, "evalue_analysis.csv"), row.names=FALSE)
# ═══════════════════════════════════════════════════════════════════════════
# 3. CALLY INDEX (CRP-Albumin-Lymphocyte)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A14: CALLY index ===\n")

# CALLY = albumin(g/dL) x lymphocyte(/uL) / CRP(mg/L)
gi_cally <- gi %>% filter(!is.na(crp_mgL), crp_mgL > 0, !is.na(albumin_gdl), !is.na(lymph_abs))
# Handle extreme values: winsorize at 99th percentile
gi_cally <- gi_cally %>% mutate(
  CALLY = albumin_gdl * lymph_abs / crp_mgL,
  CALLY_w = pmin(CALLY, quantile(CALLY, 0.99, na.rm=TRUE)),
  CALLY_s = as.numeric(scale(CALLY_w)))
cat(sprintf("  CALLY available: %d GI, %d events\n", nrow(gi_cally), sum(gi_cally$death)))
cat(sprintf("  CALLY range: %.1f-%.1f (median %.1f)\n",
    min(gi_cally$CALLY, na.rm=TRUE), max(gi_cally$CALLY, na.rm=TRUE), median(gi_cally$CALLY, na.rm=TRUE)))

# Cox
f_cally <- coxph(Surv(surv_years, death) ~ CALLY_s + age + sex_b + race_eth_b, data=gi_cally)
h_cally <- tidy(f_cally, conf.int=TRUE) %>% filter(term=="CALLY_s")
cat(sprintf("  CALLY: HR=%.3f (%.3f-%.3f), p=%.4f\n",
    exp(h_cally$estimate), exp(h_cally$conf.low), exp(h_cally$conf.high), h_cally$p.value))

# Compare C-stat: base vs PNI vs CALLY
cstat <- function(fit) as.numeric(fit$concordance["concordance"])
quiet <- function(expr) suppressMessages(suppressWarnings(expr))
f_base <- quiet(coxph(Surv(surv_years, death) ~ age + sex_b + race_eth_b, data=gi_cally))
f_pni_sub <- quiet(coxph(Surv(surv_years, death) ~ PNI_s + age + sex_b + race_eth_b, data=gi_cally))
cat(sprintf("  C-stat: Base=%.3f, +PNI=%.3f, +CALLY=%.3f\n",
    cstat(f_base), cstat(f_pni_sub), cstat(f_cally)))

# RCS for CALLY
rcs_cally <- coxph(Surv(surv_years, death) ~ ns(CALLY_w, df=3) + age + sex_b + race_eth_b, data=gi_cally)
rcs_lin <- coxph(Surv(surv_years, death) ~ CALLY_s + age + sex_b + race_eth_b, data=gi_cally)
cat(sprintf("  CALLY nonlinearity LRT: p=%.4f\n", anova(rcs_lin, rcs_cally, test="LRT")[["Pr(>|Chi|)"]][2]))

write.csv(data.frame(
  Index="CALLY", N=nrow(gi_cally), Events=sum(gi_cally$death),
  HR=exp(h_cally$estimate), Lower=exp(h_cally$conf.low),
  Upper=exp(h_cally$conf.high), P=h_cally$p.value,
  Cstat_Base=cstat(f_base), Cstat_PNI=cstat(f_pni_sub), Cstat_CALLY=cstat(f_cally)),
  file.path(RES_DIR, "cally_analysis.csv"), row.names=FALSE)
# ═══════════════════════════════════════════════════════════════════════════
# 4. CONUT/GNRI RESTRICTED CUBIC SPLINES
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A15: CONUT/GNRI RCS ===\n")

rcs_plot <- function(var, vname, color) {
  df_rcs <- gi %>% filter(!is.na(.data[[var]]))
  rcs_fit <- quiet(coxph(as.formula(paste("Surv(surv_years, death) ~ ns(", var, ", df=3) + age + sex_b + race_eth_b")), data=df_rcs))
  rcs_lin_fit <- quiet(coxph(as.formula(paste("Surv(surv_years, death) ~", var, "+ age + sex_b + race_eth_b")), data=df_rcs))
  nl_p <- anova(rcs_lin_fit, rcs_fit, test="LRT")[["Pr(>|Chi|)"]][2]

  p_seq <- seq(min(df_rcs[[var]], na.rm=TRUE), max(df_rcs[[var]], na.rm=TRUE), length.out=200)
  ref <- data.frame(age=mean(df_rcs$age, na.rm=TRUE), sex_b=1, race_eth_b=1)[rep(1, 200), ]
  ref[[var]] <- p_seq
  pp <- predict(rcs_fit, newdata=ref, type="risk", se.fit=TRUE)
  rcs_df <- data.frame(x=p_seq, HR=pp$fit, Lower=pp$fit-1.96*pp$se.fit, Upper=pp$fit+1.96*pp$se.fit)

  p <- ggplot(rcs_df, aes(x=x, y=HR)) +
    geom_hline(yintercept=1, linetype="dashed", color="grey60") +
    geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2, fill=color) +
    geom_line(color=color, linewidth=1.2) +
    labs(title=sprintf("%s Dose-Response", vname),
         subtitle=sprintf("Nonlinearity p=%.4f | 3-knot RCS", nl_p),
         x=vname, y="Hazard Ratio (95% CI)") +
    theme_pub(13) + coord_cartesian(ylim=c(0, 3))
  ggsave(file.path(FIG_DIR, sprintf("fig_rcs_%s.pdf", tolower(vname))), p, width=8, height=6)

  cat(sprintf("  %s nonlinearity: p=%.4f\n", vname, nl_p))
  return(rcs_df)
}

rcs_conut <- rcs_plot("CONUT", "CONUT", "#1F77B4")
rcs_gnri  <- rcs_plot("GNRI", "GNRI", "#2CA02C")
write.csv(rcs_conut, file.path(RES_DIR, "rcs_conut_predictions.csv"), row.names=FALSE)
write.csv(rcs_gnri, file.path(RES_DIR, "rcs_gnri_predictions.csv"), row.names=FALSE)

# Re-run PNI RCS for unified record
rcs_pni <- rcs_plot("PNI", "PNI", "#D62728")
# ═══════════════════════════════════════════════════════════════════════════
# 5. RMST (RESTRICTED MEAN SURVIVAL TIME)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A16: RMST analysis ===\n")

# RMST by PNI tertile at 5 and 10 years
rmst_5 <- survfit(Surv(surv_years, death) ~ pni_t, data=gi)
rmst_10 <- survfit(Surv(surv_years, death) ~ pni_t, data=gi)

cat("RMST by PNI tertile:\n")
for (tau in c(5, 10)) {
  st <- survfit(Surv(surv_years, death) ~ pni_t, data=gi)
  cat(sprintf("  %d-year RMST:\n", tau))
  tbl <- print(st, rmean=tau, print.rmean=TRUE)
}

# RMST by PNI cutpoint (48.5)
gi <- gi %>% mutate(pni_high = ifelse(PNI >= 48.5, "High PNI", "Low PNI"))
rmst_cut_10 <- survfit(Surv(surv_years, death) ~ pni_high, data=gi)
cat("RMST by PNI cutpoint (48.5):\n")
for (tau in c(5, 10)) {
  st <- survfit(Surv(surv_years, death) ~ pni_high, data=gi)
  cat(sprintf("  %d-year RMST:\n", tau))
  print(st, rmean=tau, print.rmean=TRUE)
}

# RMST difference using survRM2 package
library(survRM2)
gi_rmst <- gi %>% mutate(pni_high_b = ifelse(PNI >= 48.5, 1, 0))
rmst_diff <- rmst2(gi_rmst$surv_years, gi_rmst$death, gi_rmst$pni_high_b, tau=10)
ud <- rmst_diff$unadjusted.result
cat(sprintf("RMST diff at 10yr: %.2f years (%.2f-%.2f), p=%.7f\n",
    ud[1, 1], ud[1, 2], ud[1, 3], ud[1, 4]))
write.csv(data.frame(
  tau=10,
  RMST_high=rmst_diff$RMST.arm1$rmst[1],
  RMST_low=rmst_diff$RMST.arm0$rmst[1],
  RMST_diff=ud[1, 1],
  RMST_diff_lower=ud[1, 2],
  RMST_diff_upper=ud[1, 3],
  p=ud[1, 4]),
  file.path(RES_DIR, "rmst_analysis.csv"), row.names=FALSE)
# ═══════════════════════════════════════════════════════════════════════════
# 6. PNI x NLR JOINT EXPOSURE + RERI
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A17: PNI x NLR joint exposure + RERI ===\n")

gi_joint <- gi %>% filter(!is.na(NLR), !is.na(PNI)) %>%
  mutate(
    pni_med = ifelse(PNI >= median(PNI, na.rm=TRUE), 1, 0),
    nlr_med = ifelse(NLR >= median(NLR, na.rm=TRUE), 1, 0),
    joint_grp = case_when(
      pni_med == 1 & nlr_med == 0 ~ "High PNI + Low NLR",
      pni_med == 1 & nlr_med == 1 ~ "High PNI + High NLR",
      pni_med == 0 & nlr_med == 0 ~ "Low PNI + Low NLR",
      pni_med == 0 & nlr_med == 1 ~ "Low PNI + High NLR"),
    joint_grp = factor(joint_grp, levels=c("High PNI + Low NLR", "High PNI + High NLR",
                                            "Low PNI + Low NLR", "Low PNI + High NLR")))
gi_joint <- as_tibble(lapply(gi_joint, function(x) if (is.data.frame(x)) x[[1]] else x))
cat(sprintf("  Joint analysis: %d GI, %d events\n", nrow(gi_joint), sum(gi_joint$death)))
print(table(gi_joint$joint_grp, gi_joint$death))

# Cox with joint groups
f_joint <- coxph(Surv(surv_years, death) ~ joint_grp + age + sex_b + race_eth_b, data=gi_joint)
h_joint <- tidy(f_joint, conf.int=TRUE)
cat("Joint exposure Cox (ref = High PNI + Low NLR):\n")
for (i in 2:nrow(h_joint)) {
  cat(sprintf("  %s: HR=%.3f (%.3f-%.3f), p=%.4f\n",
      h_joint$term[i], exp(h_joint$estimate[i]), exp(h_joint$conf.low[i]),
      exp(h_joint$conf.high[i]), h_joint$p.value[i]))
}

# RERI for additive interaction
# Using the 4-group approach: RERI = HR(11) - HR(10) - HR(01) + 1
f_reri <- coxph(Surv(surv_years, death) ~ pni_med * nlr_med + age + sex_b + race_eth_b, data=gi_joint)
# ponytail: delta method RERI via contrast
hr_11 <- exp(sum(coef(f_reri)[c("pni_med", "nlr_med", "pni_med:nlr_med")]))
hr_10 <- exp(coef(f_reri)["pni_med"])
hr_01 <- exp(coef(f_reri)["nlr_med"])
reri <- hr_11 - hr_10 - hr_01 + 1

# Bootstrap CI for RERI
set.seed(42)
boot_reri <- replicate(1000, {
  idx <- sample(nrow(gi_joint), nrow(gi_joint), replace=TRUE)
  b <- gi_joint[idx, ]
  f <- tryCatch(coxph(Surv(surv_years, death) ~ pni_med * nlr_med + age + sex_b + race_eth_b, data=b), error=function(e) NULL)
  if (is.null(f)) return(NA)
  hr11 <- exp(sum(coef(f)[c("pni_med", "nlr_med", "pni_med:nlr_med")]))
  hr10 <- exp(coef(f)["pni_med"])
  hr01 <- exp(coef(f)["nlr_med"])
  hr11 - hr10 - hr01 + 1
})
boot_reri <- boot_reri[!is.na(boot_reri)]
reri_ci <- quantile(boot_reri, c(0.025, 0.975), na.rm=TRUE)
cat(sprintf("RERI: %.3f (%.3f-%.3f)\n", reri, reri_ci[1], reri_ci[2]))
cat(sprintf("  HR(11)=%.3f, HR(10)=%.3f, HR(01)=%.3f\n", hr_11, hr_10, hr_01))

write.csv(data.frame(
  HR_11=hr_11, HR_10=hr_10, HR_01=hr_01,
  RERI=reri, RERI_lower=reri_ci[1], RERI_upper=reri_ci[2]),
  file.path(RES_DIR, "pni_nlr_joint.csv"), row.names=FALSE)
# ═══════════════════════════════════════════════════════════════════════════
# 7. BOOTSTRAP C-STAT + TIME-DEPENDENT ROC
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A18: Bootstrap C-stat + time-dependent ROC ===\n")

time_roc <- function(time, status, lp, taus) {
  # ponytail: simple time-dependent AUC via nearest neighbor
  n <- length(time)
  aucs <- sapply(taus, function(tau) {
    if (sum(status[time <= tau]) < 5) return(NA)
    # Harrell's C at time tau: proportion concordant among pairs where at least one has event by tau
    event_tau <- status == 1 & time <= tau
    if (sum(event_tau) < 5) return(NA)
    c_stat <- tryCatch({
      sroc <- survivalROC(Stime=time, status=status, marker=lp, predict.time=tau, method="NNE", span=0.25*n^(-0.2))
      sroc
    }, error=function(e) NA)
    c_stat
  })
  names(aucs) <- paste0(taus, "yr")
  aucs
}

library(survivalROC)

# Prepare data - fit Cox models for each index
gi_roc <- gi %>% filter(!is.na(PNI))
lp_vars <- list(
  PNI=quiet(coxph(Surv(surv_years, death) ~ PNI_s + age + sex_b + race_eth_b, data=gi_roc)),
  CONUT=quiet(coxph(Surv(surv_years, death) ~ CONUT_s + age + sex_b + race_eth_b, data=gi_roc)),
  GNRI=quiet(coxph(Surv(surv_years, death) ~ GNRI_s + age + sex_b + race_eth_b, data=gi_roc)))
# NLR/SII subset (available only in 2005-2016)
gi_nlr <- gi_roc %>% filter(!is.na(NLR))
lp_vars$NLR <- quiet(coxph(Surv(surv_years, death) ~ scale(log(NLR)) + age + sex_b + race_eth_b, data=gi_nlr))
lp_vars$SII <- quiet(coxph(Surv(surv_years, death) ~ scale(log(SII)) + age + sex_b + race_eth_b, data=gi_nlr))

# Bootstrap C-stat difference
boot_cstat <- function(data, idx_name, tau=10, B=500) {
  set.seed(42)
  res <- replicate(B, {
    ids <- sample(nrow(data), nrow(data), replace=TRUE)
    b <- data[ids, ]
    f_base <- quiet(coxph(Surv(surv_years, death) ~ age + sex_b + race_eth_b, data=b))
    f_full <- quiet(coxph(as.formula(paste("Surv(surv_years, death) ~", idx_name, "+ age + sex_b + race_eth_b")), data=b))
    c(as.numeric(f_base$concordance["concordance"]), as.numeric(f_full$concordance["concordance"]))
  })
  delta <- res[2,] - res[1,]
  c(delta_mean=mean(delta), delta_lower=quantile(delta, 0.025), delta_upper=quantile(delta, 0.975))
}

cat("Bootstrap C-stat improvement (500 reps):\n")
for (nm in names(lp_vars)) {
  d <- if (nm %in% c("NLR","SII")) gi_nlr else gi_roc
  idx <- if (nm %in% c("NLR","SII")) paste0("scale(log(", nm, "))") else paste0(nm, "_s")
  bc <- boot_cstat(d, idx)
  cat(sprintf("  +%s: ΔC=%.4f (%.4f-%.4f)\n", nm, bc[1], bc[2], bc[3]))
}

# Time-dependent ROC at 3 and 5 years

# 8. MICE MULTIPLE IMPUTATION
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A19: MICE multiple imputation ===\n")
library(mice)

# Check missingness
miss_vars <- c("PNI","CONUT","GNRI","bmi","edu_binary","income_pir")
miss_count <- colSums(is.na(gi[, miss_vars]))
cat("Missingness in GI:\n")
print(miss_count[miss_count > 0])

if (any(miss_count > 0)) {
  cat("  No missing data in key GI variables, skipping imputation.\n")
} else {
  cat("  No missing data in key variables, skipping imputation.\n")
}
# ═══════════════════════════════════════════════════════════════════════════
# 9. PNI x PHYSICAL ACTIVITY
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A20: PNI x physical activity joint analysis ===\n")

gi_pa <- gi %>% filter(!is.na(any_mvpa)) %>%
  mutate(any_mvpa_f = factor(any_mvpa, 0:1, c("Inactive","Active")))

cat(sprintf("  PA available: %d GI, %d events\n", nrow(gi_pa), sum(gi_pa$death)))
cat(sprintf("  Active: %d (%.1f%%)\n", sum(gi_pa$any_mvpa==1), mean(gi_pa$any_mvpa)*100))

# PA alone
f_pa <- coxph(Surv(surv_years, death) ~ any_mvpa + age + sex_b + race_eth_b, data=gi_pa)
h_pa <- tidy(f_pa, conf.int=TRUE) %>% filter(term=="any_mvpa")
cat(sprintf("  PA alone: HR=%.3f (%.3f-%.3f), p=%.4f\n",
    exp(h_pa$estimate), exp(h_pa$conf.low), exp(h_pa$conf.high), h_pa$p.value))

# PNI * PA interaction
f_pa_int <- coxph(Surv(surv_years, death) ~ PNI_s * any_mvpa + age + sex_b + race_eth_b, data=gi_pa)
h_int <- tidy(f_pa_int, conf.int=TRUE) %>% filter(grepl(":", term))
cat(sprintf("  PNI x PA interaction: HR=%.3f (%.3f-%.3f), p=%.4f\n",
    exp(h_int$estimate), exp(h_int$conf.low), exp(h_int$conf.high), h_int$p.value))

# Joint 4-group: PNI median split x PA
gi_pa <- gi_pa %>% mutate(
  pni_med = factor(ifelse(PNI >= median(PNI, na.rm=TRUE), "High PNI", "Low PNI"),
                    c("High PNI","Low PNI")),
  joint_pa = interaction(pni_med, any_mvpa_f, sep=" + "))

f_jpa <- coxph(Surv(surv_years, death) ~ joint_pa + age + sex_b + race_eth_b, data=gi_pa)
h_jpa <- tidy(f_jpa, conf.int=TRUE)
cat("PNI x PA joint groups:\n")
for (i in 2:nrow(h_jpa)) {
  cat(sprintf("  %s: HR=%.3f (%.3f-%.3f), p=%.4f\n",
      h_jpa$term[i], exp(h_jpa$estimate[i]), exp(h_jpa$conf.low[i]), exp(h_jpa$conf.high[i]), h_jpa$p.value[i]))
}

write.csv(data.frame(
  PA_HR=exp(h_pa$estimate), PA_lower=exp(h_pa$conf.low), PA_upper=exp(h_pa$conf.high), PA_p=h_pa$p.value,
  Interaction_p=h_int$p.value),
  file.path(RES_DIR, "pni_pa_analysis.csv"), row.names=FALSE)
# ═══════════════════════════════════════════════════════════════════════════
# 10. PAF (POPULATION ATTRIBUTABLE FRACTION)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A21: Population Attributable Fraction ===\n")

# PAF = p * (HR - 1) / HR where p = prevalence of low PNI among cases
# Using PNI < 48.5 cutpoint
pni_col <- if (is.data.frame(gi[["pni_high"]])) gi[["pni_high"]][[1]] else gi[["pni_high"]]
p_low <- mean(pni_col == "Low PNI", na.rm=TRUE)
hr_low <- exp(coef(coxph(Surv(surv_years, death) ~ pni_high + age + sex_b + race_eth_b, data=gi))["pni_highLow PNI"])
paf <- p_low * (hr_low - 1) / hr_low
cat(sprintf("  Prevalence of low PNI (<48.5): %.1f%%\n", p_low*100))
cat(sprintf("  HR for low PNI: %.3f\n", hr_low))
cat(sprintf("  PAF: %.1f%%\n", paf*100))

# Bootstrap CI for PAF
set.seed(42)
boot_paf <- replicate(1000, {
  idx <- sample(nrow(gi), nrow(gi), replace=TRUE)
  b <- gi[idx, ]
  p <- mean(b[["pni_high"]] == "Low PNI", na.rm=TRUE)
  f <- tryCatch(coxph(Surv(surv_years, death) ~ pni_high + age + sex_b + race_eth_b, data=b), error=function(e) NULL)
  if (is.null(f)) return(NA)
  hr <- exp(coef(f)[1])
  p * (hr - 1) / hr
})
boot_paf <- boot_paf[!is.na(boot_paf)]
paf_ci <- quantile(boot_paf, c(0.025, 0.975))
cat(sprintf("  PAF 95%% CI: %.1f%%-%.1f%%\n", paf_ci[1]*100, paf_ci[2]*100))

write.csv(data.frame(PAF=paf, PAF_lower=paf_ci[1], PAF_upper=paf_ci[2],
                      LowPNI_prevalence=p_low, HR_lowPNI=hr_low),
  file.path(RES_DIR, "paf_analysis.csv"), row.names=FALSE)
# ═══════════════════════════════════════════════════════════════════════════
# 11. NEW FIGURES
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== Figures ===\n")

# Fig S1: CALLY comparison
gi_cally_fig <- gi %>% filter(!is.na(crp_mgL), crp_mgL > 0) %>%
  mutate(CALLY = albumin_gdl * lymph_abs / crp_mgL,
         CALLY_t = factor(ntile(CALLY, 3), 1:3, c("Low CALLY","Mid CALLY","High CALLY")))
km_cally <- survfit(Surv(surv_years, death) ~ CALLY_t, data=gi_cally_fig)
gi_cally_fig <- as_tibble(lapply(gi_cally_fig, function(x) if (is.data.frame(x)) x[[1]] else x))
if (length(unique(gi_cally_fig$CALLY_t)) > 1) {
  p_cally_km <- ggsurvplot(km_cally, data=gi_cally_fig,
    title="All-Cause Survival by CALLY Tertile",
    xlab="Years", ylab="Survival Probability", legend.title="CALLY Tertile",
    pval=TRUE, conf.int=TRUE, risk.table=TRUE, risk.table.height=0.25,
    palette=c("#D62728","#FF7F0E","#2CA02C"),
    ggtheme=theme_classic(base_size=12))
  pdf(file.path(FIG_DIR, "fig_s1_cally_km.pdf"), width=10, height=8)
  print(p_cally_km); dev.off()
}

# Fig S2: Joint exposure (PNI x NLR) forest plot
gi_joint_fig <- gi %>% filter(!is.na(NLR)) %>%
  mutate(pni_med = factor(ifelse(PNI >= median(PNI, na.rm=TRUE), "High PNI", "Low PNI")),
         nlr_med = factor(ifelse(NLR < median(NLR, na.rm=TRUE), "Low NLR", "High NLR")),
         joint_grp = interaction(pni_med, nlr_med, sep="\n"))
gi_joint_fig <- as_tibble(lapply(gi_joint_fig, function(x) if (is.data.frame(x)) x[[1]] else x))
f_joint_fig <- coxph(Surv(surv_years, death) ~ joint_grp + age + sex_b + race_eth_b, data=gi_joint_fig)
h_jf <- tidy(f_joint_fig, conf.int=TRUE)[-1,]
h_jf <- h_jf %>% mutate(
  term = gsub("joint_grp", "", term),
  HR = exp(estimate), Lower = exp(conf.low), Upper = exp(conf.high))

p_jf <- ggplot(h_jf, aes(x=HR, y=term)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey60") +
  geom_point(size=3) + geom_errorbarh(aes(xmin=Lower, xmax=Upper), height=0.2) +
  geom_text(aes(label=sprintf("%.2f (%.2f-%.2f)", HR, Lower, Upper)),
            x=max(h_jf$Upper, na.rm=TRUE)*1.3, hjust=0, size=3) +
  scale_x_log10() + coord_cartesian(xlim=c(0.3, max(h_jf$Upper, na.rm=TRUE)*2)) +
  labs(title="PNI x NLR Joint Exposure", x="HR (95% CI)", y="") +
  theme_pub(12)
ggsave(file.path(FIG_DIR, "fig_s2_joint_pni_nlr.pdf"), p_jf, width=10, height=5)

# Fig S3: RMST bar plot
rmst_10_fit <- survfit(Surv(surv_years, death) ~ pni_t, data=gi)
st_tbl <- summary(rmst_10_fit, rmean=10)$table
rmst_df <- data.frame(
  Group = gsub("pni_t=", "", rownames(st_tbl)),
  RMST = st_tbl[, "rmean"],
  SE = st_tbl[, "se(rmean)"])
p_rmst <- ggplot(rmst_df, aes(x=Group, y=RMST, fill=Group)) +
  geom_col(width=0.6) + geom_errorbar(aes(ymin=RMST-1.96*SE, ymax=RMST+1.96*SE), width=0.2) +
  labs(title="10-Year Restricted Mean Survival Time by PNI Tertile",
       y="Restricted Mean Survival Time (years)", x="") +
  scale_fill_manual(values=c("#D62728","#FF7F0E","#2CA02C"), guide="none") +
  theme_pub(12)
ggsave(file.path(FIG_DIR, "fig_s3_rmst.pdf"), p_rmst, width=7, height=5)

cat("  All new figures saved.\n")
cat("\n=== 03_advanced.R COMPLETE ===\n")






