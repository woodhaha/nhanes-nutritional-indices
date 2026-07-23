# run_gi_survival.R — NHANES GI tumor survival analysis
# Run after run_gi_data.R completes.
# Rscript run_gi_survival.R

library(dplyr)
library(survival)
library(broom)
if (!require("survRM2")) install.packages("survRM2", repos="https://cloud.r-project.org")
library(survRM2)

GI_DATA_DIR <- file.path(dirname(normalizePath(getwd())), "data", "gi_analysis")
GI_RESULTS_DIR <- file.path(dirname(normalizePath(getwd())), "results", "gi_analysis")
dir.create(GI_RESULTS_DIR, recursive=TRUE, showWarnings=FALSE)

cat("Loading nutrition data...\n")
df <- readRDS(file.path(GI_DATA_DIR, "nhanes_gi_nutrition_raw.rds"))
cat(sprintf("  N=%d, GI=%d\n", nrow(df), sum(df$gi_tumor)))
# Re-set age if dropped
if (!"age" %in% names(df)) df$age <- df$RIDAGEYR

# ── Mortality linkage ──────────────────────────────────────────────────────────
cat("Loading mortality data...\n")
mort_cache <- file.path(GI_DATA_DIR, "mort_2019.rds")
if (file.exists(mort_cache)) {
  mort <- readRDS(mort_cache)
} else {
  cat("Downloading mortality files...\n")
  base_url <- "https://ftp.cdc.gov/pub/health_statistics/nchs/datalinkage/linked_mortality"
  cycles_mort <- c("2005_2006","2007_2008","2009_2010","2011_2012","2013_2014","2015_2016")
  all_m <- list()
  for (cyc in cycles_mort) {
    fname <- sprintf("NHANES_%s_MORT_2019_PUBLIC.dat", cyc)
    tmp <- tempfile(fileext=".dat")
    cat(sprintf("  %s...", fname))
    download.file(file.path(base_url, fname), tmp, mode="wb", quiet=TRUE)
    m <- read.fwf(tmp,
      widths=c(6,8,1,1,3,1,1,21,3,3,13),
      col.names=c("SEQN","b1","eligstat","mortstat","ucod_leading",
                  "diabetes","hyperten","b2","permth_int","permth_exm","b3"),
      colClasses=c("integer","NULL","integer","integer","character",
                   "integer","integer","NULL","integer","integer","NULL"),
      na.strings=c("","."," "))
    m$ucod_leading[m$ucod_leading==""] <- NA_character_
    m$ucod_leading <- as.integer(m$ucod_leading)
    all_m[[cyc]] <- m
    cat(sprintf(" %d rows\n", nrow(m)))
    unlink(tmp)
  }
  mort <- bind_rows(all_m)
  saveRDS(mort, mort_cache)
}

# Merge mortality
df <- merge(df, mort[, c("SEQN", "eligstat", "mortstat", "ucod_leading",
                          "permth_int", "permth_exm")],
            by="SEQN", all.x=TRUE)

cat(sprintf("After merge: N=%d\n", nrow(df)))
cat(sprintf("  Eligible: %d\n", sum(df$eligstat==1, na.rm=TRUE)))

# Filter to eligible
df <- df[df$eligstat==1 & !is.na(df$mortstat), , drop=FALSE]
cat(sprintf("  Eligible + known vital status: %d\n", nrow(df)))

# Survival variables
df$surv_years <- df$permth_int / 12
df$death <- df$mortstat

# Cause of death: UCOD_LEADING=002 = malignant neoplasm (C00-C97)
# Proxy for GI cancer death: baseline GI tumor + cancer death
df$gi_cancer_death <- ifelse(df$death==1 & df$ucod_leading==2 & df$gi_tumor==1, 1L, 0L)
df$other_death <- ifelse(df$death==1 & (df$gi_cancer_death!=1 | is.na(df$gi_cancer_death)), 1L, 0L)
df$cod_status <- ifelse(df$death==0, 0L, ifelse(df$gi_cancer_death==1, 1L, 2L))

cat(sprintf("  Deaths: %d (%.1f%%)\n", sum(df$death, na.rm=TRUE), mean(df$death,na.rm=TRUE)*100))
cat(sprintf("  GI cancer deaths: %d\n", sum(df$gi_cancer_death, na.rm=TRUE)))

# Complete nutrition
df <- df[!is.na(df$PNI) & !is.na(df$CONUT) & !is.na(df$GNRI), , drop=FALSE]
cat(sprintf("  Final analytic: N=%d\n", nrow(df)))
cat(sprintf("  GI tumor patients: %d\n", sum(df$gi_tumor==1, na.rm=TRUE)))

# ═══ Table 1: Descriptive ═══════════════════════════════════════════════════════
cat("\n── Table 1 ──\n")
df$gi_status <- case_when(
  df$gi_tumor==1 ~ "GI Tumor",
  df$any_cancer==0 ~ "Non-Cancer",
  TRUE ~ "Other Cancer"
)

t1 <- df %>%
  group_by(gi_status) %>%
  summarise(
    N = n(),
    Age = sprintf("%.1f (%.1f)", mean(age, na.rm=TRUE), sd(age, na.rm=TRUE)),
    Female = sprintf("%.1f%%", mean(sex=="Female", na.rm=TRUE)*100),
    PNI = sprintf("%.1f (%.1f)", mean(PNI, na.rm=TRUE), sd(PNI, na.rm=TRUE)),
    CONUT = sprintf("%.1f (%.1f)", mean(CONUT, na.rm=TRUE), sd(CONUT, na.rm=TRUE)),
    GNRI = sprintf("%.1f (%.1f)", mean(GNRI, na.rm=TRUE), sd(GNRI, na.rm=TRUE)),
    Deaths = sprintf("%d (%.1f%%)", sum(death, na.rm=TRUE), mean(death, na.rm=TRUE)*100)
  )
print(as.data.frame(t1))
write.csv(t1, file.path(GI_RESULTS_DIR, "table1_descriptive.csv"), row.names=FALSE)

# GI site breakdown
cat("\nGI site breakdown:\n")
print(table(df$gi_site[df$gi_tumor==1]))

# ═══ Cross-sectional: nutrition × GI tumor ═════════════════════════════════════
cat("\n── Cross-sectional: Nutrition × GI Tumor (vs Non-cancer) ──\n")
df_cs <- df[df$gi_status %in% c("GI Tumor", "Non-Cancer"), ]
cs_results <- data.frame()
for (nv in c("PNI", "CONUT", "GNRI")) {
  df_cs$sex_b <- ifelse(df_cs$sex=="Female", 1L, 0L)
  df_cs$race_eth_b <- ifelse(df_cs$race_eth=="Non-Hispanic White", 1L, 0L)
  m_crude <- lm(as.formula(paste(nv, "~ gi_tumor")), data=df_cs)
  b <- tidy(m_crude) %>% filter(term=="gi_tumor")
  m_adj <- lm(as.formula(paste(nv, "~ gi_tumor + age + sex_b + race_eth_b + edu_binary")), data=df_cs)
  b2 <- tidy(m_adj) %>% filter(term=="gi_tumor")
  cs_results <- rbind(cs_results,
    data.frame(Nutrition=nv, Model="Crude", Beta=b$estimate, SE=b$std.error, P=b$p.value),
    data.frame(Nutrition=nv, Model="Adjusted", Beta=b2$estimate, SE=b2$std.error, P=b2$p.value))
}
cs_results$CI95 <- sprintf("%.3f (%.3f, %.3f)", cs_results$Beta,
                            cs_results$Beta-1.96*cs_results$SE,
                            cs_results$Beta+1.96*cs_results$SE)
print(cs_results)
write.csv(cs_results, file.path(GI_RESULTS_DIR, "cross_sectional.csv"), row.names=FALSE)

# ═══ All-cause survival (GI tumor patients) ════════════════════════════════════
cat("\n── All-cause Survival (GI tumor patients) ──\n")
df_gi <- df[df$gi_tumor==1 & df$surv_years>0, , drop=FALSE]
df_gi$sex_b <- ifelse(df_gi$sex=="Female", 1L, 0L)
df_gi$race_eth_b <- ifelse(df_gi$race_eth=="Non-Hispanic White", 1L, 0L)
n_gi <- nrow(df_gi)
n_events <- sum(df_gi$death, na.rm=TRUE)
cat(sprintf("  N=%d, events=%d\n", n_gi, n_events))

if (n_events >= 10) {
  df_gi$PNI_s <- as.numeric(scale(df_gi$PNI))
  df_gi$CONUT_s <- as.numeric(scale(-df_gi$CONUT))
  df_gi$GNRI_s <- as.numeric(scale(df_gi$GNRI))

  cox_res <- data.frame()
  for (exp in c("PNI_s", "CONUT_s", "GNRI_s")) {
    for (adj in c("Crude", "Adjusted")) {
      covs <- if (adj=="Crude") "1" else "age + sex_b + race_eth_b + edu_binary"
      f <- as.formula(paste("Surv(surv_years, death) ~", exp, "+", covs))
      fit <- tryCatch(coxph(f, data=df_gi), error=\(e)NULL)
      if (is.null(fit)) next
      hr <- tidy(fit, conf.int=TRUE) %>% filter(term==exp)
      cox_res <- rbind(cox_res, data.frame(
        Nutrition = gsub("_s","",exp), Adjustment=adj,
        HR=exp(hr$estimate), Lower=exp(hr$conf.low), Upper=exp(hr$conf.high),
        P=hr$p.value, N=n_gi, Events=n_events))
    }
  }
  cox_res$HR_CI <- sprintf("%.3f (%.3f-%.3f)", cox_res$HR, cox_res$Lower, cox_res$Upper)
  print(cox_res[, c("Nutrition","Adjustment","HR_CI","P")])
  write.csv(cox_res, file.path(GI_RESULTS_DIR, "cox_allcause.csv"), row.names=FALSE)

  # PH check
  cat("\nProportional hazards test:\n")
  for (exp in c("PNI_s", "CONUT_s", "GNRI_s")) {
    fit <- coxph(as.formula(paste("Surv(surv_years, death) ~", exp,
                                   "+ age + sex_b + race_eth_b + edu_binary")),
                 data=df_gi)
    cat(sprintf("  %s: p=%.4f\n", exp, cox.zph(fit)$table[1,3]))
  }

  # ─── Time-dependent Cox (allow HR to vary by follow-up period) ─────────────
  # PH violated for PNI and GNRI → split time at 2 years
  cat("\n── Time-dependent Cox (splitting at 2 years, PH violation fix) ──\n")
  df_gi_td <- survSplit(Surv(surv_years, death) ~ PNI_s + CONUT_s + GNRI_s +
                         age + sex_b + race_eth_b + edu_binary,
                        data=df_gi, cut=c(2, 5), episode="period",
                        id="id_td")
  df_gi_td$period_f <- factor(df_gi_td$period)

  td_results <- data.frame()
  for (exp in c("PNI_s", "CONUT_s", "GNRI_s")) {
    f <- as.formula(paste("Surv(tstart, surv_years, death) ~", exp,
                          "* period_f + age + sex_b + race_eth_b + edu_binary + cluster(id_td)"))
    fit_td <- coxph(f, data=df_gi_td)
    hr_main <- tidy(fit_td, conf.int=TRUE) %>% filter(term == exp)
    hr_int <- tidy(fit_td, conf.int=TRUE) %>%
      filter(grepl(paste0(exp, ":"), term))
    td_results <- rbind(td_results,
      data.frame(Nutrition=gsub("_s","",exp), Period="0-2 yr (main)",
                 HR=exp(hr_main$estimate), Lower=exp(hr_main$conf.low),
                 Upper=exp(hr_main$conf.high), P=hr_main$p.value))
    for (k in seq_len(nrow(hr_int))) {
      period_label <- gsub(paste0(exp, ":period_f"), "", hr_int$term[k])
      td_results <- rbind(td_results,
        data.frame(Nutrition=gsub("_s","",exp),
                   Period=paste0("vs main: period", period_label),
                   HR=exp(hr_int$estimate[k]), Lower=exp(hr_int$conf.low[k]),
                   Upper=exp(hr_int$conf.high[k]), P=hr_int$p.value[k]))
    }
  }
  td_results$HR_CI <- sprintf("%.3f (%.3f-%.3f)", td_results$HR,
                               td_results$Lower, td_results$Upper)
  print(td_results[, c("Nutrition","Period","HR_CI","P")])
  write.csv(td_results, file.path(GI_RESULTS_DIR, "cox_timedependent.csv"), row.names=FALSE)

  # ─── Landmark analysis ─────────────────────────────────────────────────────
  cat("\n── Landmark analysis (1, 3, 5 year) ──\n")
  lm_results <- data.frame()
  for (lm_time in c(1, 3, 5)) {
    df_lm <- df_gi[df_gi$surv_years >= lm_time, , drop=FALSE]
    df_lm$surv_remaining <- df_lm$surv_years - lm_time
    n_lm <- nrow(df_lm)
    e_lm <- sum(df_lm$death, na.rm=TRUE)
    if (e_lm < 10) next
    for (exp in c("PNI_s", "CONUT_s", "GNRI_s")) {
      f <- as.formula(paste("Surv(surv_remaining, death) ~", exp,
                            "+ age + sex_b + race_eth_b + edu_binary"))
      fit_lm <- coxph(f, data=df_lm)
      hr <- tidy(fit_lm, conf.int=TRUE) %>% filter(term==exp)
      lm_results <- rbind(lm_results, data.frame(
        Landmark=sprintf("%d-yr", lm_time), Nutrition=gsub("_s","",exp),
        N=n_lm, Events=e_lm,
        HR=exp(hr$estimate), Lower=exp(hr$conf.low),
        Upper=exp(hr$conf.high), P=hr$p.value))
    }
  }
  if (nrow(lm_results) > 0) {
    lm_results$HR_CI <- sprintf("%.3f (%.3f-%.3f)", lm_results$HR,
                                 lm_results$Lower, lm_results$Upper)
    print(lm_results[, c("Nutrition","Landmark","N","Events","HR_CI","P")])
    write.csv(lm_results, file.path(GI_RESULTS_DIR, "cox_landmark.csv"), row.names=FALSE)
  }

  # ─── Restricted Mean Survival Time (RMST) ──────────────────────────────────
  cat("\n── RMST at 5 and 10 years ──\n")
  rmst_results <- data.frame()
  for (rmst_t in c(5, 10)) {
    df_gi$pni_med <- ifelse(df_gi$PNI >= median(df_gi$PNI, na.rm=TRUE), "High", "Low")
    km_rmst <- survfit(Surv(surv_years, death) ~ pni_med, data=df_gi)
    # Manual RMST calculation via area under KM
    rmst_est <- tryCatch({
      s <- summary(km_rmst, time=seq(0, rmst_t, length.out=10000))
      # Split by strata
      high_t <- s$time[s$strata=="pni_med=High"]
      high_s <- s$surv[s$strata=="pni_med=High"]
      low_t <- s$time[s$strata=="pni_med=Low"]
      low_s <- s$surv[s$strata=="pni_med=Low"]
      # AUC
      rmst_high <- sum(diff(c(0, high_t)) * head(c(1, high_s), -1))
      rmst_low <- sum(diff(c(0, low_t)) * head(c(1, low_s), -1))
      c(high=rmst_high, low=rmst_low, diff=rmst_high-rmst_low)
    }, error=function(e) NULL)
    if (!is.null(rmst_est)) {
      rmst_results <- rbind(rmst_results, data.frame(
        Tau=sprintf("%d-yr", rmst_t),
        High_PNI_RMST=sprintf("%.2f yr", rmst_est["high"]),
        Low_PNI_RMST=sprintf("%.2f yr", rmst_est["low"]),
        Difference=sprintf("%.2f yr", rmst_est["diff"])))
    }
  }
  if (nrow(rmst_results) > 0) {
    print(rmst_results)
    write.csv(rmst_results, file.path(GI_RESULTS_DIR, "rmst_results.csv"), row.names=FALSE)
  }

  # ─── Era analysis: interaction of nutrition × year of diagnosis era ────────
  cat("\n── Era interaction (diagnosis period) ──\n")
  df_gi$era <- ifelse(df_gi$cycle %in% c("2005-2006","2007-2008","2009-2010"),
                      "2005-2010 (earlier)", "2011-2016 (later)")
  cat("  Era distribution:\n")
  print(table(df_gi$era))
  for (exp in c("PNI_s", "CONUT_s", "GNRI_s")) {
    f <- as.formula(paste("Surv(surv_years, death) ~", exp, "* era + age + sex_b + race_eth_b"))
    fit_era <- coxph(f, data=df_gi)
    int_term <- grep(paste0(exp, ":"), names(coef(fit_era)), value=TRUE)
    if (length(int_term) > 0) {
      hr_int <- tidy(fit_era, conf.int=TRUE) %>% filter(term %in% int_term)
      cat(sprintf("  %s × era: HR=%.3f, p=%.4f\n", exp,
                  exp(hr_int$estimate), hr_int$p.value))
    }
  }

  # KM median survival by PNI tertile
  df_gi$pni_t <- factor(ntile(df_gi$PNI, 3), 1:3, c("Low (T1)", "Mid (T2)", "High (T3)"))
  km <- survfit(Surv(surv_years, death) ~ pni_t, data=df_gi)
  cat("\nMedian survival by PNI tertile:\n")
  print(km)
  saveRDS(km, file.path(GI_RESULTS_DIR, "km_fit.rds"))

  # ─── Competing risks ──────────────────────────────────────────────────────────
  cat("\n── Competing risks (GI cancer vs other death) ──\n")
  n_gc <- sum(df_gi$gi_cancer_death, na.rm=TRUE)
  cat(sprintf("  GI cancer deaths: %d\n", n_gc))

  if (n_gc >= 5) {
    cs_cox <- data.frame()
    for (cause_nm in c("gi_cancer_death", "other_death")) {
      for (exp in c("PNI_s", "CONUT_s", "GNRI_s")) {
        f <- as.formula(paste("Surv(surv_years,", cause_nm, ") ~", exp,
                              "+ age + sex_b + race_eth_b"))
        fit <- coxph(f, data=df_gi)
        hr <- tidy(fit, conf.int=TRUE) %>% filter(term==exp)
        cs_cox <- rbind(cs_cox, data.frame(
          Nutrition=gsub("_s","",exp),
          Cause=ifelse(cause_nm=="gi_cancer_death","GI cancer death","Other death"),
          HR=exp(hr$estimate), Lower=exp(hr$conf.low), Upper=exp(hr$conf.high),
          P=hr$p.value))
      }
    }
    cs_cox$HR_CI <- sprintf("%.3f (%.3f-%.3f)", cs_cox$HR, cs_cox$Lower, cs_cox$Upper)
    print(cs_cox[, c("Nutrition","Cause","HR_CI","P")])
    write.csv(cs_cox, file.path(GI_RESULTS_DIR, "cox_causespecific.csv"), row.names=FALSE)
  } else {
    cat("  Too few GI cancer deaths for competing risks\n")
  }
} else {
  cat("  Too few events for survival analysis\n")
}

cat(sprintf("\nResults saved to: %s\n", GI_RESULTS_DIR))
cat("Done.\n")
