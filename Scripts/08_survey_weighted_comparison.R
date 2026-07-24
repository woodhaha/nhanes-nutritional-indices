# 08_survey_weighted_comparison.R — Side-by-side unweighted vs weighted Cox
# Produces: combined tables for GI + pancancer with both approaches
# Output: results/gi_analysis/cox_combined.csv, results/pancancer/cox_combined.csv

setwd("D:/Researching/NHANES")
cat("WD:", getwd(), "\n")
library(dplyr); library(survival); library(survey)
options(survey.lonely.psu = "adjust")

OUT <- normalizePath("results")

# ─── PART 1: GI Cancer ──────────────────────────────────────────────────────
cat("\n═══ GI Cancer: Unweighted vs Survey-Weighted ═══\n")
gi <- readRDS("results/gi_analysis/nhanes_clean.rds") %>%
  filter(gi_tumor == 1, surv_years > 0) %>%
  mutate(PNI_s = as.numeric(scale(PNI)),
         CONUT_s = as.numeric(scale(-CONUT)),
         GNRI_s = as.numeric(scale(GNRI)),
         age_s = as.numeric(scale(age)),
         sex_b = ifelse(sex == "Female", 1L, 0L),
         race_eth_b = ifelse(race_eth == "Non-Hispanic White", 1L, 0L),
         wt_pooled = WTMEC2YR / 2)

covars <- "age_s + sex_b + race_eth_b"

gi_combined <- data.frame()
for (idx in c("PNI_s", "CONUT_s", "GNRI_s")) {
  label <- gsub("_s", "", idx)

  # Unweighted
  f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", covars))
  ufit <- coxph(f, data = gi)
  us <- summary(ufit)
  uci <- us$conf.int[rownames(us$conf.int) == idx, ]

  # Weighted (restrict to Continuous NHANES with weights)
  gi_w <- gi %>% filter(!is.na(wt_pooled), !is.na(SDMVPSU), !is.na(SDMVSTRA))
  des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                   data = gi_w, nest = TRUE)
  wf <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", covars))
  wfit <- svycoxph(wf, design = des)
  ws <- summary(wfit)
  wci <- ws$conf.int[rownames(ws$conf.int) == idx, ]

  gi_combined <- rbind(gi_combined, data.frame(
    Index = label,
    Method = "Unweighted",
    N = nrow(gi), Events = sum(gi$death),
    HR = sprintf("%.2f", uci[[1]]),
    CI = sprintf("%.2f-%.2f", uci[[3]], uci[[4]]),
    P = sprintf("%.3f", us$coefficients[rownames(us$coefficients) == idx, 5]),
    stringsAsFactors = FALSE
  ))
  # svycoxph col 5 = z, col 6 = p
  w_p <- ws$coefficients[rownames(ws$coefficients) == idx, 6]
  gi_combined <- rbind(gi_combined, data.frame(
    Index = label,
    Method = "Survey-Weighted",
    N = nrow(gi_w), Events = sum(gi_w$death),
    HR = sprintf("%.2f", wci[[1]]),
    CI = sprintf("%.2f-%.2f", wci[[3]], wci[[4]]),
    P = sprintf("%.3f", w_p),
    stringsAsFactors = FALSE
  ))
}
write.csv(gi_combined, file.path(OUT, "gi_analysis", "cox_combined.csv"), row.names = FALSE)
print(gi_combined, row.names = FALSE)

# ─── PART 2: Pancancer Pooled ───────────────────────────────────────────────
cat("\n═══ Pancancer Pooled: Unweighted vs Survey-Weighted ═══\n")
pc <- readRDS("data/pancancer/nhanes_pancancer.rds") %>%
  filter(any_cancer == 1, surv_years > 0,
         death %in% c(0, 1), !is.na(pan_group),
         pan_group != "NonCancer") %>%
  group_by(pan_group) %>%
  mutate(PNI_z = as.numeric(scale(PNI)),
         CONUT_z = as.numeric(scale(-CONUT)),
         GNRI_z = as.numeric(scale(GNRI))) %>%
  ungroup() %>%
  mutate(wt_pooled = WTMEC2YR / 2)

pc_covars <- "age + sex + race_eth + edu_binary + income_pir"
MIN_EV <- 50
groups_keep <- pc %>% group_by(pan_group) %>%
  summarise(ev = sum(death), .groups = "drop") %>% filter(ev >= MIN_EV) %>% pull(pan_group)
pc <- pc %>% filter(pan_group %in% groups_keep)
cat(sprintf("Analytic sample: N=%d, events=%d\n", nrow(pc), sum(pc$death)))

pc_combined <- data.frame()
for (idx in c("PNI_z", "CONUT_z", "GNRI_z")) {
  label <- gsub("_z", "", idx)

  # Unweighted (stratified by cancer type — same as manuscript)
  f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+",
                         pc_covars, "+ strata(pan_group)"))
  ufit <- coxph(f, data = pc)
  us <- summary(ufit)
  uci <- us$conf.int[rownames(us$conf.int) == idx, ]

  # Weighted (restrict to observations with weights)
  pc_w <- pc %>% filter(!is.na(wt_pooled), !is.na(SDMVPSU), !is.na(SDMVSTRA))
  des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                   data = pc_w, nest = TRUE)
  wf <- as.formula(paste("Surv(surv_years, death) ~", idx, "+",
                          pc_covars, "+ pan_group"))
  wfit <- svycoxph(wf, design = des)
  ws <- summary(wfit)
  wci <- ws$conf.int[rownames(ws$conf.int) == idx, ]

  pc_combined <- rbind(pc_combined, data.frame(
    Index = label,
    Method = "Unweighted",
    N = nrow(pc), Events = sum(pc$death),
    HR = sprintf("%.2f", uci[[1]]),
    CI = sprintf("%.2f-%.2f", uci[[3]], uci[[4]]),
    P = sprintf("%.3f", us$coefficients[rownames(us$coefficients) == idx, 5]),
    stringsAsFactors = FALSE
  ))
  # svycoxph col 5 = z, col 6 = p
  w_p_pc <- ws$coefficients[rownames(ws$coefficients) == idx, 6]
  pc_combined <- rbind(pc_combined, data.frame(
    Index = label,
    Method = "Survey-Weighted",
    N = nrow(pc_w), Events = sum(pc_w$death),
    HR = sprintf("%.2f", wci[[1]]),
    CI = sprintf("%.2f-%.2f", wci[[3]], wci[[4]]),
    P = sprintf("%.3f", w_p_pc),
    stringsAsFactors = FALSE
  ))
}
dir.create(file.path(OUT, "pancancer"), showWarnings = FALSE, recursive = TRUE)
write.csv(pc_combined, file.path(OUT, "pancancer", "cox_combined.csv"), row.names = FALSE)
print(pc_combined, row.names = FALSE)

cat("\nDone. Files written:\n")
cat("  results/gi_analysis/cox_combined.csv\n")
cat("  results/pancancer/cox_combined.csv\n")
