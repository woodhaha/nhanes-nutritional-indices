# run_hei_analysis.R — HEI-2015 calculation + GI tumor survival
# Rscript run_hei_analysis.R

library(haven); library(dplyr); library(survival); library(broom)
library(nhanesA)

FPED_DIR <- "D:/Researching/NHANES_aged_GI_tumor_nutrition/data/fped"
GI_DATA_DIR <- "D:/Researching/data/gi_analysis"
GI_RESULTS_DIR <- "D:/Researching/NHANES_aged_GI_tumor_nutrition/results/gi_analysis"
dir.create(GI_RESULTS_DIR, recursive=TRUE, showWarnings=FALSE)

# ── HEI-2015 scoring standards ───────────────────────────────────────────────
# From: Krebs-Smith et al. J Acad Nutr Diet. 2018;118(9):1591-1602
# All densities per 1000 kcal (except fatty acids ratio and SFA as % energy)

# Adequacy components (higher = better): [max score, min score, max density, min density]
HEI_ADEQUACY <- list(
  total_fruit    = c(max_score=5, min_score=0, max_density=0.8, min_density=0),
  whole_fruit    = c(max_score=5, min_score=0, max_density=0.4, min_density=0),
  total_veg      = c(max_score=5, min_score=0, max_density=1.1, min_density=0),
  greens_beans   = c(max_score=5, min_score=0, max_density=0.2, min_density=0),
  whole_grains   = c(max_score=10, min_score=0, max_density=1.5, min_density=0),
  dairy          = c(max_score=10, min_score=0, max_density=1.3, min_density=0),
  total_protein  = c(max_score=5, min_score=0, max_density=2.5, min_density=0),
  seafood_plant_p = c(max_score=5, min_score=0, max_density=0.8, min_density=0),
  fatty_acids    = c(max_score=10, min_score=0, max_ratio=2.5, min_ratio=1.2)  # MUFA+PUFA/SFA
)

# Moderation components (lower = better): [max score, min score, max density, min density]
HEI_MODERATION <- list(
  refined_grains  = c(max_score=10, min_score=0, max_density=1.8, min_density=4.5),
  sodium          = c(max_score=10, min_score=0, max_density=1.1, min_density=2.0),
  added_sugars    = c(max_score=10, min_score=0, max_density=6.5, min_density=26.0),
  saturated_fats  = c(max_score=10, min_score=0, max_pct=8, min_pct=16)
)

# ── Calculate HEI-2015 component scores ─────────────────────────────────────
calc_hei2015 <- function(fped, diet) {
  # Merge FPED + diet nutrient data
  df <- merge(fped, diet[, c("SEQN", "DR1TKCAL", "DR1TSFAT", "DR1TMFAT",
                              "DR1TPFAT", "DR1TSODI")],
              by="SEQN", all=TRUE)
  df <- df[!is.na(df$DR1TKCAL) & df$DR1TKCAL > 0 & !is.na(df$DR1T_F_TOTAL), ]
  if (nrow(df) == 0) return(NULL)

  # Energy in kcal
  kcal <- df$DR1TKCAL / 1000  # per 1000 kcal

  # ── Adequacy components (per 1000 kcal) ──
  # Total Fruit (cup eq/1000kcal)
  total_fruit_d <- df$DR1T_F_TOTAL / kcal
  df$HEI_total_fruit <- pmin(5, pmax(0, 5 * (total_fruit_d - 0) / (0.8 - 0)))

  # Whole Fruit (cup eq/1000kcal)
  whole_fruit_d <- df$DR1T_F_CITMLB / kcal
  df$HEI_whole_fruit <- pmin(5, pmax(0, 5 * (whole_fruit_d - 0) / (0.4 - 0)))

  # Total Vegetables (cup eq/1000kcal)
  total_veg_d <- df$DR1T_V_TOTAL / kcal
  df$HEI_total_veg <- pmin(5, pmax(0, 5 * (total_veg_d - 0) / (1.1 - 0)))

  # Greens and Beans (cup eq/1000kcal) = dark green veg + legumes as veg
  greens_beans_d <- (df$DR1T_V_DRKGR + df$DR1T_V_LEGUMES) / kcal
  df$HEI_greens_beans <- pmin(5, pmax(0, 5 * (greens_beans_d - 0) / (0.2 - 0)))

  # Whole Grains (oz eq/1000kcal)
  whole_grains_d <- df$DR1T_G_WHOLE / kcal
  df$HEI_whole_grains <- pmin(10, pmax(0, 10 * (whole_grains_d - 0) / (1.5 - 0)))

  # Dairy (cup eq/1000kcal)
  dairy_d <- df$DR1T_D_TOTAL / kcal
  df$HEI_dairy <- pmin(10, pmax(0, 10 * (dairy_d - 0) / (1.3 - 0)))

  # Total Protein Foods (oz eq/1000kcal)
  protein_d <- df$DR1T_PF_TOTAL / kcal
  df$HEI_total_protein <- pmin(5, pmax(0, 5 * (protein_d - 0) / (2.5 - 0)))

  # Seafood & Plant Protein (oz eq/1000kcal)
  spp_d <- (df$DR1T_PF_SEAFD_HI + df$DR1T_PF_SEAFD_LOW +
            df$DR1T_PF_NUTSDS + df$DR1T_PF_SOY + df$DR1T_PF_LEGUMES) / kcal
  df$HEI_seafood_plant <- pmin(5, pmax(0, 5 * (spp_d - 0) / (0.8 - 0)))

  # Fatty Acids Ratio (MUFA+PUFA)/SFA
  fa_ratio <- (df$DR1TMFAT + df$DR1TPFAT) / df$DR1TSFAT
  fa_ratio[df$DR1TSFAT == 0] <- 0
  df$HEI_fatty_acids <- pmin(10, pmax(0, 10 * (fa_ratio - 1.2) / (2.5 - 1.2)))

  # ── Moderation components ──
  # Refined Grains (oz eq/1000kcal, lower is better)
  refined_d <- df$DR1T_G_REFINED / kcal
  df$HEI_refined_grains <- pmin(10, pmax(0, 10 * (4.5 - refined_d) / (4.5 - 1.8)))

  # Sodium (g/1000kcal, lower is better)
  sodium_d <- (df$DR1TSODI / 1000) / kcal  # g/1000kcal
  df$HEI_sodium <- pmin(10, pmax(0, 10 * (2.0 - sodium_d) / (2.0 - 1.1)))

  # Added Sugars (tsp eq/1000kcal, lower is better)
  added_s_d <- df$DR1T_ADD_SUGARS / kcal
  df$HEI_added_sugars <- pmin(10, pmax(0, 10 * (26 - added_s_d) / (26 - 6.5)))

  # Saturated Fats (% energy, lower is better)
  sfa_pct <- (df$DR1TSFAT * 9) / df$DR1TKCAL * 100
  df$HEI_saturated_fats <- pmin(10, pmax(0, 10 * (16 - sfa_pct) / (16 - 8)))

  # ── Total HEI-2015 score (0-100) ──
  hei_cols <- grep("^HEI_", names(df), value=TRUE)
  df$HEI2015_total <- rowSums(df[, hei_cols], na.rm=TRUE)
  # Flag if >2 components missing
  n_missing <- rowSums(is.na(df[, hei_cols]))
  df$HEI2015_total[n_missing > 2] <- NA_real_

  cat(sprintf("HEI-2015: N=%d, mean=%.1f (SD=%.1f)\n",
              sum(!is.na(df$HEI2015_total)),
              mean(df$HEI2015_total, na.rm=TRUE),
              sd(df$HEI2015_total, na.rm=TRUE)))

  df[, c("SEQN", hei_cols, "HEI2015_total")]
}

# ── Cycle mapping ────────────────────────────────────────────────────────────
FPED_CYCLES <- list(
  "1112" = list(fped="fped_dr1tot_1112.sas7bdat", diet="DR1TOT_G"),
  "1314" = list(fped="fped_dr1tot_1314.sas7bdat", diet="DR1TOT_H"),
  "1516" = list(fped="fped_dr1tot_1516.sas7bdat", diet="DR1TOT_I")
)

cat("=== HEI-2015 Calculation ===\n")

all_hei <- list()
for (cyc in names(FPED_CYCLES)) {
  cat(sprintf("\n--- %s ---\n", cyc))

  # Read FPED
  fp <- read_sas(file.path(FPED_DIR, paste0("extracted_", cyc), FPED_CYCLES[[cyc]]$fped))
  cat(sprintf("  FPED: %d rows\n", nrow(fp)))

  # Read DR1TOT
  diet <- nhanesA::nhanes(FPED_CYCLES[[cyc]]$diet)
  cat(sprintf("  DIET: %d rows\n", nrow(diet)))

  # Calculate HEI
  hei <- calc_hei2015(fp, diet)
  if (!is.null(hei)) {
    hei$cycle <- cyc
    all_hei[[cyc]] <- hei
  }
}

hei_all <- bind_rows(all_hei)
dir.create(GI_DATA_DIR, recursive=TRUE, showWarnings=FALSE)
saveRDS(hei_all, file.path(GI_DATA_DIR, "hei2015_scores.rds"))
cat(sprintf("\nTotal HEI records: %d\n", nrow(hei_all)))

# ── Merge with GI tumor + survival data ──────────────────────────────────────
cat("\n=== GI Tumor × HEI-2015 Analysis ===\n")

# Load existing GI data
df_gi <- readRDS(file.path(GI_RESULTS_DIR, "nhanes_combined_gi.rds"))

# Map cycle to FPED years
df_gi$fped_cycle <- case_when(
  grepl("2011-2012", df_gi$cycle) ~ "1112",
  grepl("2013-2014", df_gi$cycle) ~ "1314",
  grepl("2015-2016", df_gi$cycle) ~ "1516",
  TRUE ~ NA_character_
)

# Merge HEI
df_gi <- df_gi %>% left_join(hei_all, by=c("SEQN"="SEQN"), suffix=c("", ".hei"))
cat(sprintf("Merged: %d rows, %d with HEI\n", nrow(df_gi), sum(!is.na(df_gi$HEI2015_total))))

# ── Cross-sectional: HEI comparison ──
cat("\n── HEI by GI status ──\n")
df_gi$HEI2015_total <- as.numeric(df_gi$HEI2015_total)
t_hei <- df_gi %>% filter(!is.na(HEI2015_total)) %>%
  group_by(gi_status) %>%
  summarise(N=n(), HEI_mean=mean(HEI2015_total, na.rm=TRUE),
            HEI_sd=sd(HEI2015_total, na.rm=TRUE))
print(as.data.frame(t_hei))

# ── Cox: HEI → mortality in GI tumor patients ──
cat("\n── HEI-2015 → All-cause Mortality (GI tumor) ──\n")
df_gi_cox <- df_gi %>% filter(gi_tumor==1, surv_years>0, !is.na(HEI2015_total)) %>%
  mutate(HEI_s = as.numeric(scale(HEI2015_total)),
         sex_b = ifelse(sex=="Female", 1L, 0L),
         race_eth_b = ifelse(race_eth=="Non-Hispanic White", 1L, 0L))

n_cox <- nrow(df_gi_cox); e_cox <- sum(df_gi_cox$death)
cat(sprintf("  N=%d, events=%d\n", n_cox, e_cox))

if (e_cox >= 10) {
  # HEI alone
  for (adj in c("Crude", "Adjusted")) {
    covs <- if(adj=="Crude") "1" else "age + sex_b + race_eth_b + edu_binary"
    f <- as.formula(paste("Surv(surv_years, death) ~ HEI_s +", covs))
    fit <- coxph(f, data=df_gi_cox)
    hr <- tidy(fit, conf.int=TRUE) %>% filter(term=="HEI_s")
    cat(sprintf("  HEI-2015 (%s): HR=%.3f (%.3f-%.3f), p=%.4f\n",
                adj, exp(hr$estimate), exp(hr$conf.low), exp(hr$conf.high), hr$p.value))
  }

  # HEI + PNI together (which is stronger?)
  cat("\n── HEI vs PNI: both in model ──\n")
  df_gi_cox$PNI_s <- as.numeric(scale(df_gi_cox$PNI))
  f_both <- Surv(surv_years, death) ~ HEI_s + PNI_s + age + sex_b + race_eth_b + edu_binary
  fit_both <- coxph(f_both, data=df_gi_cox)
  both <- tidy(fit_both, conf.int=TRUE) %>% filter(term %in% c("HEI_s", "PNI_s"))
  print(both)

  # KM by HEI tertile
  df_gi_cox$hei_t <- factor(ntile(df_gi_cox$HEI2015_total, 3), 1:3,
                             c("Low HEI (T1)", "Mid HEI (T2)", "High HEI (T3)"))
  km_hei <- survfit(Surv(surv_years, death) ~ hei_t, data=df_gi_cox)
  cat("\nMedian survival by HEI tertile:\n")
  print(km_hei)
}

saveRDS(df_gi, file.path(GI_RESULTS_DIR, "nhanes_combined_gi_with_hei.rds"))
cat("\nDone.\n")
