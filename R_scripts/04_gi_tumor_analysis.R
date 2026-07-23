# =============================================================================
# 04_gi_tumor_analysis.R — NHANES GI Tumor × Nutrition × Survival
# =============================================================================
# 1. Download NHANES 2005-2018 (7 cycles) + cancer site + mortality
# 2. Classify GI tumors (esophagus, stomach, colon, rectum, liver, pancreas, biliary)
# 3. Derive PNI/CONUT/GNRI/E-DII (reusing 01_load_and_derive logic)
# 4. Cross-sectional: nutrition vs GI tumor status
# 5. Survival: Cox PH + competing risks (cancer vs non-cancer death)
# =============================================================================

source(here::here("R_scripts", "00_config.R"))
source(here::here("R_scripts", "00b_gi_cancer_config.R"))

# ── GI-specific output paths ─────────────────────────────────────────────────
GI_DATA_DIR   <- file.path(PROJ_ROOT, "data", "gi_analysis")
GI_RESULTS_DIR <- file.path(PROJ_ROOT, "results", "gi_analysis")
GI_FIG_DIR    <- file.path(PROJ_ROOT, "figures", "gi_analysis")
dir.create(GI_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GI_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GI_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION A: Download cancer site tables
# ═══════════════════════════════════════════════════════════════════════════════

#' Download cancer site data for one NHANES cycle
download_cancer_site <- function(cycle, use_cache = TRUE) {
  info <- CANCER_SITE_CYCLE_MAP[[cycle]]
  if (is.null(info)) stop("Unknown cycle: ", cycle)

  cache_file <- file.path(GI_DATA_DIR, sprintf("mcq_%s.rds", cycle))
  if (use_cache && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  cat(sprintf("  Downloading cancer site table for %s...\n", cycle))
  raw <- nhanesA::nhanes(info$table_name)

  # Rename cancer site var to unified name, keep SEQN
  if (!is.null(raw) && info$var_name %in% names(raw)) {
    result <- raw[, c("SEQN", info$var_name, "MCQ220")]
    names(result)[2] <- "cancer_site_code"
    result$cycle <- cycle
  } else {
    result <- data.frame(SEQN = raw$SEQN, cancer_site_code = NA_integer_,
                         MCQ220 = NA_real_, cycle = cycle)
  }

  saveRDS(result, cache_file)
  return(result)
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION B: Classify GI tumors
# ═══════════════════════════════════════════════════════════════════════════════

#' Classify cancer site codes into GI tumor categories
classify_gi_tumor <- function(mcq_df) {
  info <- CANCER_SITE_CYCLE_MAP[[unique(mcq_df$cycle)]]
  if (is.null(info)) {
    mcq_df$gi_tumor <- NA
    mcq_df$gi_site <- NA_character_
    return(mcq_df)
  }

  code <- as.character(mcq_df$cancer_site_code)
  mcq_df$gi_tumor <- as.integer(code %in% as.character(info$gi_codes) &
                                  !is.na(code))
  mcq_df$gi_site <- info$code_site[code]
  mcq_df$gi_site[is.na(mcq_df$gi_site)] <- NA_character_
  mcq_df$any_cancer <- as.integer(mcq_df$MCQ220 == 1 & !is.na(mcq_df$MCQ220))
  mcq_df
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION C: Derive nutritional indices (adapted from 01_load_and_derive.R)
# ═══════════════════════════════════════════════════════════════════════════════

#' Derive PNI, CONUT, GNRI, and basic biomarkers from NHANES tables
derive_nutrition <- function(demo, biopro, cbc, tchol, bmx) {
  df <- demo

  # Albumin + CRP
  if (!is.null(biopro)) {
    df <- df %>%
      left_join(biopro %>%
                  dplyr::select(SEQN, albumin_gdl = LBDSALSI,
                                crp_mgdl = LBDSCRPSI) %>%
                  mutate(across(c(albumin_gdl, crp_mgdl), as.numeric)),
                by = "SEQN")
  }

  # Lymphocyte count (from CBC)
  if (!is.null(cbc)) {
    df <- df %>%
      left_join(cbc %>%
                  dplyr::select(SEQN, wbc_1000 = LBXWBCSI, lympct = LBXLYPCT) %>%
                  mutate(across(c(wbc_1000, lympct), as.numeric),
                         lymph_abs = wbc_1000 * lympct / 100 * 1000),
                by = "SEQN")
  }

  # Total cholesterol
  if (!is.null(tchol)) {
    df <- df %>%
      left_join(tchol %>%
                  dplyr::select(SEQN, tchol_mgdl = LBXTC) %>%
                  mutate(tchol_mgdl = as.numeric(tchol_mgdl)),
                by = "SEQN")
  }

  # BMI / waist
  if (!is.null(bmx)) {
    df <- df %>%
      left_join(bmx %>%
                  dplyr::select(SEQN, bmi = BMXBMI, waist_cm = BMXWAIST) %>%
                  mutate(across(c(bmi, waist_cm), as.numeric)),
                by = "SEQN")
  }

  # ── Nutrition indices ────────────────────────────────────────────────────
  df <- df %>%
    mutate(
      PNI   = 10 * albumin_gdl + 0.005 * lymph_abs,

      CONUT_alb = case_when(
        albumin_gdl >= 3.5  ~ 0, albumin_gdl >= 3.0 ~ 2,
        albumin_gdl >= 2.5  ~ 4, !is.na(albumin_gdl) ~ 6, TRUE ~ NA_real_),
      CONUT_lymph = case_when(
        lymph_abs >= 1600   ~ 0, lymph_abs >= 1200  ~ 1,
        lymph_abs >= 800    ~ 2, !is.na(lymph_abs)  ~ 3, TRUE ~ NA_real_),
      CONUT_chol = case_when(
        tchol_mgdl >= 180   ~ 0, tchol_mgdl >= 140  ~ 1,
        tchol_mgdl >= 100   ~ 2, !is.na(tchol_mgdl) ~ 3, TRUE ~ NA_real_),
      CONUT = CONUT_alb + CONUT_lymph + CONUT_chol,

      GNRI = 14.89 * albumin_gdl + 41.7 * (bmi / 22)
    )
  df
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION D: Load mortality data
# ═══════════════════════════════════════════════════════════════════════════════

#' Load NCHS linked mortality data from CDC FTP .dat files
#'
#' Reads the public-use fixed-width mortality files directly.
#' File layout (verified from hex dump):
#'   SEQN 1-6, (7-14 blank), eligstat 15, mortstat 16,
#'   ucod_leading 17-19, diabetes 20, hyperten 21,
#'   (22-42 blank), permth_int 43-45, permth_exm 46-48
load_mortality <- function(use_cache = TRUE) {
  cache_file <- file.path(GI_DATA_DIR, "mortality_2019.rds")
  if (use_cache && file.exists(cache_file)) {
    cat("Loading mortality data from cache...\n")
    return(readRDS(cache_file))
  }

  base_url <- "https://ftp.cdc.gov/pub/health_statistics/nchs/datalinkage/linked_mortality"
  cycles_mort <- c("1999_2000", "2001_2002", "2003_2004", "2005_2006",
                   "2007_2008", "2009_2010", "2011_2012", "2013_2014",
                   "2015_2016", "2017_2018")

  all_mort <- list()
  for (cyc in cycles_mort) {
    fname <- sprintf("NHANES_%s_MORT_2019_PUBLIC.dat", cyc)
    url <- file.path(base_url, fname)
    tmp <- tempfile(fileext = ".dat")

    cat(sprintf("  Downloading %s...\n", fname))
    tryCatch({
      download.file(url, tmp, mode = "wb", quiet = TRUE)
      mort <- read.fwf(tmp,
        widths = c(6, 8, 1, 1, 3, 1, 1, 21, 3, 3, 13),
        col.names = c("SEQN", "blank1", "eligstat", "mortstat",
                      "ucod_leading", "diabetes", "hyperten",
                      "blank2", "permth_int", "permth_exm", "blank3"),
        colClasses = c("integer", "NULL", "integer", "integer",
                       "character", "integer", "integer",
                       "NULL", "integer", "integer", "NULL"),
        na.strings = c("", ".", " "))
      # ucod_leading: blank = alive/ineligible
      mort$ucod_leading[mort$ucod_leading == ""] <- NA_character_
      mort$ucod_leading <- as.integer(mort$ucod_leading)
      mort$cycle_mort <- cyc
      all_mort[[cyc]] <- mort
      cat(sprintf("    %d records\n", nrow(mort)))
    }, error = function(e) {
      warning(sprintf("Failed to download %s: %s", fname, conditionMessage(e)))
    })
    unlink(tmp)
  }

  if (length(all_mort) == 0) stop("No mortality files could be downloaded")
  mort_all <- bind_rows(all_mort)
  cat(sprintf("Total mortality records: %d\n", nrow(mort_all)))
  saveRDS(mort_all, cache_file)
  mort_all
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION E: Master loader — all cycles
# ═══════════════════════════════════════════════════════════════════════════════

cycle_table_name <- function(base, cycle) {
  suffix <- switch(cycle,
    "2005-2006" = "D", "2007-2008" = "E", "2009-2010" = "F",
    "2011-2012" = "G", "2013-2014" = "H", "2015-2016" = "I",
    "2017-2018" = "J", stop("Unknown cycle: ", cycle))
  paste0(base, "_", suffix)
}

load_gi_data <- function(use_cache = TRUE, min_age = 60) {
  all_data <- list()

  for (cyc in GI_CYCLES) {
    cat(sprintf("\n── %s ──\n", cyc))

    # Download key tables
    tbl_demo <- cycle_table_name("DEMO", cyc)
    tbl_biopro <- cycle_table_name("BIOPRO", cyc)
    tbl_cbc <- cycle_table_name("CBC", cyc)
    tbl_tchol <- cycle_table_name("TCHOL", cyc)
    tbl_bmx <- cycle_table_name("BMX", cyc)

    demo  <- nhanesA::nhanes(tbl_demo)
    biopro <- tryCatch(nhanesA::nhanes(tbl_biopro), error = \(e) NULL)
    cbc    <- tryCatch(nhanesA::nhanes(tbl_cbc), error = \(e) NULL)
    tchol  <- tryCatch(nhanesA::nhanes(tbl_tchol), error = \(e) NULL)
    bmx    <- tryCatch(nhanesA::nhanes(tbl_bmx), error = \(e) NULL)

    # Restrict to age >= min_age early for efficiency
    demo <- demo %>%
      mutate(age = RIDAGEYR) %>%
      filter(age >= min_age)

    cancer_site <- download_cancer_site(cyc, use_cache)

    # Merge cancer site
    demo <- demo %>% left_join(cancer_site, by = "SEQN", suffix = c("", ".mcq"))

    # Derive nutrition
    df <- derive_nutrition(demo, biopro, cbc, tchol, bmx)

    # Classify GI tumor
    df <- classify_gi_tumor(df)

    # Demographics recode
    df <- df %>%
      mutate(
        sex = factor(RIAGENDR, 1:2, c("Male", "Female")),
        race_eth = case_when(
          RIDRETH1 == 1 ~ "Mexican American", RIDRETH1 == 2 ~ "Other Hispanic",
          RIDRETH1 == 3 ~ "Non-Hispanic White",
          RIDRETH1 == 4 ~ "Non-Hispanic Black",
          RIDRETH1 == 5 ~ "Other Race", TRUE ~ NA_character_),
        edu_binary = ifelse(DMDEDUC2 >= 4, 1, 0),
        income_pir = INDFMPIR,
        wt_mec_2yr = WTMEC2YR
      )

    df$cycle <- cyc
    all_data[[cyc]] <- df
    cat(sprintf("  N=%d (age>=%d), GI tumors: %d\n",
                nrow(df), min_age, sum(df$gi_tumor, na.rm = TRUE)))
  }

  df_all <- bind_rows(all_data)
  cat(sprintf("\nTotal across %d cycles: N=%d\n", length(GI_CYCLES), nrow(df_all)))
  cat(sprintf("Total GI tumor cases: %d\n", sum(df_all$gi_tumor, na.rm = TRUE)))
  df_all
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION F: Execute loading
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── GI Analysis: Data Loading ──\n")
df <- load_gi_data(use_cache = TRUE, min_age = 60)

# ── NCHS mortality linkage ─────────────────────────────────────────────────
# NHANES public-use Linked Mortality Files (LMF) contain only the 10-category
# UCOD_LEADING recode (001-010). Site-specific ICD-10 underlying cause of
# death codes (e.g., C15 for esophageal cancer) and the 113-cause UCOD_113
# recode are ONLY available in restricted-use RDC files.
#
# Fix: For participants with a known GI tumor at baseline (MCQ230/MCQ240),
# a cancer death (UCOD_LEADING=002) is classified as "probable GI cancer
# death" — a reasonable proxy given the high case-fatality of GI cancers
# (esophagus ~80%, pancreas ~88%, liver ~80%, stomach ~68%, CRC ~35%).
# This approach is standard in public-use NHANES survival analyses.
# ===========================================================================

#' Classify cause of death using available public-use data
#'
#' @param ucod UCOD_LEADING value (002 = malignant neoplasm)
#' @param gi_tumor 1 if participant had GI tumor at baseline, 0 otherwise
#' @return list with cod_status (0=censor, 1=GI cancer death, 2=other death)
#'         and gi_cancer_death flag
classify_cod <- function(ucod, gi_tumor) {
  death <- !is.na(ucod) & ucod == 2
  other_death <- !is.na(ucod) & ucod != 2

  # For GI tumor patients who die of cancer, classify as GI cancer death
  gi_cancer <- death & gi_tumor == 1
  # Other cancer deaths (non-GI patients or non-cancer deaths)
  other_cancer_death <- other_death | (death & (gi_tumor == 0 | is.na(gi_tumor)))

  cod_status <- case_when(
    is.na(ucod)          ~ 0L,   # alive / censored
    gi_cancer            ~ 1L,   # probable GI cancer death
    TRUE                 ~ 2L    # non-GI-cancer death
  )

  list(cod_status = cod_status,
       gi_cancer_death = as.integer(gi_cancer),
       non_gi_cod_death = as.integer(other_cancer_death))
}

# ── Merge mortality ─────────────────────────────────────────────────────────
mortality <- load_mortality(use_cache = TRUE)
df <- df %>% left_join(mortality, by = "SEQN")

# Filter to eligible mortality participants
df <- df %>% filter(eligstat == 1, !is.na(mortstat))

# ── Classify cause of death ─────────────────────────────────────────────────
# Public-use LMF only has UCOD_LEADING (10 categories, 002=malignant neoplasm)
# Restricted-use RDC data needed for site-specific ICD-10 codes.
# Using baseline GI tumor + cancer death as proxy for GI cancer-specific death.
cod_class <- classify_cod(df$ucod_leading, df$gi_tumor)

df <- df %>%
  mutate(
    surv_years = permth_int / 12,
    death = mortstat,
    cod_status = cod_class$cod_status,         # 0=censor, 1=GI cancer, 2=other
    gi_cancer_death = cod_class$gi_cancer_death
  )

cat(sprintf("\nAfter mortality merge: N=%d\n", nrow(df)))
cat(sprintf("  Deceased: %d (%.1f%%)\n", sum(df$death, na.rm = TRUE),
            mean(df$death, na.rm = TRUE) * 100))
cat(sprintf("  GI cancer deaths (UCOD=002 + GI tumor at baseline): %d\n",
            sum(df$gi_cancer_death, na.rm = TRUE)))
cat(sprintf("  Other deaths: %d\n",
            sum(df$cod_status == 2, na.rm = TRUE)))
cat(sprintf("  Note: Public-use COD only has 10-group recode.\n"))
cat(sprintf("  Site-specific ICD-10 codes require restricted-use RDC access.\n"))

# Complete exclusion
df <- df %>%
  filter(!is.na(albumin_gdl), !is.na(lymph_abs),
         !is.na(tchol_mgdl), !is.na(bmi),
         !is.na(PNI))

cat(sprintf("Analytic sample (nutrition complete): N=%d\n", nrow(df)))
cat(sprintf("  GI tumor patients: %d\n", sum(df$gi_tumor == 1, na.rm = TRUE)))
cat(sprintf("  Non-cancer: %d\n",
            sum(df$any_cancer == 0, na.rm = TRUE)))

# Save merged dataset
saveRDS(df, file.path(GI_DATA_DIR, "nhanes_gi_merged.rds"))
cat("Merged dataset saved.\n")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION G: Survey design + Table 1
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Survey Design ──\n")

# Pooled weight: WTMEC2YR / number of cycles
n_cycles <- length(unique(df$cycle))
df <- df %>% mutate(wt_pooled = wt_mec_2yr / n_cycles)

design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  nest = TRUE, data = df
)

# ── Table 1 ──
cat("\n── Table 1: Descriptive by GI tumor status ──\n")

df <- df %>% mutate(
  gi_status = case_when(
    gi_tumor == 1 ~ "GI Tumor",
    any_cancer == 0 ~ "Non-Cancer",
    TRUE ~ "Other Cancer"
  )
)
df$gi_status <- factor(df$gi_status, levels = c("Non-Cancer", "GI Tumor", "Other Cancer"))

t1_vars <- c("age", "sex", "race_eth", "edu_binary", "income_pir",
             "PNI", "CONUT", "GNRI", "bmi", "albumin_gdl", "surv_years", "death")

t1_list <- list()
for (v in t1_vars) {
  if (is.numeric(df[[v]])) {
    means <- svyby(as.formula(paste0("~", v)), ~gi_status, design,
                   svymean, na.rm = TRUE)
    t1_list[[v]] <- data.frame(
      variable = v, level = "",
      NonCancer = sprintf("%.2f", means["Non-Cancer", 2]),
      GI_Tumor  = sprintf("%.2f", means["GI Tumor", 2]),
      OtherCancer = sprintf("%.2f", means["Other Cancer", 2])
    )
  } else {
    props <- svyby(as.formula(paste0("~", v)), ~gi_status, design,
                   svymean, na.rm = TRUE)
    # just first level
    lev <- levels(df[[v]])[1]
    t1_list[[v]] <- data.frame(
      variable = v, level = lev,
      NonCancer = sprintf("%.1f%%", props["Non-Cancer", lev] * 100),
      GI_Tumor  = sprintf("%.1f%%", props["GI Tumor", lev] * 100),
      OtherCancer = sprintf("%.1f%%", props["Other Cancer", lev] * 100)
    )
  }
}
table1 <- do.call(rbind, t1_list)
write.csv(table1, file.path(GI_RESULTS_DIR, "table1_gi_status.csv"), row.names = FALSE)
cat("Table 1 saved.\n")

# Site breakdown
site_counts <- df %>%
  filter(gi_tumor == 1) %>%
  count(gi_site, name = "N")
write.csv(site_counts, file.path(GI_RESULTS_DIR, "gi_site_counts.csv"), row.names = FALSE)
cat("Site breakdown:\n"); print(site_counts)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION H: Cross-sectional — nutrition indices × GI tumor
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Cross-sectional: Nutrition × GI Tumor ──\n")

# Restrict to GI tumor vs non-cancer for cleaner comparison
df_cs <- df %>% filter(gi_status %in% c("Non-Cancer", "GI Tumor"))
design_cs <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  nest = TRUE, data = df_cs
)

nutrition_vars <- c("PNI", "CONUT", "GNRI")
covars_cs <- c("age", "sex", "race_eth", "edu_binary", "income_pir")

cs_results <- data.frame()
for (nv in nutrition_vars) {
  # Model 1: crude
  f1 <- as.formula(paste(nv, "~ gi_tumor"))
  m1 <- svyglm(f1, design = design_cs)
  b1 <- tidy(m1) %>% filter(term == "gi_tumor")

  # Model 2: adjusted
  f2 <- as.formula(paste(nv, "~ gi_tumor +", paste(covars_cs, collapse = " + ")))
  m2 <- svyglm(f2, design = design_cs)
  b2 <- tidy(m2) %>% filter(term == "gi_tumor")

  cs_results <- rbind(cs_results,
    data.frame(nutrition = nv, model = "Crude",
               beta = b1$estimate, se = b1$std.error,
               lower = b1$conf.low, upper = b1$conf.high, p = b1$p.value),
    data.frame(nutrition = nv, model = "Adjusted",
               beta = b2$estimate, se = b2$std.error,
               lower = b2$conf.low, upper = b2$conf.high, p = b2$p.value))
}
cs_results$ci95 <- sprintf("%.3f (%.3f, %.3f)", cs_results$beta,
                            cs_results$lower, cs_results$upper)
write.csv(cs_results, file.path(GI_RESULTS_DIR, "cross_sectional_results.csv"),
          row.names = FALSE)
cat("Cross-sectional results saved.\n")
print(cs_results[, c("nutrition", "model", "ci95", "p")])

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION I: All-cause survival — Cox PH (GI tumor patients only)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── All-Cause Survival (GI tumor patients) ──\n")

df_gi <- df %>% filter(gi_tumor == 1, surv_years > 0) %>%
  mutate(
    PNI_scaled = as.numeric(scale(PNI)),
    CONUT_scaled = as.numeric(scale(-CONUT)),     # reverse: higher = better
    GNRI_scaled = as.numeric(scale(GNRI))
  )

n_gi <- nrow(df_gi)
n_events <- sum(df_gi$death, na.rm = TRUE)
cat(sprintf("GI tumor analytic sample: N=%d, events=%d\n", n_gi, n_events))

if (n_events >= 20) {
  surv_covars <- c("age", "sex", "race_eth", "edu_binary", "income_pir")

  cox_results <- data.frame()
  for (exp_var in c("PNI_scaled", "CONUT_scaled", "GNRI_scaled")) {
    for (adj in c("Crude", "Adjusted")) {
      covs <- if (adj == "Crude") "1" else paste(surv_covars, collapse = " + ")
      form <- as.formula(paste("Surv(surv_years, death) ~", exp_var, "+", covs))
      fit <- coxph(form, data = df_gi)

      hr <- tidy(fit, conf.int = TRUE) %>% filter(term == exp_var)
      cox_results <- rbind(cox_results, data.frame(
        nutrition = gsub("_scaled", "", exp_var),
        adjustment = adj,
        HR = exp(hr$estimate), se = hr$std.error,
        lower = exp(hr$conf.low), upper = exp(hr$conf.high),
        p = hr$p.value, N = n_gi, events = n_events
      ))
    }
  }

  cox_results$hr_ci <- sprintf("%.3f (%.3f, %.3f)", cox_results$HR,
                                cox_results$lower, cox_results$upper)
  write.csv(cox_results, file.path(GI_RESULTS_DIR, "cox_allcause.csv"), row.names = FALSE)
  cat("All-cause Cox results:\n")
  print(cox_results[, c("nutrition", "adjustment", "hr_ci", "p")])

  # Check PH assumption
  cat("\nProportional hazards test (Schoenfeld residuals):\n")
  for (exp_var in c("PNI_scaled", "CONUT_scaled", "GNRI_scaled")) {
    fit_adj <- coxph(as.formula(paste(
      "Surv(surv_years, death) ~", exp_var, "+",
      paste(surv_covars, collapse = " + "))), data = df_gi)
    print(cox.zph(fit_adj))
  }
} else {
  cat("Too few events for Cox regression.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION J: Competing risks — GI cancer death vs non-cancer death
# ═══════════════════════════════════════════════════════════════════════════════
# Note: Public-use NHANES LMF only provides UCOD_LEADING (10-group recode).
# GI cancer death = baseline GI tumor diagnosis + UCOD_LEADING=002 (cancer).
# This proxy has high specificity: most GI tumor patients who die of cancer
# died from their GI malignancy, especially for high-fatality sites
# (pancreas, esophagus, liver). Site-specific ICD-10 codes would require
# restricted-use RDC data.
# ===========================================================================

cat("\n── Competing Risks Analysis ──\n")

# Competing risks events among GI tumor patients: GI cancer death vs other death
if (n_events >= 30 && sum(df_gi$cod_status == 1, na.rm = TRUE) >= 10) {
  # Cause-specific Cox: GI cancer death
  cs_cox_results <- data.frame()
  for (cause_event in c("gi_cancer_death = 1", "cod_status = 2")) {
    cause_label <- c("gi_cancer_death = 1" = "GI cancer death (proxy)",
                     "cod_status = 2" = "Non-cancer death")[cause_event]
    cause_var <- ifelse(grepl("gi_cancer", cause_event), "gi_cancer_death", "other_death_flag")

    for (exp_var in c("PNI_scaled", "CONUT_scaled", "GNRI_scaled")) {
      if (cause_event == "gi_cancer_death = 1") {
        f_surv <- as.formula(paste(
          "Surv(surv_years, gi_cancer_death) ~", exp_var, "+",
          "age + sex + race_eth + edu_binary + income_pir"))
      } else {
        df_gi$other_death_flag <- ifelse(df_gi$cod_status == 2, 1L, 0L)
        f_surv <- as.formula(paste(
          "Surv(surv_years, other_death_flag) ~", exp_var, "+",
          "age + sex + race_eth + edu_binary + income_pir"))
      }

      fit <- coxph(f_surv, data = df_gi)
      hr <- tidy(fit, conf.int = TRUE) %>% filter(term == exp_var)
      cs_cox_results <- rbind(cs_cox_results, data.frame(
        nutrition = gsub("_scaled", "", exp_var),
        cause = cause_label,
        HR = exp(hr$estimate), lower = exp(hr$conf.low),
        upper = exp(hr$conf.high), p = hr$p.value
      ))
    }
  }
  cs_cox_results$hr_ci <- sprintf("%.3f (%.3f, %.3f)",
                                   cs_cox_results$HR,
                                   cs_cox_results$lower,
                                   cs_cox_results$upper)
  write.csv(cs_cox_results, file.path(GI_RESULTS_DIR, "cox_causespecific.csv"),
            row.names = FALSE)
  cat("Cause-specific Cox results (GI cancer death = baseline GI tumor + UCOD=002):\n")
  print(cs_cox_results[, c("nutrition", "cause", "hr_ci", "p")])

  # Gray's test (crude by PNI tertile)
  library(cmprsk)
  df_gi <- df_gi %>% mutate(pni_t = factor(ntile(PNI, 3), 1:3))
  cif <- cuminc(
    ftime   = df_gi$surv_years,
    fstatus = df_gi$cod_status,
    group   = df_gi$pni_t,
    cencode = 0
  )
  saveRDS(cif, file.path(GI_RESULTS_DIR, "cuminc_cif.rds"))

  # Gray's test
  cat("\nGray's test p-values:\n")
  print(cif$Tests)
} else {
  cat("Too few events for competing risks.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION K: Sensitivity analyses
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── Sensitivity Analyses ──\n")

sens_results <- data.frame()
sens_subsets <- list(
  "GI only (exclude other cancer)" = list(
    subset = expression(gi_tumor == 1),
    label = "GI tumor only"),
  "GI + non-cancer" = list(
    subset = expression(gi_tumor == 1 | any_cancer == 0),
    label = "GI vs non-cancer"),
  "Exclude 6mo deaths" = list(
    subset = expression(gi_tumor == 1 & surv_years > 0.5),
    label = "GI only, exclude deaths <6mo"),
  "Age >= 65" = list(
    subset = expression(gi_tumor == 1 & age >= 65),
    label = "GI only, age >= 65")
)

for (sn in names(sens_subsets)) {
  s <- sens_subsets[[sn]]
  sub_df <- df %>% filter(eval(s$subset), surv_years > 0)
  n_sub <- nrow(sub_df)
  e_sub <- sum(sub_df$death, na.rm = TRUE)

  if (e_sub < 15) { cat(sprintf("  Skipping [%s]: %d events\n", s$label, e_sub)); next }

  sub_df <- sub_df %>% mutate(PNI_scaled = as.numeric(scale(PNI)))
  fit <- coxph(Surv(surv_years, death) ~ PNI_scaled + age + sex + edu_binary,
               data = sub_df)
  hr <- tidy(fit, conf.int = TRUE) %>% filter(term == "PNI_scaled")
  sens_results <- rbind(sens_results, data.frame(
    analysis = s$label, N = n_sub, events = e_sub,
    HR = exp(hr$estimate), lower = exp(hr$conf.low),
    upper = exp(hr$conf.high), p = hr$p.value
  ))
}
if (nrow(sens_results) > 0) {
  write.csv(sens_results, file.path(GI_RESULTS_DIR, "sensitivity.csv"), row.names = FALSE)
  cat("Sensitivity results:\n"); print(sens_results)
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION L: Save workspace
# ═══════════════════════════════════════════════════════════════════════════════

save.image(file.path(GI_RESULTS_DIR, "gi_analysis_workspace.RData"))
cat(sprintf("\n── 04_gi_tumor_analysis.R complete. Results in %s ──\n",
            GI_RESULTS_DIR))
