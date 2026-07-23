# ==============================================================================
# 01_load_and_derive.R — NHANES data loading & variable derivation
# Loads NHANES 2011-2014 data, merges tables, derives all analysis variables
# ==============================================================================

source(here::here("NHANES_cognition_nutrition", "R_scripts", "00_config.R"))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION A: Download or load NHANES data
# ═══════════════════════════════════════════════════════════════════════════════

#' Download a single NHANES table for one cycle
#' @param cycle NHANES cycle string, e.g., "2011-2012"
#' @param table_name NHANES table name, e.g., "DEMO_G"
#' @param vars Character vector of variable names. If NULL, downloads all.
download_nhanes_table <- function(cycle, table_name, vars = NULL) {
  cycle_letter <- switch(cycle,
    "2011-2012" = "G",
    "2013-2014" = "H",
    stop("Unknown cycle: ", cycle)
  )
  cat(sprintf("  Downloading %s from %s...\n", table_name, cycle))

  result <- tryCatch({
    if (is.null(vars)) {
      nhanesA::nhanes(table_name)
    } else {
      nhanesA::nhanes(table_name)
    }
  }, error = function(e) {
    warning("Failed to download ", table_name, ": ", e$message)
    return(NULL)
  })
  return(result)
}

#' Load or download NHANES data with caching
load_nhanes_cycle <- function(cycle, use_cache = TRUE) {
  cycle_letter <- switch(cycle,
    "2011-2012" = "G",
    "2013-2014" = "H"
  )
  cache_file <- file.path(DATA_DIR, sprintf("nhanes_%s_raw.RData", cycle))

  if (use_cache && file.exists(cache_file)) {
    cat(sprintf("[%s] Loading from cache...\n", cycle))
    load(cache_file, envir = .GlobalEnv)
    return(TRUE)
  }

  cat(sprintf("[%s] Downloading NHANES data...\n", cycle))

  # Build table names for this cycle
  tables <- list(
    DEMO   = list(name = paste0("DEMO_", cycle_letter), vars = DEMO_VARS),
    DR1TOT = list(name = paste0("DR1TOT_", cycle_letter), vars = DIET_VARS),
    DPQ    = list(name = paste0("DPQ_", cycle_letter), vars = PHQ9_VARS),
    CFQ    = list(name = paste0("CFQ_", cycle_letter), vars = CFQ_VARS),
    BIOPRO = list(name = paste0("BIOPRO_", cycle_letter), vars = BIOPRO_VARS),
    CBC    = list(name = paste0("CBC_", cycle_letter), vars = CBC_VARS),
    TCHOL  = list(name = paste0("TCHOL_", cycle_letter), vars = TCHOL_VARS),
    BMX    = list(name = paste0("BMX_", cycle_letter), vars = BMX_VARS),
    BPX    = list(name = paste0("BPX_", cycle_letter), vars = BPX_VARS),
    SMQ    = list(name = paste0("SMQ_", cycle_letter), vars = SMQ_VARS),
    ALQ    = list(name = paste0("ALQ_", cycle_letter), vars = ALQ_VARS),
    DIQ    = list(name = paste0("DIQ_", cycle_letter), vars = DIQ_VARS),
    MCQ    = list(name = paste0("MCQ_", cycle_letter), vars = MCQ_VARS),
    FSQ    = list(name = paste0("FSQ_", cycle_letter), vars = FSQ_VARS),
    PAQ    = list(name = paste0("PAQ_", cycle_letter), vars = PAQ_VARS)
  )

  # Download all tables
  raw_data <- list()
  for (tbl_key in names(tables)) {
    tbl <- tables[[tbl_key]]
    raw_data[[tbl_key]] <- download_nhanes_table(cycle, tbl$name, tbl$vars)
  }

  # Cache
  save(raw_data, file = cache_file)
  assign(paste0("raw_", cycle_letter), raw_data, envir = .GlobalEnv)
  cat(sprintf("[%s] Data saved to %s\n", cycle, cache_file))
  return(TRUE)
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION B: Variable derivation
# ═══════════════════════════════════════════════════════════════════════════════

#' Derive all analysis variables from raw NHANES data for a single cycle
derive_cycle_variables <- function(raw_list, cycle) {
  cycle_letter <- switch(cycle,
    "2011-2012" = "G",
    "2013-2014" = "H"
  )

  cat(sprintf("[%s] Deriving variables...\n", cycle))

  # Start with DEMO
  df <- raw_list$DEMO %>%
    dplyr::select(SEQN, SDDSRVYR, RIDAGEYR, RIAGENDR, RIDRETH1,
                  DMDEDUC2, DMDMARTL, INDFMPIR, SDMVPSU, SDMVSTRA,
                  WTMEC2YR, WTINT2YR)

  # ── Demographics ────────────────────────────────────────────────────────
  df <- df %>%
    mutate(
      age        = RIDAGEYR,
      age_group  = cut(age, breaks = c(60, 70, 80, Inf),
                       labels = c("60-69", "70-79", "80+"),
                       right = FALSE),
      sex        = factor(RIAGENDR, levels = c(1, 2),
                          labels = c("Male", "Female")),
      sex_binary = ifelse(RIAGENDR == 2, 1, 0),  # 0=Male, 1=Female

      race_eth = case_when(
        RIDRETH1 == 1 ~ "Mexican American",
        RIDRETH1 == 2 ~ "Other Hispanic",
        RIDRETH1 == 3 ~ "Non-Hispanic White",
        RIDRETH1 == 4 ~ "Non-Hispanic Black",
        RIDRETH1 == 5 ~ "Other Race",
        TRUE ~ NA_character_
      ),

      edu = case_when(
        DMDEDUC2 %in% c(1, 2) ~ "Less than HS",
        DMDEDUC2 == 3         ~ "High School",
        DMDEDUC2 == 4         ~ "Some College",
        DMDEDUC2 == 5         ~ "College Graduate",
        TRUE ~ NA_character_
      ),
      edu_binary = ifelse(DMDEDUC2 >= 4, 1, 0),

      income_pir = INDFMPIR,
      income_low = ifelse(INDFMPIR <= 1.30, 1, 0)  # ≤130% FPL
    )

  # ── PHQ-9 Depression (DPQ) ──────────────────────────────────────────────
  # DPQ010-DPQ090: "Not at all"=0 to "Nearly every day"=3
  if (!is.null(raw_list$DPQ)) {
    phq9_cols <- paste0("DPQ0", sprintf("%02d", 10:19))
    phq9_available <- intersect(phq9_cols, names(raw_list$DPQ))

    if (length(phq9_available) >= 8) {
      # Convert to numeric (NHANES data may have labels or codes)
      phq9_df <- raw_list$DPQ %>%
        dplyr::select(SEQN, all_of(phq9_available)) %>%
        mutate(across(-SEQN, ~ as.numeric(as.character(.))))

      # Recode: "Not at all"=0, "Several days"=1, "More than half"=2, "Nearly every day"=3
      # NHANES stores these as 0-3 directly
      phq9_items <- phq9_df %>%
        dplyr::select(-SEQN) %>%
        mutate(across(everything(), ~ ifelse(. > 3 | . < 0, NA, .)))

      phq9_df$phq9_total <- rowSums(phq9_items, na.rm = FALSE)
      # Only keep if >= 8 items are non-missing
      phq9_df$phq9_missing <- rowSums(is.na(phq9_items))
      phq9_df$phq9_total[phq9_df$phq9_missing > 1] <- NA_real_
      phq9_df$phq9_depression <- ifelse(phq9_df$phq9_total >= 10, 1, 0)

      df <- df %>%
        left_join(phq9_df %>% dplyr::select(SEQN, phq9_total, phq9_depression),
                  by = "SEQN")
    } else {
      warning("[", cycle, "] PHQ-9 data incomplete, skipping depression variables")
      df$phq9_total <- NA_real_
      df$phq9_depression <- NA_integer_
    }
  }

  # ── Cognitive Function (CFQ) ─────────────────────────────────────────────
  if (!is.null(raw_list$CFQ)) {
    cfq_df <- raw_list$CFQ %>%
      dplyr::select(SEQN, CFD_WL_IMM, CFD_WL_DEL, CFD_AFT, CFD_DSST) %>%
      mutate(
        cerad_imm  = as.numeric(CFD_WL_IMM),    # immediate recall, 0-30
        cerad_del  = as.numeric(CFD_WL_DEL),    # delayed recall, 0-10
        animal_flu = as.numeric(CFD_AFT),       # animal fluency
        dsst       = as.numeric(CFD_DSST)       # digit symbol substitution
      ) %>%
      dplyr::select(SEQN, cerad_imm, cerad_del, animal_flu, dsst)

    # Compute composite Z-score (internal standardization within cycle)
    cfq_df <- cfq_df %>%
      mutate(
        z_cerad_imm  = (cerad_imm - mean(cerad_imm, na.rm = TRUE)) /
                        sd(cerad_imm, na.rm = TRUE),
        z_cerad_del  = (cerad_del - mean(cerad_del, na.rm = TRUE)) /
                        sd(cerad_del, na.rm = TRUE),
        z_animal_flu = (animal_flu - mean(animal_flu, na.rm = TRUE)) /
                        sd(animal_flu, na.rm = TRUE),
        z_dsst       = (dsst - mean(dsst, na.rm = TRUE)) /
                        sd(dsst, na.rm = TRUE)
      )

    cfq_df$cog_z_composite <- rowMeans(
      cfq_df[, c("z_cerad_imm", "z_cerad_del", "z_animal_flu", "z_dsst")],
      na.rm = FALSE
    )

    # Probable MCI: composite z < -1 (similar to Du 2024)
    cfq_df$probable_mci <- ifelse(cfq_df$cog_z_composite < -1, 1, 0)

    df <- df %>% left_join(cfq_df, by = "SEQN")
  } else {
    warning("[", cycle, "] Cognitive data not available")
    df$cog_z_composite <- NA_real_
  }

  # ── Blood Biomarkers ────────────────────────────────────────────────────
  # Albumin (BIOPRO)
  if (!is.null(raw_list$BIOPRO)) {
    df <- df %>%
      left_join(raw_list$BIOPRO %>%
                dplyr::select(SEQN, albumin_gdl = LBDSALSI, crp_mgdl = LBDSCRPSI) %>%
                mutate(
                  albumin_gdl = as.numeric(albumin_gdl),
                  crp_mgdl    = as.numeric(crp_mgdl)
                ), by = "SEQN")
  }

  # CBC (WBC + Lymphocyte %)
  if (!is.null(raw_list$CBC)) {
    df <- df %>%
      left_join(raw_list$CBC %>%
                dplyr::select(SEQN, wbc_1000 = LBXWBCSI, lympct = LBXLYPCT) %>%
                mutate(
                  wbc_1000     = as.numeric(wbc_1000),
                  lympct       = as.numeric(lympct),
                  lymph_abs    = wbc_1000 * lympct / 100 * 1000  # cells/µL
                ), by = "SEQN")
  }

  # Total Cholesterol (TCHOL)
  if (!is.null(raw_list$TCHOL)) {
    df <- df %>%
      left_join(raw_list$TCHOL %>%
                dplyr::select(SEQN, tchol_mgdl = LBXTC) %>%
                mutate(tchol_mgdl = as.numeric(tchol_mgdl)),
                by = "SEQN")
  }

  # ── Anthropometrics (BMX) ───────────────────────────────────────────────
  if (!is.null(raw_list$BMX)) {
    df <- df %>%
      left_join(raw_list$BMX %>%
                dplyr::select(SEQN, bmi = BMXBMI, waist_cm = BMXWAIST) %>%
                mutate(bmi = as.numeric(bmi), waist_cm = as.numeric(waist_cm)),
                by = "SEQN")
  }

  # ── Blood Pressure (BPX) ────────────────────────────────────────────────
  if (!is.null(raw_list$BPX)) {
    df <- df %>%
      left_join(raw_list$BPX %>%
                dplyr::select(SEQN, starts_with("BPX")) %>%
                mutate(
                  sbp_avg = rowMeans(dplyr::select(., starts_with("BPXSY")),
                                     na.rm = TRUE),
                  dbp_avg = rowMeans(dplyr::select(., starts_with("BPXDI")),
                                     na.rm = TRUE)
                ) %>%
                dplyr::select(SEQN, sbp_avg, dbp_avg),
                by = "SEQN")
  }

  # ── Smoking (SMQ) ───────────────────────────────────────────────────────
  if (!is.null(raw_list$SMQ)) {
    df <- df %>%
      left_join(raw_list$SMQ %>%
                dplyr::select(SEQN, SMQ020, SMQ040) %>%
                mutate(
                  smk_100cigs = as.numeric(SMQ020),
                  smk_now     = as.numeric(SMQ040),
                  smoking = case_when(
                    SMQ020 == 2  ~ "Never",
                    SMQ020 == 1 & SMQ040 == 3 ~ "Former",
                    SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current",
                    TRUE ~ NA_character_
                  )
                ) %>%
                dplyr::select(SEQN, smoking),
                by = "SEQN")
  }

  # ── Alcohol (ALQ) ───────────────────────────────────────────────────────
  if (!is.null(raw_list$ALQ)) {
    df <- df %>%
      left_join(raw_list$ALQ %>%
                dplyr::select(SEQN, ALQ111, ALQ121, ALQ130) %>%
                mutate(
                  alc_ever    = as.numeric(ALQ111),
                  alc_freq_yr = as.numeric(ALQ121),
                  alc_avg_drinks = as.numeric(ALQ130),
                  alcohol = case_when(
                    ALQ111 == 2 ~ "Never",
                    ALQ121 >= 2 & ALQ121 <= 6 ~ "Current",
                    ALQ121 == 1 ~ "Former",
                    TRUE ~ NA_character_
                  )
                ) %>%
                dplyr::select(SEQN, alcohol),
                by = "SEQN")
  }

  # ── Medical History (MCQ + DIQ) ─────────────────────────────────────────
  # Diabetes
  if (!is.null(raw_list$DIQ)) {
    df <- df %>%
      left_join(raw_list$DIQ %>%
                dplyr::select(SEQN, DIQ010) %>%
                mutate(dm_diagnosed = ifelse(DIQ010 == 1, 1, 0)),
                by = "SEQN")
  }

  # CVD, Asthma, Cancer
  if (!is.null(raw_list$MCQ)) {
    df <- df %>%
      left_join(raw_list$MCQ %>%
                dplyr::select(SEQN, MCQ160B, MCQ160C, MCQ160D, MCQ160E,
                              MCQ160F, MCQ010, MCQ220) %>%
                mutate(
                  chf     = ifelse(MCQ160B == 1, 1, 0),
                  chd     = ifelse(MCQ160C == 1, 1, 0),
                  angina  = ifelse(MCQ160D == 1, 1, 0),
                  mi      = ifelse(MCQ160E == 1, 1, 0),
                  stroke  = ifelse(MCQ160F == 1, 1, 0),
                  cvd_any = ifelse(chf == 1 | chd == 1 | angina == 1 |
                                   mi == 1 | stroke == 1, 1, 0),
                  asthma  = ifelse(MCQ010 == 1, 1, 0),
                  cancer_hx = ifelse(MCQ220 == 1, 1, 0)
                ) %>%
                dplyr::select(SEQN, chf, chd, angina, mi, stroke, cvd_any,
                              asthma, cancer_hx),
                by = "SEQN")
  }

  # ── Food Security (FSQ) ─────────────────────────────────────────────────
  if (!is.null(raw_list$FSQ)) {
    df <- df %>%
      left_join(raw_list$FSQ %>%
                dplyr::select(SEQN, FSDHH) %>%
                mutate(
                  food_security = factor(FSDHH, levels = c(1, 2, 3, 4),
                    labels = c("Full", "Marginal", "Low", "Very Low")),
                  food_insecure = ifelse(FSDHH %in% c(3, 4), 1, 0)
                ) %>%
                dplyr::select(SEQN, food_security, food_insecure),
                by = "SEQN")
  }

  # ── Physical Activity (PAQ) ─────────────────────────────────────────────
  if (!is.null(raw_list$PAQ)) {
    paq_df <- raw_list$PAQ %>%
      dplyr::select(SEQN, PAQ605, PAQ610, PAD615, PAQ620, PAQ625,
                    PAD630, PAQ650, PAQ655, PAD660, PAQ665, PAQ670, PAD675) %>%
      mutate(
        across(-SEQN, ~ as.numeric(as.character(.))),
        # Vigorous work MET-min/week
        vig_work_mets = ifelse(PAQ605 == 1, PAD615 * PAQ610 * 8.0, 0),
        # Moderate work MET-min/week
        mod_work_mets = ifelse(PAQ620 == 1, PAD630 * PAQ625 * 4.0, 0),
        # Vigorous recreation MET-min/week
        vig_rec_mets  = ifelse(PAQ665 == 1, PAD675 * PAQ670 * 8.0, 0),
        # Moderate recreation MET-min/week
        mod_rec_mets  = ifelse(PAQ650 == 1, PAD660 * PAQ655 * 4.0, 0),
        # Total MET-min/week
        pa_total_mets = vig_work_mets + mod_work_mets +
                        vig_rec_mets + mod_rec_mets,
        pa_active = ifelse(pa_total_mets >= 500, 1, 0)
      )
    df <- df %>% left_join(paq_df %>% dplyr::select(SEQN, pa_total_mets, pa_active),
                          by = "SEQN")
  }

  # ── Nutritional Indices ─────────────────────────────────────────────────
  df <- df %>%
    mutate(
      # PNI (Prognostic Nutritional Index): Onodera 1984
      # PNI = 10 × albumin(g/dL) + 0.005 × lymphocyte(/µL)
      PNI = 10 * albumin_gdl + 0.005 * lymph_abs,

      # CONUT score (Controlling Nutritional Status): Ignacio de Ulíbarri 2005
      # Albumin: ≥3.5=0, 3.0-3.49=2, 2.5-2.99=4, <2.5=6
      conut_alb = case_when(
        albumin_gdl >= 3.5         ~ 0,
        albumin_gdl >= 3.0         ~ 2,
        albumin_gdl >= 2.5         ~ 4,
        !is.na(albumin_gdl)        ~ 6,
        TRUE ~ NA_real_
      ),
      # Lymphocyte: ≥1600=0, 1200-1599=1, 800-1199=2, <800=3
      conut_lymph = case_when(
        lymph_abs >= 1600          ~ 0,
        lymph_abs >= 1200          ~ 1,
        lymph_abs >= 800           ~ 2,
        !is.na(lymph_abs)          ~ 3,
        TRUE ~ NA_real_
      ),
      # Cholesterol: ≥180=0, 140-179=1, 100-139=2, <100=3
      conut_chol = case_when(
        tchol_mgdl >= 180          ~ 0,
        tchol_mgdl >= 140          ~ 1,
        tchol_mgdl >= 100          ~ 2,
        !is.na(tchol_mgdl)         ~ 3,
        TRUE ~ NA_real_
      ),
      CONUT = conut_alb + conut_lymph + conut_chol,

      # CONUT categories
      conut_cat = case_when(
        CONUT <= 1 ~ "Normal (0-1)",
        CONUT <= 4 ~ "Mild (2-4)",
        CONUT >= 5 ~ "Moderate-Severe (5-12)",
        TRUE ~ NA_character_
      ),

      # GNRI (Geriatric Nutritional Risk Index)
      # GNRI = 1.489 × albumin(g/L) + 41.7 × (actual/ideal weight)
      # Simplified: GNRI = 14.89 × albumin(g/dL) + 41.7 × (BMI/22)
      GNRI = 14.89 * albumin_gdl + 41.7 * (bmi / 22),

      # Obesity (WHO: ≥30 for US population)
      obesity = ifelse(bmi >= 30, 1, 0),

      # Hypertension: SBP≥140 or DBP≥90 or on treatment
      htn = ifelse(sbp_avg >= 140 | dbp_avg >= 90, 1, 0),

      # CKD: eGFR (CKD-EPI 2021) — requires creatinine which isn't in BIOPRO
      # Placeholder: will be added if creatinine data is available
      # For now, use CKD self-report from MCQ
    )

  # ── Comorbidity count ────────────────────────────────────────────────────
  df <- df %>%
    mutate(
      comorb_count = rowSums(
        dplyr::select(., dm_diagnosed, htn, cvd_any, obesity, asthma),
        na.rm = TRUE
      )
    )

  # ── Cycle indicator ──────────────────────────────────────────────────────
  df$cycle <- cycle

  cat(sprintf("  [%s] N = %d before exclusion\n", cycle, nrow(df)))
  return(df)
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION C: DII (Dietary Inflammatory Index) calculation
# ═══════════════════════════════════════════════════════════════════════════════

#' Calculate Energy-Adjusted Dietary Inflammatory Index (E-DII)
#'
#' Implements the Shivappa et al. (2014) method:
#' 1. Energy-adjust each nutrient via residual method
#' 2. Z-score using global reference means and SDs
#' 3. Center to percentile (z → percentile → centered percentile × 2 − 1)
#' 4. Multiply by inflammatory effect score
#' 5. Sum all components
#'
#' @param diet_df Data frame with dietary variables (per day 1 24h recall)
#' @return Data frame with SEQN and E_DII score
calculate_dii <- function(diet_df) {
  cat("  Calculating Energy-Adjusted DII...\n")

  # ── Step 1: Energy adjustment (residual method) ──────────────────────────
  # Regress each nutrient on total energy, use residual
  nutrient_vars <- c(
    "DR1TPROT", "DR1TCARB", "DR1TFIBE", "DR1TCHOL",
    "DR1TSFAT", "DR1TMFAT", "DR1TPFAT",
    "DR1TVB1", "DR1TVB2", "DR1TNIAC", "DR1TVB6", "DR1TFOLA",
    "DR1TVB12", "DR1TVC", "DR1TVD", "DR1TVE",
    "DR1TCALC", "DR1TIRON", "DR1TMAGN", "DR1TZINC", "DR1TSELE",
    "DR1TSODI", "DR1TPOTA"
  )

  energy <- diet_df$DR1TKCAL

  dii_df <- data.frame(SEQN = diet_df$SEQN)

  for (nv in nutrient_vars) {
    if (!nv %in% names(diet_df)) next

    nutrient_val <- as.numeric(as.character(diet_df[[nv]]))
    nutrient_val[nutrient_val < 0] <- NA_real_

    # Residual method: nutrient ~ energy
    valid <- !is.na(nutrient_val) & !is.na(energy) & energy > 0
    if (sum(valid) < 50) {
      dii_df[[nv]] <- NA_real_
      next
    }

    lm_fit <- lm(nutrient_val[valid] ~ energy[valid])
    residual_vals <- rep(NA_real_, length(nutrient_val))
    residual_vals[valid] <- resid(lm_fit)

    # Store energy-adjusted value
    dii_df[[nv]] <- residual_vals
  }

  # ── Step 2: Global reference values (Shivappa et al. 2014, Table 1) ──────
  # Global mean, SD, and inflammatory effect score
  dii_ref <- list(
    DR1TPROT = c(mean = 76.4, sd = 26.8, score = 0.02),
    DR1TCARB = c(mean = 274.8, sd = 64.0, score = 0.10),
    DR1TFIBE = c(mean = 18.5, sd = 5.9, score = -0.66),
    DR1TCHOL = c(mean = 278.1, sd = 75.0, score = 0.11),
    DR1TSFAT = c(mean = 24.7, sd = 6.3, score = 0.37),
    DR1TMFAT = c(mean = 27.0, sd = 5.2, score = -0.03),
    DR1TPFAT = c(mean = 13.3, sd = 5.2, score = -0.34),
    DR1TVB1  = c(mean = 1.7, sd = 0.65, score = -0.35),
    DR1TVB2  = c(mean = 1.8, sd = 0.61, score = -0.40),
    DR1TNIAC = c(mean = 26.1, sd = 7.6, score = -0.25),
    DR1TVB6  = c(mean = 1.7, sd = 0.55, score = -0.08),
    DR1TFOLA = c(mean = 282.2, sd = 89.0, score = -0.21),
    DR1TVB12 = c(mean = 5.1, sd = 2.3, score = -0.14),
    DR1TVC   = c(mean = 113.0, sd = 48.0, score = -0.42),
    DR1TVD   = c(mean = 5.6, sd = 3.0, score = -0.42),
    DR1TVE   = c(mean = 8.7, sd = 3.1, score = -0.42),
    DR1TCALC = c(mean = 842.6, sd = 265.0, score = -0.46),
    DR1TIRON = c(mean = 14.1, sd = 4.0, score = -0.16),
    DR1TMAGN = c(mean = 284.2, sd = 74.0, score = -0.48),
    DR1TZINC = c(mean = 8.7, sd = 2.0, score = -0.41),
    DR1TSELE = c(mean = 79.0, sd = 23.0, score = -0.19),
    DR1TSODI = c(mean = 3446.0, sd = 622.0, score = 0.02),
    DR1TPOTA = c(mean = 2711.0, sd = 597.0, score = -0.36)
  )

  # ── Step 3: Z-score → centered percentile → DII contribution ─────────────
  # For each nutrient: Z = (value - global_mean) / global_sd
  # percentile = pnorm(Z) → centered = 2*percentile − 1
  # contribution = centered × inflammatory_score

  dii_components <- matrix(NA_real_, nrow = nrow(dii_df),
                           ncol = length(dii_ref))
  colnames(dii_components) <- names(dii_ref)

  for (nv in names(dii_ref)) {
    if (!nv %in% names(dii_df)) next

    ref <- dii_ref[[nv]]
    z_score <- (dii_df[[nv]] - ref["mean"]) / ref["sd"]
    percentile <- pnorm(z_score)
    centered <- 2 * percentile - 1
    dii_components[, nv] <- centered * ref["score"]
  }

  dii_df$E_DII <- rowSums(dii_components, na.rm = FALSE)
  # Flag those with many missing nutrients
  dii_df$E_DII_n_components <- rowSums(!is.na(dii_components))

  # Only compute DII if ≥20 of 23 components are available
  dii_df$E_DII[dii_df$E_DII_n_components < 20] <- NA_real_

  cat(sprintf("  E-DII computed for %d participants (mean %.2f, SD %.2f)\n",
              sum(!is.na(dii_df$E_DII)),
              mean(dii_df$E_DII, na.rm = TRUE),
              sd(dii_df$E_DII, na.rm = TRUE)))

  return(dii_df)
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION D: Master loading function
# ═══════════════════════════════════════════════════════════════════════════════

#' Load, merge, and derive all NHANES variables for specified cycles
load_and_derive_nhanes <- function(cycles = CYCLES, use_cache = TRUE,
                                   min_age = 60) {
  all_cycles <- list()

  for (cyc in cycles) {
    # Download/load
    load_nhanes_cycle(cyc, use_cache = use_cache)

    # Get the cached raw data
    cycle_letter <- switch(cyc, "2011-2012" = "G", "2013-2014" = "H")
    raw_name <- paste0("raw_", cycle_letter)
    raw_list <- get(raw_name, envir = .GlobalEnv)

    # Derive variables
    df <- derive_cycle_variables(raw_list, cyc)
    all_cycles[[cyc]] <- df
  }

  # Merge cycles
  df_all <- bind_rows(all_cycles)

  # ── Apply exclusion criteria ───────────────────────────────────────────────
  df_all <- df_all %>%
    filter(
      age >= min_age,                                   # Age ≥ 60
      !is.na(cog_z_composite),                          # Cognitive data complete
      !is.na(phq9_total),                               # PHQ-9 complete
      !is.na(albumin_gdl),                              # Has blood work
      !is.na(bmi)                                       # Has anthropometrics
    )

  # Add exclusion flags for sensitivity analyses
  df_all <- df_all %>%
    mutate(
      excl_stroke  = ifelse(stroke == 1, 1, 0),
      excl_cancer  = ifelse(cancer_hx == 1, 1, 0),
      excl_low_alb = ifelse(albumin_gdl < 3.0, 1, 0)
    )

  cat(sprintf("\nFinal analytic sample: N = %d\n", nrow(df_all)))
  cat(sprintf("  Cycle 2011-2012: %d\n",
              sum(df_all$cycle == "2011-2012")))
  cat(sprintf("  Cycle 2013-2014: %d\n",
              sum(df_all$cycle == "2013-2014")))
  cat(sprintf("  Age range: %.0f-%.0f\n",
              min(df_all$age), max(df_all$age)))
  cat(sprintf("  Female: %.1f%%\n",
              mean(df_all$sex_binary == 1, na.rm = TRUE) * 100))

  return(df_all)
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION E: Execute
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── NHANES Data Loading & Variable Derivation ──\n")

# Load & derive
df <- load_and_derive_nhanes(cycles = CYCLES, use_cache = TRUE, min_age = 60)

# Calculate DII (requires diet data)
# Merge diet data back for DII calculation
dii_dfs <- list()
for (cyc in CYCLES) {
  cycle_letter <- switch(cyc, "2011-2012" = "G", "2013-2014" = "H")
  raw_name <- paste0("raw_", cycle_letter)
  raw_list <- get(raw_name, envir = .GlobalEnv)

  if (!is.null(raw_list$DR1TOT)) {
    dii_res <- calculate_dii(raw_list$DR1TOT)
    dii_res$cycle <- cyc
    dii_dfs[[cyc]] <- dii_res
  }
}

if (length(dii_dfs) > 0) {
  dii_all <- bind_rows(dii_dfs)
  df <- df %>% left_join(dii_all %>% dplyr::select(SEQN, E_DII), by = "SEQN")

  # Categorize E-DII
  df <- df %>%
    mutate(
      dii_quartile = factor(ntile(E_DII, 4), levels = 1:4,
                            labels = c("Q1 (Anti-inflam)", "Q2", "Q3",
                                       "Q4 (Pro-inflam)")),
      dii_tertile  = factor(ntile(E_DII, 3), levels = 1:3,
                            labels = c("T1 (Low)", "T2 (Mid)", "T3 (High)"))
    )

  cat(sprintf("\nE-DII merged. Mean: %.2f (SD %.2f), N with DII: %d\n",
              mean(df$E_DII, na.rm = TRUE),
              sd(df$E_DII, na.rm = TRUE),
              sum(!is.na(df$E_DII))))
}

# Save derived dataset
saveRDS(df, file.path(DATA_DIR, "nhanes_2011_2014_derived.rds"))
cat(sprintf("\nDerived dataset saved to %s\n",
            file.path(DATA_DIR, "nhanes_2011_2014_derived.rds")))

cat("── 01_load_and_derive.R complete ──\n")
