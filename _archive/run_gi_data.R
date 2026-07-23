# run_gi_data.R — NHANES GI data download (2005-2016, 6 cycles)
# Run: Rscript run_gi_data.R
#
# GI cancer identification: MCQ240[letter] age-at-diagnosis variables
# Non-missing = participant had that specific cancer
# GI sites: G=Colon, H=Esophageal, L=Leukemia(not GI), M=Liver,
#           T=Pancreatic, V=Rectal, Z=Stomach, I=Gallbladder

if (!require("nhanesA")) install.packages("nhanesA", repos="https://cloud.r-project.org")
if (!require("dplyr")) install.packages("dplyr", repos="https://cloud.r-project.org")
library(nhanesA)
library(dplyr)

# Convert haven_labelled to plain vectors
unlabel <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "haven_labelled")) as.numeric(x) else as.numeric(x)
}

GI_DATA_DIR <- file.path(dirname(normalizePath(getwd())), "data", "gi_analysis")
dir.create(GI_DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# GI cancer variables: MCQ240[letter] = age at first dx for that site
GI_LETTERS <- c(G="Colon", H="Esophageal", I="Gallbladder",
                M="Liver", T="Pancreatic", V="Rectal", Z="Stomach")

# Use 2005-2016 (exclude 2017-2018: different variable structure, short FU)
GI_CYCLES <- c("2005-2006"="D","2007-2008"="E","2009-2010"="F",
               "2011-2012"="G","2013-2014"="H","2015-2016"="I")

all_data <- list()
for (cyc in names(GI_CYCLES)) {
  sfx <- GI_CYCLES[cyc]
  cat(sprintf("\n=== %s ===\n", cyc)); flush.console()

  cat("  DEMO..."); flush.console()
  demo <- nhanesA::nhanes(paste0("DEMO_", sfx))
  cat(sprintf(" %d (all ages)\n", nrow(demo))); flush.console()

  cat("  MCQ..."); flush.console()
  mcq <- nhanesA::nhanes(paste0("MCQ_", sfx))

  # GI tumor: non-missing age-at-dx for any GI site variable
  gi_vars <- paste0("MCQ240", names(GI_LETTERS))
  gi_sites_present <- gi_vars[gi_vars %in% names(mcq)]

  # Initialize
  mcq$gi_tumor <- 0L
  mcq$gi_site <- NA_character_
  mcq$gi_count <- 0L

  for (gv in gi_sites_present) {
    ages <- unlabel(mcq[[gv]])
    is_gi <- !is.na(ages) & ages > 0
    mcq$gi_tumor[is_gi] <- 1L
    letter <- sub("MCQ240", "", gv)
    site_name <- GI_LETTERS[letter]
    mcq$gi_site[is_gi & is.na(mcq$gi_site)] <- site_name
    mcq$gi_count[is_gi] <- mcq$gi_count[is_gi] + 1L
  }

  mcq$any_cancer <- ifelse(is.na(unlabel(mcq$MCQ220)), 0L,
                           as.integer(unlabel(mcq$MCQ220) == 1))
  cat(" merged")

  df <- merge(demo, mcq[, c("SEQN", "gi_tumor", "gi_site", "gi_count", "any_cancer")],
              by = "SEQN", all.x = TRUE)

  cat(" BIOPRO..."); flush.console()
  bp <- tryCatch(nhanesA::nhanes(paste0("BIOPRO_", sfx)), error=\(e)NULL)
  if (!is.null(bp)) {
    idx <- match(df$SEQN, bp$SEQN)
    df$albumin_gdl <- unlabel(bp[["LBDSALSI"]])[idx]
    if ("LBDSCRSI" %in% names(bp)) df$crp_mgdl <- unlabel(bp[["LBDSCRSI"]])[idx]
    else df$crp_mgdl <- NA_real_
  }

  cat(" CBC..."); flush.console()
  cb <- tryCatch(nhanesA::nhanes(paste0("CBC_", sfx)), error=\(e)NULL)
  if (!is.null(cb)) {
    idx <- match(df$SEQN, cb$SEQN)
    wbc <- unlabel(cb[["LBXWBCSI"]])[idx]
    lpct <- unlabel(cb[["LBXLYPCT"]])[idx]
    df$lymph_abs <- wbc * lpct / 100 * 1000
  }

  cat(" TCHOL..."); flush.console()
  tc <- tryCatch(nhanesA::nhanes(paste0("TCHOL_", sfx)), error=\(e)NULL)
  if (!is.null(tc)) {
    idx <- match(df$SEQN, tc$SEQN)
    df$tchol_mgdl <- unlabel(tc[["LBXTC"]])[idx]
  }

  cat(" BMX..."); flush.console()
  bx <- tryCatch(nhanesA::nhanes(paste0("BMX_", sfx)), error=\(e)NULL)
  if (!is.null(bx)) {
    idx <- match(df$SEQN, bx$SEQN)
    df$bmi <- unlabel(bx[["BMXBMI"]])[idx]
  }

  n0 <- nrow(df)
  df <- df[!is.na(df$albumin_gdl) & !is.na(df$lymph_abs) &
           !is.na(df$tchol_mgdl) & !is.na(df$bmi), , drop = FALSE]
  cat(sprintf(" complete: %d/%d\n", nrow(df), n0)); flush.console()
  if (nrow(df) == 0) next

  # Nutrition indices
  df$PNI <- 10 * df$albumin_gdl + 0.005 * df$lymph_abs

  conut_a <- rep(6, nrow(df))
  conut_a[df$albumin_gdl >= 2.5] <- 4
  conut_a[df$albumin_gdl >= 3.0] <- 2
  conut_a[df$albumin_gdl >= 3.5] <- 0

  conut_l <- rep(3, nrow(df))
  conut_l[df$lymph_abs >= 800] <- 2
  conut_l[df$lymph_abs >= 1200] <- 1
  conut_l[df$lymph_abs >= 1600] <- 0

  conut_c <- rep(3, nrow(df))
  conut_c[df$tchol_mgdl >= 100] <- 2
  conut_c[df$tchol_mgdl >= 140] <- 1
  conut_c[df$tchol_mgdl >= 180] <- 0

  df$CONUT <- conut_a + conut_l + conut_c
  df$GNRI <- 14.89 * df$albumin_gdl + 41.7 * (df$bmi / 22)

  df$sex <- factor(ifelse(df$RIAGENDR == 2, "Female", "Male"))
  reth1 <- unlabel(df$RIDRETH1)
  df$race_eth <- case_when(
    reth1 %in% c(1,2) ~ "Hispanic",
    reth1 == 3 ~ "Non-Hispanic White",
    reth1 == 4 ~ "Non-Hispanic Black",
    TRUE ~ "Other")
  # DMDEDUC2 might be factor in some cycles
  educ_val <- unlabel(df$DMDEDUC2)
  df$edu_binary <- ifelse(is.na(educ_val), 0L, as.integer(educ_val >= 4))
  df$income_pir <- df$INDFMPIR
  df$wt_mec_2yr <- df$WTMEC2YR
  df$psu <- df$SDMVPSU
  df$strata <- df$SDMVSTRA
  df$cycle <- cyc

  # Keep only needed columns for bind_rows compatibility
  df$age <- df$RIDAGEYR
  keep_cols <- c("SEQN", "gi_tumor", "gi_site", "gi_count", "any_cancer",
                 "albumin_gdl", "crp_mgdl", "lymph_abs", "tchol_mgdl", "bmi",
                 "age", "PNI", "CONUT", "GNRI", "sex", "race_eth", "edu_binary",
                 "income_pir", "wt_mec_2yr", "psu", "strata", "cycle")
  all_data[[cyc]] <- df[, intersect(keep_cols, names(df))]
  cat(sprintf("  GI tumors: %d\n", sum(df$gi_tumor, na.rm=TRUE))); flush.console()
}

df_all <- bind_rows(all_data)
cat(sprintf("\n=== TOTAL: N=%d, GI tumors=%d ===\n",
            nrow(df_all), sum(df_all$gi_tumor, na.rm=TRUE)))

# GI site breakdown
gi_df <- df_all[df_all$gi_tumor == 1, ]
cat("GI sites:\n")
print(table(gi_df$gi_site, useNA = "ifany"))

saveRDS(df_all, file.path(GI_DATA_DIR, "nhanes_gi_nutrition_raw.rds"))
cat("Saved.\n")
