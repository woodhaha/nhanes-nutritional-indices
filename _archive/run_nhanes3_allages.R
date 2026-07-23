# run_nhanes3_allages.R — NHANES III extraction WITHOUT age filter
# Rscript run_nhanes3_allages.R

library(dplyr)

DATA_DIR <- "D:/Researching/NHANES_aged_GI_tumor_nutrition/data/nhanes3"
GI_DATA_DIR <- "D:/Researching/NHANES_aged_GI_tumor_nutrition/data/gi_analysis"

# ── Read columns from .dat files (same as run_nhanes3_v2.R) ──────────────
read_columns <- function(file, col_spec) {
  con <- file(file, "r", blocking=TRUE)
  on.exit(close(con))
  lines <- readLines(con, warn=FALSE)
  n <- length(lines)
  result <- list()
  for (spec in col_spec) result[[spec$name]] <- vector("character", n)
  for (i in seq_len(n)) {
    line <- lines[i]
    for (spec in col_spec) {
      result[[spec$name]][i] <- substr(line, spec$start, spec$end)
    }
  }
  df <- as.data.frame(result, stringsAsFactors=FALSE)
  for (nm in names(df)) df[[nm]] <- suppressWarnings(as.numeric(trimws(df[[nm]])))
  df
}

adult_cols <- list(
  list(name="SEQN", start=1, end=5),
  list(name="HSSEX", start=15, end=15),
  list(name="HSAGEIR", start=18, end=19),
  list(name="DMAETHNR", start=12, end=12),
  list(name="DMARACER", start=13, end=13),
  list(name="DMPPIR", start=36, end=41),
  list(name="WTPFEX6", start=61, end=69),
  list(name="SDPPSU6", start=43, end=43),
  list(name="SDPSTRA6", start=44, end=45),
  list(name="HFHEDUCR", start=1423, end=1424),
  list(name="HFF1", start=1356, end=1356),
  list(name="HFF20R", start=1411, end=1412),
  list(name="HFF19R", start=1409, end=1410),
  list(name="HAC1N", start=1478, end=1478),
  list(name="HAC1O", start=1479, end=1479),
  list(name="HAC3OS", start=1527, end=1528)
)

lab_cols <- list(
  list(name="SEQN", start=1, end=5),
  list(name="AMP", start=1846, end=1848),
  list(name="WCP", start=1273, end=1277),
  list(name="LMPPCNT", start=1283, end=1287),
  list(name="LMP", start=1298, end=1302),
  list(name="TCP", start=1598, end=1600)
)

exam_cols <- list(
  list(name="SEQN", start=1, end=5),
  list(name="BMPBMI", start=1524, end=1527),
  list(name="BMPWT", start=1508, end=1513),
  list(name="BMPHT", start=1528, end=1532)
)

cat("=== NHANES III all ages ===\n")
adult <- read_columns(file.path(DATA_DIR, "adult.dat"), adult_cols)
lab <- read_columns(file.path(DATA_DIR, "lab.dat"), lab_cols)
exam <- read_columns(file.path(DATA_DIR, "exam.dat"), exam_cols)

df <- merge(adult, lab, by="SEQN", all.x=TRUE)
df <- merge(df, exam, by="SEQN", all.x=TRUE)
cat(sprintf("Raw: %d rows\n", nrow(df)))

# Plausibility (same as before, no age filter)
df <- df %>% filter(
  AMP >= 1.0 & AMP <= 7.0,
  LMP >= 0.1 & LMP <= 20,
  TCP >= 50 & TCP <= 500,
  BMPBMI >= 10 & BMPBMI <= 80
)
cat(sprintf("After plausibility: %d\n", nrow(df)))

# Cancer classification
GI_CODES <- c(6, 12, 13, 14, 15, 25)
df <- df %>% mutate(
  any_cancer = ifelse((!is.na(HAC1N) & HAC1N==1) | (!is.na(HAC1O) & HAC1O==1), 1L, 0L),
  gi_tumor = ifelse(!is.na(HAC3OS) & HAC3OS %in% GI_CODES, 1L, 0L),
  gi_site = case_when(
    HAC3OS == 6 ~ "Colon", HAC3OS == 12 ~ "Stomach",
    HAC3OS == 13 ~ "Esophagus", HAC3OS == 14 ~ "Pancreas",
    HAC3OS == 15 ~ "Liver", HAC3OS == 25 ~ "Gallbladder",
    TRUE ~ NA_character_),
  age = HSAGEIR,
  sex = ifelse(HSSEX==1, "Male", "Female"),
  race_eth = case_when(
    DMAETHNR %in% c(1,2) ~ "Hispanic",
    DMARACER == 1 ~ "Non-Hispanic White",
    DMARACER == 2 ~ "Non-Hispanic Black",
    TRUE ~ "Other"),
  edu_binary = ifelse(HFHEDUCR >= 4, 1, 0),
  income_pir = DMPPIR,
  albumin_gdl = AMP,
  lymph_abs = LMP * 1000,
  tchol_mgdl = TCP,
  bmi = BMPBMI,
  PNI = 10*AMP + 0.005*(LMP*1000),
  CONUT = case_when(AMP>=3.5~0, AMP>=3.0~2, AMP>=2.5~4, TRUE~6) +
          case_when(LMP*1000>=1600~0, LMP*1000>=1200~1, LMP*1000>=800~2, TRUE~3) +
          case_when(TCP>=180~0, TCP>=140~1, TCP>=100~2, TRUE~3),
  GNRI = 14.89*AMP + 41.7*(bmi/22),
  cycle = "NHANES_III_1988-1994"
)

cat(sprintf("GI tumors (all ages): %d\n", sum(df$gi_tumor, na.rm=TRUE)))
cat(sprintf("Age range: %d-%d\n", min(df$age, na.rm=TRUE), max(df$age, na.rm=TRUE)))

# Mortality
mort <- readRDS(file.path(GI_DATA_DIR, "mort_nhanes3.rds"))
df <- merge(df, mort[, c("SEQN","eligstat","mortstat","ucod_leading","permth_int")],
            by="SEQN", all.x=TRUE)
df <- df[df$eligstat==1 & !is.na(df$mortstat), , drop=FALSE]
df$surv_years <- df$permth_int/12
df$death <- df$mortstat
df$gi_cancer_death <- ifelse(df$death==1 & df$ucod_leading==2 & df$gi_tumor==1, 1L, 0L)

saveRDS(df, file.path(GI_DATA_DIR, "nhanes3_gi_allages.rds"))
cat(sprintf("Saved: N=%d, GI tumors=%d, age range %d-%d\n",
    nrow(df), sum(df$gi_tumor, na.rm=TRUE),
    min(df$age, na.rm=TRUE), max(df$age, na.rm=TRUE)))
cat("Done.\n")
