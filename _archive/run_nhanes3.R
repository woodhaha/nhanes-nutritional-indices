# run_nhanes3.R — NHANES III (1988-1994) data extraction + GI analysis
# Rscript run_nhanes3.R

library(dplyr); library(survival); library(broom)

DATA_DIR <- "D:/Researching/NHANES_aged_GI_tumor_nutrition/data/nhanes3"
GI_DATA_DIR <- "D:/Researching/data/gi_analysis"
GI_RESULTS_DIR <- "D:/Researching/results/gi_analysis"
dir.create(GI_RESULTS_DIR, recursive=TRUE, showWarnings=FALSE)

# ── Parse SAS INPUT format → data frame of variable positions ────────────────
parse_sas_input <- function(sas_file) {
  sas <- readLines(sas_file, warn=FALSE)
  idx <- grep("^\\s+INPUT", sas, ignore.case=TRUE)[1]
  if (is.na(idx)) stop("No INPUT found")
  lines <- sas[idx:length(sas)]
  # Stop at LABEL or first blank after a variable with =
  end <- grep("^\\s+LABEL|^\\s*;|^$", lines)[1]
  lines <- lines[2:(end-1)]  # skip INPUT, stop before LABEL

  vars <- data.frame(name=character(), start=integer(), end=integer(),
                     type=character(), stringsAsFactors=FALSE)
  for (ln in lines) {
    ln <- gsub("\\s+", " ", trimws(ln))
    if (nchar(ln) == 0 || grepl("^[*/]", ln)) next
    parts <- strsplit(ln, "\\s+")[[1]]
    name <- parts[1]
    pos <- parts[2]
    if (grepl("^\\$", pos)) {  # character var
      type <- "character"
      pos <- sub("^\\$", "", pos)
    } else {
      type <- "numeric"
    }
    if (grepl("-", pos)) {
      range <- as.integer(strsplit(pos, "-")[[1]])
      vars <- rbind(vars, data.frame(name, start=range[1], end=range[2], type, stringsAsFactors=FALSE))
    } else {
      pos_val <- suppressWarnings(as.integer(pos))
      if (!is.na(pos_val)) {
        vars <- rbind(vars, data.frame(name, start=pos_val, end=pos_val, type, stringsAsFactors=FALSE))
      }
    }
  }
  vars
}

# ── Parse WIDTHS from LENGTH section (for complete column widths) ────────────
parse_sas_length <- function(sas_file) {
  sas <- readLines(sas_file, warn=FALSE)
  idx <- grep("^\\s+LENGTH", sas, ignore.case=TRUE)[1]
  if (is.na(idx)) return(NULL)
  lines <- sas[idx:length(sas)]
  end <- grep("^\\s*;", lines)[1]
  lines <- lines[2:(end-1)]

  widths <- list()
  for (ln in lines) {
    ln <- gsub("\\s+", " ", trimws(ln))
    if (nchar(ln) == 0) next
    parts <- strsplit(ln, "\\s+")[[1]]
    if (length(parts) >= 2) {
      w <- suppressWarnings(as.integer(parts[2]))
      if (!is.na(w)) widths[[parts[1]]] <- w
    }
  }
  widths
}

# ── Read .dat file by parsing SAS format ─────────────────────────────────────
read_nhanes3_dat <- function(dat_path, sas_path, keep_vars=NULL) {
  vars <- parse_sas_input(sas_path)
  widths <- parse_sas_length(sas_path)

  if (!is.null(keep_vars)) {
    vars <- vars[vars$name %in% keep_vars, , drop=FALSE]
  }
  # Calculate column widths
  sorted <- vars[order(vars$start),]
  widths_list <- sorted$end - sorted$start + 1
  # Also need to handle gaps between variables - add dummy widths
  rec_len <- 3348  # from adult.sas LRECL

  cat(sprintf("  Reading %d variables from %s\n", nrow(sorted), basename(dat_path)))
  df <- read.fwf(dat_path,
    widths=widths_list,
    col.names=sorted$name,
    na.strings=c(".", "", " ", "  ", "   ", "    "),
    buffersize=50000)
  # Convert character-like numeric columns
  for (nm in names(df)) {
    if (sorted$type[sorted$name==nm] == "character") {
      df[[nm]] <- trimws(as.character(df[[nm]]))
      df[[nm]][df[[nm]] == ""] <- NA_character_
    } else {
      df[[nm]] <- suppressWarnings(as.numeric(as.character(df[[nm]])))
    }
  }
  df
}

# ── Step 1: Parse variable positions ─────────────────────────────────────────
cat("=== Parsing SAS format files ===\n")

# Key variables from adult: SEQN, age, sex, race, education, PIR, cancer
# Key variables from lab: AMP (albumin), TCP/CHP (cholesterol),
#   WCP (WBC), LMPPCNT (lymphocyte pct), LMP (lymphocyte #)
# Key variables from exam: BMI (calculated from weight/height)

# We need full variable lists to ensure correct widths
adult_vars_all <- parse_sas_input("data/nhanes3/adult.sas")
lab_vars_all <- parse_sas_input("data/nhanes3/lab.sas")
exam_vars_all <- parse_sas_input("data/nhanes3/exam.sas")
cat(sprintf("Adult vars: %d, Lab vars: %d, Exam vars: %d\n",
    nrow(adult_vars_all), nrow(lab_vars_all), nrow(exam_vars_all)))

# ── Step 2: Read the data files (all vars needed for correct widths) ─────────
cat("\n=== Reading data files ===\n")

cat("Reading adult...\n")
# Actually let me use the simple approach - read all with correct widths
adult_widths <- adult_vars_all$end - adult_vars_all$start + 1
adult <- read.fwf("data/nhanes3/adult.dat",
  widths=adult_widths, col.names=adult_vars_all$name,
  na.strings=c(".",""," "), buffersize=50000)
# Convert
for (nm in names(adult)) adult[[nm]] <- suppressWarnings(as.numeric(as.character(adult[[nm]])))
cat(sprintf("  Adult: %d rows, %d cols\n", nrow(adult), ncol(adult)))

cat("Reading lab...\n")
lab_widths <- lab_vars_all$end - lab_vars_all$start + 1
lab <- read.fwf("data/nhanes3/lab.dat",
  widths=lab_widths, col.names=lab_vars_all$name,
  na.strings=c(".",""," "), buffersize=50000)
for (nm in names(lab)) lab[[nm]] <- suppressWarnings(as.numeric(as.character(lab[[nm]])))
cat(sprintf("  Lab: %d rows, %d cols\n", nrow(lab), ncol(lab)))

cat("Reading exam...\n")
exam_widths <- exam_vars_all$end - exam_vars_all$start + 1
exam <- read.fwf("data/nhanes3/exam.dat",
  widths=exam_widths, col.names=exam_vars_all$name,
  na.strings=c(".",""," "), buffersize=50000)
for (nm in names(exam)) exam[[nm]] <- suppressWarnings(as.numeric(as.character(exam[[nm]])))
cat(sprintf("  Exam: %d rows, %d cols\n", nrow(exam), ncol(exam)))

# ── Step 3: HAC3OS codes for cancer site ─────────────────────────────────────
# From NHANES III documentation:
# HAC1N = skin cancer (yes/no)
# HAC1O = other cancer (yes/no)
# HAC3OS = cancer site code where first told (used for "other" cancer site)
# Codes for GI sites:
# Need to check the codebook - typical values would be:
# 10-14 = GI sites (colon, rectum, stomach, esophagus, pancreas, liver)

cat("\n=== Cancer variables check ===\n")
cat("HAC1N (skin cancer):", table(adult$HAC1N, useNA="ifany"), "\n")
cat("HAC1O (other cancer):", table(adult$HAC1O, useNA="ifany"), "\n")
cat("HAC3OS unique values:", sort(unique(adult$HAC3OS[!is.na(adult$HAC3OS)])), "\n")

# ── Step 4: Merge datasets ───────────────────────────────────────────────────
cat("\n=== Merging ===\n")

# Keep all age >= 20 for now (filter to >=60 later)
df <- adult

# Merge lab
df <- merge(df, lab[, c("SEQN", "AMP", "AMPSI", "TCP", "TCPSI",
                         "WCP", "WCPSI", "LMPPCNT", "LMP", "CRP",
                         "TPP", "GBP", "CHP", "CHPSI", "HGP", "HGPSI")],
            by="SEQN", all.x=TRUE, suffixes=c("",".lab"))

# Merge exam
exam_vars <- names(exam)
exam_keep <- intersect(exam_vars, c("SEQN", "BMPBMI", "BMPWT", "BMPHT",
                                     "BMPWST", "BMARECUM"))
df <- merge(df, exam[, exam_keep], by="SEQN", all.x=TRUE, suffixes=c("",".ex"))
cat(sprintf("Merged: %d rows\n", nrow(df)))

# ── Step 5: Derive variables ─────────────────────────────────────────────────
cat("\n=== Deriving variables ===\n")

# Age filter
df <- df %>% filter(HSAGEIR >= 60)
cat(sprintf("Age >= 60: %d rows\n", nrow(df)))

# Cancer classification
# HAC1N=1 skin cancer, HAC1O=1 other cancer
# HAC3OS codes for cancer site (need to map from codebook)
# Based on NHANES III codebook, HAC3OS values:
# 1 = Oral, 2 = Skin (non-melanoma), 3 = Skin (melanoma),
# 4 = Breast, 5 = Lung, 6 = Colon/Rectal, 7 = Prostate,
# 8 = Bladder, 9 = Uterus, 10 = Cervix, 11 = Ovary,
# 12 = Stomach, 13 = Esophagus, 14 = Pancreas, 15 = Liver,
# 16 = Kidney, 17 = Thyroid, 18 = Brain, 19 = Lymphoma,
# 20 = Leukemia, 21 = Other, 22 = Bone, 23 = Testicular,
# 24 = Blood, 25 = Gallbladder, 26 = Larynx
GI_CODES_N3 <- c(6, 12, 13, 14, 15, 25)  # colon/rectal, stomach, esophagus, pancreas, liver, gallbladder

df <- df %>% mutate(
  # Has any cancer (HAC1N=1 or HAC1O=1)
  any_cancer = ifelse(!is.na(HAC1N) & HAC1N==1, 1L,
                      ifelse(!is.na(HAC1O) & HAC1O==1, 1L, 0L)),
  # GI tumor: site code in GI codes
  gi_tumor = ifelse(!is.na(HAC3OS) & HAC3OS %in% GI_CODES_N3, 1L, 0L),
  # Site name mapping
  gi_site = case_when(
    HAC3OS == 6 ~ "Colon",
    HAC3OS == 12 ~ "Stomach",
    HAC3OS == 13 ~ "Esophagus",
    HAC3OS == 14 ~ "Pancreas",
    HAC3OS == 15 ~ "Liver",
    HAC3OS == 25 ~ "Gallbladder",
    TRUE ~ NA_character_
  ),
  # Demographics
  age = HSAGEIR,
  sex = ifelse(HSSEX==1, "Male", "Female"),
  race_eth = case_when(
    DMAETHNR %in% c(1,2) ~ "Hispanic",
    DMARACER == 1 ~ "Non-Hispanic White",
    DMARACER == 2 ~ "Non-Hispanic Black",
    TRUE ~ "Other"
  ),
  # Nutrition indices
  albumin_gdl = AMP,  # serum albumin in g/dL
  lymph_abs = LMP,    # lymphocyte number from Coulter
  tchol_mgdl = TCP,   # total cholesterol mg/dL
  bmi = BMPBMI,       # calculated BMI from exam
  # Nutrition indices
  PNI = 10 * AMP + 0.005 * LMP,
  CONUT = case_when(AMP >= 3.5 ~ 0, AMP >= 3.0 ~ 2, AMP >= 2.5 ~ 4, TRUE ~ 6) +
          case_when(LMP >= 1600 ~ 0, LMP >= 1200 ~ 1, LMP >= 800 ~ 2, TRUE ~ 3) +
          case_when(TCP >= 180 ~ 0, TCP >= 140 ~ 1, TCP >= 100 ~ 2, TRUE ~ 3),
  GNRI = 14.89*AMP + 41.7*(bmi/22)
)

cat(sprintf("GI tumors in NHANES III: %d\n", sum(df$gi_tumor, na.rm=TRUE)))
if (sum(df$gi_tumor, na.rm=TRUE) > 0) {
  cat("By site:\n")
  print(table(df$gi_site[df$gi_tumor==1], useNA="ifany"))
}

# Complete cases
df <- df %>% filter(!is.na(PNI), !is.na(CONUT), !is.na(GNRI))
cat(sprintf("Complete nutrition: %d (GI tumors: %d)\n",
            nrow(df), sum(df$gi_tumor, na.rm=TRUE)))

# ── Step 6: Load + merge mortality ───────────────────────────────────────────
cat("\n=== Mortality merge ===\n")
mort_path <- file.path(GI_DATA_DIR, "mort_nhanes3.rds")
if (file.exists(mort_path)) {
  mort <- readRDS(mort_path)
} else {
  base_url <- "https://ftp.cdc.gov/pub/health_statistics/nchs/datalinkage/linked_mortality"
  tmp <- tempfile(fileext=".dat")
  cat("Downloading NHANES III mortality file...\n")
  download.file(file.path(base_url, "NHANES_III_MORT_2019_PUBLIC.dat"),
                tmp, mode="wb", quiet=TRUE)
  # NHANES III mortality file has a different record layout!
  # From the SAS read-in program, NHANES III uses same layout as continuous
  mort <- read.fwf(tmp,
    widths=c(6,8,1,1,3,1,1,21,3,3,13),
    col.names=c("SEQN","b1","eligstat","mortstat","ucod_leading",
                "diabetes","hyperten","b2","permth_int","permth_exm","b3"),
    colClasses=c("integer","NULL","integer","integer","character",
                 "integer","integer","NULL","integer","integer","NULL"),
    na.strings=c("","."," "))
  mort$ucod_leading[mort$ucod_leading==""] <- NA_character_
  mort$ucod_leading <- as.integer(mort$ucod_leading)
  unlink(tmp)
  saveRDS(mort, mort_path)
}

df <- merge(df, mort[, c("SEQN", "eligstat", "mortstat", "ucod_leading",
                          "permth_int", "permth_exm")],
            by="SEQN", all.x=TRUE)
cat(sprintf("Mortality merged: %d rows\n", nrow(df)))

# Filter to eligible + known vital status
df <- df[df$eligstat==1 & !is.na(df$mortstat), , drop=FALSE]
cat(sprintf("Follow-up eligible: %d\n", nrow(df)))

# Survival vars
df$surv_years <- df$permth_int / 12
df$death <- df$mortstat
df$gi_cancer_death <- ifelse(df$death==1 & df$ucod_leading==2 & df$gi_tumor==1, 1L, 0L)

cat(sprintf("  Deaths: %d (%.1f%%)\n", sum(df$death, na.rm=TRUE),
            mean(df$death, na.rm=TRUE)*100))
cat(sprintf("  GI tumor: %d\n", sum(df$gi_tumor==1, na.rm=TRUE)))
cat(sprintf("  GI cancer deaths: %d\n", sum(df$gi_cancer_death, na.rm=TRUE)))

# Save NHANES III sample
saveRDS(df, file.path(GI_DATA_DIR, "nhanes3_gi_merged.rds"))
cat("Saved.\n")

# ── Step 7: Combine with existing NHANES 2005-2016 data ─────────────────────
cat("\n=== Combining with NHANES 2005-2016 ===\n")
df_existing <- readRDS(file.path(GI_DATA_DIR, "nhanes_gi_nutrition_raw.rds"))
mort_existing <- readRDS(file.path(GI_DATA_DIR, "mort_2019.rds"))
df_existing <- merge(df_existing, mort_existing[,c("SEQN","eligstat","mortstat","ucod_leading","permth_int")],
                     by="SEQN", all.x=TRUE)
df_existing <- df_existing[df_existing$eligstat==1 & !is.na(df_existing$mortstat),]
df_existing$surv_years <- df_existing$permth_int/12
df_existing$death <- df_existing$mortstat
df_existing$gi_cancer_death <- ifelse(df_existing$death==1 & df_existing$ucod_leading==2 & df_existing$gi_tumor==1, 1L, 0L)
df_existing$cycle <- paste0("NHANES_", df_existing$cycle)
df_existing$gi_status <- case_when(df_existing$gi_tumor==1 ~ "GI Tumor",
                                    df_existing$any_cancer==0 ~ "Non-Cancer",
                                    TRUE ~ "Other Cancer")

# NHANES III
df3 <- df
df3$race_eth <- as.character(df3$race_eth)
df3$sex <- as.character(df3$sex)
df3$cycle <- "NHANES_III_1988-1994"
df3$gi_status <- case_when(df3$gi_tumor==1 ~ "GI Tumor",
                            df3$any_cancer==0 ~ "Non-Cancer",
                            TRUE ~ "Other Cancer")

# Combine
common_cols <- intersect(names(df_existing), names(df3))
df_combined <- bind_rows(
  df_existing[, common_cols],
  df3[, common_cols]
)
cat(sprintf("\n=== COMBINED: N=%d, GI tumors=%d ===\n",
            nrow(df_combined), sum(df_combined$gi_tumor, na.rm=TRUE)))
cat(sprintf("  NHANES III GI tumors: %d\n",
            sum(df3$gi_tumor, na.rm=TRUE)))
cat(sprintf("  2005-2016 GI tumors: %d\n",
            sum(df_existing$gi_tumor, na.rm=TRUE)))
cat(sprintf("  Total GI cancer deaths: %d\n",
            sum(df_combined$gi_cancer_death, na.rm=TRUE)))

# ── Step 8: Combined survival analysis ───────────────────────────────────────
cat("\n── Combined: All-cause Survival ──\n")
df_gi <- df_combined %>% filter(gi_tumor==1, surv_years>0) %>%
  mutate(PNI_s=as.numeric(scale(PNI)),
         CONUT_s=as.numeric(scale(-CONUT)),
         GNRI_s=as.numeric(scale(GNRI)))

n_gi <- nrow(df_gi)
n_events <- sum(df_gi$death, na.rm=TRUE)
cat(sprintf("  N=%d, events=%d\n", n_gi, n_events))

if (n_events >= 10) {
  cox_res <- data.frame()
  for (exp in c("PNI_s", "CONUT_s", "GNRI_s")) {
    f <- as.formula(paste("Surv(surv_years, death) ~", exp, "+ age + sex + race_eth"))
    fit <- coxph(f, data=df_gi)
    hr <- tidy(fit, conf.int=TRUE) %>% filter(term==exp)
    cox_res <- rbind(cox_res, data.frame(
      Nutrition=gsub("_s","",exp),
      HR=exp(hr$estimate), Lower=exp(hr$conf.low),
      Upper=exp(hr$conf.high), P=hr$p.value,
      N=n_gi, Events=n_events))
  }
  cox_res$HR_CI <- sprintf("%.3f (%.3f-%.3f)", cox_res$HR, cox_res$Lower, cox_res$Upper)
  print(cox_res[, c("Nutrition","HR_CI","P")])
  write.csv(cox_res, file.path(GI_RESULTS_DIR, "cox_allcause_combined.csv"), row.names=FALSE)

  # KM by PNI tertile
  df_gi$pni_t <- factor(ntile(df_gi$PNI, 3), 1:3, c("Low","Mid","High"))
  km <- survfit(Surv(surv_years, death) ~ pni_t, data=df_gi)
  cat("\nMedian survival by PNI tertile:\n")
  print(km)
}

cat("\n── Combined: Competing risks ──\n")
n_gc <- sum(df_gi$gi_cancer_death, na.rm=TRUE)
cat(sprintf("  GI cancer deaths: %d\n", n_gc))

cat(sprintf("\nAll results in: %s\n", GI_RESULTS_DIR))
cat("Done.\n")
