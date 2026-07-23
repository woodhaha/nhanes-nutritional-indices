# run_nhanes3_fast.R — NHANES III fast extraction (only needed vars)
# Rscript run_nhanes3_fast.R

library(dplyr); library(survival); library(broom)

DATA_DIR <- "D:/Researching/NHANES_aged_GI_tumor_nutrition/data/nhanes3"
GI_DATA_DIR <- "D:/Researching/data/gi_analysis"
GI_RESULTS_DIR <- "D:/Researching/results/gi_analysis"
dir.create(GI_RESULTS_DIR, recursive=TRUE, showWarnings=FALSE)

# ── Parse SAS INPUT to get variable positions ────────────────────────────────
get_var_pos <- function(sas_file, target_vars) {
  sas <- readLines(sas_file, warn=FALSE)
  idx <- grep("^\\s+INPUT", sas, ignore.case=TRUE)[1]
  lines <- sas[idx:length(sas)]
  end <- grep("^\\s+LABEL|^\\s*;|^$", lines)[1]
  lines <- lines[2:(end-1)]

  result <- list()
  for (ln in lines) {
    ln <- gsub("\\s+", " ", trimws(ln))
    parts <- strsplit(ln, " ")[[1]]
    if (length(parts) < 2) next
    name <- parts[1]
    if (!(name %in% target_vars)) next
    pos <- gsub("\\$", "", parts[2])
    if (grepl("-", pos)) {
      r <- as.integer(strsplit(pos, "-")[[1]])
      result[[name]] <- list(start=r[1], end=r[2])
    } else {
      p <- suppressWarnings(as.integer(pos))
      if (!is.na(p)) result[[name]] <- list(start=p, end=p)
    }
  }
  result
}

# ── Fast read with only needed columns ───────────────────────────────────────
fast_read <- function(dat_path, var_positions) {
  # Sort by start position
  sorted <- var_positions[order(sapply(var_positions, `[[`, "start"))]

  # Calculate widths including gaps
  positions <- list()
  prev_end <- 0
  for (nm in names(sorted)) {
    s <- sorted[[nm]]$start
    e <- sorted[[nm]]$end
    if (s > prev_end + 1) {
      # GAP - add dummy
      positions[[length(positions)+1]] <- list(name=NULL, width=s-prev_end-1, skip=TRUE)
    }
    positions[[length(positions)+1]] <- list(name=nm, width=e-s+1, skip=FALSE)
    prev_end <- e
  }

  widths <- sapply(positions, `[[`, "width")
  is_skip <- sapply(positions, `[[`, "skip")
  col_names <- sapply(positions, function(p) if(is.null(p$name)) NA else p$name)

  # Read with only needed columns
  cat(sprintf("  Reading %s...", basename(dat_path)))
  df <- read.fwf(dat_path, widths=widths,
                 col.names=paste0("V", seq_along(widths)),
                 na.strings=c(".",""," ","  "),
                 buffersize=100000)
  cat(sprintf(" %d rows\n", nrow(df)))

  # Extract only non-skip columns
  for (i in seq_along(positions)) {
    if (!is_skip[i]) {
      nm <- col_names[i]
      # For sequential processing
    }
  }

  # Build result
  result <- data.frame(SEQN=df$V1)
  for (i in seq_along(positions)) {
    if (!is_skip[i]) {
      nm <- col_names[i]
      vals <- suppressWarnings(as.numeric(as.character(df[[i]])))
      result[[nm]] <- vals
    }
  }
  result
}

# ── Adult variables ──────────────────────────────────────────────────────────
cat("=== Adult file ===\n")
adult_vars <- c("SEQN", "HSSEX", "HSAGEIR", "DMAETHNR", "DMARACER",
                "DMPPIR", "HAC1N", "HAC1O", "HAC3OS",
                "HFHEDUCR", "WTPFEX6", "SDPPSU6", "SDPSTRA6",
                "HFF1", "HFF20R", "HFF19R")
adult_pos <- get_var_pos("data/nhanes3/adult.sas", adult_vars)
cat("Positions found:", length(adult_pos), "\n")
adult <- fast_read("data/nhanes3/adult.dat", adult_pos)
cat(names(adult), "\n")

# ── Lab variables ────────────────────────────────────────────────────────────
cat("\n=== Lab file ===\n")
lab_vars <- c("SEQN", "AMP", "AMPSI", "TCP", "TCPSI",
              "WCP", "WCPSI", "LMPPCNT", "LMP",
              "CRP", "TPP", "GBP", "CHP", "CHPSI",
              "HGP", "HGPSI", "MCPSI", "MVPSI")
lab_pos <- get_var_pos("data/nhanes3/lab.sas", lab_vars)
cat("Positions found:", length(lab_pos), "\n")
lab <- fast_read("data/nhanes3/lab.dat", lab_pos)
cat(names(lab), "\n")

# ── Exam variables ───────────────────────────────────────────────────────────
cat("\n=== Exam file ===\n")
exam_vars <- c("SEQN", "BMPBMI", "BMPWT", "BMPHT")
exam_pos <- get_var_pos("data/nhanes3/exam.sas", exam_vars)
cat("Positions found:", length(exam_pos), "\n")
exam <- fast_read("data/nhanes3/exam.dat", exam_pos)
cat(names(exam), "\n")

# ── Merge ────────────────────────────────────────────────────────────────────
cat("\n=== Merging ===\n")
df <- adult
df <- merge(df, lab, by="SEQN", all.x=TRUE, suffixes=c("",".lab"))
df <- merge(df, exam, by="SEQN", all.x=TRUE, suffixes=c("",".ex"))

# Age filter
df <- df %>% filter(HSAGEIR >= 60)
cat(sprintf("Age >= 60: %d rows\n", nrow(df)))

# Cancer classification
GI_CODES_N3 <- c(6, 12, 13, 14, 15, 25)  # colon, stomach, esophagus, pancreas, liver, gallbladder
df <- df %>% mutate(
  any_cancer = ifelse((!is.na(HAC1N) & HAC1N==1) | (!is.na(HAC1O) & HAC1O==1), 1L, 0L),
  gi_tumor = ifelse(!is.na(HAC3OS) & HAC3OS %in% GI_CODES_N3, 1L, 0L),
  # Demographics
  age = HSAGEIR,
  sex = ifelse(HSSEX==1, "Male", "Female"),
  race_eth = case_when(
    DMAETHNR %in% c(1,2) ~ "Hispanic",
    DMARACER == 1 ~ "Non-Hispanic White",
    DMARACER == 2 ~ "Non-Hispanic Black",
    TRUE ~ "Other"),
  edu_binary = ifelse(HFHEDUCR >= 4, 1, 0),
  income_pir = DMPPIR,
  # Nutrition
  albumin_gdl = AMP,
  lymph_abs = LMP,
  tchol_mgdl = ifelse(is.na(TCP), CHP, TCP),
  bmi = BMPBMI,
  # Indices
  PNI = 10*AMP + 0.005*LMP,
  CONUT = case_when(AMP>=3.5~0, AMP>=3.0~2, AMP>=2.5~4, TRUE~6) +
          case_when(LMP>=1600~0, LMP>=1200~1, LMP>=800~2, TRUE~3) +
          case_when(tchol_mgdl>=180~0, tchol_mgdl>=140~1, tchol_mgdl>=100~2, TRUE~3),
  GNRI = 14.89*AMP + 41.7*(bmi/22),
  wt_mec_2yr = WTPFEX6,
  psu = SDPPSU6,
  strata = SDPSTRA6,
  cycle = "NHANES_III_1988-1994"
)

cat(sprintf("GI tumors: %d\n", sum(df$gi_tumor, na.rm=TRUE)))
gi_df <- df[df$gi_tumor==1 & !is.na(df$gi_tumor),]
if (nrow(gi_df) > 0) {
  cat("By site:\n")
  print(table(HAC3OS=gi_df$HAC3OS, useNA="ifany"))
}

# Complete nutrition
df <- df %>% filter(!is.na(PNI), !is.na(CONUT), !is.na(GNRI))
cat(sprintf("Complete nutrition: %d (GI tumors: %d)\n",
            nrow(df), sum(df$gi_tumor, na.rm=TRUE)))

# ── Mortality ────────────────────────────────────────────────────────────────
cat("\n=== Mortality ===\n")
mort_path <- file.path(GI_DATA_DIR, "mort_nhanes3.rds")
if (!file.exists(mort_path)) {
  tmp <- tempfile(fileext=".dat")
  download.file(
    "https://ftp.cdc.gov/pub/health_statistics/nchs/datalinkage/linked_mortality/NHANES_III_MORT_2019_PUBLIC.dat",
    tmp, mode="wb", quiet=TRUE)
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
} else {
  mort <- readRDS(mort_path)
}
cat(sprintf("Mortality: %d records\n", nrow(mort)))

df <- merge(df, mort[, c("SEQN", "eligstat", "mortstat", "ucod_leading", "permth_int")],
            by="SEQN", all.x=TRUE)
df <- df[df$eligstat==1 & !is.na(df$mortstat), , drop=FALSE]
df$surv_years <- df$permth_int/12
df$death <- df$mortstat
df$gi_cancer_death <- ifelse(df$death==1 & df$ucod_leading==2 & df$gi_tumor==1, 1L, 0L)

cat(sprintf("Follow-up eligible: %d\n", nrow(df)))
cat(sprintf("  Deaths: %d (%.1f%%)\n", sum(df$death), mean(df$death)*100))
cat(sprintf("  GI cancer deaths: %d\n", sum(df$gi_cancer_death, na.rm=TRUE)))

# Save
saveRDS(df, file.path(GI_DATA_DIR, "nhanes3_gi_merged.rds"))
cat("Saved.\n")

# ── Combine with NHANES 2005-2016 ───────────────────────────────────────────
cat("\n=== Combining ===\n")
df_old <- readRDS(file.path(GI_DATA_DIR, "nhanes_gi_nutrition_raw.rds"))
mort_old <- readRDS(file.path(GI_DATA_DIR, "mort_2019.rds"))
df_old <- merge(df_old, mort_old[,c("SEQN","eligstat","mortstat","ucod_leading","permth_int")],
                by="SEQN", all.x=TRUE)
df_old <- df_old[df_old$eligstat==1 & !is.na(df_old$mortstat),]
df_old$surv_years <- df_old$permth_int/12
df_old$death <- df_old$mortstat
df_old$gi_cancer_death <- ifelse(df_old$death==1 & df_old$ucod_leading==2 & df_old$gi_tumor==1, 1L, 0L)
df_old$cycle <- paste0("NHANES_", df_old$cycle)

# Common columns
common <- intersect(names(df_old), names(df))
df3 <- df[, common]
df_old <- df_old[, common]
df_all <- bind_rows(df_old, df3)

cat(sprintf("COMBINED: N=%d, GI tumors=%d\n",
            nrow(df_all), sum(df_all$gi_tumor, na.rm=TRUE)))
cat(sprintf("  NHANES III: %d GI tumors\n", sum(df$gi_tumor, na.rm=TRUE)))
cat(sprintf("  2005-2016:  %d GI tumors\n", sum(df_old$gi_tumor, na.rm=TRUE)))

# ── Survival analysis ────────────────────────────────────────────────────────
cat("\n── Combined analysis ──\n")
df_gi <- df_all %>% filter(gi_tumor==1, surv_years>0)
df_gi$sex_b <- ifelse(df_gi$sex=="Female", 1L, 0L)
df_gi$race_eth_b <- ifelse(df_gi$race_eth=="Non-Hispanic White", 1L, 0L)
df_gi$PNI_s <- as.numeric(scale(df_gi$PNI))
df_gi$CONUT_s <- as.numeric(scale(-df_gi$CONUT))
df_gi$GNRI_s <- as.numeric(scale(df_gi$GNRI))

n_gi <- nrow(df_gi)
n_ev <- sum(df_gi$death, na.rm=TRUE)
cat(sprintf("  N=%d, events=%d\n", n_gi, n_ev))

if (n_ev >= 10) {
  cox_res <- data.frame()
  for (exp in c("PNI_s", "CONUT_s", "GNRI_s")) {
    for (adj in c("Crude", "Adjusted")) {
      covs <- if (adj=="Crude") "1" else "age + sex_b + race_eth_b"
      f <- as.formula(paste("Surv(surv_years, death) ~", exp, "+", covs))
      fit <- coxph(f, data=df_gi)
      hr <- tidy(fit, conf.int=TRUE) %>% filter(term==exp)
      cox_res <- rbind(cox_res, data.frame(
        Nutrition=gsub("_s","",exp), Adjustment=adj,
        HR=exp(hr$estimate), Lower=exp(hr$conf.low),
        Upper=exp(hr$conf.high), P=hr$p.value, N=n_gi, Events=n_ev))
    }
  }
  cox_res$HR_CI <- sprintf("%.3f (%.3f-%.3f)", cox_res$HR, cox_res$Lower, cox_res$Upper)
  cat("\nAll-cause Cox:\n")
  print(cox_res[, c("Nutrition","Adjustment","HR_CI","P")])
  write.csv(cox_res, file.path(GI_RESULTS_DIR, "cox_allcause_combined.csv"), row.names=FALSE)

  # Median survival
  df_gi$pni_t <- factor(ntile(df_gi$PNI, 3), 1:3, c("Low","Mid","High"))
  km <- survfit(Surv(surv_years, death) ~ pni_t, data=df_gi)
  cat("\nMedian survival:\n")
  print(km)
  saveRDS(km, file.path(GI_RESULTS_DIR, "km_fit_combined.rds"))
}

cat("\nDone.\n")
