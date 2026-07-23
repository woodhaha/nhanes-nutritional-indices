# rebuild_clean_data.R — Fix nested data.frame columns and rebuild combined
# Rscript rebuild_clean_data.R
# Then: Rscript run_manuscript_tables.R

library(dplyr)

PROJ <- "D:/Researching/NHANES_aged_GI_tumor_nutrition"
DATA_DIR <- file.path(PROJ, "data/gi_analysis")
RES_DIR <- file.path(PROJ, "results/gi_analysis")

# Each column is an embedded data.frame; extract the matching inner column
unnest_df <- function(df) {
  cols <- lapply(names(df), function(nm) {
    v <- df[[nm]]
    if (is.data.frame(v)) {
      if (nm %in% names(v)) x <- v[[nm]]
      else x <- v[[1]]
    } else if (is.list(v)) {
      x <- unlist(v)
    } else {
      x <- v
    }
    if (is.character(x)) return(x)
    if (is.factor(x)) return(as.character(x))
    as.numeric(x)
  })
  names(cols) <- names(df)
  as_tibble(cols)
}

cat("=== Rebuilding clean combined data ===\n")

# ── NHANES III ──────────────────────────────────────────────────────────────
cat("NHANES III (all ages)...\n")
df3 <- readRDS(file.path(DATA_DIR, "nhanes3_gi_allages.rds"))
cat(sprintf("  Raw: %d x %d\n", nrow(df3), ncol(df3)))

# Unnest
df3 <- unnest_df(df3)
cat(sprintf("  Clean: %d x %d\n", nrow(df3), ncol(df3)))
cat(sprintf("  PNI range: %.1f-%.1f\n", min(df3$PNI, na.rm=TRUE), max(df3$PNI, na.rm=TRUE)))
cat(sprintf("  GI tumors: %d\n", sum(df3$gi_tumor==1, na.rm=TRUE)))

# ── NHANES 2005-2016 ────────────────────────────────────────────────────────
cat("\nNHANES 2005-2016...\n")
df5 <- readRDS(file.path(DATA_DIR, "nhanes_gi_nutrition_raw.rds"))
df5 <- unnest_df(df5)
cat(sprintf("  Clean: %d x %d\n", nrow(df5), ncol(df5)))

# BUGFIX: albumin in 2005-2016 NHANES is stored in g/L (SI units),
# but PNI formula requires g/dL. Convert: g/L / 10 = g/dL
# This also fixes PNI/CONUT/GNRI which were computed on g/L values
cat("\n  Fixing albumin units (g/L → g/dL) and recalculating indices...\n")
df5$albumin_gdl <- df5$albumin_gdl / 10
df5$PNI <- 10 * df5$albumin_gdl + 0.005 * df5$lymph_abs
df5$GNRI <- 14.89 * df5$albumin_gdl + 41.7 * (df5$bmi / 22)
df5$CONUT <- case_when(df5$albumin_gdl>=3.5~0, df5$albumin_gdl>=3.0~2,
                       df5$albumin_gdl>=2.5~4, TRUE~6) +
             case_when(df5$lymph_abs>=1600~0, df5$lymph_abs>=1200~1,
                       df5$lymph_abs>=800~2, TRUE~3) +
             case_when(df5$tchol_mgdl>=180~0, df5$tchol_mgdl>=140~1,
                       df5$tchol_mgdl>=100~2, TRUE~3)
cat(sprintf("  Fixed PNI range: %.1f-%.1f\n", min(df5$PNI, na.rm=TRUE), max(df5$PNI, na.rm=TRUE)))
cat(sprintf("  Fixed PNI mean: %.1f\n", mean(df5$PNI, na.rm=TRUE)))

# Load mortality for 2005-2016
mort5 <- readRDS(file.path(DATA_DIR, "mort_2019.rds"))
mort5 <- unnest_df(mort5)
df5 <- merge(df5, mort5[,c("SEQN","eligstat","mortstat","ucod_leading","permth_int")],
             by="SEQN", all.x=TRUE)

df5 <- df5[!is.na(df5$eligstat) & df5$eligstat==1 & !is.na(df5$mortstat), , drop=FALSE]
df5$surv_years <- df5$permth_int / 12
df5$death <- df5$mortstat
df5 <- df5[!is.na(df5$PNI) & !is.na(df5$CONUT) & !is.na(df5$GNRI), , drop=FALSE]
df5$gi_cancer_death <- ifelse(df5$death==1 & df5$ucod_leading==2 & df5$gi_tumor==1, 1L, 0L)
df5$other_death <- ifelse(df5$death==1 & (df5$gi_cancer_death!=1 | is.na(df5$gi_cancer_death)), 1L, 0L)
df5$cycle <- paste0("NHANES_", df5$cycle)
cat(sprintf("  Final: %d (GI: %d)\n", nrow(df5), sum(df5$gi_tumor==1, na.rm=TRUE)))

# ── NHANES III mortality filter ─────────────────────────────────────────────
cat("\nNHANES III filtering...\n")
df3 <- df3[!is.na(df3$eligstat) & df3$eligstat==1 & !is.na(df3$mortstat), , drop=FALSE]
df3 <- df3[!is.na(df3$PNI) & !is.na(df3$CONUT) & !is.na(df3$GNRI), , drop=FALSE]
df3$gi_cancer_death <- ifelse(df3$death==1 & df3$ucod_leading==2 & df3$gi_tumor==1, 1L, 0L)
df3$other_death <- ifelse(df3$death==1 & (df3$gi_cancer_death!=1 | is.na(df3$gi_cancer_death)), 1L, 0L)
cat(sprintf("  Final: %d (GI: %d)\n", nrow(df3), sum(df3$gi_tumor==1, na.rm=TRUE)))

# ── Combine ─────────────────────────────────────────────────────────────────
cat("\nCombining...\n")
common <- intersect(names(df5), names(df3))
cat(sprintf("  Common columns: %d\n", length(common)))

df_all <- rbind(df5[, common, drop=FALSE], df3[, common, drop=FALSE])
cat(sprintf("  Combined: N=%d, GI=%d\n", nrow(df_all), sum(df_all$gi_tumor==1, na.rm=TRUE)))

# Key numeric conversions
for (v in c("age","PNI","CONUT","GNRI","albumin_gdl","lymph_abs","tchol_mgdl","bmi","surv_years")) {
  if (v %in% names(df_all)) df_all[[v]] <- as.numeric(df_all[[v]])
}

# Plausibility filters (physiologically impossible values)
df_all <- df_all[which(df_all$albumin_gdl >= 1.5 & df_all$albumin_gdl <= 6.0), , drop=FALSE]
df_all <- df_all[which(df_all$lymph_abs >= 200 & df_all$lymph_abs <= 15000), , drop=FALSE]
df_all <- df_all[which(df_all$tchol_mgdl >= 50 & df_all$tchol_mgdl <= 500), , drop=FALSE]
df_all <- df_all[which(df_all$bmi >= 12 & df_all$bmi <= 80), , drop=FALSE]
df_all <- df_all[which(df_all$age >= 18 & df_all$age <= 90), , drop=FALSE]
df_all <- df_all[which(df_all$income_pir < 888888), , drop=FALSE]
cat(sprintf("  After plausibility filter: N=%d, GI=%d\n", nrow(df_all), sum(df_all$gi_tumor, na.rm=TRUE)))

# Add gi_status for plotting
df_all$gi_status <- case_when(
  df_all$gi_tumor==1 ~ "GI Tumor",
  df_all$any_cancer==0 ~ "Non-Cancer",
  TRUE ~ "Other Cancer"
)

saveRDS(df_all, file.path(RES_DIR, "nhanes_combined_gi_fixed.rds"))
cat("  Saved: nhanes_combined_gi_fixed.rds\n")

# Quick verification
gi <- df_all[df_all$gi_tumor==1 & df_all$surv_years>0, ]
cat(sprintf("\nVerification: GI N=%d, deaths=%d, PNI=%.1f±%.1f\n",
    nrow(gi), sum(gi$death, na.rm=TRUE),
    mean(gi$PNI, na.rm=TRUE), sd(gi$PNI, na.rm=TRUE)))

cat("\nDone.\n")
