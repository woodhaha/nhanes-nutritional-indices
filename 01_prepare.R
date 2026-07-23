# ─── 01_prepare.R — Wrangle existing data to clean format ────────────────
# Rscript 01_prepare.R
# Loads existing cleaned RDS → unnests → saves as tibble

library(dplyr)

PROJ <- normalizePath(".")
DATA_DIR <- file.path(PROJ, "data/gi_analysis")
RES_DIR <- file.path(PROJ, "results/gi_analysis")

# Load the previously cleaned combined data (with fixed albumin, all ages)
df <- readRDS(file.path(RES_DIR, "nhanes_combined_gi_fixed.rds"))

# Unnest any remaining embedded columns (R 4.6.0 workaround)
unnest_col <- function(x) {
  if (is.data.frame(x) && ncol(x) > 1) as_tibble(x)[[1]]  # embedded df → first col
  else if (is.data.frame(x) && ncol(x) == 1) x[[1]]
  else x
}
df <- as_tibble(lapply(df, unnest_col))

# Merge CRP from raw 2005-2016 data (stored as LBDSCRSI mg/L → rename)
df5_raw <- tryCatch(readRDS(file.path(DATA_DIR, "nhanes_gi_nutrition_raw.rds")), error=\(e)NULL)
if (!is.null(df5_raw)) {
  crp <- df5_raw[["crp_mgdl"]]
  if (is.data.frame(crp)) crp <- crp[["crp_mgdl"]]
  crp_df <- data.frame(SEQN = as.numeric(df5_raw[["SEQN"]]), crp_mgL = as.numeric(crp))
  df <- merge(df, crp_df, by="SEQN", all.x=TRUE)
  rm(crp, crp_df, df5_raw)
}
if (!"crp_mgL" %in% names(df)) df$crp_mgL <- NA_real_

# Recompute gi_status
df <- df %>% mutate(
  gi_status = case_when(gi_tumor==1 ~ "GI Tumor",
                        any_cancer==0 ~ "Non-Cancer",
                        TRUE ~ "Other Cancer"),
  gi_cancer_death = ifelse(is.na(gi_cancer_death), 0L, gi_cancer_death),
  other_death = ifelse(is.na(other_death), 0L, other_death)
)

saveRDS(df, file.path(RES_DIR, "nhanes_clean.rds"))
cat(sprintf("Saved: N=%d, GI=%d, deaths=%d\n",
    nrow(df), sum(df$gi_tumor, na.rm=TRUE), sum(df$death, na.rm=TRUE)))
