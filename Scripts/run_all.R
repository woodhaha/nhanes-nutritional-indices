# ─── run_all.R — Full GI tumor × nutrition pipeline ────────────────────────
# Rscript run_all.R

cat("========================================\n")
cat("NHANES GI Tumor Analysis Pipeline\n")
cat("========================================\n\n")

for (step in c("01_prepare.R", "02_analyze.R", "03_advanced.R",
                "nlr_sii_extract.R", "hei_merge.R")) {
  if (!file.exists(step)) { cat(sprintf("Skipping %s (not found)\n", step)); next }
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), step))
  system(paste("Rscript", step))
}

cat("\nDone. Results in results/gi_analysis/\n")
cat("Figures in figures/gi_analysis/\n")
