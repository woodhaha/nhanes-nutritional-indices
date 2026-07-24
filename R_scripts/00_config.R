# ==============================================================================
# 00_config.R — NHANES Nutrition + Cognition + Depression Study
# Configuration, paths, and package loading
# Study: Multi-index nutritional status → cognitive function ↔ depression
# Population: US adults ≥60, NHANES 2011-2014
# ==============================================================================

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJ_ROOT  <- here::here()
DATA_DIR   <- file.path(PROJ_ROOT, "data")
RESULTS_DIR <- file.path(PROJ_ROOT, "results")
FIG_DIR    <- file.path(PROJ_ROOT, "figures")

# ── Package management ────────────────────────────────────────────────────────
required_packages <- c(
  "survey", "dplyr", "tidyr", "tibble", "purrr", "stringr", "forcats",
  "ggplot2", "cowplot", "ggforestplot", "gtsummary", "broom",
  "marginaleffects", "splines", "rms", "Hmisc",
  "lavaan", "mice", "VIM", "naniar",
  "nhanesA", "here", "data.table", "patchwork"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    try(install.packages(pkg, repos = "https://cloud.r-project.org"), silent = TRUE)
  }
  try(library(pkg, character.only = TRUE), silent = TRUE)
}

# ── NHANES cycle selection ────────────────────────────────────────────────────
# Cognitive function data (CERAD/AFT/DSST) ONLY available in 2011-2014
CYCLES <- c("2011-2012", "2013-2014")

# ── Key variable lists (for nhanesA selective download) ───────────────────────
# These are the exact NHANES variable names

DEMO_VARS <- c("SEQN", "SDDSRVYR", "RIDAGEYR", "RIAGENDR", "RIDRETH1",
               "DMDEDUC2", "DMDMARTL", "INDFMPIR", "SDMVPSU", "SDMVSTRA",
               "WTMEC2YR", "WTINT2YR")

DIET_VARS <- c("SEQN", "DR1TKCAL", "DR1TPROT", "DR1TCARB", "DR1TTFAT",
               "DR1TSFAT", "DR1TMFAT", "DR1TPFAT", "DR1TFIBE", "DR1TCHOL",
               "DR1TVB1", "DR1TVB2", "DR1TNIAC", "DR1TVB6", "DR1TFOLA",
               "DR1TVB12", "DR1TVC", "DR1TVD", "DR1TVE", "DR1TVK",
               "DR1TCALC", "DR1TIRON", "DR1TMAGN", "DR1TZINC", "DR1TSELE",
               "DR1TSODI", "DR1TPOTA", "DR1TCAFF", "DR1TALCO",
               "DR1TN3", "DR1TN6", "DR1TVARA", "DR1TBETA", "DR1TLYCO",
               "DR1TSUGR", "DR1TFA")

PHQ9_VARS <- c("SEQN", paste0("DPQ0", sprintf("%02d", 10:19)))

CFQ_VARS <- c("SEQN", "CFD_WL_IMM", "CFD_WL_DEL", "CFD_AFT", "CFD_DSST")

BIOPRO_VARS <- c("SEQN", "LBDSALSI", "LBDSCRPSI")

CBC_VARS <- c("SEQN", "LBXWBCSI", "LBXLYPCT", "LBXNEPCT", "LBXRBCSI")

TCHOL_VARS <- c("SEQN", "LBXTC")

BMX_VARS <- c("SEQN", "BMXBMI", "BMXWAIST", "BMXWT", "BMXHT")

BPX_VARS <- c("SEQN", "BPXSY1", "BPXSY2", "BPXSY3", "BPXDI1", "BPXDI2", "BPXDI3")

SMQ_VARS <- c("SEQN", "SMQ020", "SMQ040")

ALQ_VARS <- c("SEQN", "ALQ111", "ALQ121", "ALQ130")

DIQ_VARS <- c("SEQN", "DIQ010")

MCQ_VARS <- c("SEQN", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F",
              "MCQ010", "MCQ220")

FSQ_VARS <- c("SEQN", "FSDHH", "FSDAD", "FSD151")

PAQ_VARS <- c("SEQN", "PAQ605", "PAQ610", "PAD615", "PAQ620", "PAQ625",
              "PAD630", "PAQ635", "PAQ640", "PAD645", "PAQ650", "PAQ655",
              "PAD660", "PAQ665", "PAQ670", "PAD675")

BIOPRO_F_VARS <- c("SEQN", "LBDSCRPSI")

# ── ggplot2 theme ──────────────────────────────────────────────────────────────
theme_nhanes <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      plot.title    = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(color = "grey40"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = base_size - 1),
      legend.position = "bottom",
      axis.title   = element_text(size = base_size),
      plot.caption = element_text(color = "grey50", size = base_size - 3,
                                  hjust = 0, face = "italic")
    )
}

cat("[00_config.R] Configuration loaded.\n")
cat(sprintf("  Project: %s\n", PROJ_ROOT))
cat(sprintf("  NHANES cycles: %s\n", paste(CYCLES, collapse = ", ")))
