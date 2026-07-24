# ==============================================================================
# 01a_load_local_xpt.R — Load NHANES XPT files from local cache, then derive
# Bypasses nhanesA timeouts by reading XPT files downloaded via curl
# ==============================================================================
library(haven)
library(dplyr)
library(tidyr)
library(stringr)

source(here::here("R_scripts", "00_config.R"))

XPT_DIR <- file.path(DATA_DIR, "nhanes_xpt")

#' Load a single NHANES XPT table from local cache
load_xpt <- function(tbl_name, cycle_letter) {
  file_path <- file.path(XPT_DIR, sprintf("%s_%s.XPT", tbl_name, cycle_letter))
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(NULL)
  }
  haven::read_xpt(file_path)
}

load_nhanes_local <- function(cycle) {
  cycle_letter <- switch(cycle,
    "2011-2012" = "G",
    "2013-2014" = "H"
  )

  cat(sprintf("[%s] Loading from local XPT cache...\n", cycle))

  tables <- list(
    DEMO   = load_xpt("DEMO", cycle_letter),
    DR1TOT = load_xpt("DR1TOT", cycle_letter),
    DPQ    = load_xpt("DPQ", cycle_letter),
    CFQ    = load_xpt("CFQ", cycle_letter),
    BIOPRO = load_xpt("BIOPRO", cycle_letter),
    CBC    = load_xpt("CBC", cycle_letter),
    TCHOL  = load_xpt("TCHOL", cycle_letter),
    BMX    = load_xpt("BMX", cycle_letter),
    BPX    = load_xpt("BPX", cycle_letter),
    SMQ    = load_xpt("SMQ", cycle_letter),
    ALQ    = load_xpt("ALQ", cycle_letter),
    DIQ    = load_xpt("DIQ", cycle_letter),
    MCQ    = load_xpt("MCQ", cycle_letter),
    FSQ    = load_xpt("FSQ", cycle_letter),
    PAQ    = load_xpt("PAQ", cycle_letter)
  )

  assign(paste0("raw_", cycle_letter), tables, envir = .GlobalEnv)
  cat(sprintf("[%s] All tables loaded (%d tables)\n", cycle, length(tables)))
  TRUE
}

# ── Execute ──────────────────────────────────────────────────────────────────
cat("── NHANES Local XPT Loading ──\n")

# Load both cycles
for (cyc in CYCLES) {
  load_nhanes_local(cyc)
}

# Now source the derivation script
source(here::here("R_scripts", "01_load_and_derive.R"))
