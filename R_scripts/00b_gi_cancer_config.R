# =============================================================================
# 00b_gi_cancer_config.R — GI Tumor + Mortality Linkage Configuration
# Adds GI cancer site mapping and NCHS mortality linkage to NHANES
# =============================================================================

# ── Additional packages ──────────────────────────────────────────────────────
gi_required <- c("survival", "cmprsk", "nhanesdata")
for (pkg in gi_required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# ── GI Anatomic Sites (ICD-10) ───────────────────────────────────────────────
# Digestive system malignant neoplasms C15-C25
GI_SITES <- list(
  esophagus        = c(icd = "C15",  label = "Esophagus"),
  stomach          = c(icd = "C16",  label = "Stomach"),
  small_intestine  = c(icd = "C17",  label = "Small Intestine"),
  colon            = c(icd = "C18",  label = "Colon"),
  rectosigmoid     = c(icd = "C19",  label = "Rectosigmoid"),
  rectum           = c(icd = "C20",  label = "Rectum"),
  anus             = c(icd = "C21",  label = "Anus/Anal Canal"),
  liver            = c(icd = "C22",  label = "Liver/Intrahepatic Bile Duct"),
  gallbladder      = c(icd = "C23",  label = "Gallbladder"),
  biliary_other    = c(icd = "C24",  label = "Other Biliary"),
  pancreas         = c(icd = "C25",  label = "Pancreas")
)
GI_ICD_CODES <- paste0("C", 15:25)  # C15-C25 inclusive

# ── NHANES cancer site variable mapping (MCQ230/MCQ240) ─────────────────────
# Different NHANES cycles use MCQ230 or MCQ240 for "what kind of cancer"
# Values are NHANES codebook response codes for GI sites
CANCER_SITE_CYCLE_MAP <- list(
  "2005-2006" = list(
    var_name = "MCQ230",
    gi_codes = c(7, 8, 9, 13, 22, 24, 25),     # colon, esoph, gb, liver, pancr, rectum, stomach
    code_site = c("7" = "Colon", "8" = "Esophagus", "9" = "Gallbladder",
                  "13" = "Liver", "22" = "Pancreas",
                  "24" = "Rectum", "25" = "Stomach"),
    table_name = "MCQ_D"   # 2005-2006 table suffix
  ),
  "2007-2008" = list(
    var_name = "MCQ240",
    gi_codes = c(7, 8, 9, 13, 22, 24, 25),
    code_site = c("7" = "Colon", "8" = "Esophagus", "9" = "Gallbladder",
                  "13" = "Liver", "22" = "Pancreas",
                  "24" = "Rectum", "25" = "Stomach"),
    table_name = "MCQ_E"
  ),
  "2009-2010" = list(
    var_name = "MCQ240",
    gi_codes = c(7, 8, 9, 13, 22, 24, 25),
    code_site = c("7" = "Colon", "8" = "Esophagus", "9" = "Gallbladder",
                  "13" = "Liver", "22" = "Pancreas",
                  "24" = "Rectum", "25" = "Stomach"),
    table_name = "MCQ_F"
  ),
  "2011-2012" = list(
    var_name = "MCQ240",
    gi_codes = c(7, 8, 9, 13, 22, 24, 25),
    code_site = c("7" = "Colon", "8" = "Esophagus", "9" = "Gallbladder",
                  "13" = "Liver", "22" = "Pancreas",
                  "24" = "Rectum", "25" = "Stomach"),
    table_name = "MCQ_G"
  ),
  "2013-2014" = list(
    var_name = "MCQ240",
    gi_codes = c(7, 8, 9, 13, 22, 24, 25),
    code_site = c("7" = "Colon", "8" = "Esophagus", "9" = "Gallbladder",
                  "13" = "Liver", "22" = "Pancreas",
                  "24" = "Rectum", "25" = "Stomach"),
    table_name = "MCQ_H"
  ),
  "2015-2016" = list(
    var_name = "MCQ240",
    gi_codes = c(7, 8, 9, 13, 22, 24, 25),
    code_site = c("7" = "Colon", "8" = "Esophagus", "9" = "Gallbladder",
                  "13" = "Liver", "22" = "Pancreas",
                  "24" = "Rectum", "25" = "Stomach"),
    table_name = "MCQ_I"
  ),
  "2017-2018" = list(
    var_name = "MCQ240",
    gi_codes = c(7, 8, 9, 13, 22, 24, 25),
    code_site = c("7" = "Colon", "8" = "Esophagus", "9" = "Gallbladder",
                  "13" = "Liver", "22" = "Pancreas",
                  "24" = "Rectum", "25" = "Stomach"),
    table_name = "MCQ_J"
  )
)

# ── Mortality linkage ───────────────────────────────────────────────────────
# NCHS 2019 Public-Use Linked Mortality File
# Variables: mortstat (0/1), permth_int (months), ucod_leading (COD recode)
#
# LIMITATION: Public-use UCOD_LEADING only has 10 broad categories.
#   code 002 = "Malignant neoplasms (C00-C97)" — ALL cancers combined.
#   Site-specific ICD-10 codes (e.g. C15 for esophageal cancer) and the
#   detailed UCOD_113 recode (113 causes, including codes 019-023 for GI
#   sites) are ONLY available in NCHS restricted-use RDC files.
#   Reference: https://www.cdc.gov/nchs/data-linkage/mortality-public.htm
#
# FIX: In survival analysis, participants with a GI tumor at baseline who
# die with UCOD_LEADING=002 are classified as "probable GI cancer death."
# This proxy has high specificity for high-mortality GI cancers (pancreas,
# esophagus, liver, stomach) where 5-year survival is <35%.
# For colorectal cancer (5-yr ~65%), specificity is lower but still
# plausible since CRC is the most common GI cancer and recurrent/metastatic
# disease drives most cancer deaths among CRC patients.
CANCER_COD_CODE <- 2  # ucod_leading value for malignant neoplasms

# ── Expand NHANES cycles for GI analysis ─────────────────────────────────────
GI_CYCLES <- c("2005-2006", "2007-2008", "2009-2010", "2011-2012",
               "2013-2014", "2015-2016", "2017-2018")

# ── Extended demo variables (include cancer site vars) ───────────────────────
GI_DEMO_VARS <- c("SEQN", "SDDSRVYR", "RIDAGEYR", "RIAGENDR", "RIDRETH1",
                  "DMDEDUC2", "DMDMARTL", "INDFMPIR", "SDMVPSU", "SDMVSTRA",
                  "WTMEC2YR", "WTINT2YR")
GI_MCQ_VARS <- c("SEQN", "MCQ010", "MCQ160B", "MCQ160C", "MCQ160D",
                  "MCQ160E", "MCQ160F", "MCQ220")

# ── ggplot2 theme (reuse from 00_config or redefine) ─────────────────────────
theme_gi <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      plot.title    = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(color = "grey40"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      axis.title   = element_text(size = base_size),
      plot.caption = element_text(color = "grey50", size = base_size - 3,
                                  hjust = 0, face = "italic")
    )
}

cat("[00b_gi_cancer_config.R] GI cancer + mortality config loaded.\n")
cat(sprintf("  Cycles: %s\n", paste(GI_CYCLES, collapse = ", ")))
cat(sprintf("  GI sites: %s\n", paste(sapply(GI_SITES, `[[`, "label"),
                                       collapse = ", ")))
