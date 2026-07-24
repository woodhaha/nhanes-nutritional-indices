# 06_pancancer_prep.R — Classify all NHANES cancer patients into pancancer groups
# Extends GI analysis to cover breast, lung, female repro, and other sites
# ==============================================================================

library(dplyr)
library(nhanesA)
library(data.table)

PROJ <- normalizePath(".")
DATA_DIR <- file.path(PROJ, "data/gi_analysis")
RES_DIR <- file.path(PROJ, "results/gi_analysis")
OUTPUT <- file.path(PROJ, "data/pancancer")
dir.create(OUTPUT, showWarnings = FALSE, recursive = TRUE)

# ── Harmonized pancancer group mapping (text labels from nhanesA) ────────────
# nhanesA returns MCQ230A as text labels, not numeric codes

PAN_GROUPS <- list(
  GI = c("Colon", "Rectum (rectal)", "Esophagus (esophageal)", "Stomach",
         "Liver", "Pancreas", "Gallbladder", "Gallbladder (gallbladder)"),
  Breast = c("Breast"),
  Lung = c("Lung", "Lung (lung)"),
  FemaleRepro = c("Cervix (cervical)", "Cervix", "Ovary (ovarian)", "Ovary",
                   "Uterus (uterine)", "Uterine"),
  ProstateUrinary = c("Prostate", "Bladder", "Kidney", "Testis (testicular)", "Testicular"),
  Hematologic = c("Blood", "Leukemia", "Lymphoma/Hodgkin's disease", "Lymphoma",
                  "Lymphoma/Hodgkin", "Lymphoma/Non-Hodgkin", "Lymphoma/NonHodgkin"),
  OtherSolid = c("Melanoma", "Thyroid", "Brain", "Bone", "Mouth/tongue/lip",
                 "Larynx/ windpipe", "Larynx/Trachea", "Soft tissue (muscle or fat)",
                 "Nervous System", "Skin (non-melanoma)", "Skin (nonmelanoma)",
                 "Skin (don't know what kind)", "Skin_NonMelanoma", "Skin_Unknown",
                 "Other")
)

# NHANES III HAC3OS map (codes 1-22 → group)
HAC3OS_TO_GROUP <- list(
  "1"="Bladder", "2"="Bone", "3"="Brain", "4"="Breast",
  "5"="Cervix", "6"="Colon", "7"="Esophagus", "8"="Gallbladder",
  "9"="Kidney", "10"="Larynx", "11"="Leukemia", "12"="Liver",
  "13"="Lung", "14"="Lymphoma", "15"="Melanoma", "16"="Other",
  "17"="Ovary", "18"="Pancreas", "19"="Prostate", "20"="Rectum",
  "21"="Stomach", "22"="Thyroid"
)

# Map from harmonized site name → pancancer group
build_label_to_group <- function() {
  m <- c()
  for (grp in names(PAN_GROUPS)) {
    for (label in PAN_GROUPS[[grp]]) {
      m[label] <- grp
    }
  }
  # Also HAC3OS-derived names
  hac_to_group <- list(
    "Bladder"="ProstateUrinary", "Bone"="OtherSolid", "Brain"="OtherSolid",
    "Breast"="Breast", "Cervix"="FemaleRepro", "Colon"="GI",
    "Esophagus"="GI", "Gallbladder"="GI", "Kidney"="ProstateUrinary",
    "Larynx"="OtherSolid", "Leukemia"="Hematologic", "Liver"="GI",
    "Lung"="Lung", "Lymphoma"="Hematologic", "Melanoma"="OtherSolid",
    "Other"="OtherSolid", "Ovary"="FemaleRepro", "Pancreas"="GI",
    "Prostate"="ProstateUrinary", "Rectum"="GI", "Stomach"="GI",
    "Thyroid"="OtherSolid"
  )
  for (label in names(hac_to_group)) {
    m[label] <- hac_to_group[[label]]
  }
  m
}
label_to_group <- build_label_to_group()

classify_pan_group <- function(site_label) {
  if (is.na(site_label) || nchar(trimws(site_label)) == 0) return(NA_character_)
  sl <- trimws(as.character(site_label))
  if (sl %in% c("Don't know", "Refused", "88", "99")) return(NA_character_)
  if (sl %in% names(label_to_group)) return(label_to_group[[sl]])
  # Try normalization
  sl_lower <- tolower(sl)
  for (grp in names(PAN_GROUPS)) {
    for (label in PAN_GROUPS[[grp]]) {
      if (tolower(label) == sl_lower) return(grp)
    }
  }
  return("OtherSolid")
}

# ── Step 1: Load clean data ──────────────────────────────────────────────────
cat("\n── Step 1: Loading clean data ──\n")
clean <- readRDS(file.path(RES_DIR, "nhanes_clean.rds"))
cat(sprintf("Total N: %d\n", nrow(clean)))
cat(sprintf("Cancer patients (any_cancer=1): %d\n", sum(clean$any_cancer == 1, na.rm=TRUE)))
cat(sprintf("Cancer patients needing site classification: %d\n",
            sum(clean$any_cancer == 1 & clean$gi_tumor != 1, na.rm=TRUE)))

# ── Step 2: NHANES III site codes ────────────────────────────────────────────
cat("\n── Step 2: NHANES III cancer sites (HAC3OS) ──\n")
n3 <- readRDS("data/gi_analysis/nhanes3_gi_allages.rds")
n3_sites <- n3[!is.na(n3$HAC3OS) & !n3$HAC3OS %in% c(88, 99), c("SEQN", "HAC3OS")]
n3_sites$site_label <- unlist(HAC3OS_TO_GROUP[as.character(n3_sites$HAC3OS)])
n3_sites$pan_group <- sapply(n3_sites$site_label, function(x) {
  g <- label_to_group[x]
  if (is.null(g) || is.na(g)) "OtherSolid" else g
}, USE.NAMES = FALSE)

cat(sprintf("NHANES III cancer patients with site: %d\n", nrow(n3_sites)))

# ── Step 3: Continuous NHANES site codes ─────────────────────────────────────
cat("\n── Step 3: Continuous NHANES cancer sites (MCQ230A) ──\n")

mcq_cycles <- list(
  "NHANES_2005-2006" = "MCQ_D",
  "NHANES_2007-2008" = "MCQ_E",
  "NHANES_2009-2010" = "MCQ_F",
  "NHANES_2011-2012" = "MCQ_G",
  "NHANES_2013-2014" = "MCQ_H",
  "NHANES_2015-2016" = "MCQ_I"
)

all_cancer_sites <- list()
for (cycle_name in names(mcq_cycles)) {
  tbl <- mcq_cycles[[cycle_name]]
  raw <- nhanesA::nhanes(tbl)

  if (!"MCQ230A" %in% names(raw)) {
    cat(sprintf("  %s: no MCQ230A column, skipping\n", cycle_name))
    next
  }

  sites <- data.frame(
    SEQN = raw$SEQN,
    site_label = as.character(raw$MCQ230A),
    stringsAsFactors = FALSE
  )
  # Filter valid cancer site labels
  sites <- sites[!is.na(sites$site_label) & sites$site_label != "" &
                   !sites$site_label %in% c("Don't know", "Refused", "9"), ]
  # Remove purely numeric labels
  sites <- sites[grepl("[A-Za-z]", sites$site_label), ]

  sites$pan_group <- sapply(sites$site_label, classify_pan_group)
  sites$cycle_src <- cycle_name
  sites$SEQN <- as.numeric(sites$SEQN)

  all_cancer_sites[[cycle_name]] <- sites
  cat(sprintf("  %s: %d cancer site records\n", cycle_name, nrow(sites)))
}

cont_sites <- bind_rows(all_cancer_sites)
cat(sprintf("Total continuous NHANES: %d\n", nrow(cont_sites)))

# ── Step 4: Merge ────────────────────────────────────────────────────────────
cat("\n── Step 4: Merging site codes with clean data ──\n")

# Combine site data
all_site_data <- bind_rows(
  n3_sites %>% dplyr::select(SEQN, site_label, pan_group),
  cont_sites %>% dplyr::select(SEQN, site_label, pan_group)
)
# Deduplicate (keep first match per SEQN)
all_site_data <- all_site_data[!duplicated(all_site_data$SEQN), ]
cat(sprintf("Total unique SEQNs with site: %d\n", nrow(all_site_data)))

cat("\nPancancer group distribution:\n")
print(table(all_site_data$pan_group, useNA="ifany"))

# Merge into clean data
clean$SEQN <- as.numeric(clean$SEQN)
clean <- clean %>%
  left_join(all_site_data %>% dplyr::select(SEQN, site_label, pan_group), by = "SEQN")

# For GI tumor patients (existing classification), override with "GI"
gi_idx <- which(clean$gi_tumor == 1 & !is.na(clean$gi_site))
clean$pan_group[gi_idx] <- "GI"
clean$site_label[gi_idx] <- as.character(clean$gi_site[gi_idx])

# Non-cancer patients
clean$pan_group[clean$any_cancer == 0 & is.na(clean$pan_group)] <- "NonCancer"

# Death type classification
clean <- clean %>%
  mutate(
    ucod_num = as.numeric(as.character(ucod_leading)),
    cancer_death = as.integer(!is.na(ucod_num) & ucod_num == 2),
    noncancer_death = as.integer(!is.na(ucod_num) & ucod_num != 2 & ucod_num > 0),
    cancer_type_death = as.integer(cancer_death == 1 & !is.na(pan_group) & pan_group != "NonCancer"),
    cod_status2 = case_when(
      is.na(ucod_num) | ucod_num <= 0 ~ 0L,
      cancer_type_death == 1 ~ 1L,
      TRUE ~ 2L
    )
  )

# Summary
cat(sprintf("\nFinal sample: N=%d\n", nrow(clean)))
pan_tab <- table(clean$pan_group, useNA="ifany")
print(pan_tab)

# Deaths by group
for (g in names(pan_tab)[pan_tab > 0]) {
  sub <- clean[clean$pan_group == g, ]
  cat(sprintf("  %-20s: N=%-5d deaths=%d (%.1f%%)\n",
              g, nrow(sub), sum(sub$death, na.rm=TRUE),
              mean(sub$death, na.rm=TRUE)*100))
}

saveRDS(clean, file.path(OUTPUT, "nhanes_pancancer.rds"))
cat("\nSaved: data/pancancer/nhanes_pancancer.rds\n")
