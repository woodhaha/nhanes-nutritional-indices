# nlr_sii_extract.R — Extract NLR (neutrophil/lymphocyte) and SII (platelet×neutrophil/lymphocyte)
# Only available for NHANES 2005-2016 (not NHANES III)
library(dplyr); library(broom); library(survival); library(nhanesA)

PROJ <- normalizePath(".")
RES_DIR <- file.path(PROJ, "results/gi_analysis")

# Step 1: Download raw CBC data → extract neutrophil + platelet
unlabel <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "haven_labelled")) as.numeric(x) else as.numeric(x)
}

CYCLES <- c("D","E","F","G","H","I")
cbc_all <- list()
for (s in CYCLES) {
  cb <- nhanesA::nhanes(paste0("CBC_", s))
  cbc_all[[s]] <- data.frame(
    SEQN = as.numeric(unlabel(cb$SEQN)),
    neut_abs = unlabel(cb[["LBDNENO"]]) * 1000,      # ×10³/µL → cells/µL
    platelet = unlabel(cb[["LBXPLTSI"]]) * 1000,      # ×10³/µL → cells/µL
    wbc = unlabel(cb[["LBXWBCSI"]]),
    lymph_pct = unlabel(cb[["LBXLYPCT"]]),
    stringsAsFactors = FALSE)
  cat(sprintf("CBC_%s: %d\n", s, nrow(cbc_all[[s]])))
}
cbc_all <- bind_rows(cbc_all)
saveRDS(cbc_all, file.path(RES_DIR, "cbc_extracted.rds"))
cat(sprintf("Total CBC records: %d\n", nrow(cbc_all)))

# Step 2: Merge into clean data (2005-2016 subset)
df <- readRDS(file.path(RES_DIR, "nhanes_clean.rds"))
df <- merge(df, cbc_all[, c("SEQN","neut_abs","platelet")], by="SEQN", all.x=TRUE)
cat(sprintf("Merged: %d with CBC, %d missing\n",
    sum(!is.na(df$neut_abs)), sum(is.na(df$neut_abs))))

# Compute NLR and SII
df <- df %>% mutate(
  NLR = neut_abs / lymph_abs,          # neutrophil count / lymphocyte count
  SII = (platelet * neut_abs) / lymph_abs  # platelet × neutrophil / lymphocyte
)

# Flag NLR outliers for capping
cat(sprintf("NLR range: %.1f-%.1f (median %.1f)\n",
    min(df$NLR, na.rm=TRUE), max(df$NLR, na.rm=TRUE),
    median(df$NLR, na.rm=TRUE)))
cat(sprintf("SII range: %.0f-%.0f (median %.0f)\n",
    min(df$SII, na.rm=TRUE), max(df$SII, na.rm=TRUE),
    median(df$SII, na.rm=TRUE)))

saveRDS(df, file.path(RES_DIR, "nhanes_clean.rds"))

# Step 3: Analyze in GI tumor subset
gi <- df %>% filter(gi_tumor==1, surv_years>0, !is.na(NLR), !is.na(SII)) %>%
  mutate(
    NLR_s = as.numeric(scale(NLR)),     # z-score for comparability
    logNLR_s = as.numeric(scale(log(NLR))),
    SII_s = as.numeric(scale(SII)),
    logSII_s = as.numeric(scale(log(SII))),
    PNI_s = as.numeric(scale(PNI)),
    sex_b = ifelse(sex=="Female", 1L, 0L),
    race_eth_b = ifelse(race_eth=="Non-Hispanic White", 1L, 0L))

cat(sprintf("\nNLR/SII GI subset: N=%d, events=%d\n", nrow(gi), sum(gi[["death"]])))

adj <- "age + sex_b + race_eth_b"

res <- bind_rows(lapply(c("PNI_s", "NLR_s", "logNLR_s", "SII_s", "logSII_s"), function(idx) {
  f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", adj))
  h <- tidy(coxph(f, data=gi), conf.int=TRUE) %>% filter(term==idx)
  data.frame(Index = gsub("_s","",idx), N=nrow(gi), Events=sum(gi[["death"]]),
             HR=exp(h[["estimate"]]), Lower=exp(h[["conf.low"]]),
             Upper=exp(h[["conf.high"]]), P=h[["p.value"]])
}))
res <- res %>% mutate(HR_CI = sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
cat("\n=== PNI vs NLR vs SII ===\n")
print(res[, c("Index","N","Events","HR_CI","P")])
write.csv(res, file.path(RES_DIR, "nlr_sii_comparison.csv"), row.names=FALSE)

# C-stat comparison
cstat <- function(fit) fit$concordance[["concordance"]]
f_base <- coxph(Surv(surv_years, death) ~ age + sex_b + race_eth_b, data=gi)
f_pni <- coxph(Surv(surv_years, death) ~ PNI_s + age + sex_b + race_eth_b, data=gi)
f_nlr <- coxph(Surv(surv_years, death) ~ logNLR_s + age + sex_b + race_eth_b, data=gi)
f_sii <- coxph(Surv(surv_years, death) ~ logSII_s + age + sex_b + race_eth_b, data=gi)
cat(sprintf("\nC-stat: Base=%.3f, +PNI=%.3f, +logNLR=%.3f, +logSII=%.3f\n",
    cstat(f_base), cstat(f_pni), cstat(f_nlr), cstat(f_sii)))

# Mutual adjustment
cat("\n=== Mutual adjustment ===\n")
f <- coxph(Surv(surv_years, death) ~ PNI_s + logNLR_s + age + sex_b + race_eth_b, data=gi)
for (i in c("PNI_s", "logNLR_s")) {
  h <- tidy(f, conf.int=TRUE) %>% filter(term==i)
  cat(sprintf("  %s: HR=%.3f (%.3f-%.3f), p=%.4f\n",
      gsub("_s","",i), exp(h[["estimate"]]), exp(h[["conf.low"]]), exp(h[["conf.high"]]), h[["p.value"]]))
}

cat("\nDone.\n")
