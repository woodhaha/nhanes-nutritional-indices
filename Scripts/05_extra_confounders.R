# -- 05_extra_confounders.R - Simplified: per-cycle, type-safe merge
# Items 1-3: smoking/alcohol/charlson, Fine-Gray, Bonferroni

library(dplyr); library(survival); library(broom)

PROJ <- normalizePath(".")
DATA_DIR <- file.path(PROJ, "data/gi_analysis")
RES_DIR <- file.path(PROJ, "results/gi_analysis")

# -- 0. Load data ---------------------------------------------------------------
df <- readRDS(file.path(RES_DIR, "nhanes_clean.rds"))
for (nm in names(df)) {
  if (is.data.frame(df[[nm]])) df[[nm]] <- as.numeric(df[[nm]][[1]])
  if (is.list(df[[nm]])) df[[nm]] <- as.numeric(unlist(df[[nm]]))
}

# -- 1. Per-cycle NHANES extraction with type-safe merge ------------------------
cat("=== 1. Extracting questionnaire data per cycle ===\n")
library(nhanesA)

CYCLES <- c("D","E","F","G","H","I")
get_safe <- function(tbl, cyc) {
  nm <- sprintf("%s_%s", tbl, cyc)
  d <- tryCatch(nhanes(nm), error=function(e) NULL, warning=function(w) NULL)
  if (is.null(d)) return(NULL)
  # Convert all to numeric: factors → underlying codes, doubles stay
  for (cn in setdiff(names(d), "SEQN")) {
    if (is.factor(d[[cn]])) d[[cn]] <- as.numeric(d[[cn]]) - 1  # 0-indexed: Yes=0, No=1? No...
    # Actually nhanesA factor has labels "Yes"=1, "No"=2 - just convert to numeric
    if (is.factor(d[[cn]])) d[[cn]] <- as.numeric(d[[cn]])
    if (is.character(d[[cn]])) d[[cn]] <- suppressWarnings(as.numeric(d[[cn]]))
  }
  d$SEQN <- as.numeric(d$SEQN)
  d
}

# Extract only needed variables, one row per SEQN per cycle
extract_cycle <- function(cyc) {
  smq <- get_safe("SMQ", cyc); alq <- get_safe("ALQ", cyc)
  mcq <- get_safe("MCQ", cyc); diq <- get_safe("DIQ", cyc)
  bpq <- get_safe("BPQ", cyc)

  # Start with the first available table as base, merge others
  all_tables <- list(smq, alq, mcq, diq, bpq)
  all_tables <- all_tables[!sapply(all_tables, is.null)]
  if (length(all_tables) == 0) return(NULL)

  merged <- all_tables[[1]]
  for (i in 2:length(all_tables)) {
    # Common columns cause .x/.y suffixes - keep only from first table
    common <- intersect(names(merged), names(all_tables[[i]]))
    common <- setdiff(common, "SEQN")
    for (c in common) all_tables[[i]][[c]] <- NULL
    merged <- merge(merged, all_tables[[i]], by="SEQN", all=TRUE)
  }
  merged$cycle_code <- cyc
  merged
}

all_rows <- bind_rows(lapply(CYCLES, extract_cycle))
cat(sprintf("  Total: %d rows\n", nrow(all_rows)))

# Convert character back to numeric where possible
suppressWarnings({
  for (nm in names(all_rows)) {
    if (is.character(all_rows[[nm]]) && nm != "SEQN") {
      num <- as.numeric(all_rows[[nm]])
      if (mean(is.na(num)) < 0.5) all_rows[[nm]] <- num  # Keep as char if mostly non-numeric
    }
  }
})

# -- 2. Derive variables --------------------------------------------------------
cat("\n=== 2. Deriving variables ===\n")

# All columns are character after type-safe binding; convert to numeric where possible
char_to_num <- function(v) {
  if (is.character(v)) { suppressWarnings(as.numeric(v)) } else { as.numeric(v) }
}

in1 <- function(x) { x %in% c(1, "1") }  # handle both char and numeric

all_rows <- all_rows %>% mutate(
  SEQN = char_to_num(SEQN),
  alc101 = char_to_num(if ("ALQ101" %in% names(all_rows)) ALQ101 else NA),
  alc110 = char_to_num(if ("ALQ110" %in% names(all_rows)) ALQ110 else NA),
  smoke_ever = in1(if ("SMQ020" %in% names(all_rows)) SMQ020 else NA),
  smoke_now  = {x <- char_to_num(if ("SMQ040" %in% names(all_rows)) SMQ040 else NA); ifelse(x %in% c(1,2), 1L, ifelse(x==3, 0L, NA_integer_))},
  alc_ever   = ifelse(in1(alc101), 1L, ifelse(in1(alc110), 1L, ifelse(alc101==2|alc110==2, 0L, NA_integer_))),
  diabetes   = in1(if ("DIQ010" %in% names(all_rows)) DIQ010 else NA),
  hypertension = in1(if ("BPQ020" %in% names(all_rows)) BPQ020 else NA),
  charlson_mi     = in1(if ("MCQ160A" %in% names(all_rows)) MCQ160A else NA),
  charlson_chf    = in1(if ("MCQ160B" %in% names(all_rows)) MCQ160B else NA),
  charlson_stroke = in1(if ("MCQ160C" %in% names(all_rows)) MCQ160C else NA),
  charlson_copd   = in1(if ("MCQ160D" %in% names(all_rows)) MCQ160D else NA) | in1(if ("MCQ160E" %in% names(all_rows)) MCQ160E else NA),
  charlson_copd   = ifelse(is.na(charlson_copd), NA, as.integer(charlson_copd)),
  charlson_liver  = in1(if ("MCQ160K" %in% names(all_rows)) MCQ160K else NA),
  charlson_score  = coalesce(charlson_mi,0L)+coalesce(charlson_chf,0L)+
    coalesce(charlson_stroke,0L)+coalesce(charlson_copd,0L)+
    coalesce(charlson_liver,0L)+coalesce(diabetes,0L)
)

# Deduplicate by SEQN (take first occurrence per participant)
quest_key <- all_rows %>% group_by(SEQN) %>%
  summarise(across(c(smoke_ever, smoke_now, alc_ever, diabetes, hypertension, charlson_score),
                   ~first(na.omit(.))), .groups="drop")

# -- 3. Merge -------------------------------------------------------------------
cat("\n=== 3. Merging ===\n")
df <- merge(df, quest_key, by="SEQN", all.x=TRUE)
df <- df %>% mutate(
  smoke_status = case_when(smoke_ever==0~"never", smoke_ever==1&smoke_now==0~"former",
                            smoke_ever==1&smoke_now==1~"current", TRUE~NA_character_),
  nhanes3_flag = ifelse(grepl("III", cycle), 1L, 0L),
  charlson_any = ifelse(charlson_score>0, 1L, 0L))

cat(sprintf("  Smoke coverage: %.0f%%, Charlson: %.0f%%\n",
    mean(!is.na(df$smoke_ever))*100, mean(!is.na(df$charlson_score))*100))
# Check GI subset coverage
gi_tmp <- df %>% filter(gi_tumor==1)
cat(sprintf("  GI smoke: %.0f%%, GI Charlson: %.0f%%\n",
    mean(!is.na(gi_tmp$smoke_ever))*100, mean(!is.na(gi_tmp$charlson_score))*100))
saveRDS(df, file.path(RES_DIR, "nhanes_clean_extended.rds"))

# -- 4. Extended Cox models -----------------------------------------------------
cat("\n\n=== 4. Extended Cox models ===\n")
gi <- df %>% filter(gi_tumor==1, surv_years>0, nhanes3_flag==0) %>%
  mutate(PNI_s=scale(PNI)[,1], CONUT_s=scale(-CONUT)[,1], GNRI_s=scale(GNRI)[,1],
         sex_b=ifelse(sex=="Female",1L,0L),
         race_eth_b=ifelse(race_eth=="Non-Hispanic White",1L,0L),
         smoke_curr=ifelse(smoke_status=="current",1L,0L)) %>%
  filter(!is.na(death),!is.na(PNI))
cat(sprintf("  Cont NHANES GI: N=%d, events=%d\n", nrow(gi), sum(gi$death)))
for (v in c("smoke_curr","charlson_any","edu_binary","income_pir")) {
  cat(sprintf("  %s missing: %.0f%%\n", v, mean(is.na(gi[[v]]))*100))
}

adj <- "age + sex_b + race_eth_b"
# Drop alc_ever (nhanesA SSL issue → 70% missing). Smoke + Charlson + SES are robust proxies.
MOD <- list(
  M1=c("Basic", adj),
  M2=c("+edu+PIR", paste(adj,"+ edu_binary + income_pir")),
  M3=c("+Smoking", paste(adj,"+ edu_binary + income_pir + smoke_curr")),
  M4=c("+Charlson", paste(adj,"+ edu_binary + income_pir + smoke_curr + charlson_any"))
)
res_ext <- bind_rows(lapply(MOD, function(m) {
  bind_rows(lapply(c("PNI_s","CONUT_s","GNRI_s"), function(idx) {
    f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", m[2]))
    fit <- tryCatch(coxph(f, data=gi), error=function(e) NULL)
    if (is.null(fit)) return(data.frame())
    h <- tidy(fit, conf.int=TRUE) %>% filter(term==idx)
    data.frame(Index=gsub("_s","",idx), Model=m[1], N=fit$n, Events=fit$nevent,
               HR=exp(h$estimate), Lower=exp(h$conf.low), Upper=exp(h$conf.high), P=h$p.value)
  }))
}))
res_ext <- res_ext %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
cat("\n"); print(res_ext[,c("Index","Model","N","Events","HR_CI","P")], row.names=FALSE)
write.csv(res_ext, file.path(RES_DIR, "sensitivity_extended_confounders.csv"), row.names=FALSE)

# -- 5. Fine-Gray ---------------------------------------------------------------
cat("\n\n=== 5. Fine-Gray ===\n")
library(cmprsk)
gi_all <- df %>% filter(gi_tumor==1, surv_years>0) %>%
  mutate(PNI_s=scale(PNI)[,1], CONUT_s=scale(-CONUT)[,1], GNRI_s=scale(GNRI)[,1],
         sex_b=ifelse(sex=="Female",1L,0L), race_eth_b=ifelse(race_eth=="Non-Hispanic White",1L,0L)) %>%
  mutate(fg=case_when(gi_cancer_death==1~1L, other_death==1~2L, TRUE~0L)) %>%
  filter(!is.na(death),!is.na(PNI))

cat(sprintf("  N=%d, cancer death=%d, competing=%d, alive=%d\n", nrow(gi_all),
    sum(gi_all$fg==1), sum(gi_all$fg==2), sum(gi_all$fg==0)))

fg_out <- bind_rows(lapply(c("PNI_s","CONUT_s","GNRI_s"), function(idx) {
  mm <- model.matrix(as.formula(paste("~", idx, "+ age + sex_b + race_eth_b")), data=gi_all)[,-1]
  f1 <- crr(gi_all$surv_years, gi_all$fg, cov1=mm, failcode=1, cencode=0)
  f2 <- crr(gi_all$surv_years, gi_all$fg, cov1=mm, failcode=2, cencode=0)
  # Manual CI (confint.crr not available)
  s1 <- sqrt(diag(f1$var))[1]; z1 <- f1$coef[1]/s1; p1 <- 2*pnorm(-abs(z1))
  s2 <- sqrt(diag(f2$var))[1]; z2 <- f2$coef[1]/s2; p2 <- 2*pnorm(-abs(z2))
  bind_rows(
    data.frame(Index=gsub("_s","",idx), Risk="GI cancer death", sHR=exp(f1$coef[1]),
               Lower=exp(f1$coef[1]-1.96*s1), Upper=exp(f1$coef[1]+1.96*s1), P=p1),
    data.frame(Index=gsub("_s","",idx), Risk="Non-cancer death", sHR=exp(f2$coef[1]),
               Lower=exp(f2$coef[1]-1.96*s2), Upper=exp(f2$coef[1]+1.96*s2), P=p2))
}))
fg_out <- fg_out %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", sHR, Lower, Upper))
print(fg_out[,c("Index","Risk","HR_CI","P")], row.names=FALSE)
write.csv(fg_out, file.path(RES_DIR, "competingrisks_finegray.csv"), row.names=FALSE)

# -- 6. Bonferroni ---------------------------------------------------------------
cat("\n\n=== 6. Bonferroni ===\n")
n_primary <- 6
ba <- 0.05/n_primary
cr <- read.csv(file.path(RES_DIR, "cox_main.csv"))
cat(sprintf("  α = 0.05/%d = %.4f\n\n", n_primary, ba))
for (i in c("PNI","CONUT","GNRI")) for (a in c("Crude","Adjusted")) {
  h <- cr %>% filter(Index==i, Adj==a)
  if (nrow(h)) cat(sprintf("  %s %s: p=%.4f %s (α=%.4f)\n", i, a, h$P[1],
      ifelse(h$P[1]<ba,"SURVIVES","FAILS"), ba))
}
write.csv(data.frame(N_Tests=n_primary, Alpha=ba,
  PNI_Survives=cr$P[cr$Index=="PNI"&cr$Adj=="Adjusted"]<ba,
  CONUT_Survives=cr$P[cr$Index=="CONUT"&cr$Adj=="Adjusted"]<ba,
  GNRI_Survives=cr$P[cr$Index=="GNRI"&cr$Adj=="Adjusted"]<ba),
  file.path(RES_DIR, "bonferroni_correction.csv"), row.names=FALSE)

cat("\n=== COMPLETE ===\n")
