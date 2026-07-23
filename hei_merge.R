# hei_merge.R — Quick HEI re-merge with clean data
library(dplyr); library(broom); library(survival)

RES_DIR <- "results/gi_analysis"
DATA_DIR <- "data/gi_analysis"

df <- readRDS(file.path(RES_DIR, "nhanes_clean.rds"))
hei <- readRDS(file.path(DATA_DIR, "hei2015_scores.rds"))

df <- df %>% mutate(fped_cycle = case_when(
  grepl("2011-2012", cycle) ~ "1112",
  grepl("2013-2014", cycle) ~ "1314",
  grepl("2015-2016", cycle) ~ "1516",
  TRUE ~ NA_character_))
df <- left_join(df, hei, by="SEQN")
saveRDS(df, file.path(RES_DIR, "nhanes_clean_with_hei.rds"))

gi_hei <- df %>% filter(gi_tumor==1, surv_years>0, !is.na(HEI2015_total))
cat(sprintf("HEI GI: N=%d, deaths=%d\n", nrow(gi_hei), sum(gi_hei[["death"]])))

if (sum(gi_hei[["death"]]) >= 10) {
  gi_hei <- gi_hei %>% mutate(
    HEI_s = as.numeric(scale(as.numeric(HEI2015_total))),
    PNI_s = as.numeric(scale(PNI)),
    sex_b = ifelse(sex=="Female", 1L, 0L),
    race_eth_b = ifelse(race_eth=="Non-Hispanic White", 1L, 0L))

  for (i in c("HEI_s", "PNI_s")) {
    f <- as.formula(paste("Surv(surv_years, death) ~", i, "+ age + sex_b + race_eth_b"))
    h <- tidy(coxph(f, data=gi_hei), conf.int=TRUE) %>% filter(term==i)
    cat(sprintf("%s alone: HR=%.3f (%.3f-%.3f), p=%.4f\n",
        gsub("_s","",i), exp(h[["estimate"]]), exp(h[["conf.low"]]), exp(h[["conf.high"]]), h[["p.value"]]))
  }

  f2 <- coxph(Surv(surv_years, death) ~ PNI_s + HEI_s + age + sex_b + race_eth_b, data=gi_hei)
  for (i in c("PNI_s", "HEI_s")) {
    h2 <- tidy(f2, conf.int=TRUE) %>% filter(term==i)
    cat(sprintf("%s mutual: HR=%.3f (%.3f-%.3f), p=%.4f\n",
        gsub("_s","",i), exp(h2[["estimate"]]), exp(h2[["conf.low"]]), exp(h2[["conf.high"]]), h2[["p.value"]]))
  }

  cat(sprintf("HEI-PNI r: %.3f\n", with(gi_hei, cor(as.numeric(HEI2015_total), PNI, use="comp"))))
}
