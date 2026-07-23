# run_allages_analysis.R — Full analysis WITHOUT age filter
# Rscript run_allages_analysis.R

library(dplyr); library(survival); library(broom); library(ggplot2)

PROJ <- "D:/Researching/NHANES_aged_GI_tumor_nutrition"
RES_DIR <- file.path(PROJ, "results/gi_analysis")
FIG_DIR <- file.path(PROJ, "figures/gi_analysis")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

theme_pub <- function(base_size=12) {
  theme_classic(base_size=base_size) +
    theme(panel.grid.major.y=element_line(color="grey90", linewidth=0.3),
          plot.title=element_text(face="bold", size=base_size+2),
          plot.subtitle=element_text(color="grey40", size=base_size-1),
          legend.position="bottom")
}

df <- readRDS(file.path(RES_DIR, "nhanes_combined_gi_fixed.rds"))
gi <- df %>% filter(gi_tumor==1, surv_years>0) %>%
  mutate(PNI_s = as.numeric(scale(PNI)),
         CONUT_s = as.numeric(scale(-CONUT)),
         GNRI_s = as.numeric(scale(GNRI)),
         sex_b = ifelse(sex == "Female", 1L, 0L),
         race_eth_b = ifelse(race_eth == "Non-Hispanic White", 1L, 0L),
         pni_t = factor(ntile(PNI, 3), 1:3, c("Low (T1)","Mid (T2)","High (T3)")))

cat(sprintf("GI N=%d, events=%d\n", nrow(gi), sum(gi[["death"]])))
cat(sprintf("Age range: %.0f-%.0f (median %.0f)\n",
    min(gi[["age"]], na.rm=TRUE), max(gi[["age"]], na.rm=TRUE),
    median(gi[["age"]], na.rm=TRUE)))

# ── Table 1: Demographics ──────────────────────────────────────────────────
cat("\n=== Table 1: Demographics ===\n")
gi_status <- df %>% filter(!is.na(gi_tumor)) %>%
  mutate(gi_grp = case_when(gi_tumor==1 ~ "GI Tumor", any_cancer==0 ~ "Non-Cancer",
                             TRUE ~ "Other Cancer"))
t1 <- gi_status %>% group_by(gi_grp) %>%
  summarise(N=n(), Age=sprintf("%.1f (%.1f)", mean(age, na.rm=TRUE), sd(age, na.rm=TRUE)),
            Female=sprintf("%.1f%%", mean(sex=="Female", na.rm=TRUE)*100),
            PNI=sprintf("%.1f (%.1f)", mean(PNI, na.rm=TRUE), sd(PNI, na.rm=TRUE)),
            CONUT=sprintf("%.1f (%.1f)", mean(CONUT, na.rm=TRUE), sd(CONUT, na.rm=TRUE)),
            GNRI=sprintf("%.1f (%.1f)", mean(GNRI, na.rm=TRUE), sd(GNRI, na.rm=TRUE)),
            Deaths=sprintf("%.1f%%", mean(death, na.rm=TRUE)*100),
            .groups="drop")
print(as.data.frame(t1))
write.csv(t1, file.path(RES_DIR, "table1_allages.csv"), row.names=FALSE)

# ── Main Cox ────────────────────────────────────────────────────────────────
cat("\n=== All-cause Cox ===\n")
res <- data.frame()
for (idx in c("PNI_s", "CONUT_s", "GNRI_s")) {
  for (adj in c("Crude", "Adjusted")) {
    covs <- if (adj=="Crude") "1" else "age + sex_b + race_eth_b"
    f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", covs))
    fit <- coxph(f, data=gi)
    hr <- tidy(fit, conf.int=TRUE) %>% filter(term==idx)
    res <- rbind(res, data.frame(
      Index=gsub("_s","",idx), Adj=adj,
      HR=exp(hr[["estimate"]]), Lower=exp(hr[["conf.low"]]), Upper=exp(hr[["conf.high"]]),
      P=hr[["p.value"]]))
  }
}
res <- res %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
print(res[, c("Index","Adj","HR_CI","P")])
write.csv(res, file.path(RES_DIR, "cox_allcause_allages.csv"), row.names=FALSE)

# PH check
cat("\nProportional hazards test:\n")
for (idx in c("PNI_s", "CONUT_s", "GNRI_s")) {
  fit <- coxph(as.formula(paste("Surv(surv_years, death) ~", idx, "+ age + sex_b + race_eth_b")), data=gi)
  cat(sprintf("  %s: p=%.4f\n", idx, cox.zph(fit)$table[1,3]))
}

# ── Time-dependent Cox ────────────────────────────────────────────────────
cat("\n=== Time-dependent Cox ===\n")
gi_td <- survSplit(Surv(surv_years, death) ~ PNI_s + CONUT_s + GNRI_s +
                     age + sex_b + race_eth_b,
                    data=gi, cut=c(2,5), episode="period", id="id_td")
gi_td <- gi_td %>% mutate(period_f = factor(period))

td_res <- data.frame()
for (idx in c("PNI_s", "CONUT_s", "GNRI_s")) {
  f <- as.formula(paste("Surv(tstart, surv_years, death) ~", idx,
                        "* period_f + age + sex_b + race_eth_b + cluster(id_td)"))
  fit <- coxph(f, data=gi_td)
  hr_main <- tidy(fit, conf.int=TRUE) %>% filter(term==idx)
  hr_int <- tidy(fit, conf.int=TRUE) %>% filter(grepl(paste0(idx, ":"), term))
  td_res <- rbind(td_res, data.frame(
    Index=gsub("_s","",idx), Period="0-2 yr (main)",
    HR=exp(hr_main[["estimate"]]), Lower=exp(hr_main[["conf.low"]]),
    Upper=exp(hr_main[["conf.high"]]), P=hr_main[["p.value"]]))
  for (k in seq_len(nrow(hr_int))) {
    period_label <- gsub(paste0(idx, ":period_f"), "", hr_int[["term"]][k])
    td_res <- rbind(td_res, data.frame(
      Index=gsub("_s","",idx), Period=paste0("Interaction: period", period_label),
      HR=exp(hr_int[["estimate"]][k]), Lower=exp(hr_int[["conf.low"]][k]),
      Upper=exp(hr_int[["conf.high"]][k]), P=hr_int[["p.value"]][k]))
  }
}
td_res <- td_res %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
print(td_res[, c("Index","Period","HR_CI","P")])
write.csv(td_res, file.path(RES_DIR, "cox_timedependent_allages.csv"), row.names=FALSE)

# ── Landmark ────────────────────────────────────────────────────────────────
cat("\n=== Landmark ===\n")
lm_res <- data.frame()
for (lm_t in c(1,3,5)) {
  lm <- gi %>% filter(surv_years >= lm_t) %>% mutate(surv_rem = surv_years - lm_t)
  for (idx in c("PNI_s", "CONUT_s", "GNRI_s")) {
    f <- as.formula(paste("Surv(surv_rem, death) ~", idx, "+ age + sex_b + race_eth_b"))
    fit <- coxph(f, data=lm)
    hr <- tidy(fit, conf.int=TRUE) %>% filter(term==idx)
    lm_res <- rbind(lm_res, data.frame(
      Index=gsub("_s","",idx), Landmark=paste0(lm_t,"-yr"),
      N=nrow(lm), Events=sum(lm[["death"]]),
      HR=exp(hr[["estimate"]]), Lower=exp(hr[["conf.low"]]),
      Upper=exp(hr[["conf.high"]]), P=hr[["p.value"]]))
  }
}
lm_res <- lm_res %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
print(lm_res[, c("Index","Landmark","N","Events","HR_CI","P")])
write.csv(lm_res, file.path(RES_DIR, "cox_landmark_allages.csv"), row.names=FALSE)

# ── Competing risks ────────────────────────────────────────────────────────
cat("\n=== Competing risks ===\n")
cr_res <- data.frame()
gi <- gi %>% mutate(gi_cancer_death = ifelse(is.na(gi_cancer_death), 0L, gi_cancer_death),
                     other_death = ifelse(is.na(other_death), 0L, other_death))
for (cause_nm in c("gi_cancer_death", "other_death")) {
  for (idx in c("PNI_s", "CONUT_s", "GNRI_s")) {
    f <- as.formula(paste("Surv(surv_years,", cause_nm, ") ~", idx, "+ age + sex_b + race_eth_b"))
    fit <- coxph(f, data=gi)
    hr <- tidy(fit, conf.int=TRUE) %>% filter(term==idx)
    cr_res <- rbind(cr_res, data.frame(
      Index=gsub("_s","",idx),
      Cause=ifelse(cause_nm=="gi_cancer_death","GI cancer death","Other death"),
      HR=exp(hr[["estimate"]]), Lower=exp(hr[["conf.low"]]),
      Upper=exp(hr[["conf.high"]]), P=hr[["p.value"]]))
  }
}
cr_res <- cr_res %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
print(cr_res[, c("Index","Cause","HR_CI","P")])
write.csv(cr_res, file.path(RES_DIR, "cox_causespecific_allages.csv"), row.names=FALSE)

# ── KM ─────────────────────────────────────────────────────────────────────
cat("\n=== KM ===\n")
km <- survfit(Surv(surv_years, death) ~ pni_t, data=gi)
cat("Median survival:\n")
print(km)
saveRDS(km, file.path(RES_DIR, "km_fit_allages.rds"))

# ── Table 2: Baseline by PNI tertile ──────────────────────────────────────
cat("\n=== Table 2: Baseline by PNI tertile ===\n")
t2 <- gi %>% group_by(pni_t) %>%
  summarise(N=n(),
    Age=sprintf("%.1f (%.1f)", mean(age), sd(age)),
    Female=sprintf("%d (%.1f%%)", sum(sex=="Female"), mean(sex=="Female")*100),
    PNI=sprintf("%.1f (%.1f)", mean(PNI), sd(PNI)),
    CONUT=sprintf("%.1f (%.1f)", mean(CONUT), sd(CONUT)),
    GNRI=sprintf("%.1f (%.1f)", mean(GNRI), sd(GNRI)),
    Albumin=sprintf("%.2f (%.2f)", mean(albumin_gdl), sd(albumin_gdl)),
    Lymphocyte=sprintf("%.0f (%.0f)", mean(lymph_abs), sd(lymph_abs)),
    BMI=sprintf("%.1f (%.1f)", mean(bmi), sd(bmi)),
    Deaths=sprintf("%d (%.1f%%)", sum(death), mean(death)*100),
    Med_Surv=sprintf("%.1f yr", median(surv_years)),
    .groups="drop")
print(as.data.frame(t2))
write.csv(t2, file.path(RES_DIR, "table2_allages.csv"), row.names=FALSE)

# ── Figures ────────────────────────────────────────────────────────────────
# Fig 1: KM
library(survminer)
p_km <- ggsurvplot(km, data=gi,
  title="All-Cause Survival by PNI Tertile (All Ages, NHANES 1988-2016)",
  xlab="Years", ylab="Survival Probability",
  legend.title="PNI Tertile", pval=TRUE, pval.coord=c(0, 0.15),
  conf.int=TRUE, risk.table=TRUE, risk.table.height=0.25,
  palette=c("#D62728","#FF7F0E","#2CA02C"),
  surv.median.line="hv", ggtheme=theme_classic(base_size=12))
pdf(file.path(FIG_DIR, "fig1_km_allages.pdf"), width=10, height=8)
print(p_km); dev.off()
png(file.path(FIG_DIR, "fig1_km_allages.png"), width=10, height=8, units="in", res=300)
print(p_km); dev.off()
cat("Fig1 saved.\n")

# Fig 2: Forest
cox <- read.csv(file.path(RES_DIR, "cox_allcause_allages.csv"))
cox_adj <- cox %>% filter(Adj == "Adjusted") %>%
  mutate(Nutrition = factor(Index, levels=rev(c("PNI","CONUT","GNRI"))),
         sig = P < 0.05)
p_f <- ggplot(cox_adj, aes(x=HR, y=Nutrition)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey60") +
  geom_point(aes(color=sig), size=4) +
  geom_errorbarh(aes(xmin=Lower, xmax=Upper, color=sig), height=0.15, linewidth=1.3) +
  scale_color_manual(values=c("TRUE"="#D62728","FALSE"="#1F77B4"), guide="none") +
  scale_x_log10(breaks=c(0.4, 0.6, 0.8, 1.0)) +
  labs(title="Nutrition Indices and All-Cause Mortality (All Ages)",
       subtitle="Adjusted HR per 1-SD | GI Tumor Patients",
       x="Hazard Ratio (95% CI)", y="") +
  geom_text(aes(label=sprintf("HR=%.2f (%.2f-%.2f), p=%.4f", HR, Lower, Upper, P)),
            x=1.4, hjust=0, size=3.2, color="grey30") +
  coord_cartesian(xlim=c(0.4, 2.0)) + theme_classic(base_size=12)
ggsave(file.path(FIG_DIR, "fig2_cox_forest_allages.pdf"), p_f, width=9, height=4)
ggsave(file.path(FIG_DIR, "fig2_cox_forest_allages.png"), p_f, width=9, height=4, dpi=300)
cat("Fig2 saved.\n")

# Fig 3: Time-dependent
td <- read.csv(file.path(RES_DIR, "cox_timedependent_allages.csv"))
td_plot <- data.frame()
for (nut in c("PNI","CONUT","GNRI")) {
  sub <- td[grep(nut, td$Index), ]
  td_plot <- rbind(td_plot, data.frame(
    Nutrition=nut,
    Period=factor(c("0-2 yr","2-5 yr","5+ yr"), levels=c("0-2 yr","2-5 yr","5+ yr")),
    HR=c(sub$HR[1], sub$HR[1]*sub$HR[2], sub$HR[1]*sub$HR[3])))
}
td_plot$Nutrition <- factor(td_plot$Nutrition, levels=c("PNI","CONUT","GNRI"))
p_td <- ggplot(td_plot, aes(x=Period, y=HR, color=Nutrition, group=Nutrition)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey60") +
  geom_line(linewidth=1, position=position_dodge(0.2)) +
  geom_point(size=3.5, position=position_dodge(0.2)) +
  scale_color_manual(values=c("PNI"="#D62728","CONUT"="#FF7F0E","GNRI"="#2CA02C")) +
  labs(title="Time-Varying Effect of Nutrition on Mortality (All Ages)",
       x="Follow-up Period", y="Hazard Ratio", color="Index") +
  geom_text(data=subset(td_plot, Nutrition=="PNI"),
            aes(label=sprintf("HR=%.2f", HR)), vjust=-1.2, size=3) +
  theme_pub(13)
ggsave(file.path(FIG_DIR, "fig3_timedependent_allages.pdf"), p_td, width=9, height=5.5)
ggsave(file.path(FIG_DIR, "fig3_timedependent_allages.png"), p_td, width=9, height=5.5, dpi=300)

cat("All done. Results saved to:", RES_DIR, "\n")
