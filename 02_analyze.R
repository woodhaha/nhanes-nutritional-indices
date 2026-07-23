# ─── 02_analyze.R — All analyses + figures ─────────────────────────────────
# Rscript 02_analyze.R
# Produces: Table 1-2, Cox, time-dependent, landmark, competing risks,
#            cutpoint, RCS, subgroup forest, PNI decomposition, NLR/SII

library(dplyr); library(tidyr); library(survival); library(broom)
library(ggplot2); library(survminer); library(splines)

PROJ <- normalizePath(".")
RES_DIR <- file.path(PROJ, "results/gi_analysis")
FIG_DIR <- file.path(PROJ, "figures/gi_analysis")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

theme_pub <- function(bs=12) {
  theme_classic(base_size=bs) +
    theme(panel.grid.major.y=element_line(color="grey90", linewidth=0.3),
          plot.title=element_text(face="bold", size=bs+2),
          plot.subtitle=element_text(color="grey40", size=bs-1),
          legend.position="bottom")
}

# ── Load data ─────────────────────────────────────────────────────────────
cat("Loading...\n")
df <- readRDS(file.path(RES_DIR, "nhanes_clean.rds"))
gi <- df %>% filter(gi_tumor==1, surv_years>0) %>%
  mutate(PNI_s = as.numeric(scale(PNI)),
         CONUT_s = as.numeric(scale(-CONUT)),
         GNRI_s = as.numeric(scale(GNRI)),
         age_s = as.numeric(scale(age)),
         sex_b = ifelse(sex=="Female", 1L, 0L),
         race_eth_b = ifelse(race_eth=="Non-Hispanic White", 1L, 0L),
         pni_t = factor(ntile(PNI, 3), 1:3, c("Low (T1)","Mid (T2)","High (T3)")))

cat(sprintf("GI N=%d, events=%d\n", nrow(gi), sum(gi[["death"]])))

# ═══════════════════════════════════════════════════════════════════════════
# TABLE 1 — Comprehensive
# ═══════════════════════════════════════════════════════════════════════════
cat("=== Table 1 ===\n")
fmt_m <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
fmt_n <- function(x) sprintf("%d (%.1f%%)", sum(x, na.rm=TRUE), mean(x, na.rm=TRUE)*100)

t1 <- df %>% filter(!is.na(gi_tumor), surv_years>0) %>%
  mutate(grp = case_when(gi_tumor==1~"GI_Tumor", any_cancer==0~"Non_Cancer", TRUE~"Other_Cancer")) %>%
  group_by(grp) %>%
  summarise(N=n(), Age=fmt_m(age), Female=fmt_n(sex=="Female"),
            Race_NH_White=fmt_n(race_eth=="Non-Hispanic White"),
            Race_NH_Black=fmt_n(race_eth=="Non-Hispanic Black"),
            Race_Hispanic=fmt_n(race_eth=="Hispanic"),
            Edu_AboveHS=fmt_n(edu_binary==1),
            BMI=fmt_m(bmi), Albumin=fmt_m(albumin_gdl),
            Lymphocyte=fmt_m(lymph_abs), Chol=fmt_m(tchol_mgdl),
            PNI=fmt_m(PNI), CONUT=fmt_m(CONUT), GNRI=fmt_m(GNRI),
            Deaths=fmt_n(death==1), Surv=fmt_m(surv_years), .groups="drop")
write.csv(t1, file.path(RES_DIR, "table1.csv"), row.names=FALSE)
print(as.data.frame(t1), row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════
# TABLE 2 — Baseline by PNI tertile (GI only)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== Table 2 ===\n")
t2 <- gi %>% group_by(pni_t) %>%
  summarise(N=n(), Age=fmt_m(age), Female=fmt_n(sex=="Female"),
            NH_White=fmt_n(race_eth=="Non-Hispanic White"),
            Edu_AboveHS=fmt_n(edu_binary==1),
            PNI=fmt_m(PNI), CONUT=fmt_m(CONUT), GNRI=fmt_m(GNRI),
            Albumin=fmt_m(albumin_gdl), Lymphocyte=fmt_m(lymph_abs),
            Chol=fmt_m(tchol_mgdl), BMI=fmt_m(bmi),
            Deaths=fmt_n(death==1), Med_Surv=sprintf("%.1f yr", median(surv_years)),
            .groups="drop")
write.csv(t2, file.path(RES_DIR, "table2.csv"), row.names=FALSE)
print(as.data.frame(t2), row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 1 — Main Cox (3 indices)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A1: Main Cox ===\n")
run_cox <- function(data, idx, covs) {
  f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", covs))
  fit <- coxph(f, data=data)
  hr <- tidy(fit, conf.int=TRUE) %>% filter(term==idx)
  data.frame(Index=gsub("_s","",idx), HR=exp(hr$estimate),
             Lower=exp(hr$conf.low), Upper=exp(hr$conf.high), P=hr$p.value)
}

adj_cov <- "age + sex_b + race_eth_b"
cox1 <- bind_rows(lapply(c("PNI_s","CONUT_s","GNRI_s"),
                          function(i) bind_rows(run_cox(gi,i,"1"),
                                                run_cox(gi,i,adj_cov))))
cox1$Adj <- rep(c("Crude","Adjusted"), 3)
cox1 <- cox1 %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
write.csv(cox1, file.path(RES_DIR, "cox_main.csv"), row.names=FALSE)
print(cox1[, c("Index","Adj","HR_CI","P")])

# PH test
cat("\nPH test:\n")
for (i in c("PNI_s","CONUT_s","GNRI_s")) {
  f <- as.formula(paste("Surv(surv_years, death) ~", i, "+", adj_cov))
  cat(sprintf("  %s: p=%.4f\n", i, cox.zph(coxph(f, data=gi))$table[1,3]))
}

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 2 — Time-dependent Cox
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A2: Time-dependent Cox ===\n")
gi_td <- survSplit(Surv(surv_years, death) ~ PNI_s + CONUT_s + GNRI_s + age + sex_b + race_eth_b,
                    data=gi, cut=c(2,5), episode="period", id="id_td") %>%
  mutate(period_f = factor(period))

td_res <- bind_rows(lapply(c("PNI_s","CONUT_s","GNRI_s"), function(idx) {
  f <- as.formula(paste("Surv(tstart, surv_years, death) ~", idx,
                        "* period_f + age + sex_b + race_eth_b + cluster(id_td)"))
  fit <- coxph(f, data=gi_td)
  h <- tidy(fit, conf.int=TRUE); m <- h %>% filter(term==idx)
  int <- h %>% filter(grepl(paste0(idx,":"), term))
  rbind(data.frame(Index=gsub("_s","",idx), Period="0-2 yr",
                   HR=exp(m$estimate), Lower=exp(m$conf.low), Upper=exp(m$conf.high), P=m$p.value),
        data.frame(Index=gsub("_s","",idx),
                   Period=paste0("vs: period", gsub(paste0(idx,":period_f"),"",int$term)),
                   HR=exp(int$estimate), Lower=exp(int$conf.low),
                   Upper=exp(int$conf.high), P=int$p.value))
}))
td_res <- td_res %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
write.csv(td_res, file.path(RES_DIR, "cox_timedependent.csv"), row.names=FALSE)
print(td_res[, c("Index","Period","HR_CI","P")])

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 3 — Landmark
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A3: Landmark ===\n")
lm_res <- bind_rows(lapply(c(1,3,5), function(lm_t) {
  lm <- gi %>% filter(surv_years >= lm_t) %>% mutate(surv_rem=surv_years-lm_t)
  bind_rows(lapply(c("PNI_s","CONUT_s","GNRI_s"), function(idx) {
    f <- as.formula(paste("Surv(surv_rem, death) ~", idx, "+", adj_cov))
    hr <- tidy(coxph(f, data=lm), conf.int=TRUE) %>% filter(term==idx)
    data.frame(Index=gsub("_s","",idx), Landmark=paste0(lm_t,"-yr"),
               N=nrow(lm), Events=sum(lm$death),
               HR=exp(hr$estimate), Lower=exp(hr$conf.low),
               Upper=exp(hr$conf.high), P=hr$p.value)
  }))
}))
lm_res <- lm_res %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
write.csv(lm_res, file.path(RES_DIR, "cox_landmark.csv"), row.names=FALSE)
print(lm_res[, c("Index","Landmark","N","Events","HR_CI","P")])

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 4 — Competing risks
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A4: Competing risks ===\n")
cr_res <- bind_rows(lapply(c("gi_cancer_death","other_death"), function(cause) {
  bind_rows(lapply(c("PNI_s","CONUT_s","GNRI_s"), function(idx) {
    f <- as.formula(paste("Surv(surv_years,",cause,") ~", idx, "+", adj_cov))
    hr <- tidy(coxph(f, data=gi), conf.int=TRUE) %>% filter(term==idx)
    data.frame(Index=gsub("_s","",idx),
               Cause=ifelse(grepl("gi",cause),"GI cancer death","Other death"),
               HR=exp(hr$estimate), Lower=exp(hr$conf.low),
               Upper=exp(hr$conf.high), P=hr$p.value)
  }))
}))
cr_res <- cr_res %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
write.csv(cr_res, file.path(RES_DIR, "cox_causespecific.csv"), row.names=FALSE)
print(cr_res[, c("Index","Cause","HR_CI","P")])

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 5 — PNI cutpoint (optimal threshold)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A5: PNI cutpoint ===\n")
# Contal & O'Quigley method via maximally selected log-rank statistic
suppressPackageStartupMessages(if (!require("maxstat")) install.packages("maxstat", repos="https://cloud.r-project.org"))
library(maxstat)
mt <- maxstat.test(Surv(surv_years, death) ~ PNI, data=gi, smethod="LogRank",
                    pmethod="Lau92", minprop=0.2, maxprop=0.8)
cut_pni <- as.numeric(mt$estimate)
cat(sprintf("  Optimal PNI cutpoint: %.1f (p=%.4f)\n", cut_pni, mt$p.value))

gi <- gi %>% mutate(pni_high = ifelse(PNI >= cut_pni, "High", "Low"))
km_cut <- survfit(Surv(surv_years, death) ~ pni_high, data=gi)
cat("  Median survival by cutpoint:\n")
print(km_cut)

p_cut <- ggsurvplot(km_cut, data=gi,
  title=sprintf("PNI Cutpoint Analysis: %s (optimal threshold = %.1f)", cut_pni, cut_pni),
  subtitle=sprintf("Log-rank p=%.4f | Maximal selected test", mt$p.value),
  xlab="Years", ylab="Survival Probability", legend.title="PNI Group",
  pval=TRUE, pval.coord=c(0,0.15), risk.table=TRUE, risk.table.height=0.25,
  palette=c("#D62728","#2CA02C"), surv.median.line="hv",
  ggtheme=theme_classic(base_size=12))
pdf(file.path(FIG_DIR, "fig_cutpoint.pdf"), width=10, height=8)
print(p_cut); dev.off()

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 6 — RCS dose-response
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A6: RCS ===\n")
# 3-knot restricted cubic spline
rcs <- coxph(Surv(surv_years, death) ~ ns(PNI, df=3) + age + sex_b + race_eth_b, data=gi)
rcs_pred <- predict(rcs, type="terms", se=TRUE)
rcs_idx <- grep("PNI", colnames(rcs_pred$fit))
pni_seq <- seq(min(gi$PNI, na.rm=TRUE), max(gi$PNI, na.rm=TRUE), length.out=200)

# Predict at reference values
ref <- data.frame(PNI=pni_seq, age=mean(gi$age), sex_b=1, race_eth_b=1)
pp <- predict(rcs, newdata=ref, type="risk", se.fit=TRUE)
rcs_df <- data.frame(PNI=pni_seq, HR=pp$fit, Lower=pp$fit-1.96*pp$se.fit, Upper=pp$fit+1.96*pp$se.fit)
rcs_df$HR[rcs_df$HR > 5] <- NA; rcs_df$Upper[rcs_df$Upper > 10] <- 10

p_rcs <- ggplot(rcs_df, aes(x=PNI, y=HR)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey60") +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2, fill="#D62728") +
  geom_line(color="#D62728", linewidth=1.2) +
  geom_rug(data=gi, aes(x=PNI, y=1), sides="b", alpha=0.1) +
  labs(title="PNI Dose-Response: Restricted Cubic Spline",
       subtitle="Adjusted for age, sex, race | 3 knots | NHANES 1988-2016",
       x="PNI", y="Hazard Ratio (95% CI)") +
  theme_pub(13) + coord_cartesian(ylim=c(0, 3))
ggsave(file.path(FIG_DIR, "fig_rcs.pdf"), p_rcs, width=8, height=6)

# Test nonlinearity
rcs_lin <- coxph(Surv(surv_years, death) ~ PNI_s + age + sex_b + race_eth_b, data=gi)
cat(sprintf("  Nonlinearity LRT: p=%.4f\n", anova(rcs_lin, rcs, test="LRT")$`Pr(>|Chi|)`[2]))
write.csv(rcs_df, file.path(RES_DIR, "rcs_predictions.csv"), row.names=FALSE)

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 7 — Subgroup forest
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A7: Subgroup forest ===\n")
sgrp <- function(label, subset, adj=adj_cov) {
  s <- gi[subset, ]
  e <- sum(s$death)
  if (e < 10) return(data.frame(Subgroup=label, N=nrow(s), Events=e, HR=NA, Lower=NA, Upper=NA, P=NA))
  f <- as.formula(paste("Surv(surv_years, death) ~ PNI_s +", adj))
  hr <- tidy(coxph(f, data=s), conf.int=TRUE) %>% filter(term=="PNI_s")
  data.frame(Subgroup=label, N=nrow(s), Events=e,
             HR=exp(hr$estimate), Lower=exp(hr$conf.low),
             Upper=exp(hr$conf.high), P=hr$p.value)
}

sg <- bind_rows(list(
  sgrp("All patients", TRUE),
  sgrp("Age < 60", gi$age < 60),
  sgrp("Age 60-74", gi$age >= 60 & gi$age < 75),
  sgrp("Age >= 75", gi$age >= 75),
  sgrp("Female", gi$sex == "Female"),
  sgrp("Male", gi$sex == "Male"),
  sgrp("Non-Hispanic White", gi$race_eth == "Non-Hispanic White"),
  sgrp("Non-Hispanic Black", gi$race_eth == "Non-Hispanic Black"),
  sgrp("Hispanic", gi$race_eth == "Hispanic"),
  sgrp("BMI < 30", gi$bmi < 30),
  sgrp("BMI >= 30", gi$bmi >= 30),
  sgrp("CRC", gi$gi_site %in% c("Colon","Rectum")),
  sgrp("Non-CRC GI", !gi$gi_site %in% c("Colon","Rectum"))
))
sg <- sg %>% mutate(HR_CI=sprintf("%.2f (%.2f-%.2f)", HR, Lower, Upper),
                     Subgroup=factor(Subgroup, levels=rev(sg$Subgroup)))
write.csv(sg, file.path(RES_DIR, "subgroup_forest.csv"), row.names=FALSE)

p_sg <- ggplot(sg, aes(x=HR, y=Subgroup)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey60") +
  geom_point(aes(color=P<0.05), size=3.5) +
  geom_errorbarh(aes(xmin=Lower, xmax=Upper, color=P<0.05), height=0.15, linewidth=1.2) +
  geom_text(aes(label=HR_CI), x=max(sg$Upper)*1.3, hjust=0, size=3, color="grey40") +
  scale_color_manual(values=c("TRUE"="#D62728","FALSE"="#1F77B4"), guide="none") +
  scale_x_log10() + labs(title="PNI Subgroup Analysis", x="HR per 1-SD (95% CI)", y="") +
  theme_pub(12) + coord_cartesian(xlim=c(0.2, max(sg$Upper)*1.8))
ggsave(file.path(FIG_DIR, "fig_subgroup_forest.pdf"), p_sg, width=10, height=6)
ggsave(file.path(FIG_DIR, "fig_subgroup_forest.png"), p_sg, width=10, height=6, dpi=300)
cat(sprintf("  %d subgroups\n", nrow(sg)))

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 8 — PNI decomposition (component vs composite)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A8: PNI decomposition ===\n")
gi <- gi %>% mutate(alb_s = as.numeric(scale(albumin_gdl)),
                     lym_s = as.numeric(scale(lymph_abs)))

comp <- bind_rows(lapply(c("alb_s","lym_s","PNI_s"), function(idx) {
  run_cox(gi, idx, adj_cov)
}))
comp <- comp %>% mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, Lower, Upper))
write.csv(comp, file.path(RES_DIR, "pni_decomposition.csv"), row.names=FALSE)
print(comp)

# C-statistic comparison
cstat <- function(fit) fit$concordance[["concordance"]]
f_alb <- coxph(Surv(surv_years, death) ~ alb_s + age + sex_b + race_eth_b, data=gi)
f_lym <- coxph(Surv(surv_years, death) ~ lym_s + age + sex_b + race_eth_b, data=gi)
f_pni <- coxph(Surv(surv_years, death) ~ PNI_s + age + sex_b + race_eth_b, data=gi)
f_base <- coxph(Surv(surv_years, death) ~ age + sex_b + race_eth_b, data=gi)
cat(sprintf("  C-stat: Base=%.3f, +Albumin=%.3f, +Lymph=%.3f, +PNI=%.3f\n",
    cstat(f_base), cstat(f_alb), cstat(f_lym), cstat(f_pni)))

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 9 — CRP × PNI interaction
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A9: CRP interaction ===\n")
if ("crp_mgL" %in% names(gi)) {
  gi_crp <- gi %>% filter(!is.na(crp_mgL))
  cat(sprintf("  CRP avail: %d GI, %d events\n", nrow(gi_crp), sum(gi_crp$death)))
  if (sum(gi_crp$death) >= 10) {
    gi_crp <- gi_crp %>% mutate(crp_s = as.numeric(scale(log(crp_mgL + 1))))
    f_crp <- coxph(Surv(surv_years, death) ~ PNI_s * crp_s + age + sex_b + race_eth_b, data=gi_crp)
    crp_res <- tidy(f_crp, conf.int=TRUE) %>% mutate(HR=exp(estimate), HR_CI=sprintf("%.3f (%.3f-%.3f)", HR, exp(conf.low), exp(conf.high)))
    write.csv(crp_res, file.path(RES_DIR, "crp_interaction.csv"), row.names=FALSE)
    print(crp_res[, c("term","HR_CI","p.value")])
  }
} else { cat("  CRP not in data\n") }

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 10 — Age × PNI interaction
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== A10: Age interaction ===\n")
f_age_int <- coxph(Surv(surv_years, death) ~ PNI_s * age_s + sex_b + race_eth_b, data=gi)
hr_age <- tidy(f_age_int, conf.int=TRUE) %>% filter(grepl(":", term)) %>%
  mutate(HR_CI=sprintf("%.3f (%.3f-%.3f)", exp(estimate), exp(conf.low), exp(conf.high)))
cat(sprintf("  PNI × age interaction: %s, p=%.4f\n", hr_age$HR_CI, hr_age$p.value))

# ═══════════════════════════════════════════════════════════════════════════
# FIGURES — KM + forest (main)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== Figures ===\n")
km <- survfit(Surv(surv_years, death) ~ pni_t, data=gi)

# Fig 1: KM
p_km <- ggsurvplot(km, data=gi,
  title="All-Cause Survival by PNI Tertile (GI Tumor, NHANES 1988-2016)",
  xlab="Years", ylab="Survival Probability", legend.title="PNI Tertile",
  pval=TRUE, pval.coord=c(0,0.15), conf.int=TRUE, risk.table=TRUE,
  risk.table.height=0.25, palette=c("#D62728","#FF7F0E","#2CA02C"),
  surv.median.line="hv", ggtheme=theme_classic(base_size=12))
pdf(file.path(FIG_DIR, "fig1_km.pdf"), width=10, height=8); print(p_km); dev.off()
png(file.path(FIG_DIR, "fig1_km.png"), width=10, height=8, units="in", res=300); print(p_km); dev.off()

# Fig 2: Forest
cox_adj <- cox1 %>% filter(Adj=="Adjusted") %>%
  mutate(Nutrition=factor(Index, levels=rev(c("PNI","CONUT","GNRI"))))
p_f <- ggplot(cox_adj, aes(x=HR, y=Nutrition)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey60") +
  geom_point(aes(color=P<0.05), size=4) +
  geom_errorbarh(aes(xmin=Lower, xmax=Upper, color=P<0.05), height=0.15, linewidth=1.3) +
  scale_color_manual(values=c("TRUE"="#D62728","FALSE"="#1F77B4"), guide="none") +
  scale_x_log10(breaks=c(0.4,0.6,0.8,1.0)) +
  labs(title="Nutrition Indices and All-Cause Mortality",
       subtitle="Adjusted HR per 1-SD | GI Tumor Patients",
       x="Hazard Ratio (95% CI)", y="") +
  geom_text(aes(label=HR_CI), x=1.4, hjust=0, size=3.2, color="grey30") +
  coord_cartesian(xlim=c(0.4, 2.0)) + theme_classic(base_size=12)
ggsave(file.path(FIG_DIR, "fig2_forest.pdf"), p_f, width=9, height=4)
ggsave(file.path(FIG_DIR, "fig2_forest.png"), p_f, width=9, height=4, dpi=300)

cat("All done.\n")
