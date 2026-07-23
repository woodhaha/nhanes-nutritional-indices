# run_manuscript_tables.R — Final tables & figures for manuscript
# Rscript run_manuscript_tables.R

library(dplyr); library(ggplot2); library(survival); library(broom)
library(tidyr); library(survminer)

PROJ <- "D:/Researching/NHANES_aged_GI_tumor_nutrition"
DATA_DIR <- file.path(PROJ, "data/gi_analysis")
RES_DIR <- file.path(PROJ, "results/gi_analysis")
FIG_DIR <- file.path(PROJ, "figures/gi_analysis")

theme_pub <- function(base_size=12) {
  theme_classic(base_size=base_size) +
    theme(panel.grid.major.y=element_line(color="grey90", linewidth=0.3),
          plot.title=element_text(face="bold", size=base_size+2),
          plot.subtitle=element_text(color="grey40", size=base_size-1),
          legend.position="bottom", axis.title=element_text(size=base_size))
}

# ── Load combined data (NHANES III + 2005-2016) ──────────────────────────
cat("=== Loading data ===\n")
df_comb <- readRDS(file.path(RES_DIR, "nhanes_combined_gi_fixed.rds"))
cat(sprintf("Combined: N=%d, GI=%d\n", nrow(df_comb), sum(df_comb[["gi_tumor"]], na.rm=TRUE)))

# Also load HEI data for panels that need it
df_hei <- readRDS(file.path(RES_DIR, "nhanes_combined_gi_with_hei.rds"))

# ═══════════════════════════════════════════════════════════════════════════
# TABLE 2: Baseline characteristics by PNI tertile (GI tumor patients only)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== Table 2: Baseline by PNI Tertile ===\n")
df_gi_comb <- df_comb %>% filter(gi_tumor==1, surv_years>0)

df_gi_comb$pni_t <- factor(ntile(df_gi_comb$PNI, 3), 1:3,
                            c("Low (T1)", "Mid (T2)", "High (T3)"))

# Function: mean (SD) or n (%)
fmt_mean <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
fmt_npct <- function(x) sprintf("%d (%.1f%%)", sum(x, na.rm=TRUE), mean(x, na.rm=TRUE)*100)

t2 <- df_gi_comb %>%
  group_by(pni_t) %>%
  summarise(
    N = n(),
    Age = fmt_mean(age),
    Female = fmt_npct(sex=="Female"),
    NH_White = fmt_npct(race_eth=="Non-Hispanic White"),
    Edu_AboveHigh = fmt_npct(edu_binary==1),
    PNI = fmt_mean(PNI),
    CONUT = fmt_mean(CONUT),
    GNRI = fmt_mean(GNRI),
    Albumin_gdl = fmt_mean(albumin_gdl),
    Lymphocyte = fmt_mean(lymph_abs),
    BMI = fmt_mean(bmi),
    Deaths = fmt_npct(death==1),
    Surv_Years = fmt_mean(surv_years),
    .groups="drop"
  )
cat("--- Table 2 ---\n")
print(as.data.frame(t2))
write.csv(t2, file.path(RES_DIR, "table2_baseline_pni_tertile.csv"), row.names=FALSE)

# Also compute overall column
t2_overall <- df_gi_comb %>%
  summarise(
    N = n(), Age = fmt_mean(age), Female = fmt_npct(sex=="Female"),
    NH_White = fmt_npct(race_eth=="Non-Hispanic White"),
    Edu_AboveHigh = fmt_npct(edu_binary==1),
    PNI = fmt_mean(PNI), CONUT = fmt_mean(CONUT), GNRI = fmt_mean(GNRI),
    Albumin_gdl = fmt_mean(albumin_gdl), Lymphocyte = fmt_mean(lymph_abs),
    BMI = fmt_mean(bmi), Deaths = fmt_npct(death==1), Surv_Years = fmt_mean(surv_years))
t2_overall$pni_t <- "Overall"
t2_overall <- t2_overall[, names(t2)]
t2_all <- bind_rows(t2, t2_overall)
write.csv(t2_all, file.path(RES_DIR, "table2_baseline_pni_tertile_with_overall.csv"), row.names=FALSE)
cat("  Saved to table2_baseline_pni_tertile.csv\n")

# ═══════════════════════════════════════════════════════════════════════════
# SENSITIVITY: Excluding CRC
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== Sensitivity: Excluding CRC ===\n")

CRC_SITES <- c("Colon", "Rectum", "Colon/Rectum")
df_nocrc <- df_comb %>% filter(gi_tumor==1, surv_years>0,
                                !(gi_site %in% CRC_SITES))
n_nocrc <- nrow(df_nocrc)
e_nocrc <- sum(df_nocrc$death, na.rm=TRUE)
cat(sprintf("Non-CRC GI tumors: N=%d, events=%d\n", n_nocrc, e_nocrc))

if (e_nocrc >= 5) {
  df_nocrc <- df_nocrc %>% mutate(
    PNI_s = as.numeric(scale(PNI)),
    CONUT_s = as.numeric(scale(-CONUT)),
    GNRI_s = as.numeric(scale(GNRI)),
    sex_b = ifelse(sex=="Female", 1L, 0L),
    race_eth_b = ifelse(race_eth=="Non-Hispanic White", 1L, 0L))

  sens_res <- data.frame()
  for (exp in c("PNI_s", "CONUT_s", "GNRI_s")) {
    for (adj in c("Crude", "Adjusted")) {
      covs <- if (adj=="Crude") "1" else "age + sex_b + race_eth_b + edu_binary"
      f <- as.formula(paste("Surv(surv_years, death) ~", exp, "+", covs))
      fit <- tryCatch(coxph(f, data=df_nocrc), error=function(e) NULL)
      if (is.null(fit)) next
      hr <- tidy(fit, conf.int=TRUE) %>% filter(term==exp)
      sens_res <- rbind(sens_res, data.frame(
        Nutrition = gsub("_s","",exp), Adjustment=adj,
        N=n_nocrc, Events=e_nocrc,
        HR=exp(hr$estimate), Lower=exp(hr$conf.low),
        Upper=exp(hr$conf.high), P=hr$p.value))
    }
  }
  sens_res$HR_CI <- sprintf("%.3f (%.3f-%.3f)", sens_res$HR, sens_res$Lower, sens_res$Upper)
  cat("--- Non-CRC Sensitivity ---\n")
  print(sens_res[, c("Nutrition","Adjustment","N","Events","HR_CI","P")])
  write.csv(sens_res, file.path(RES_DIR, "sensitivity_nocrc.csv"), row.names=FALSE)

  # Table: GI site breakdown of the excluded CRC vs remaining
  cat("\nSite breakdown of excluded CRC:\n")
  site_tab <- df_comb %>% filter(gi_tumor==1, surv_years>0) %>%
    mutate(group=ifelse(gi_site %in% CRC_SITES, "CRC", "Non-CRC")) %>%
    group_by(group, gi_site) %>% summarise(N=n(), .groups="drop") %>%
    pivot_wider(names_from=group, values_from=N, values_fill=0)
  print(as.data.frame(site_tab))
  write.csv(site_tab, file.path(RES_DIR, "sensitivity_nocrc_site_breakdown.csv"), row.names=FALSE)
} else {
  cat("  Too few events for analysis.\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE: HEI-2015 vs PNI comparison (side-by-side forest)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== Figure: HEI vs PNI comparison ===\n")

df_hei_cox <- df_hei %>% filter(gi_tumor==1, surv_years>0, !is.na(HEI2015_total)) %>%
  mutate(HEI_s = as.numeric(scale(HEI2015_total)),
         PNI_s = as.numeric(scale(PNI)),
         sex_b = ifelse(sex=="Female", 1L, 0L),
         race_eth_b = ifelse(race_eth=="Non-Hispanic White", 1L, 0L))

n_hei <- nrow(df_hei_cox)
e_hei <- sum(df_hei_cox$death, na.rm=TRUE)
cat(sprintf("HEI subset: N=%d, events=%d\n", n_hei, e_hei))

if (e_hei >= 10) {
  # Run models for each index separately + together
  comp_res <- data.frame()
  for (idx in c("PNI_s", "HEI_s")) {
    f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+ age + sex_b + race_eth_b"))
    fit <- coxph(f, data=df_hei_cox)
    hr <- tidy(fit, conf.int=TRUE) %>% filter(term==idx)
    comp_res <- rbind(comp_res, data.frame(
      Index=gsub("_s","",idx), Model="Separate",
      HR=exp(hr$estimate), Lower=exp(hr$conf.low),
      Upper=exp(hr$conf.high), P=hr$p.value))
  }
  # Both together
  f_both <- Surv(surv_years, death) ~ PNI_s + HEI_s + age + sex_b + race_eth_b
  fit_both <- coxph(f_both, data=df_hei_cox)
  for (idx in c("PNI_s", "HEI_s")) {
    hr <- tidy(fit_both, conf.int=TRUE) %>% filter(term==idx)
    comp_res <- rbind(comp_res, data.frame(
      Index=gsub("_s","",idx), Model="Mutually Adjusted",
      HR=exp(hr$estimate), Lower=exp(hr$conf.low),
      Upper=exp(hr$conf.high), P=hr$p.value))
  }
  comp_res$HR_CI <- sprintf("%.3f (%.3f-%.3f)", comp_res$HR, comp_res$Lower, comp_res$Upper)
  cat("--- HEI vs PNI Cox ---\n")
  print(comp_res[, c("Index","Model","HR_CI","P")])
  write.csv(comp_res, file.path(RES_DIR, "hei_vs_pni_cox.csv"), row.names=FALSE)

  # Plot: side-by-side forest
  comp_plot <- comp_res
  comp_plot$Label <- paste0(comp_plot$Index, " (", comp_plot$Model, ")")
  comp_plot$Label <- factor(comp_plot$Label, levels=rev(unique(comp_plot$Label)))
  comp_plot$sig <- comp_plot$P < 0.05

  p_comp <- ggplot(comp_plot, aes(x=HR, y=Label)) +
    geom_vline(xintercept=1, linetype="dashed", color="grey60", linewidth=0.5) +
    geom_point(aes(color=sig), size=4) +
    geom_errorbarh(aes(xmin=Lower, xmax=Upper, color=sig), height=0.15, linewidth=1.3) +
    geom_text(aes(label=sprintf("HR=%.2f (%.2f-%.2f), p=%.3f", HR, Lower, Upper, P)),
              x=2.2, hjust=0, size=3.2, color="grey30") +
    scale_color_manual(values=c("TRUE"="#D62728","FALSE"="#1F77B4"), name="") +
    scale_x_log10(breaks=c(0.25, 0.50, 0.75, 1.0, 1.5)) +
    coord_cartesian(xlim=c(0.2, 3.0)) +
    labs(title="Dietary vs Biochemical Nutrition in GI Cancer Prognosis",
         subtitle="HEI-2015 (diet quality) vs PNI (albumin+lymphocyte) | Adjusted HR per 1-SD",
         x="Hazard Ratio (95% CI)", y="") +
    theme_pub(12) + theme(legend.position="none")
  ggsave(file.path(FIG_DIR, "fig_hei_vs_pni_forest.pdf"), p_comp, width=9, height=4)
  ggsave(file.path(FIG_DIR, "fig_hei_vs_pni_forest.png"), p_comp, width=9, height=4, dpi=300)
  cat("  Saved HEI vs PNI forest plot.\n")

  # KM by HEI tertile (for comparison with PNI KM)
  df_hei_cox$hei_t <- factor(ntile(df_hei_cox$HEI2015_total, 3), 1:3,
                              c("Low HEI (T1)", "Mid HEI (T2)", "High HEI (T3)"))
  km_hei <- survfit(Surv(surv_years, death) ~ hei_t, data=df_hei_cox)

  p_km_hei <- ggsurvplot(km_hei, data=df_hei_cox,
    title="All-Cause Survival by HEI-2015 Tertile (GI Tumor Patients)",
    subtitle="NHANES 2011-2016 | p=0.81 (adjusted Cox)",
    xlab="Years", ylab="Survival Probability",
    legend.title="HEI-2015 Tertile",
    pval=TRUE, pval.coord=c(0, 0.15), conf.int=TRUE,
    risk.table=TRUE, risk.table.height=0.25,
    palette=c("#D62728","#FF7F0E","#2CA02C"),
    surv.median.line="hv",
    ggtheme=theme_classic(base_size=12))
  pdf(file.path(FIG_DIR, "fig_hei_km.pdf"), width=10, height=8)
  print(p_km_hei); dev.off()
  png(file.path(FIG_DIR, "fig_hei_km.png"), width=10, height=8, units="in", res=300)
  print(p_km_hei); dev.off()
  cat("  Saved HEI KM curves.\n")
} else {
  cat("  Insufficient events for HEI Cox.\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE: HEI-PNI correlation scatter plot
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== Figure: HEI-PNI correlation ===\n")

df_cor <- df_hei %>% filter(!is.na(HEI2015_total), !is.na(PNI))
cat(sprintf("N for correlation: %d\n", nrow(df_cor)))

cor_val <- with(df_cor, cor(HEI2015_total, PNI, use="complete.obs"))
cat(sprintf("Pearson r = %.3f\n", cor_val))

# Also by GI status
cor_gi <- df_cor %>% filter(gi_tumor==1) %>% summarise(r=cor(HEI2015_total, PNI, use="comp")) %>% pull(r)
cor_nc <- df_cor %>% filter(gi_tumor==0) %>% summarise(r=cor(HEI2015_total, PNI, use="comp")) %>% pull(r)
cat(sprintf("  GI tumor: r=%.3f\n  Non-cancer: r=%.3f\n", cor_gi, cor_nc))

p_scatter <- ggplot(df_cor, aes(x=HEI2015_total, y=PNI, color=gi_status)) +
  geom_point(alpha=0.3, size=1.5) +
  geom_smooth(method="lm", se=TRUE, linewidth=1) +
  scale_color_manual(values=c("GI Tumor"="#D62728", "Non-Cancer"="#2CA02C",
                               "Other Cancer"="#1F77B4"),
                     name="") +
  annotate("text", x=max(df_cor$HEI2015_total, na.rm=TRUE)*0.7,
           y=max(df_cor$PNI, na.rm=TRUE)*0.95,
           label=sprintf("Overall r = %.3f", cor_val),
           size=4, fontface="italic") +
  labs(title="Correlation: Dietary Quality (HEI-2015) vs Biochemical Nutrition (PNI)",
       subtitle=sprintf("NHANES 2011-2016, Age ≥ 60 | GI tumor r=%.3f, Non-cancer r=%.3f", cor_gi, cor_nc),
       x="Healthy Eating Index 2015", y="Prognostic Nutritional Index") +
  theme_pub(13)
ggsave(file.path(FIG_DIR, "fig_hei_pni_correlation.pdf"), p_scatter, width=8, height=6)
ggsave(file.path(FIG_DIR, "fig_hei_pni_correlation.png"), p_scatter, width=8, height=6, dpi=300)
cat("  Saved HEI-PNI correlation scatter.\n")

# ═══════════════════════════════════════════════════════════════════════════
# SUPPLEMENTARY: Site-specific breakdown table
# ═══════════════════════════════════════════════════════════════════════════
cat("\n=== Supplementary: GI site breakdown ===\n")
site_summary <- df_comb %>% filter(gi_tumor==1, surv_years>0) %>%
  group_by(gi_site) %>%
  summarise(N=n(), Deaths=sum(death, na.rm=TRUE),
            PNI_mean=sprintf("%.1f (%.1f)", mean(PNI, na.rm=TRUE), sd(PNI, na.rm=TRUE)),
            Surv_median=sprintf("%.1f yr", median(surv_years, na.rm=TRUE)),
            .groups="drop")
print(as.data.frame(site_summary))
write.csv(site_summary, file.path(RES_DIR, "supplementary_site_breakdown.csv"), row.names=FALSE)

cat("\n=== All manuscript tables/figures complete ===\n")
