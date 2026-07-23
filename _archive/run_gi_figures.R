# run_gi_figures.R — GI Tumor survival figures (self-contained)
# Rscript run_gi_figures.R

library(ggplot2); library(survminer); library(survival); library(dplyr)

GI_RESULTS_DIR <- "D:/Researching/results/gi_analysis"
GI_FIG_DIR <- "D:/Researching/figures/gi_analysis"
dir.create(GI_FIG_DIR, recursive=TRUE, showWarnings=FALSE)

theme_gi <- function(base_size=12) {
  theme_classic(base_size=base_size) +
    theme(panel.grid.major.y=element_line(color="grey90", linewidth=0.3),
          plot.title=element_text(face="bold", size=base_size+2),
          plot.subtitle=element_text(color="grey40"),
          legend.position="bottom", axis.title=element_text(size=base_size))
}

# ── Load data ────────────────────────────────────────────────────────────────
cat("Loading data...\n")
df <- readRDS("D:/Researching/data/gi_analysis/nhanes_gi_nutrition_raw.rds")
mort <- readRDS("D:/Researching/data/gi_analysis/mort_2019.rds")
df <- merge(df, mort[,c("SEQN","eligstat","mortstat","permth_int")], by="SEQN", all.x=TRUE)
df <- df[df$eligstat==1 & !is.na(df$mortstat) & !is.na(df$PNI),]
df$surv_years <- df$permth_int/12; df$death <- df$mortstat
df_gi <- df[df$gi_tumor==1 & df$surv_years>0,]
df_gi$pni_t <- factor(ntile(df_gi$PNI, 3), 1:3, c("Low (T1)","Mid (T2)","High (T3)"))
df_gi$sex_b <- ifelse(df_gi$sex=="Female",1L,0L)
df_gi$race_eth_b <- ifelse(df_gi$race_eth=="Non-Hispanic White",1L,0L)
df_gi$PNI_s <- as.numeric(scale(df_gi$PNI))
df_gi$CONUT_s <- as.numeric(scale(-df_gi$CONUT))
df_gi$GNRI_s <- as.numeric(scale(df_gi$GNRI))
cat(sprintf("GI tumor N=%d, events=%d\n", nrow(df_gi), sum(df_gi$death)))

# ── 1. KM curves by PNI tertile ──────────────────────────────────────────────
cat("Figure 1: KM curves...\n")
km <- survfit(Surv(surv_years, death) ~ pni_t, data=df_gi)

p_km <- ggsurvplot(km, data=df_gi,
  title="All-Cause Survival by PNI Tertile (GI Tumor Patients, NHANES 2005-2016)",
  xlab="Years", ylab="Survival Probability",
  legend.title="PNI Tertile",
  pval=TRUE, pval.coord=c(0, 0.15), conf.int=TRUE,
  risk.table=TRUE, risk.table.height=0.25,
  palette=c("#D62728","#FF7F0E","#2CA02C"),
  surv.median.line="hv",
  ggtheme=theme_classic(base_size=12))
pdf(file.path(GI_FIG_DIR, "fig1_km_pni.pdf"), width=10, height=8); print(p_km); dev.off()
png(file.path(GI_FIG_DIR, "fig1_km_pni.png"), width=10, height=8, units="in", res=300); print(p_km); dev.off()
saveRDS(list(km=km, data=df_gi), file.path(GI_RESULTS_DIR, "km_fit.rds"))

# ── 2. Cox forest plot ──────────────────────────────────────────────────────
cat("Figure 2: Cox forest plot...\n")
cox <- read.csv(file.path(GI_RESULTS_DIR, "cox_allcause.csv"))
cox_plot <- cox[cox$Adjustment=="Adjusted",]
cox_plot$Nutrition <- factor(cox_plot$Nutrition, levels=rev(c("PNI","CONUT","GNRI")))
cox_plot$sig <- cox_plot$P < 0.05

p_forest <- ggplot(cox_plot, aes(x=HR, y=Nutrition)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey60", linewidth=0.5) +
  geom_point(aes(color=sig), size=4) +
  geom_errorbarh(aes(xmin=Lower, xmax=Upper, color=sig), height=0.15, linewidth=1.3) +
  scale_color_manual(values=c("TRUE"="#D62728","FALSE"="#1F77B4"), name="") +
  scale_x_log10(breaks=c(0.4, 0.6, 0.8, 1.0)) +
  labs(title="Nutrition Indices & All-Cause Mortality (GI Tumor Patients)",
       subtitle="Adjusted HR per 1-SD increase | NHANES 2005-2016, Age ≥ 60",
       x="Hazard Ratio (95% CI)", y="") +
  geom_text(aes(label=sprintf("HR=%.2f (%.2f-%.2f), p=%.4f", HR, Lower, Upper, P)),
            x=1.4, hjust=0, size=3.2, color="grey30") +
  coord_cartesian(xlim=c(0.4, 2.5)) +
  theme_gi(13) + theme(legend.position="none")
ggsave(file.path(GI_FIG_DIR, "fig2_cox_forest.pdf"), p_forest, width=10, height=5)
ggsave(file.path(GI_FIG_DIR, "fig2_cox_forest.png"), p_forest, width=10, height=5, dpi=300)

# ── 3. Time-dependent HR ────────────────────────────────────────────────────
cat("Figure 3: Time-dependent HR...\n")
td <- read.csv(file.path(GI_RESULTS_DIR, "cox_timedependent.csv"))
td_plot <- data.frame()
for (nut in c("PNI","CONUT","GNRI")) {
  sub <- td[grep(nut, td$Nutrition),]
  td_plot <- rbind(td_plot, data.frame(
    Nutrition=nut,
    Period=factor(c("0-2 yr","2-5 yr","5+ yr"), levels=c("0-2 yr","2-5 yr","5+ yr")),
    HR=c(sub$HR[1], sub$HR[1]*sub$HR[2], sub$HR[1]*sub$HR[3]),
    Lower=c(sub$Lower[1], NA, NA), Upper=c(sub$Upper[1], NA, NA)))
}
td_plot$Nutrition <- factor(td_plot$Nutrition, levels=c("PNI","CONUT","GNRI"))

p_td <- ggplot(td_plot, aes(x=Period, y=HR, color=Nutrition, group=Nutrition)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey60") +
  geom_line(linewidth=1, position=position_dodge(0.2)) +
  geom_point(size=3.5, position=position_dodge(0.2)) +
  scale_color_manual(values=c("PNI"="#D62728","CONUT"="#FF7F0E","GNRI"="#2CA02C")) +
  labs(title="Time-Varying Effect of Nutrition on Mortality",
       subtitle="HR per 1-SD increase, by follow-up period | GI tumor patients",
       x="Follow-up Period", y="Hazard Ratio (95% CI)", color="Index") +
  geom_text(data=subset(td_plot, Nutrition=="PNI"),
            aes(label=sprintf("HR=%.2f", HR)), vjust=-1.2, size=3) +
  theme_gi(13)
ggsave(file.path(GI_FIG_DIR, "fig3_timedependent.pdf"), p_td, width=9, height=5.5)
ggsave(file.path(GI_FIG_DIR, "fig3_timedependent.png"), p_td, width=9, height=5.5, dpi=300)

# ── 4. Landmark ─────────────────────────────────────────────────────────────
cat("Figure 4: Landmark...\n")
lm <- read.csv(file.path(GI_RESULTS_DIR, "cox_landmark.csv"))
lm$Landmark <- factor(lm$Landmark, levels=c("1-yr","3-yr","5-yr"))
lm$Nutrition <- factor(lm$Nutrition, levels=c("PNI","CONUT","GNRI"))
lm$sig <- lm$P < 0.05

p_lm <- ggplot(lm, aes(x=Landmark, y=HR, color=Nutrition, group=Nutrition)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey60") +
  geom_line(linewidth=1, position=position_dodge(0.3)) +
  geom_point(aes(shape=sig), size=3.5, position=position_dodge(0.3)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.15, position=position_dodge(0.3)) +
  scale_color_manual(values=c("PNI"="#D62728","CONUT"="#FF7F0E","GNRI"="#2CA02C")) +
  scale_shape_manual(values=c("TRUE"=16, "FALSE"=1), guide="none") +
  labs(title="Landmark Analysis: Conditional Survival",
       subtitle="GI tumor patients surviving to landmark | Adjusted HR per 1-SD",
       x="Landmark Time", y="Hazard Ratio (95% CI)", color="Index") +
  theme_gi(13)
ggsave(file.path(GI_FIG_DIR, "fig4_landmark.pdf"), p_lm, width=9, height=5.5)
ggsave(file.path(GI_FIG_DIR, "fig4_landmark.png"), p_lm, width=9, height=5.5, dpi=300)

# ── 5. Competing risks ──────────────────────────────────────────────────────
cat("Figure 5: Competing risks...\n")
cr <- read.csv(file.path(GI_RESULTS_DIR, "cox_causespecific.csv"))
cr$Label <- paste0(cr$Nutrition, ": ", cr$Cause)
cr$Label <- factor(cr$Label, levels=rev(unique(cr$Label)))
cr$sig <- cr$P < 0.05

p_cr <- ggplot(cr, aes(x=HR, y=Label)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey60") +
  geom_point(aes(color=sig), size=3.5) +
  geom_errorbarh(aes(xmin=Lower, xmax=Upper, color=sig), height=0.15, linewidth=1.2) +
  geom_text(aes(label=sprintf("HR=%.2f, p=%.3f", HR, P)), x=2.2, hjust=0, size=3) +
  scale_color_manual(values=c("TRUE"="#D62728","FALSE"="#1F77B4"), name="") +
  scale_x_log10(limits=c(0.4, 3.5)) +
  labs(title="Competing Risks: GI Cancer Death vs Non-Cancer Death",
       subtitle="GI tumor patients | Per 1-SD increase in nutrition index",
       x="Hazard Ratio (95% CI)", y="") +
  theme_gi(12) + theme(legend.position="none")
ggsave(file.path(GI_FIG_DIR, "fig5_competing_risks.pdf"), p_cr, width=10, height=5)
ggsave(file.path(GI_FIG_DIR, "fig5_competing_risks.png"), p_cr, width=10, height=5, dpi=300)

# ── 6. Cross-sectional boxplot (bonus) ──────────────────────────────────────
cat("Figure 6: PNI boxplot by GI status...\n")
df_plot <- df[df$gi_status %in% c("GI Tumor","Non-Cancer"),]
df_plot$gi_status <- factor(df_plot$gi_status, levels=c("Non-Cancer","GI Tumor"))

p_box <- ggplot(df_plot, aes(x=gi_status, y=PNI, fill=gi_status)) +
  geom_violin(alpha=0.5, trim=FALSE) +
  geom_boxplot(width=0.2, fill="white", outlier.size=0.5) +
  scale_fill_manual(values=c("Non-Cancer"="#2CA02C","GI Tumor"="#D62728"), guide="none") +
  labs(title="PNI Distribution: GI Tumor vs Non-Cancer",
       subtitle="NHANES 2005-2016, Age ≥ 60",
       x="", y="Prognostic Nutritional Index") +
  theme_gi(13) +
  stat_summary(geom="text", fun=median, aes(label=sprintf("Median=%.1f", ..y..)),
               vjust=-1.5, size=3.2, color="grey30")
ggsave(file.path(GI_FIG_DIR, "fig6_pni_boxplot.pdf"), p_box, width=7, height=6)
ggsave(file.path(GI_FIG_DIR, "fig6_pni_boxplot.png"), p_box, width=7, height=6, dpi=300)

cat(sprintf("\nAll 6 figures saved to: %s\n", GI_FIG_DIR))
cat(sprintf("  %d PDF + %d PNG\n",
    length(list.files(GI_FIG_DIR, pattern="\\.pdf$")),
    length(list.files(GI_FIG_DIR, pattern="\\.png$"))))
cat("Done.\n")
