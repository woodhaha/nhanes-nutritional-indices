# -- fig9_fig10_generate.R
# Figure 9: PNI vs NLR vs SII (ΔC-stat)
# Figure 10: PNI decomposition — Albumin vs PNI vs Lymphocyte (ΔC-stat)
# Values from bootstrap C-stat (2000 reps, bootstrap CI) in 04_revision_sensitivity.R

library(ggplot2)
FIG_DIR <- normalizePath("figures/gi_analysis")

# -- Figure 9 data (from nlr_sii_extract.R + main analysis) --
# Manuscript: +0.017 for PNI, +0.017 for logNLR, +0.014 for logSII
fig9 <- data.frame(
  Index = c("PNI", "logNLR", "logSII"),
  DeltaCstat = c(0.017, 0.017, 0.014)
)

# -- Figure 10 data (from sensitivity_pni_decomposition.csv bootstrap) --
fig10 <- data.frame(
  Index = c("Albumin", "PNI", "Lymphocyte"),
  DeltaCstat = c(0.032, 0.022, 0.002)
)

# === FIGURE 9 ===
p9 <- ggplot(fig9, aes(x = reorder(Index, DeltaCstat), y = DeltaCstat, fill = Index)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_text(aes(label = sprintf("+%.3f", DeltaCstat)), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("PNI" = "#1F77B4", "logNLR" = "#D62728", "logSII" = "#2CA02C"), guide = "none") +
  scale_y_continuous(limits = c(0, 0.035), expand = c(0, 0)) +
  labs(title = expression(Delta*"C-statistic: PNI vs Inflammatory Markers"),
       subtitle = "NHANES 2005-2016, GI cancer (n=258, events=125)",
       x = "", y = expression(Delta*"C-statistic (vs baseline)")) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        axis.text = element_text(size = 13))

ggsave(file.path(FIG_DIR, "fig9_inflammatory_comparison.pdf"), p9, width = 6, height = 5)
ggsave(file.path(FIG_DIR, "fig9_inflammatory_comparison.png"), p9, width = 6, height = 5, dpi = 300)
cat("Figure 9 saved.\n")

# === FIGURE 10 ===
p10 <- ggplot(fig10, aes(x = factor(Index, levels = c("Lymphocyte", "PNI", "Albumin")),
                          y = DeltaCstat, fill = Index)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_text(aes(label = sprintf("+%.3f", DeltaCstat)), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Albumin" = "#2CA02C", "PNI" = "#1F77B4", "Lymphocyte" = "#D62728"), guide = "none") +
  scale_y_continuous(limits = c(0, 0.05), expand = c(0, 0)) +
  labs(title = expression(Delta*"C-statistic: PNI Decomposition"),
       subtitle = "NHANES 2005-2016, GI cancer (n=313, events=169)",
       x = "", y = expression(Delta*"C-statistic (vs baseline)")) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        axis.text = element_text(size = 13))

ggsave(file.path(FIG_DIR, "fig10_decomposition.pdf"), p10, width = 6, height = 5)
ggsave(file.path(FIG_DIR, "fig10_decomposition.png"), p10, width = 6, height = 5, dpi = 300)
cat("Figure 10 saved.\n")
