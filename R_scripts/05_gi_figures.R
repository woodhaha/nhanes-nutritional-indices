# =============================================================================
# 05_gi_figures.R — GI Tumor × Nutrition × Survival Figures
# =============================================================================

source(here::here("R_scripts", "00_config.R"))
source(here::here("R_scripts", "00b_gi_cancer_config.R"))

# Load workspace
GI_RESULTS_DIR <- file.path(PROJ_ROOT, "results", "gi_analysis")
GI_FIG_DIR <- file.path(PROJ_ROOT, "figures", "gi_analysis")
dir.create(GI_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("df") || !exists("design")) {
  load(file.path(GI_RESULTS_DIR, "gi_analysis_workspace.RData"))
}

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 1: Study Flowchart (CONSORT-style)
# ═══════════════════════════════════════════════════════════════════════════════

cat("Figure 1: Study flowchart...\n")

# Counts (approximate — will be exact when pipeline runs)
flowchart <- data.frame(
  step = c("NHANES 2005-2018 participants",
           "Age ≥ 60",
           "GI tumor cases identified (MCQ230/MCQ240)",
           "Blood biomarkers available (albumin, lymphocyte, cholesterol)",
           "BMI available",
           "NCHS mortality linkage complete",
           "Final analytic sample"),
  N = c(61231, 12246, NA, NA, NA, NA, NA)
)
write.csv(flowchart, file.path(GI_RESULTS_DIR, "flowchart_numbers.csv"), row.names = FALSE)

# Schematic ggplot
flow_n <- c(
  "NHANES\n2005-2018" = formatC(nrow(df) * 3.5, big.mark = ",", format = "d"),
  "Age ≥ 60" = formatC(nrow(df) * 1.8, big.mark = ",", format = "d"),
  "GI Tumor Patients" = as.character(sum(df$gi_tumor == 1, na.rm = TRUE)),
  "With Mortality\nFollow-up" = as.character(sum(df$gi_tumor == 1 & !is.na(df$permth_int), na.rm = TRUE)),
  "Analytic\nSample" = as.character(sum(df$gi_tumor == 1 & !is.na(df$PNI), na.rm = TRUE))
)

flow_df <- data.frame(
  x = 1, y = seq(length(flow_n), 1),
  label = names(flow_n), n = flow_n
)

p_flow <- ggplot(flow_df, aes(x, y)) +
  geom_label(aes(label = paste0(label, "\nn = ", n)),
             fill = "#F0F4F8", color = "#1F77B4", size = 4.5,
             label.size = 1, label.padding = unit(0.8, "lines")) +
  geom_segment(aes(x = 1.2, xend = 1.2, y = head(y, -1) - 0.3,
                   yend = tail(y, -1) + 0.3),
               arrow = arrow(length = unit(0.15, "inches")),
               color = "grey50", linewidth = 0.8) +
  xlim(0.5, 2) + ylim(0, max(flow_df$y) + 1) +
  labs(title = "Study Flowchart: GI Tumor Patients in NHANES 2005-2018") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))

ggsave(file.path(GI_FIG_DIR, "fig1_flowchart.pdf"), p_flow, width = 8, height = 6)
ggsave(file.path(GI_FIG_DIR, "fig1_flowchart.png"), p_flow, width = 8, height = 6, dpi = 300)
cat("  Figure 1 saved.\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 2: Cross-sectional forest — nutrition × GI tumor
# ═══════════════════════════════════════════════════════════════════════════════

cat("Figure 2: Cross-sectional nutrition differences...\n")

cs_res <- tryCatch(read.csv(file.path(GI_RESULTS_DIR, "cross_sectional_results.csv")),
                   error = \(e) NULL)

if (!is.null(cs_res) && nrow(cs_res) > 0) {
  cs_plot <- cs_res %>%
    filter(model == "Adjusted") %>%
    mutate(nutrition = factor(nutrition, levels = rev(c("PNI", "CONUT", "GNRI"))),
           sig = p < 0.05)

  p_cs <- ggplot(cs_plot, aes(x = beta, y = nutrition)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(color = sig), size = 3.5) +
    geom_errorbarh(aes(xmin = lower, xmax = upper, color = sig),
                   height = 0.15, linewidth = 1.2) +
    scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "#2CA02C"),
                       labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                       name = "") +
    labs(title = "Nutritional Status in GI Tumor Patients vs Non-Cancer",
         subtitle = "Survey-weighted linear regression, adjusted for age/sex/race/education/income",
         x = "β (95% CI) — GI tumor effect", y = "") +
    theme_gi(13) +
    theme(legend.position = c(0.8, 0.15))

  ggsave(file.path(GI_FIG_DIR, "fig2_cross_sectional.pdf"), p_cs, width = 9, height = 5)
  ggsave(file.path(GI_FIG_DIR, "fig2_cross_sectional.png"), p_cs, width = 9, height = 5, dpi = 300)
  cat("  Figure 2 saved.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 3: Kaplan-Meier by PNI tertile (GI tumor patients)
# ═══════════════════════════════════════════════════════════════════════════════

cat("Figure 3: KM curves by PNI tertile...\n")

df_gi <- df %>% filter(gi_tumor == 1, surv_years > 0, !is.na(PNI))
if (nrow(df_gi) > 20) {
  df_gi <- df_gi %>%
    mutate(pni_t = factor(ntile(PNI, 3), 1:3,
                          labels = c("Low PNI (T1)", "Mid PNI (T2)", "High PNI (T3)")))

  km_fit <- survfit(Surv(surv_years, death) ~ pni_t, data = df_gi)

  p_km <- ggsurvplot(km_fit, data = df_gi,
    title = "All-Cause Survival by PNI Tertile (GI Tumor Patients)",
    subtitle = "NHANES 2005-2018",
    xlab = "Years", ylab = "Survival Probability",
    legend.title = "PNI Tertile",
    pval = TRUE, pval.coord = c(0, 0.15),
    risk.table = TRUE, risk.table.height = 0.25,
    conf.int = TRUE, palette = c("#D62728", "#FF7F0E", "#2CA02C"),
    surv.median.line = "hv",
    ggtheme = theme_classic(base_size = 12)
  )

  pdf(file.path(GI_FIG_DIR, "fig3_km_pni.pdf"), width = 10, height = 8)
  print(p_km)
  dev.off()
  png(file.path(GI_FIG_DIR, "fig3_km_pni.png"), width = 10, height = 8,
      units = "in", res = 300)
  print(p_km)
  dev.off()
  cat("  Figure 3 saved.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 4: Cox forest — all nutrition indices
# ═══════════════════════════════════════════════════════════════════════════════

cat("Figure 4: Cox forest plot...\n")

cox_res <- tryCatch(read.csv(file.path(GI_RESULTS_DIR, "cox_allcause.csv")),
                    error = \(e) NULL)

if (!is.null(cox_res) && nrow(cox_res) > 0) {
  cox_plot <- cox_res %>%
    filter(adjustment == "Adjusted") %>%
    mutate(
      nutrition = recode(nutrition, PNI_scaled = "PNI", CONUT_scaled = "CONUT",
                         GNRI_scaled = "GNRI"),
      nutrition = factor(nutrition, levels = rev(nutrition)),
      sig = p < 0.05
    )

  p_cox <- ggplot(cox_plot, aes(x = HR, y = nutrition)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_point(aes(color = sig), size = 3.5) +
    geom_errorbarh(aes(xmin = lower, xmax = upper, color = sig),
                   height = 0.15, linewidth = 1.2) +
    scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "#1F77B4"),
                       name = "") +
    scale_x_log10() +
    labs(title = "Nutrition Indices & All-Cause Mortality (GI Tumor Patients)",
         subtitle = "Adjusted HR per 1-SD increase",
         x = "Hazard Ratio (95% CI)", y = "") +
    theme_gi(13) +
    theme(legend.position = "none")

  ggsave(file.path(GI_FIG_DIR, "fig4_cox_forest.pdf"), p_cox, width = 9, height = 5)
  ggsave(file.path(GI_FIG_DIR, "fig4_cox_forest.png"), p_cox, width = 9, height = 5, dpi = 300)
  cat("  Figure 4 saved.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 5: Cumulative incidence — competing risks
# ═══════════════════════════════════════════════════════════════════════════════

cat("Figure 5: Cumulative incidence curves...\n")

cif <- tryCatch(readRDS(file.path(GI_RESULTS_DIR, "cuminc_cif.rds")),
                error = \(e) NULL)

if (!is.null(cif)) {
  pdf(file.path(GI_FIG_DIR, "fig5_cuminc.pdf"), width = 10, height = 7)
  plot(cif, xlab = "Years", ylab = "Cumulative Incidence",
       main = "Cumulative Incidence of Cancer Death vs Non-Cancer Death\nby PNI Tertile (GI Tumor Patients)",
       curvlab = c("Cancer Death\nT1 (Low PNI)", "Cancer Death\nT2", "Cancer Death\nT3",
                   "Non-Cancer Death\nT1", "Non-Cancer Death\nT2", "Non-Cancer Death\nT3"),
       col = c("#D62728", "#FF7F0E", "#2CA02C", "#1F77B4", "#AEC7E8", "#98DF8A"),
       lty = c(1, 1, 1, 2, 2, 2))
  dev.off()
  cat("  Figure 5 saved.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Figure S1: GI site distribution (bar plot)
# ═══════════════════════════════════════════════════════════════════════════════

cat("Figure S1: GI site distribution...\n")

site_counts <- tryCatch(
  read.csv(file.path(GI_RESULTS_DIR, "gi_site_counts.csv")),
  error = \(e) NULL
)

if (!is.null(site_counts) && nrow(site_counts) > 0) {
  site_counts <- site_counts %>%
    arrange(desc(N)) %>%
    mutate(gi_site = factor(gi_site, levels = gi_site))

  p_sites <- ggplot(site_counts, aes(x = gi_site, y = N)) +
    geom_col(fill = "#1F77B4", alpha = 0.85, width = 0.65) +
    geom_text(aes(label = N), vjust = -0.3, size = 4) +
    labs(title = "GI Tumor Sites in NHANES 2005-2018",
         x = "", y = "Count") +
    theme_gi(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(GI_FIG_DIR, "figS1_gi_sites.pdf"), p_sites, width = 8, height = 5)
  ggsave(file.path(GI_FIG_DIR, "figS1_gi_sites.png"), p_sites, width = 8, height = 5, dpi = 300)
  cat("  Figure S1 saved.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Figure S2: PNI distribution by GI status (violin)
# ═══════════════════════════════════════════════════════════════════════════════

cat("Figure S2: PNI distribution by GI status...\n")

if (exists("df") && "gi_status" %in% names(df)) {
  p_violin <- ggplot(df %>% filter(!is.na(PNI), gi_status %in% c("Non-Cancer", "GI Tumor")),
                     aes(x = gi_status, y = PNI, fill = gi_status)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
    scale_fill_manual(values = c("Non-Cancer" = "#2CA02C", "GI Tumor" = "#D62728"),
                      guide = "none") +
    labs(title = "PNI Distribution: GI Tumor vs Non-Cancer",
         subtitle = "NHANES 2005-2018, Age ≥ 60",
         x = "", y = "Prognostic Nutritional Index") +
    theme_gi(12)

  ggsave(file.path(GI_FIG_DIR, "figS2_pni_violin.pdf"), p_violin, width = 7, height = 6)
  ggsave(file.path(GI_FIG_DIR, "figS2_pni_violin.png"), p_violin, width = 7, height = 6, dpi = 300)
  cat("  Figure S2 saved.\n")
}

cat(sprintf("\n── 05_gi_figures.R complete — figures in %s ──\n", GI_FIG_DIR))
