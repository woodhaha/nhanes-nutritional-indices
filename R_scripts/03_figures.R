# ==============================================================================
# 03_figures.R — Publication-quality figures
# ==============================================================================

source(here::here("R_scripts", "00_config.R"))

# Load results if not already in environment
if (!exists("main_results")) {
  load(file.path(RESULTS_DIR, "analysis_workspace.RData"))
}

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 1: Study Flowchart (CONSORT-style)
# ═══════════════════════════════════════════════════════════════════════════════

cat("Generating Figure 1: Study flowchart...\n")

# Will be created as a separate PDF using DiagrammeR or manual ggplot
# For now, save the numbers
flowchart_numbers <- data.frame(
  step = c("NHANES 2011-2014 participants",
           "Age ≥ 60",
           "Cognitive assessment complete (CERAD/AFT/DSST)",
           "PHQ-9 depression screener complete",
           "Dietary recall (Day 1) valid",
           "Blood biomarkers (albumin/lymphocyte/cholesterol) available",
           "Anthropometrics (BMI) available",
           "Final analytic sample"),
  N = c(19931, 3921, 3428, 3389, 3321, 2842, 2815, 2815)
)
write.csv(flowchart_numbers, file.path(RESULTS_DIR, "flowchart_numbers.csv"),
          row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 2: Head-to-head forest plot — Model 3 standardized betas
# ═══════════════════════════════════════════════════════════════════════════════

cat("Generating Figure 2: Head-to-head forest plot...\n")

fig2_data <- head2head %>%
  mutate(
    exposure_label = factor(exposure_label,
      levels = rev(c("E-DII (anti-inflam)", "PNI", "CONUT (reversed)", "GNRI")))
  )

p_head2head <- ggplot(fig2_data, aes(x = beta, y = exposure_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_point(aes(color = p_value < 0.001), size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = p_value < 0.001),
                height = 0.2, linewidth = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "#D62728", "FALSE" = "#2CA02C"),
    labels = c("TRUE" = "p < 0.001", "FALSE" = "p ≥ 0.001"),
    name = ""
  ) +
  labs(
    title = "Association of Nutritional Indices with Cognitive Function",
    subtitle = "Standardized β (95% CI) | NHANES 2011-2014, Age ≥ 60",
    x = "Standardized β (per 1-SD increase in better nutrition)",
    y = "",
    caption = "Model adjusted for: age, sex, race/ethnicity, education, income, BMI,
    smoking, alcohol, physical activity, comorbidity count, CRP"
  ) +
  theme_nhanes(13) +
  theme(
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = "white", color = "grey80"),
    panel.grid.major.y = element_blank()
  )

# Add exact beta annotations
p_head2head <- p_head2head +
  geom_text(aes(
    x = upper + 0.03,
    label = sprintf("β=%.3f\n(%.3f, %.3f)", beta, lower, upper)
  ), hjust = 0, size = 3.2, color = "grey30")

ggsave(file.path(FIG_DIR, "fig2_head2head_forest.pdf"),
       p_head2head, width = 10, height = 6.5, device = "pdf")
ggsave(file.path(FIG_DIR, "fig2_head2head_forest.png"),
       p_head2head, width = 10, height = 6.5, dpi = 300)

cat("  Figure 2 saved.\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 3: RCS dose-response curve
# ═══════════════════════════════════════════════════════════════════════════════

cat("Generating Figure 3: RCS dose-response...\n")

if (exists("pred_df") && nrow(pred_df) > 0) {

  # Rug plot data
  rug_data <- df_cc %>%
    filter(!is.na(E_DII_raw)) %>%
    dplyr::select(E_DII_raw, cog_z_composite)

  p_rcs <- ggplot() +
    # Confidence band
    geom_ribbon(data = pred_df,
                aes(x = E_DII_raw, ymin = y_lower, ymax = y_upper),
                fill = "#1F77B4", alpha = 0.15) +
    # Main curve
    geom_line(data = pred_df,
              aes(x = E_DII_raw, y = y_hat),
              color = "#1F77B4", linewidth = 1.3) +
    # Reference line
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.5) +
    # Ruq (sample distribution)
    geom_rug(data = rug_data, aes(x = E_DII_raw),
             sides = "b", alpha = 0.1, color = "grey40") +
    # Knot locations
    geom_vline(xintercept = knots, linetype = "dotted", color = "#D62728",
               alpha = 0.5, linewidth = 0.5) +
    labs(
      title = "Dose-Response: Dietary Inflammatory Index & Cognitive Function",
      subtitle = sprintf("Restricted cubic spline (3 knots at 10th, 50th, 90th %ile)
      Non-linearity p = %.4f", nonlin_p),
      x = "Energy-Adjusted Dietary Inflammatory Index (E-DII)",
      y = "Predicted Cognitive Composite Z-score",
      caption = "Adjusted for: age, sex, race/ethnicity, education, income, BMI,
      smoking, alcohol, physical activity, comorbidities, CRP\n
      Red dotted lines = knot locations | Shaded band = 95% CI\n
      More negative DII = anti-inflammatory diet | More positive = pro-inflammatory"
    ) +
    annotate("text", x = -2, y = max(pred_df$y_upper, na.rm = TRUE) * 0.9,
             label = "Anti-inflammatory\ndiet →", hjust = 0.5, size = 3.5,
             color = "grey40") +
    annotate("text", x = 3, y = max(pred_df$y_upper, na.rm = TRUE) * 0.9,
             label = "← Pro-inflammatory\ndiet", hjust = 0.5, size = 3.5,
             color = "grey40") +
    theme_nhanes(13)

  ggsave(file.path(FIG_DIR, "fig3_rcs_doseresponse.pdf"),
         p_rcs, width = 10, height = 7, device = "pdf")
  ggsave(file.path(FIG_DIR, "fig3_rcs_doseresponse.png"),
         p_rcs, width = 10, height = 7, dpi = 300)

  cat("  Figure 3 saved.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 4: Subgroup forest plot
# ═══════════════════════════════════════════════════════════════════════════════

cat("Generating Figure 4: Subgroup forest plot...\n")

if (exists("subgroup_results") && nrow(subgroup_results) > 0) {

  # Add interaction p-values
  subgroup_results$interaction_p <- NA_real_
  if (exists("interaction_results") && nrow(interaction_results) > 0) {
    for (i in 1:nrow(subgroup_results)) {
      match_row <- interaction_results[
        interaction_results$subgroup_var == subgroup_results$subgroup_var[i], ]
      if (nrow(match_row) > 0) {
        subgroup_results$interaction_p[i] <- match_row$p_interaction[1]
      }
    }
  }

  # Create display labels
  subgroup_plot_data <- subgroup_results %>%
    mutate(
      display_label = paste(subgroup_var, subgroup_level, sep = ": "),
      display_label = factor(display_label, levels = rev(unique(display_label))),
      sig = ifelse(p_value < 0.05, "p<0.05", "p≥0.05")
    )

  p_subgroup <- ggplot(subgroup_plot_data,
                       aes(x = beta, y = display_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
    geom_point(aes(color = sig, size = N), alpha = 0.85) +
    geom_errorbarh(aes(xmin = lower, xmax = upper, color = sig),
                  height = 0.2, linewidth = 0.9) +
    scale_color_manual(
      values = c("p<0.05" = "#D62728", "p≥0.05" = "#7F7F7F"),
      name = ""
    ) +
    scale_size_continuous(range = c(2, 6), guide = "none") +
    labs(
      title = "Subgroup Analysis: E-DII & Cognitive Function",
      subtitle = "Per 1-SD increase in anti-inflammatory diet | Model 3 covariates",
      x = "Standardized β (95% CI)",
      y = "",
      caption = "Point size proportional to subgroup N. Red = p < 0.05"
    ) +
    theme_nhanes(12) +
    theme(panel.grid.major.y = element_blank())

  ggsave(file.path(FIG_DIR, "fig4_subgroup_forest.pdf"),
         p_subgroup, width = 11, height = 8, device = "pdf")
  ggsave(file.path(FIG_DIR, "fig4_subgroup_forest.png"),
         p_subgroup, width = 11, height = 8, dpi = 300)

  cat("  Figure 4 saved.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 5: Sensitivity analysis tornado plot
# ═══════════════════════════════════════════════════════════════════════════════

cat("Generating Figure 5: Sensitivity analysis...\n")

if (exists("sensitivity_results") && nrow(sensitivity_results) > 0) {

  sens_plot_data <- sensitivity_results %>%
    mutate(
      analysis = factor(analysis, levels = rev(unique(analysis))),
      analysis = forcats::fct_rev(analysis)
    )

  # Add the primary analysis as reference line
  primary_beta <- sens_plot_data$beta[sens_plot_data$analysis == "Primary analysis"]

  p_sens <- ggplot(sens_plot_data, aes(x = beta, y = analysis)) +
    geom_vline(xintercept = primary_beta, linetype = "dashed",
               color = "#1F77B4", linewidth = 0.7, alpha = 0.6) +
    geom_point(aes(size = N), color = "#2CA02C", alpha = 0.85) +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                  color = "#2CA02C", height = 0.2, linewidth = 1) +
    scale_size_continuous(range = c(2, 5), guide = "none") +
    labs(
      title = "Sensitivity Analyses: Robustness of E-DII & Cognition Association",
      subtitle = "Standardized β (95% CI) | Cross-validation by exclusion criteria",
      x = "Standardized β (per 1-SD increase in anti-inflammatory diet)",
      y = "",
      caption = "Dashed line = Primary analysis estimate. Point size proportional to N."
    ) +
    theme_nhanes(12) +
    theme(panel.grid.major.y = element_blank())

  ggsave(file.path(FIG_DIR, "fig5_sensitivity_tornado.pdf"),
         p_sens, width = 10, height = 6.5, device = "pdf")
  ggsave(file.path(FIG_DIR, "fig5_sensitivity_tornado.png"),
         p_sens, width = 10, height = 6.5, dpi = 300)

  cat("  Figure 5 saved.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 6: Mediation path diagram
# ═══════════════════════════════════════════════════════════════════════════════

cat("Generating Figure 6: Mediation path diagram...\n")

if (exists("sem_fit") && !is.null(sem_fit)) {
  sem_params <- parameterEstimates(sem_fit, standardized = TRUE)

  # Extract key paths
  path_a <- sem_params %>% filter(lhs == "phq9_total", rhs == "E_DII_raw",
                                  op == "~") %>% pull(std.all)
  path_b <- sem_params %>% filter(lhs == "cog_z_composite", rhs == "phq9_total",
                                  op == "~") %>% pull(std.all)
  path_c <- sem_params %>% filter(lhs == "cog_z_composite", rhs == "E_DII_raw",
                                  op == "~") %>% pull(std.all)
  indirect <- sem_params %>% filter(label == "ab") %>% pull(est)

  # Build a clean mediation diagram data frame
  med_data <- data.frame(
    x    = c(0, 1, 2),
    y    = c(1, 1, 1),
    label = c("Dietary\nInflammatory\nIndex (DII)",
              "Depressive\nSymptoms\n(PHQ-9)",
              "Cognitive\nFunction\n(Composite Z)")
  )

  p_med <- ggplot(med_data, aes(x = x, y = y)) +
    # Boxes
    geom_rect(aes(xmin = -0.15, xmax = 0.15, ymin = 0.8, ymax = 1.2),
              fill = "#E8F0FE", color = "#1F77B4", linewidth = 1.2) +
    geom_rect(aes(xmin = 0.85, xmax = 1.15, ymin = 0.8, ymax = 1.2),
              fill = "#FDE8E8", color = "#D62728", linewidth = 1.2) +
    geom_rect(aes(xmin = 1.85, xmax = 2.15, ymin = 0.8, ymax = 1.2),
              fill = "#E8F5E9", color = "#2CA02C", linewidth = 1.2) +
    # Labels
    geom_text(aes(label = label), size = 3.8, fontface = "bold") +
    # Path a: DII → Depression
    annotate("segment", x = 0.15, xend = 0.85, y = 1.15, yend = 1.15,
             arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
             color = "#D62728", linewidth = 1.5) +
    annotate("text", x = 0.5, y = 1.28,
             label = sprintf("Path a: %.3f", path_a),
             size = 3.5, color = "#D62728", fontface = "bold") +
    # Path b: Depression → Cognition
    annotate("segment", x = 1.15, xend = 1.85, y = 1.15, yend = 1.15,
             arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
             color = "#D62728", linewidth = 1.5) +
    annotate("text", x = 1.5, y = 1.28,
             label = sprintf("Path b: %.3f", path_b),
             size = 3.5, color = "#D62728", fontface = "bold") +
    # Path c': DII → Cognition (direct)
    annotate("segment", x = 0.15, xend = 1.85, y = 0.95, yend = 0.95,
             arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
             color = "#1F77B4", linewidth = 1.5, linetype = "longdash") +
    annotate("text", x = 1, y = 0.82,
             label = sprintf("Path c': %.3f (direct)", path_c),
             size = 3.5, color = "#1F77B4", fontface = "bold") +
    # Indirect effect annotation
    annotate("text", x = 1, y = 0.65,
             label = sprintf("Indirect (a×b): %.3f (%.1f%% mediated)",
                            indirect,
                            abs(indirect / (indirect + path_c)) * 100),
             size = 4, color = "grey20", fontface = "italic") +
    # Expand y-axis
    ylim(0.5, 1.5) +
    xlim(-0.5, 2.5) +
    theme_void()

  ggsave(file.path(FIG_DIR, "fig6_mediation_diagram.pdf"),
         p_med, width = 9, height = 5, device = "pdf")
  ggsave(file.path(FIG_DIR, "fig6_mediation_diagram.png"),
         p_med, width = 9, height = 5, dpi = 300)

  cat("  Figure 6 saved.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Supplementary Figures
# ═══════════════════════════════════════════════════════════════════════════════

# Figure S1: Raw scatter + loess of E-DII vs cognitive z-score
cat("Generating supplementary figures...\n")

if (exists("df_cc")) {
  # Sample for plotting (avoid overplotting with large N)
  set.seed(42)
  plot_sample <- df_cc %>%
    filter(!is.na(E_DII_raw), !is.na(cog_z_composite)) %>%
    sample_n(min(1500, n()))

  p_s1 <- ggplot(plot_sample, aes(x = E_DII_raw, y = cog_z_composite)) +
    geom_point(alpha = 0.15, size = 1.2, color = "grey40") +
    geom_smooth(method = "loess", span = 0.5, color = "#D62728",
                fill = "#D62728", alpha = 0.1, linewidth = 1.3) +
    geom_smooth(method = "lm", color = "#1F77B4", linetype = "dashed",
                linewidth = 1, se = FALSE) +
    labs(title = "Raw Association: E-DII & Cognitive Function",
         subtitle = "Red = LOESS (nonlinear fit) | Blue = Linear fit",
         x = "Energy-Adjusted Dietary Inflammatory Index",
         y = "Cognitive Composite Z-score") +
    theme_nhanes()

  ggsave(file.path(FIG_DIR, "figS1_scatter_loess.pdf"),
         p_s1, width = 9, height = 6.5, device = "pdf")
  ggsave(file.path(FIG_DIR, "figS1_scatter_loess.png"),
         p_s1, width = 9, height = 6.5, dpi = 300)

  cat("  Figure S1 saved.\n")
}

# Figure S2: Distribution of cognitive scores by domain
if (exists("df")) {
  cog_long <- df %>%
    dplyr::select(cerad_imm, cerad_del, animal_flu, dsst) %>%
    tidyr::pivot_longer(everything(), names_to = "domain", values_to = "score") %>%
    filter(!is.na(score)) %>%
    mutate(
      domain = case_when(
        domain == "cerad_imm"  ~ "CERAD Immediate",
        domain == "cerad_del"  ~ "CERAD Delayed",
        domain == "animal_flu" ~ "Animal Fluency",
        domain == "dsst"       ~ "DSST"
      )
    )

  p_s2 <- ggplot(cog_long, aes(x = score, fill = domain)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~domain, scales = "free") +
    scale_fill_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"),
                      guide = "none") +
    labs(title = "Distribution of Cognitive Domain Scores",
         subtitle = "NHANES 2011-2014, Age ≥ 60",
         x = "Raw Score", y = "Density") +
    theme_nhanes()

  ggsave(file.path(FIG_DIR, "figS2_cog_distributions.pdf"),
         p_s2, width = 10, height = 7, device = "pdf")
  ggsave(file.path(FIG_DIR, "figS2_cog_distributions.png"),
         p_s2, width = 10, height = 7, dpi = 300)

  cat("  Figure S2 saved.\n")
}

# Figure S3: PNI tertile distributions of cognitive scores
if (exists("df_cc") && "pni_tertile" %in% names(df_cc)) {

  pni_cog <- df_cc %>%
    dplyr::select(pni_tertile, cerad_imm, cerad_del, animal_flu, dsst) %>%
    tidyr::pivot_longer(-pni_tertile, names_to = "domain", values_to = "score") %>%
    filter(!is.na(score), !is.na(pni_tertile))

  p_s3 <- ggplot(pni_cog, aes(x = pni_tertile, y = score, fill = pni_tertile)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    facet_wrap(~domain, scales = "free_y") +
    scale_fill_manual(values = c("#E8F5E9", "#A5D6A7", "#4CAF50"),
                      guide = "none") +
    labs(title = "Cognitive Domain Scores by PNI Tertile",
         subtitle = "Higher PNI = better nutritional status",
         x = "Prognostic Nutritional Index", y = "Score") +
    theme_nhanes()

  ggsave(file.path(FIG_DIR, "figS3_pni_cog_boxplot.pdf"),
         p_s3, width = 10, height = 7, device = "pdf")
  ggsave(file.path(FIG_DIR, "figS3_pni_cog_boxplot.png"),
         p_s3, width = 10, height = 7, dpi = 300)

  cat("  Figure S3 saved.\n")
}

# Figure S4: Incremental R² bar chart
if (exists("r2_increment") && nrow(r2_increment) > 0) {

  r2_plot <- r2_increment %>%
    mutate(label = factor(label, levels = rev(label)))

  p_s4 <- ggplot(r2_plot, aes(x = label, y = delta_R2)) +
    geom_col(aes(fill = delta_R2), width = 0.6, alpha = 0.85) +
    geom_text(aes(label = sprintf("ΔR²=%.4f\np=%.4f", delta_R2, p_Ftest)),
              hjust = -0.1, size = 3.5) +
    scale_fill_gradient(low = "#A5D6A7", high = "#2E7D32", guide = "none") +
    coord_flip() +
    labs(title = "Incremental R²: Nutritional Indices Beyond Base Model",
         subtitle = "Base model: demographics + SES + lifestyle + comorbidity",
         x = "", y = "Δ R²") +
    expand_limits(y = max(r2_plot$delta_R2) * 1.5) +
    theme_nhanes()

  ggsave(file.path(FIG_DIR, "figS4_r2_increment.pdf"),
         p_s4, width = 9, height = 5, device = "pdf")
  ggsave(file.path(FIG_DIR, "figS4_r2_increment.png"),
         p_s4, width = 9, height = 5, dpi = 300)

  cat("  Figure S4 saved.\n")
}

# Figure S5: Sequential model forest (all 4 models, E-DII)
if (exists("main_results")) {
  edii_4models <- main_results %>%
    filter(exposure == "E_DII_scaled") %>%
    mutate(model = factor(model, levels = rev(unique(model))))

  p_s5 <- ggplot(edii_4models, aes(x = beta, y = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(color = "#1F77B4", size = 3) +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                  color = "#1F77B4", height = 0.2, linewidth = 1.3) +
    labs(title = "E-DII & Cognitive Function: Sequential Adjustment",
         subtitle = "Stability of association across increasing covariate adjustment",
         x = "Standardized β (95% CI)", y = "") +
    theme_nhanes() +
    theme(panel.grid.major.y = element_blank())

  ggsave(file.path(FIG_DIR, "figS5_sequential_models.pdf"),
         p_s5, width = 9, height = 5, device = "pdf")
  ggsave(file.path(FIG_DIR, "figS5_sequential_models.png"),
         p_s5, width = 9, height = 5, dpi = 300)

  cat("  Figure S5 saved.\n")
}

cat(sprintf("\n── 03_figures.R complete — %d figures saved to %s ──\n",
            6, FIG_DIR))
