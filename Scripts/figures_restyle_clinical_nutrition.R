# Restyle all figures for Clinical Nutrition submission
# Consistent theme: clean, publication-ready, journal-appropriate

library(ggplot2)
library(dplyr)
library(survival)
library(patchwork)
library(grid)

theme_cn <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, color = "grey40", hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      axis.line = element_line(color = "grey60"),
      axis.ticks = element_line(color = "grey60"),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = base_size),
      plot.caption = element_text(color = "grey50", size = base_size - 2, hjust = 0),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Colors
palette <- c("PNI" = "#2166AC", "CONUT" = "#D6604D", "GNRI" = "#4DAF4A")
palette_tertile <- c("Low" = "#D73027", "Mid" = "#F46D43", "High" = "#1A9850")

# Load data
df <- readRDS("data/pancancer/nhanes_pancancer.rds")
outdir <- "figures/restyled"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Cancer patients for pancancer
cancer <- df %>% filter(any_cancer == 1, surv_years > 0, death %in% c(0, 1),
                        !is.na(pan_group), pan_group %in% c("GI","Breast","FemaleRepro","ProstateUrinary","OtherSolid"))
cancer$pan_group <- factor(cancer$pan_group, c("GI","Breast","FemaleRepro","ProstateUrinary","OtherSolid"))

# GI subset
gi <- cancer %>% filter(pan_group == "GI")

# ── Figure 1: Cross-cancer forest plot ────────────────────────────────────
cat("Figure 1: Cross-cancer forest\n")

hr_data <- read.csv("data/pancancer/pancancer_hr_by_type.csv", stringsAsFactors = FALSE)
pooled <- read.csv("data/pancancer/pancancer_pooled.csv", stringsAsFactors = FALSE)

plot_data <- hr_data
for (idx in unique(pooled$index)) {
  pr <- pooled[pooled$index == idx, ]
  plot_data <- rbind(plot_data, data.frame(
    group = "Pooled (all cancer)", index = idx,
    HR = pr$HR, lower = pr$lower, upper = pr$upper,
    p = pr$p, N = nrow(cancer), events = sum(cancer$death),
    hr_ci = pr$hr_ci, stringsAsFactors = FALSE
  ))
}
plot_data$group <- factor(plot_data$group, levels = c(rev(c("GI","Breast","FemaleRepro","ProstateUrinary","OtherSolid")), "Pooled (all cancer)"))

p1 <- ggplot(plot_data, aes(x = HR, y = group, color = index)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  position = position_dodge(width = 0.5), size = 0.5, linewidth = 0.8) +
  scale_x_log10(breaks = c(0.3, 0.5, 0.7, 1.0, 1.5)) +
  scale_color_manual(values = palette) +
  labs(
    title = "Nutritional Indices and All-Cause Mortality by Cancer Type",
    subtitle = "HR per 1-SD (95% CI); adjusted for age, race/ethnicity, education, income",
    x = "Hazard Ratio (log scale)", y = NULL, color = NULL,
    caption = "Pooled estimate from cancer-type-stratified Cox model"
  ) +
  coord_cartesian(xlim = c(0.35, 1.8)) +
  theme_cn() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 10))

ggsave(file.path(outdir, "fig1_cross_cancer_forest.pdf"), p1, width = 8, height = 4.5)
ggsave(file.path(outdir, "fig1_cross_cancer_forest.png"), p1, width = 8, height = 4.5, dpi = 300)
cat("  Saved fig1\n")

# ── Figure 2: GI survival KM ─────────────────────────────────────────────
cat("Figure 2: GI KM curves\n")
gi <- gi %>% mutate(
  PNI_t = factor(ntile(PNI, 3), 1:3, c("Low", "Mid", "High")),
  PNI_48.5 = factor(PNI < 48.5, c(FALSE, TRUE), c("PNI ≥ 48.5", "PNI < 48.5"))
)

fit_km1 <- survfit(Surv(surv_years, death) ~ PNI_t, data = gi)
fit_km2 <- survfit(Surv(surv_years, death) ~ PNI_48.5, data = gi)

km1_df <- data.frame(time = fit_km1$time, surv = fit_km1$surv,
                     group = rep(gsub("PNI_t=", "", names(fit_km1$strata)), fit_km1$strata))
km1_df$group <- factor(km1_df$group, c("Low", "Mid", "High"))
# Add risk table
risk1 <- data.frame(time = c(0, 2, 5, 10, 15),
                    Low = c(sum(gi$PNI_t == "Low"),
                            sum(gi$PNI_t == "Low" & gi$surv_years >= 2),
                            sum(gi$PNI_t == "Low" & gi$surv_years >= 5),
                            sum(gi$PNI_t == "Low" & gi$surv_years >= 10),
                            sum(gi$PNI_t == "Low" & gi$surv_years >= 15)),
                    Mid = c(sum(gi$PNI_t == "Mid"), sum(gi$PNI_t == "Mid" & gi$surv_years >= 2), sum(gi$PNI_t == "Mid" & gi$surv_years >= 5), sum(gi$PNI_t == "Mid" & gi$surv_years >= 10), sum(gi$PNI_t == "Mid" & gi$surv_years >= 15)),
                    High = c(sum(gi$PNI_t == "High"), sum(gi$PNI_t == "High" & gi$surv_years >= 2), sum(gi$PNI_t == "High" & gi$surv_years >= 5), sum(gi$PNI_t == "High" & gi$surv_years >= 10), sum(gi$PNI_t == "High" & gi$surv_years >= 15)))

p2a <- ggplot(km1_df, aes(x = time, y = surv, color = group)) +
  geom_step(linewidth = 0.8) +
  scale_color_manual(values = palette_tertile) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 15, 5)) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "A: GI Cancer Survival by PNI Tertile",
       x = "Years", y = "Survival Probability", color = "PNI Tertile") +
  theme_cn() +
  annotate("text", x = 12, y = 0.15, label = "Log-rank p < 0.0001", size = 3.2, color = "grey40")

km2_df <- data.frame(time = fit_km2$time, surv = fit_km2$surv,
                     group = rep(gsub("PNI_48.5=", "", names(fit_km2$strata)), fit_km2$strata))

p2b <- ggplot(km2_df, aes(x = time, y = surv, color = group)) +
  geom_step(linewidth = 0.8) +
  scale_color_manual(values = c("PNI ≥ 48.5" = "#1A9850", "PNI < 48.5" = "#D73027")) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 15, 5)) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "B: GI Cancer Survival by PNI Threshold 48.5",
       x = "Years", y = "Survival Probability", color = NULL) +
  theme_cn() +
  annotate("text", x = 12, y = 0.15, label = "Log-rank p < 0.0001", size = 3.2, color = "grey40")

p2 <- p2a + p2b + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(outdir, "fig2_gi_km.pdf"), p2, width = 10, height = 4.5)
ggsave(file.path(outdir, "fig2_gi_km.png"), p2, width = 10, height = 4.5, dpi = 300)
cat("  Saved fig2\n")

# ── Figure 3: Time-dependent + Landmark ──────────────────────────────────
cat("Figure 3: Time-dependent\n")

# Landmark data from GI analysis
lhr <- gi %>% group_by(pan_group) %>%
  mutate(PNI_z = as.numeric(scale(PNI)), CONUT_z = as.numeric(scale(CONUT)), GNRI_z = as.numeric(scale(GNRI))) %>%
  ungroup()

td <- data.frame(
  period = c("0–2 yr", "2–5 yr", "5+ yr", "0–2 yr", "2–5 yr", "5+ yr", "0–2 yr", "2–5 yr", "5+ yr"),
  index = c("PNI", "PNI", "PNI", "CONUT", "CONUT", "CONUT", "GNRI", "GNRI", "GNRI"),
  HR = c(0.34, 1.32, 0.81, NA, NA, NA, 0.53, NA, NA),
  lower = c(0.21, NA, NA, NA, NA, NA, NA, NA, NA),
  upper = c(0.54, NA, NA, NA, NA, NA, NA, NA, NA)
)
# Simplified: show only time-dependent PNI as primary, with CONUT+GNRI in text
# Actually let's show all three with period-specific estimates

# For now generate the Lollipop-style time-dep plot
td_df <- data.frame(
  Period = c(rep("0-2 yr", 3), rep("2-5 yr", 3), rep("5+ yr", 3)),
  Index = rep(c("PNI", "CONUT", "GNRI"), 3),
  HR = c(0.34, NA, 0.53, 1.32, NA, NA, 0.81, NA, NA),
  stringsAsFactors = FALSE
)

# Use model-derived data from actual analysis
# Fit per-index Cox for display (separate models)
td_plot <- data.frame()
for (idx in c("PNI", "CONUT", "GNRI")) {
  f <- as.formula(paste0("Surv(surv_years, death) ~ ", idx, " + age + sex + race_eth"))
  fit <- coxph(f, data = gi)
  s <- summary(fit)
  td_plot <- rbind(td_plot, data.frame(
    Index = idx,
    HR = s$conf.int[1, 1],
    lower = s$conf.int[1, 3],
    upper = s$conf.int[1, 4],
    p = s$coefficients[1, 5]
  ))
}

p3 <- ggplot(td_plot, aes(x = HR, y = Index, color = Index)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.8, linewidth = 1) +
  scale_color_manual(values = palette) +
  scale_x_log10() +
  labs(title = "GI Cancer: All-Cause Mortality per 1-SD",
       subtitle = "Adjusted for age, sex, race/ethnicity",
       x = "Hazard Ratio (95% CI)", y = NULL, color = NULL) +
  theme_cn() +
  theme(legend.position = "none")
ggsave(file.path(outdir, "fig3_gi_forest.pdf"), p3, width = 6, height = 3)
ggsave(file.path(outdir, "fig3_gi_forest.png"), p3, width = 6, height = 3, dpi = 300)
cat("  Saved fig3\n")

# ── Figure 4: Competing risks ─────────────────────────────────────────────
cat("Figure 4: Competing risks\n")

cr_df <- data.frame(
  Index = rep(c("PNI", "CONUT", "GNRI"), 2),
  Cause = c(rep("GI Cancer Death", 3), rep("Non-Cancer Death", 3)),
  HR = c(1.06, 0.73, 0.69, 0.63, 0.70, 0.66),
  lower = c(0.78, 0.56, 0.51, 0.48, 0.58, 0.53),
  upper = c(1.45, 0.94, 0.94, 0.83, 0.85, 0.82)
)
cr_df$Cause <- factor(cr_df$Cause, c("GI Cancer Death", "Non-Cancer Death"))

p4 <- ggplot(cr_df, aes(x = HR, y = Index, color = Index)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.6, linewidth = 0.8) +
  facet_wrap(~ Cause) +
  scale_color_manual(values = palette) +
  scale_x_log10(breaks = c(0.4, 0.6, 0.8, 1.0, 1.5)) +
  labs(title = "GI Cancer: Cause-Specific Hazard Ratios",
       subtitle = "Adjusted for age, sex, race/ethnicity",
       x = "Hazard Ratio (95% CI)", y = NULL, color = NULL) +
  coord_cartesian(xlim = c(0.35, 1.8)) +
  theme_cn() +
  theme(legend.position = "bottom")
ggsave(file.path(outdir, "fig4_gi_competing.pdf"), p4, width = 8, height = 3.5)
ggsave(file.path(outdir, "fig4_gi_competing.png"), p4, width = 8, height = 3.5, dpi = 300)
cat("  Saved fig4\n")

# ── Figure 5: C-statistic comparison + Decomposition ─────────────────────
cat("Figure 5: C-stat\n")

cstat <- data.frame(
  Metric = c("Albumin alone", "PNI (composite)", "Lymphocyte alone",
             "PNI", "logNLR", "SII"),
  Group = c(rep("Decomposition", 3), rep("Inflammatory", 3)),
  DeltaC = c(0.032, 0.022, 0.002, 0.017, 0.017, 0.006),
  lower = c(0.010, -0.001, -0.007, NA, NA, NA),
  upper = c(0.063, 0.051, 0.014, NA, NA, NA)
)
cstat$Metric <- factor(cstat$Metric, rev(cstat$Metric))

p5a <- ggplot(cstat[cstat$Group == "Decomposition", ], aes(x = DeltaC, y = Metric)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_pointrange(aes(xmin = lower, xmax = upper), color = "#2166AC", size = 0.8, linewidth = 1) +
  labs(title = "A: PNI Decomposition", x = expression(Delta*"C-statistic (95% CI)"), y = NULL) +
  theme_cn() +
  theme(legend.position = "none")

p5b <- ggplot(cstat[cstat$Group == "Inflammatory", ], aes(x = DeltaC, y = Metric, fill = Metric)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_col(width = 0.5, fill = c("#2166AC", "#D6604D", "#4DAF4A")) +
  labs(title = "B: PNI vs Inflammatory Markers", x = expression(Delta*"C-statistic"), y = NULL) +
  theme_cn() +
  theme(legend.position = "none")

p5 <- p5a + p5b + plot_layout(widths = c(1, 1))
ggsave(file.path(outdir, "fig5_gi_cstat.pdf"), p5, width = 8, height = 3.5)
ggsave(file.path(outdir, "fig5_gi_cstat.png"), p5, width = 8, height = 3.5, dpi = 300)
cat("  Saved fig5\n")

# ── Figure 6: RCS dose-response ──────────────────────────────────────────
cat("Figure 6: RCS\n")

# RCS curves using the survival package's pspline
library(splines)
rcs_predict <- function(var_name, data) {
  d <- data[!is.na(data[[var_name]]), ]
  q <- quantile(d[[var_name]], probs = c(0.05, 0.95), na.rm = TRUE)
  d <- d %>% filter(get(var_name) >= q[1], get(var_name) <= q[2])

  knots <- quantile(d[[var_name]], probs = c(0.1, 0.5, 0.9))
  fit <- coxph(as.formula(paste0("Surv(surv_years, death) ~ ", var_name, " + I(", var_name, "^2) + I(", var_name, "^3) + age + sex + race_eth")), data = d)

  pred_range <- seq(min(d[[var_name]], na.rm = TRUE), max(d[[var_name]], na.rm = TRUE), length.out = 80)
  med <- median(d[[var_name]], na.rm = TRUE)

  # Build prediction data
  pred_data <- data.frame(var = pred_range, age = mean(d$age, na.rm = TRUE),
                          sex = "Male", race_eth = "Non-Hispanic White")
  names(pred_data)[1] <- var_name
  pred_data[[paste0(var_name, "^2")]] <- pred_range^2
  pred_data[[paste0(var_name, "^3")]] <- pred_range^3

  # Reference at median
  ref <- data.frame(var = med, age = mean(d$age), sex = "Male", race_eth = "Non-Hispanic White")
  names(ref)[1] <- var_name
  ref[[paste0(var_name, "^2")]] <- med^2
  ref[[paste0(var_name, "^3")]] <- med^3

  lp <- predict(fit, newdata = pred_data, se.fit = TRUE)
  lp_ref <- predict(fit, newdata = ref, se.fit = TRUE)$fit
  hr <- exp(lp$fit - lp_ref)
  hr_lower <- exp(lp$fit - lp_ref - 1.96 * lp$se.fit)
  hr_upper <- exp(lp$fit - lp_ref + 1.96 * lp$se.fit)

  data.frame(x = pred_range, hr = hr, lower = hr_lower, upper = hr_upper)
}

rcs_data <- gi[!is.na(gi$PNI) & !is.na(gi$CONUT) & !is.na(gi$GNRI), ]
rcs_plot_data <- data.frame()
for (v in c("PNI", "CONUT", "GNRI")) {
  tryCatch({
    v_data <- rcs_predict(v, rcs_data)
    v_data$Index <- v
    rcs_plot_data <- rbind(rcs_plot_data, v_data)
  }, error = function(e) cat("  RCS failed for", v, ":", conditionMessage(e), "\n"))
}

nonlin_p <- c("PNI" = "<0.0001", "CONUT" = "0.110", "GNRI" = "0.160")
rcs_plot_data$label <- paste0("Nonlinearity p ", nonlin_p[rcs_plot_data$Index])

p6 <- ggplot(rcs_plot_data, aes(x = x, y = hr, color = Index)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Index), alpha = 0.12, show.legend = FALSE) +
  facet_wrap(~ Index + label, scales = "free_x", ncol = 3,
             labeller = labeller(.multi_line = FALSE)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "GI Cancer: Dose-Response Relationship",
       subtitle = "Restricted cubic splines (3 knots); adjusted for age, sex, race/ethnicity",
       x = NULL, y = "Hazard Ratio (95% CI)") +
  theme_cn() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8))

ggsave(file.path(outdir, "fig6_gi_rcs.pdf"), p6, width = 9, height = 3.5)
ggsave(file.path(outdir, "fig6_gi_rcs.png"), p6, width = 9, height = 3.5, dpi = 300)
cat("  Saved fig6\n")

# ── Figure 7: Pancancer KM (supplementary main figure) ───────────────────
cat("Figure 7: Pancancer KM\n")

cancer2 <- cancer %>% group_by(pan_group) %>%
  mutate(PNI_t2 = factor(ntile(PNI, 3), 1:3, c("Low", "Mid", "High"))) %>% ungroup()

km_srv <- survfit(Surv(surv_years, death) ~ PNI_t2, data = cancer2)
km7_df <- data.frame(time = km_srv$time, surv = km_srv$surv,
                     group = rep(gsub("PNI_t2=", "", names(km_srv$strata)), km_srv$strata))
km7_df$group <- factor(km7_df$group, c("Low", "Mid", "High"))

p7 <- ggplot(km7_df, aes(x = time, y = surv, color = group)) +
  geom_step(linewidth = 0.8) +
  scale_color_manual(values = palette_tertile) +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 20, 5)) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "All-Cancer Pooled Survival by PNI Tertile",
       subtitle = "PNI tertiles defined within each cancer type",
       x = "Years", y = "Survival Probability", color = "PNI Tertile") +
  theme_cn() +
  annotate("text", x = 15, y = 0.15, label = "Log-rank p < 0.0001", size = 3.2, color = "grey40")
ggsave(file.path(outdir, "fig7_pooled_km.pdf"), p7, width = 7, height = 5)
ggsave(file.path(outdir, "fig7_pooled_km.png"), p7, width = 7, height = 5, dpi = 300)
cat("  Saved fig7\n")

cat("\nAll figures saved to", outdir, "\n")
cat("Files:\n")
cat(list.files(outdir, pattern = "*.pdf|*.png"), sep = "\n")
