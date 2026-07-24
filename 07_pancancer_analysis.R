# 07_pancancer_analysis.R — Pancancer and across-cancer nutritional index analysis
# 1. Cancer-type-specific HRs (forest plot)
# 2. Across-cancer pooled analysis (cancer type as strata)
# 3. Cancer type × PNI interaction
# ==============================================================================

library(dplyr)
library(survival)
library(ggplot2)
library(patchwork)

OUTPUT <- normalizePath("data/pancancer")
FIG_DIR <- file.path(getwd(), "figures/pancancer")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load data ────────────────────────────────────────────────────────────────
df <- readRDS(file.path(OUTPUT, "nhanes_pancancer.rds"))
cat(sprintf("Loaded: N=%d, cancer patients=%d\n",
            nrow(df), sum(df$any_cancer == 1, na.rm=TRUE)))

# Filter to cancer patients with complete survival data
cancer <- df %>% filter(
  any_cancer == 1,
  surv_years > 0,
  death %in% c(0, 1),
  !is.na(pan_group),
  pan_group != "NonCancer"
)
cat(sprintf("Analytic cancer sample: N=%d, events=%d\n", nrow(cancer), sum(cancer$death)))

# Focus on groups with ≥50 events for meaningful Cox
MIN_EVENTS <- 50
groups_keep <- c()
for (g in unique(cancer$pan_group)) {
  e <- sum(cancer$death[cancer$pan_group == g], na.rm=TRUE)
  if (e >= MIN_EVENTS) {
    groups_keep <- c(groups_keep, g)
    cat(sprintf("  %-20s: N=%-5d events=%-4d\n", g, sum(cancer$pan_group == g), e))
  } else {
    cat(sprintf("  %-20s: N=%-5d events=%-4d SKIP (<%d events)\n", g, sum(cancer$pan_group == g), e, MIN_EVENTS))
  }
}
cancer <- cancer %>% filter(pan_group %in% groups_keep)
cancer$pan_group <- droplevels(factor(cancer$pan_group))

# Standardize nutritional indices per group for comparability
cancer <- cancer %>%
  group_by(pan_group) %>%
  mutate(
    PNI_z = as.numeric(scale(PNI)),
    CONUT_z = as.numeric(scale(-CONUT)),     # reverse: higher = better nutrition
    GNRI_z = as.numeric(scale(GNRI))
  ) %>%
  ungroup()

# ── 1. Cancer-type-specific HRs ──────────────────────────────────────────────
cat("\n══ Cancer-type-specific HRs ══\n")

covars <- c("age", "sex", "race_eth", "edu_binary", "income_pir")
indexes <- c("PNI_z", "CONUT_z", "GNRI_z")
index_labels <- c("PNI_z"="PNI", "CONUT_z"="CONUT", "GNRI_z"="GNRI")

type_results <- data.frame()
for (g in groups_keep) {
  sub <- cancer %>% filter(pan_group == g)
  for (idx in indexes) {
    f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+", paste(covars, collapse=" + ")))
    fit <- tryCatch(coxph(f, data=sub), error=function(e) NULL)
    if (is.null(fit)) next

    s <- summary(fit)
    ci <- s$conf.int
    cc <- s$coefficients
    hr_row <- ci[rownames(ci) == idx, , drop=FALSE]
    p_row <- cc[rownames(cc) == idx, , drop=FALSE]
    if (nrow(hr_row) == 0) { cat(sprintf("  %s %s: failed\n", g, idx)); next }

    type_results <- rbind(type_results, data.frame(
      group = g,
      index = index_labels[idx],
      HR = hr_row[1, 1],
      lower = hr_row[1, 3],
      upper = hr_row[1, 4],
      p = p_row[1, 5],
      N = nrow(sub),
      events = sum(sub$death)
    ))
  }
}

type_results$hr_ci <- with(type_results, sprintf("%.2f (%.2f-%.2f)", HR, lower, upper))
cat("\nPancancer-specific Cox results (adjusted HR per 1-SD):\n")
print(type_results[, c("group", "index", "hr_ci", "p", "N", "events")])

write.csv(type_results, file.path(OUTPUT, "pancancer_hr_by_type.csv"), row.names=FALSE)

# ── 2. Across-cancer pooled analysis ─────────────────────────────────────────
cat("\n══ Across-cancer pooled analysis (stratified by cancer type) ══\n")

pooled_results <- data.frame()
for (idx in indexes) {
  # Stratified Cox: different baseline hazard per cancer type
  f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+",
                         paste(covars, collapse=" + "), "+ strata(pan_group)"))
  fit <- coxph(f, data=cancer)
  s <- summary(fit)
  ci <- s$conf.int[rownames(s$conf.int) == idx, ]
  cc <- s$coefficients[rownames(s$coefficients) == idx, ]

  pooled_results <- rbind(pooled_results, data.frame(
    index = index_labels[idx],
    HR = ci[[1]], lower = ci[[3]],
    upper = ci[[4]], p = cc[[5]]
  ))
}
pooled_results$hr_ci <- with(pooled_results, sprintf("%.2f (%.2f-%.2f)", HR, lower, upper))
cat("Pooled (cancer-type stratified):\n")
print(pooled_results[, c("index", "hr_ci", "p")])
write.csv(pooled_results, file.path(OUTPUT, "pancancer_pooled.csv"), row.names=FALSE)

# Also try with cancer type as covariate (not strata) — fully adjusted
cat("\n-- Pooled with cancer type as covariate --\n")
pooled_results2 <- data.frame()
for (idx in indexes) {
  f <- as.formula(paste("Surv(surv_years, death) ~", idx, "+",
                         paste(covars, collapse=" + "), "+ pan_group"))
  fit <- coxph(f, data=cancer)
  s <- summary(fit)
  ci <- s$conf.int[rownames(s$conf.int) == idx, ]
  cc <- s$coefficients[rownames(s$coefficients) == idx, ]
  pooled_results2 <- rbind(pooled_results2, data.frame(
    index = index_labels[idx],
    HR = ci[[1]], lower = ci[[3]],
    upper = ci[[4]], p = cc[[5]]
  ))
}
pooled_results2$hr_ci <- with(pooled_results2, sprintf("%.2f (%.2f-%.2f)", HR, lower, upper))
print(pooled_results2[, c("index", "hr_ci", "p")])

# ── 3. Cancer type × PNI interaction ─────────────────────────────────────────
cat("\n══ Cancer type × Nutrition index interaction ══\n")

interaction_results <- data.frame()
for (idx in indexes) {
  f <- as.formula(paste("Surv(surv_years, death) ~", idx, "* pan_group +",
                         paste(covars, collapse=" + ")))
  fit <- coxph(f, data=cancer)
  # Extract interaction terms
  cc <- summary(fit)$coefficients
  interxs <- cc[grepl(":", rownames(cc)), , drop=FALSE]
  # LRT for interaction
  f_no_interx <- as.formula(paste("Surv(surv_years, death) ~", idx, "+ pan_group +",
                                   paste(covars, collapse=" + ")))
  fit_null <- coxph(f_no_interx, data=cancer)
  lrt <- anova(fit_null, fit)

  cat(sprintf("\n%s × pan_group:\n", index_labels[idx]))
  if (nrow(interxs) > 0) print(interxs[, c(1, 2, 5), drop=FALSE])
  cat(sprintf("  LRT p = %.4f\n", lrt$`Pr(>|Chi|)`[2]))

  interaction_results <- rbind(interaction_results, data.frame(
    index = index_labels[idx],
    lrt_p = lrt$`Pr(>|Chi|)`[2]
  ))
}

write.csv(interaction_results, file.path(OUTPUT, "pancancer_interaction.csv"), row.names=FALSE)

# ── Figure: Forest plot ──────────────────────────────────────────────────────
cat("\n══ Generating forest plot ══\n")

plot_data <- type_results
# Add pooled results
for (idx in unique(plot_data$index)) {
  pr <- pooled_results[pooled_results$index == idx, ]
  plot_data <- rbind(plot_data, data.frame(
    group = "Pooled (all cancer)", index = idx,
    HR = pr$HR, lower = pr$lower, upper = pr$upper,
    p = pr$p, N = sum(cancer$death), events = nrow(cancer),
    hr_ci = pr$hr_ci
  ))
}
plot_data$group <- factor(plot_data$group, levels = c(rev(groups_keep), "Pooled (all cancer)"))

p <- ggplot(plot_data, aes(x = HR, y = group, color = index)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  position = position_dodge(width = 0.5), size = 0.6) +
  scale_x_log10(breaks = c(0.3, 0.5, 0.7, 1.0, 1.5)) +
  scale_color_manual(values = c("PNI" = "#E74C3C", "CONUT" = "#3498DB", "GNRI" = "#2ECC71")) +
  labs(
    title = "Nutritional Indices and All-Cause Mortality by Cancer Type",
    subtitle = "NHANES III + Continuous NHANES (1988-2019)",
    x = "HR per 1-SD (95% CI)", y = NULL, color = NULL,
    caption = "Adjusted for age, sex, race/ethnicity, education, income"
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  ) +
  coord_cartesian(xlim = c(0.3, 1.8))

ggsave(file.path(FIG_DIR, "fig_pancancer_forest.pdf"), p, width = 9, height = 5)
ggsave(file.path(FIG_DIR, "fig_pancancer_forest.png"), p, width = 9, height = 5, dpi = 300)
cat("Saved: fig_pancancer_forest.{pdf,png}\n")

# ── Figure: Across-cancer KM by PNI tertile ─────────────────────────────────
cat("\n-- KM by PNI tertile (pooled cancer) --\n")
cancer <- cancer %>%
  group_by(pan_group) %>%
  mutate(PNI_t = factor(ntile(PNI, 3), 1:3, c("Low", "Mid", "High"))) %>%
  ungroup()

# Weighted KM fit
fit_km <- survfit(Surv(surv_years, death) ~ PNI_t, data = cancer)
km_df <- data.frame(
  time = fit_km$time,
  surv = fit_km$surv,
  group = rep(gsub("PNI_t=", "", names(fit_km$strata)), fit_km$strata)
)

p_km <- ggplot(km_df, aes(x = time, y = surv, color = group)) +
  geom_step(linewidth = 1) +
  scale_color_manual(values = c("Low" = "#E74C3C", "Mid" = "#F39C12", "High" = "#27AE60")) +
  labs(
    title = "All-Cancer Pooled: PNI Tertile and Survival",
    x = "Years", y = "Survival Probability", color = "PNI Tertile"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "fig_pancancer_km.pdf"), p_km, width = 7, height = 5)
cat("Saved: fig_pancancer_km.{pdf,png}\n")

cat(sprintf("\n── All results in %s ──\n", OUTPUT))
cat(sprintf("── Figures in %s ──\n", FIG_DIR))
