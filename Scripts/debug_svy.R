library(survey)
library(dplyr)
options(survey.lonely.psu = "adjust")

d <- readRDS("results/gi_analysis/nhanes_clean.rds")
gi <- filter(d, gi_tumor == 1, surv_years > 0)
gi <- mutate(gi, PNI_s = as.numeric(scale(PNI)),
             wt_pooled = WTMEC2YR / 2,
             age_s = as.numeric(scale(age)),
             sex_b = ifelse(sex == "Female", 1L, 0L),
             race_eth_b = ifelse(race_eth == "Non-Hispanic White", 1L, 0L))

w <- filter(gi, !is.na(wt_pooled), !is.na(SDMVPSU), !is.na(SDMVSTRA))
des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled, data = w, nest = TRUE)
fit <- svycoxph(Surv(surv_years, death) ~ PNI_s + age_s + sex_b + race_eth_b, design = des)
s <- summary(fit)

cat("=== coefficients matrix ===\n")
cat("cols:", paste(colnames(s$coefficients), collapse = ", "), "\n")
print(s$coefficients)
cat("\n=== conf.int matrix ===\n")
cat("cols:", paste(colnames(s$conf.int), collapse = ", "), "\n")
print(s$conf.int)

# Extract PNI_s row
cat("\n=== PNI_s row ===\n")
cat("coef row:", paste(s$coefficients["PNI_s", ], collapse = " | "), "\n")
cat("conf.int row:", paste(s$conf.int["PNI_s", ], collapse = " | "), "\n")

# Manual p-value from z
z_val <- s$coefficients["PNI_s", 3]
p_manual <- 2 * pnorm(abs(z_val), lower.tail = FALSE)
cat(sprintf("z = %.4f, p (from z) = %.4f\n", z_val, p_manual))
