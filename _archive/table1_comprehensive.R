# table1_comprehensive.R — Publication-quality Table 1 (all ages)
# Rscript table1_comprehensive.R

library(dplyr)

PROJ <- "D:/Researching/NHANES_aged_GI_tumor_nutrition"
RES_DIR <- file.path(PROJ, "results/gi_analysis")

df <- readRDS(file.path(RES_DIR, "nhanes_combined_gi_fixed.rds"))
df <- df %>% filter(!is.na(gi_tumor), surv_years > 0)

df$grp <- case_when(df$gi_tumor==1 ~ "GI Tumor",
                     df$any_cancer==0 ~ "Non-Cancer",
                     TRUE ~ "Other Cancer")

fmt_m <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
fmt_n <- function(x) sprintf("%d (%.1f%%)", sum(x, na.rm=TRUE), mean(x, na.rm=TRUE)*100)
fmt_m2 <- function(x) sprintf("%.2f (%.2f)", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))

rows <- list()
for (g in c("GI Tumor", "Non-Cancer", "Other Cancer")) {
  s <- df %>% filter(grp == g)
  n <- nrow(s)

  rows[[g]] <- data.frame(
    Variable = c("N", "Age, years", "Sex, Female", "Sex, Male",
                 "Race: Non-Hispanic White", "Race: Non-Hispanic Black",
                 "Race: Hispanic", "Race: Other",
                 "Education > High school", "Education <= High school",
                 "BMI, kg/m2", "Albumin, g/dL", "Lymphocyte, /uL",
                 "Total cholesterol, mg/dL",
                 "PNI", "CONUT", "GNRI",
                 "Deaths", "Survival, years"),
    Value = c(
      as.character(n),
      fmt_m(s$age),
      fmt_n(s$sex=="Female"),
      fmt_n(s$sex=="Male"),
      fmt_n(s$race_eth=="Non-Hispanic White"),
      fmt_n(s$race_eth=="Non-Hispanic Black"),
      fmt_n(s$race_eth=="Hispanic"),
      fmt_n(!s$race_eth %in% c("Non-Hispanic White","Non-Hispanic Black","Hispanic")),
      fmt_n(s$edu_binary==1),
      fmt_n(s$edu_binary==0),
      fmt_m(s$bmi),
      fmt_m2(s$albumin_gdl),
      sprintf("%.0f (%.0f)", mean(s$lymph_abs, na.rm=TRUE), sd(s$lymph_abs, na.rm=TRUE)),
      fmt_m(s$tchol_mgdl),
      fmt_m(s$PNI),
      fmt_m(s$CONUT),
      fmt_m(s$GNRI),
      fmt_n(s$death==1),
      fmt_m(s$surv_years)
    ),
    stringsAsFactors=FALSE
  )
}

t1 <- bind_cols(Variable=rows[["GI Tumor"]]$Variable,
                GI_Tumor=rows[["GI Tumor"]]$Value,
                Non_Cancer=rows[["Non-Cancer"]]$Value,
                Other_Cancer=rows[["Other Cancer"]]$Value)
names(t1) <- c("Variable", "GI Tumor (n=313)", "Non-Cancer (n=40,480)", "Other Cancer (n=3,273)")

write.csv(t1, file.path(RES_DIR, "table1_comprehensive.csv"), row.names=FALSE)
cat("=== Table 1: Comprehensive ===\n")
print(as.data.frame(t1), row.names=FALSE)
cat(sprintf("\nSaved to: %s\n", file.path(RES_DIR, "table1_comprehensive.csv")))
