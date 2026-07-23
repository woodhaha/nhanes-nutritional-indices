# fix_race_eth.R — Fix race/ethnicity for continuous NHANES by re-downloading DEMO
library(dplyr); library(nhanesA)

PROJ <- "D:/Researching/NHANES_aged_GI_tumor_nutrition"
DATA_DIR <- file.path(PROJ, "data/gi_analysis")
RES_DIR <- file.path(PROJ, "results/gi_analysis")

df <- readRDS(file.path(DATA_DIR, "nhanes_gi_nutrition_raw.rds"))
cat(sprintf("Loaded: %d rows\n", nrow(df)))

# Fix: re-download DEMO tables and extract correct RIDRETH1 -> race_eth mapping
CYCLES <- c("2005-2006"="D","2007-2008"="E","2009-2010"="F",
            "2011-2012"="G","2013-2014"="H","2015-2016"="I")

reth1_map <- data.frame()
for (cyc in names(CYCLES)) {
  sfx <- CYCLES[cyc]
  cat(sprintf("  DEMO_%s...", sfx))
  demo <- tryCatch(nhanesA::nhanes(paste0("DEMO_", sfx)), error=function(e) NULL)
  if (is.null(demo)) { cat("FAILED\n"); next }

  # Extract SEQN + RIDRETH1
  reth1 <- if (inherits(demo$RIDRETH1, "haven_labelled")) as.numeric(demo$RIDRETH1) else as.numeric(demo$RIDRETH1)
  map <- data.frame(SEQN = as.numeric(demo$SEQN), RIDRETH1 = reth1, stringsAsFactors=FALSE)
  reth1_map <- bind_rows(reth1_map, map)
  cat(sprintf(" %d\n", nrow(map)))
}

cat(sprintf("Total DEMO records: %d\n", nrow(reth1_map)))

# Merge race_eth back to data
reth1_map$race_eth_fixed <- case_when(
  reth1_map$RIDRETH1 %in% c(1,2) ~ "Hispanic",
  reth1_map$RIDRETH1 == 3 ~ "Non-Hispanic White",
  reth1_map$RIDRETH1 == 4 ~ "Non-Hispanic Black",
  TRUE ~ "Other")

df <- merge(df, reth1_map[, c("SEQN","race_eth_fixed")], by="SEQN", all.x=TRUE)
cat(sprintf("Merged: %d matched, %d NA\n", sum(!is.na(df$race_eth_fixed)), sum(is.na(df$race_eth_fixed))))

# Replace race_eth
df$race_eth <- df$race_eth_fixed
df$race_eth_fixed <- NULL

cat("Corrected race_eth distribution:\n")
print(table(df$race_eth, useNA="ifany"))

saveRDS(df, file.path(DATA_DIR, "nhanes_gi_nutrition_raw.rds"))
cat("Saved.\n")
