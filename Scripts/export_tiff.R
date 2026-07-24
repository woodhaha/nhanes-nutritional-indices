# Export all figures as TIFF (LZW compressed) for Clinical Nutrition
library(ggplot2)
library(png)

# Read existing R objects
source("Scripts/figures_restyle_clinical_nutrition.R", local=TRUE)

outdir <- "figures/restyled"

tiff_save <- function(plot, filename, w, h, dpi=300) {
  fp <- file.path(outdir, filename)
  tiff(fp, width=w, height=h, units="in", res=dpi, compression="lzw")
  print(plot)
  dev.off()
  cat(paste0("  ", filename, ": ", w, "x", h, " in, ", dpi, " DPI, LZW\n"))
}

# Regenerate all plots
# Fig 1: Cross-cancer forest (p1)
if (exists("p1")) tiff_save(p1, "fig1_cross_cancer_forest.tiff", 8, 4.5)

# Fig 2: GI KM (p2)
if (exists("p2")) tiff_save(p2, "fig2_gi_km.tiff", 10, 4.5)

# Fig 3: Time-dependent forest (p3)
if (exists("p3")) tiff_save(p3, "fig3_gi_forest.tiff", 6, 3)

# Fig 4: Competing risks (p4)
if (exists("p4")) tiff_save(p4, "fig4_gi_competing.tiff", 8, 3.5)

# Fig 5: C-stat (p5)
if (exists("p5")) tiff_save(p5, "fig5_gi_cstat.tiff", 8, 3.5)

# Fig 6: RCS (p6)
if (exists("p6")) tiff_save(p6, "fig6_gi_rcs.tiff", 9, 3.5)

# Fig 7: Pooled KM (p7)
if (exists("p7")) tiff_save(p7, "fig7_pooled_km.tiff", 7, 5)

cat("\nAll TIFF files saved to", outdir, "\n")
