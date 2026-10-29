# =============================================================================
# HLA-II PTM Analysis Pipeline - Master Script
# =============================================================================
# Run from the HLA-II/ directory:
#   Rscript run_all_figures_hla2.R
#
# PREREQUISITES:
#   1. R packages: readxl, dplyr, tidyr, ggplot2, circlize, pheatmap, cowplot,
#                  ggridges, scales, ggrepel, geomtextpath
#   2. Data file: HLA-II_JY_DDA_PTMs_ResultsSummary_05012026_4.xlsx (this dir)
#
# All figures land in ./figure_panels/ as PNG, PDF and (manuscript pairs) TIF.
# The HLA-I pipeline mirror lives at the repo root (run_all_figures.R).
# =============================================================================

xlsx_file <- "HLA-II_JY_DDA_PTMs_ResultsSummary_05012026_4.xlsx"
if (!file.exists(xlsx_file)) {
  stop("ERROR: Data file not found: ", xlsx_file,
       "\nRun this script from the HLA-II/ directory.")
}
dir.create("figure_panels", showWarnings = FALSE)

pipeline_scripts <- c(
  "scripts/data_loader_hla2.R",
  "scripts/Figure1_Composition.R",
  "scripts/Figure2_LengthDistribution.R",
  "scripts/Figure4A_Circos.R",
  "scripts/Figure4B_Position_Heatmap.R",
  "scripts/Figure4B2_Position_Heatmap_ByPTM.R",
  "scripts/Figure4C_PositionDensity.R",
  "scripts/Figure4C2_PositionDensity_Core.R",
  "scripts/Figure4X_Position_P1P25.R",
  "scripts/Figure4X3_Position_Core_P1P9.R",
  "scripts/Figure5A_Binders.R",
  "scripts/Figure5B_ELRank.R",
  "scripts/Figure6ABCD_SpecificPTMs.R",
  "scripts/Figure6E_Glycan.R"
)

results <- data.frame(Script = character(), Status = character(), Time = numeric())
for (script in pipeline_scripts) {
  cat("\n=== Running:", script, "===\n")
  if (!file.exists(script)) {
    results <- rbind(results, data.frame(Script = script, Status = "SKIPPED", Time = NA))
    next
  }
  start_time <- Sys.time()
  tryCatch({
    source(script, local = new.env())
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    results <- rbind(results, data.frame(Script = script, Status = "SUCCESS", Time = elapsed))
  }, error = function(e) {
    cat("\nERROR:", conditionMessage(e), "\n")
    results <<- rbind(results, data.frame(Script = script, Status = "FAILED",
      Time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
  })
}

cat("\n=== Pipeline summary ===\n")
for (i in seq_len(nrow(results))) {
  cat(sprintf("  [%s] %s (%.1fs)\n", results$Status[i], results$Script[i],
              ifelse(is.na(results$Time[i]), 0, results$Time[i])))
}
if (any(results$Status == "FAILED")) stop("One or more scripts failed")
