# =============================================================================
# HLA-I PTM Analysis Pipeline - Master Script
# =============================================================================
# This script runs all figure generation scripts in the correct order.
#
# USAGE:
#   Rscript run_all_figures.R
#
# PREREQUISITES:
#   1. R packages: readxl, dplyr, tidyr, ggplot2, circlize, pheatmap, cowplot,
#                  ggridges, scales
#   2. Data file: HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx
#
# OUTPUT:
#   All figures are saved to ./figure_panels/ directory as PNG and PDF
#   Data tables are saved as CSV files
#
# PIPELINE ORDER:
#   1. data_loader.R        - Extract and prepare data from Excel
#   2. Figure1_Composition.R - PTM composition donut chart
#   3. Figure2_LengthDistribution.R - Peptide length ridgeline
#   4. Figure4A_Circos.R    - Circos plot (PTM-residue-allele)
#   5. Figure4B_Position_Heatmap.R - Position enrichment (by residue)
#   6. Figure4B2_Position_Heatmap_ByPTM.R - Position enrichment (by PTM)
#   7. Figure4C_PositionDensity.R - PTM position density curves
#   8. Figure5A_StrongBinders.R - Binder heatmaps
#   9. Figure5B_ELRank_Distribution.R - EL Rank distributions
#  10. Figure6_SpecificPTMs.R - Specific PTM analyses (6A-6C)
#  11. Figure6DE_Improved.R - Improved 6D and 6E analyses
#
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("HLA-I PTM Analysis Pipeline\n")
cat("=============================================================================\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
xlsx_file <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"
output_dir <- "figure_panels"

# Check for required data file
if (!file.exists(xlsx_file)) {
  stop("ERROR: Data file not found: ", xlsx_file, "\n",
       "Please ensure the Excel file is in the working directory.")
}

# Create output directory
dir.create(output_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Check/Install Required Packages
# -----------------------------------------------------------------------------
required_packages <- c("readxl", "dplyr", "tidyr", "ggplot2", "circlize",
                       "pheatmap", "cowplot", "ggridges", "scales")

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  cat("Missing packages:", paste(missing_packages, collapse = ", "), "\n")
  cat("Installing missing packages...\n")
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

# -----------------------------------------------------------------------------
# Define Pipeline Scripts
# -----------------------------------------------------------------------------
pipeline_scripts <- c(
  "data_loader.R",
  "Figure1_Composition.R",
  "Figure2_LengthDistribution.R",
  "Figure4A_Circos.R",
  "Figure4B_Position_Heatmap.R",
  "Figure4B2_Position_Heatmap_ByPTM.R",
  "Figure4C_PositionDensity.R",
  "Figure5A_StrongBinders.R",
  "Figure5B_ELRank_Distribution.R",
  "Figure6_SpecificPTMs.R",
  "Figure6DE_Improved.R"
)

# -----------------------------------------------------------------------------
# Run Pipeline
# -----------------------------------------------------------------------------
results <- data.frame(
  Script = character(),
  Status = character(),
  Time = numeric(),
  stringsAsFactors = FALSE
)

for (script in pipeline_scripts) {
  cat("\n")
  cat("=============================================================================\n")
  cat("Running:", script, "\n")
  cat("=============================================================================\n")

  if (!file.exists(script)) {
    cat("WARNING: Script not found, skipping.\n")
    results <- rbind(results, data.frame(Script = script, Status = "SKIPPED", Time = NA))
    next
  }

  start_time <- Sys.time()
  tryCatch({
    source(script, local = new.env())
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat("\nCompleted in", round(elapsed, 1), "seconds\n")
    results <- rbind(results, data.frame(Script = script, Status = "SUCCESS", Time = elapsed))
  }, error = function(e) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat("\nERROR:", conditionMessage(e), "\n")
    results <- rbind(results, data.frame(Script = script, Status = "FAILED", Time = elapsed))
  })
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n")
cat("=============================================================================\n")
cat("Pipeline Complete\n")
cat("=============================================================================\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Script Results:\n")
for (i in 1:nrow(results)) {
  status_symbol <- ifelse(results$Status[i] == "SUCCESS", "[OK]",
                          ifelse(results$Status[i] == "FAILED", "[FAIL]", "[SKIP]"))
  time_str <- ifelse(is.na(results$Time[i]), "", paste0(" (", round(results$Time[i], 1), "s)"))
  cat(sprintf("  %s %s%s\n", status_symbol, results$Script[i], time_str))
}

# List generated figures
cat("\nGenerated Figures:\n")
figure_files <- list.files(output_dir, pattern = "\\.png$", full.names = FALSE)
for (f in sort(figure_files)) {
  cat("  -", f, "\n")
}

cat("\nGenerated Data Files:\n")
data_files <- list.files(output_dir, pattern = "\\.csv$", full.names = FALSE)
for (f in sort(data_files)) {
  cat("  -", f, "\n")
}

cat("\nAll outputs saved to:", output_dir, "/\n")
