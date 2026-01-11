# =============================================================================
# FIGURE 4C: PTM Site Distribution vs Background
# =============================================================================
# Density plots comparing modified AA position distribution vs background
# Matching style of ptm_site_distrubution_vs_background.png example

library(readxl)
library(ggplot2)
library(dplyr)
library(cowplot)

# =============================================================================
# CONFIGURATION
# =============================================================================
xlsx_file <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"
positions <- 1:15  # Focus on positions 1-15 (match circos plot)

# PTM+Residue combinations to plot - selected for distinct patterns
# Sorted by correlation (most different from background first)
ptms_to_plot <- c(
  "T ubiquitination",     # Lowest corr (-0.04) - very different from background
  "K acetylation",        # Low corr (0.27)
  "S ubiquitination",     # Low corr (0.31), n=163
  "S phosphorylation",    # Most common phospho, position 4 enrichment
  "C cysteinylation",     # Large n, position 1 enrichment
  "M oxidation",          # Common oxidation
  "R methylation",        # Methylation example
  "N deamidation"         # Deamidation
)

# Panel labels
panel_labels <- c("a", "b", "c", "d", "e", "f", "g", "h")

# =============================================================================
# LOAD DATA
# =============================================================================
cat("=== Figure 4C: PTM Site Distribution vs Background ===\n\n")

# Load PTM data
ptm_data <- read.csv("figure_panels/data_figure4A_circos.csv", stringsAsFactors = FALSE)
ptm_data$PTM_Residue <- paste(ptm_data$Residue, tolower(ptm_data$PTM))
cat("Loaded", nrow(ptm_data), "PTM records\n")

# Load background peptides
cat("Loading background peptides...\n")
suppressMessages({
  bg <- read_excel(xlsx_file, sheet = "Background (wo PTMs)")
})
bg_peptides <- bg$Peptide[!is.na(bg$Peptide)]
cat("Background peptides:", length(bg_peptides), "\n\n")

# =============================================================================
# CALCULATE BACKGROUND POSITION DISTRIBUTION FOR EACH AMINO ACID
# =============================================================================
aa_list <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

# For each AA, get its position distribution in background peptides
bg_position_data <- list()

for (aa in aa_list) {
  positions_found <- c()
  for (pep in bg_peptides) {
    # Find all positions of this AA in the peptide
    pep_chars <- strsplit(pep, "")[[1]]
    aa_positions <- which(pep_chars == aa)
    # Only keep positions 1-15
    aa_positions <- aa_positions[aa_positions <= 15]
    positions_found <- c(positions_found, aa_positions)
  }
  bg_position_data[[aa]] <- positions_found
}

cat("Background position distributions calculated\n\n")

# =============================================================================
# CREATE DENSITY PLOTS - MATCHING EXAMPLE STYLE
# =============================================================================
plot_list <- list()
panel_idx <- 1

for (ptm_res in ptms_to_plot) {
  # Get PTM data for this combination
  ptm_subset <- ptm_data[ptm_data$PTM_Residue == ptm_res, ]

  if (nrow(ptm_subset) < 10) {
    cat("Skipping", ptm_res, "- insufficient data\n")
    next
  }

  target_aa <- ptm_subset$Residue[1]
  ptm_name <- ptm_subset$PTM[1]

  # Get modified positions (filter to 1-15)
  mod_positions <- ptm_subset$Site[ptm_subset$Site >= 1 & ptm_subset$Site <= 15]

  # Get background positions for this AA
  bg_positions <- bg_position_data[[target_aa]]

  if (length(bg_positions) == 0) {
    cat("Skipping", ptm_res, "- no background data for", target_aa, "\n")
    next
  }

  # Calculate counts for legend
  n_mod <- length(mod_positions)
  n_bg <- length(bg_positions)

  # Create histogram data - normalize to density
  bg_hist <- hist(bg_positions, breaks = seq(0.5, 15.5, by = 1), plot = FALSE)
  mod_hist <- hist(mod_positions, breaks = seq(0.5, 15.5, by = 1), plot = FALSE)

  # Convert counts to density (proportion)
  bg_density <- bg_hist$counts / sum(bg_hist$counts)
  mod_density <- mod_hist$counts / sum(mod_hist$counts)

  # Create data frame for plotting
  hist_df <- data.frame(
    Position = 1:15,
    Background = bg_density,
    Modified = mod_density
  )

  # Determine y-axis max
  y_max <- max(c(bg_density, mod_density)) * 1.15

  # Create plot matching example style - STEP HISTOGRAMS
  # Legend labels
  bg_label <- paste0("Background, n = ", format(n_bg, big.mark = ","))
  mod_label <- paste0("Modified, n = ", format(n_mod, big.mark = ","))

  # Create step data for proper stair-step histogram appearance
  # Each bar needs: left edge at x-0.5, right edge at x+0.5
  step_df <- data.frame(
    x = c(0.5, rep(1:15, each = 2) + rep(c(-0.5, 0.5), 15), 15.5),
    bg_y = c(0, rep(bg_density, each = 2), 0),
    mod_y = c(0, rep(mod_density, each = 2), 0)
  )
  # Fix the step pattern - need proper staircase
  step_x <- numeric(0)
  step_bg <- numeric(0)
  step_mod <- numeric(0)
  for (i in 1:15) {
    step_x <- c(step_x, i - 0.5, i + 0.5)
    step_bg <- c(step_bg, bg_density[i], bg_density[i])
    step_mod <- c(step_mod, mod_density[i], mod_density[i])
  }
  # Add endpoints
  step_x <- c(0.5, step_x, 15.5)
  step_bg <- c(0, step_bg, 0)
  step_mod <- c(0, step_mod, 0)

  step_df <- data.frame(x = step_x, bg_y = step_bg, mod_y = step_mod)

  p <- ggplot(step_df) +
    # Background - gray filled area with black outline
    geom_polygon(aes(x = x, y = bg_y), fill = "gray80", color = NA) +
    geom_step(aes(x = x, y = bg_y), color = "black", linewidth = 0.5, direction = "mid") +
    # Modified - red filled area with red outline (semi-transparent)
    geom_polygon(aes(x = x, y = mod_y), fill = "#E53935", alpha = 0.5, color = NA) +
    geom_step(aes(x = x, y = mod_y), color = "#B71C1C", linewidth = 0.6, direction = "mid") +
    # Title
    ggtitle(paste0(target_aa, " ", tolower(ptm_name))) +
    # Axis labels
    labs(x = "Position", y = "Density") +
    # X-axis breaks
    scale_x_continuous(breaks = seq(1, 15, by = 2), limits = c(0.5, 15.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    # Theme matching example
    theme_classic() +
    theme(
      plot.title = element_text(size = 11, face = "plain", hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      axis.line = element_line(linewidth = 0.5),
      plot.margin = margin(5, 10, 5, 5)
    ) +
    # Panel label in top-left
    annotate("text", x = 1.2, y = y_max * 0.97,
             label = panel_labels[panel_idx], fontface = "bold", size = 5, hjust = 0) +
    # Legend box
    annotate("rect", xmin = 6.5, xmax = 15.3, ymin = y_max * 0.73, ymax = y_max * 0.98,
             fill = "white", color = "gray50", linewidth = 0.3) +
    # Background legend entry
    annotate("segment", x = 7, xend = 8.2, y = y_max * 0.91, yend = y_max * 0.91,
             color = "black", linewidth = 0.8) +
    annotate("text", x = 8.5, y = y_max * 0.91, label = bg_label,
             size = 2.5, hjust = 0) +
    # Modified legend entry
    annotate("segment", x = 7, xend = 8.2, y = y_max * 0.79, yend = y_max * 0.79,
             color = "#B71C1C", linewidth = 0.8) +
    annotate("text", x = 8.5, y = y_max * 0.79, label = mod_label,
             size = 2.5, hjust = 0)

  plot_list[[ptm_res]] <- p
  panel_idx <- panel_idx + 1
  cat("Created plot for", ptm_res, "(n =", n_mod, ")\n")
}

# =============================================================================
# ARRANGE AND SAVE PLOTS
# =============================================================================
cat("\nArranging plots...\n")

# Arrange in 2x4 grid
n_plots <- length(plot_list)
n_cols <- 4
n_rows <- ceiling(n_plots / n_cols)

# Combine plots in grid
combined_plot <- plot_grid(plotlist = plot_list, nrow = n_rows, ncol = n_cols, align = "hv")

# PNG
ggsave("figure_panels/Figure4C_Position_Density.png", combined_plot,
       width = 3.2 * n_cols, height = 3.2 * n_rows, units = "in", dpi = 300)
cat("Saved: figure_panels/Figure4C_Position_Density.png\n")

# PDF
ggsave("figure_panels/Figure4C_Position_Density.pdf", combined_plot,
       width = 3.2 * n_cols, height = 3.2 * n_rows, units = "in")
cat("Saved: figure_panels/Figure4C_Position_Density.pdf\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=== Figure 4C Summary ===\n")
cat("PTM types plotted:", length(plot_list), "\n")
cat("Position range: 1-15\n")
cat("Style: Matching ptm_site_distrubution_vs_background.png example\n")
cat("Comparison: Modified site position distribution vs background AA position distribution\n")
