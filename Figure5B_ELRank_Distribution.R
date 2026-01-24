# =============================================================================
# Figure 5B: EL Rank Distribution by PTM Type
# =============================================================================
# Shows EL Rank distribution (0-2 range) comparing modified vs background
# Data source: data_binding.csv, data_background.csv (from data_loader.R)
# =============================================================================

library(ggplot2)
library(dplyr)
library(cowplot)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
binding_data <- read.csv("figure_panels/data_binding.csv", stringsAsFactors = FALSE)
background <- read.csv("figure_panels/data_background.csv", stringsAsFactors = FALSE)

output_dir <- "figure_panels"
dir.create(output_dir, showWarnings = FALSE)

cat("=== Figure 5B: EL Rank Distribution ===\n\n")

# Filter to length 8-14
binding_data <- binding_data %>% filter(Length >= 8 & Length <= 14)
background <- background %>% filter(Length >= 8 & Length <= 14)

# PTM colors (circos palette)
ptm_colors <- c(
  "Cysteinylation"  = "#2E7D32",
  "Deamidation"     = "#1565C0",
  "Phosphorylation" = "#6A1B9A",
  "Acetylation"     = "#E65100",
  "Ubiquitination"  = "#C62828",
  "Citrullination"  = "#AD1457",
  "Dimethylation"   = "#4527A0",
  "Methylation"     = "#558B2F",
  "Oxidation"       = "#0277BD",
  "N-Glycosylation" = "#00695C"
)

# =============================================================================
# Prepare background EL Rank data
# =============================================================================

# Background EL Rank (limit to 0-2)
bg_el <- background$EL_Rank[!is.na(background$EL_Rank)]
bg_el <- bg_el[bg_el <= 2]
n_bg <- length(bg_el)

cat("Background peptides with EL_Rank <= 2:", n_bg, "\n\n")

# =============================================================================
# Create histograms for each PTM
# =============================================================================

# Get PTMs sorted by count
ptm_counts <- binding_data %>%
  filter(!is.na(EL_Rank)) %>%
  group_by(PTM) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

# Select PTMs with enough data
ptms_to_plot <- ptm_counts %>% filter(n >= 50) %>% pull(PTM)

cat("PTMs to plot:\n")
print(ptm_counts %>% filter(PTM %in% ptms_to_plot))

plot_list <- list()

for (ptm in ptms_to_plot) {
  # Get PTM EL_Rank (limit to 0-2)
  ptm_el <- binding_data %>%
    filter(PTM == ptm, !is.na(EL_Rank), EL_Rank <= 2) %>%
    pull(EL_Rank)

  n_ptm <- length(ptm_el)

  if (n_ptm < 20) next

  # Create combined data for plotting
  plot_df <- data.frame(
    EL_Rank = c(bg_el, ptm_el),
    Type = c(rep("Background", length(bg_el)), rep(ptm, length(ptm_el)))
  )

  # Get color for this PTM
  ptm_color <- ifelse(ptm %in% names(ptm_colors), ptm_colors[ptm], "#666666")

  p <- ggplot() +
    # Background histogram (gray filled)
    geom_histogram(data = data.frame(EL_Rank = bg_el),
                   aes(x = EL_Rank, y = after_stat(density)),
                   bins = 40, fill = "gray70", color = "gray50",
                   linewidth = 0.2, alpha = 0.7) +
    # PTM histogram (colored outline only)
    geom_histogram(data = data.frame(EL_Rank = ptm_el),
                   aes(x = EL_Rank, y = after_stat(density)),
                   bins = 40, fill = NA, color = ptm_color,
                   linewidth = 0.8) +
    # Threshold lines
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "darkgreen", linewidth = 0.5) +
    geom_vline(xintercept = 2.0, linetype = "dashed", color = "darkorange", linewidth = 0.5) +
    # Labels
    labs(
      title = ptm,
      x = "EL Rank",
      y = "Density"
    ) +
    # Add legend text
    annotate("text", x = 1.5, y = Inf, vjust = 2, hjust = 0.5, size = 2.5,
             label = paste0("Background, n = ", format(n_bg, big.mark = ",")), color = "gray40") +
    annotate("text", x = 1.5, y = Inf, vjust = 3.5, hjust = 0.5, size = 2.5,
             label = paste0(ptm, ", n = ", format(n_ptm, big.mark = ",")), color = ptm_color) +
    scale_x_continuous(limits = c(0, 2), breaks = c(0, 0.5, 1, 1.5, 2)) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5, color = ptm_color),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      plot.margin = margin(5, 10, 5, 5)
    )

  plot_list[[ptm]] <- p
  cat("Created histogram for", ptm, "(n =", n_ptm, ")\n")
}

# =============================================================================
# Combine plots
# =============================================================================

n_plots <- length(plot_list)
n_cols <- 3
n_rows <- ceiling(n_plots / n_cols)

combined_plot <- plot_grid(
  plotlist = plot_list,
  ncol = n_cols,
  align = "hv"
)

# Add title
title_gg <- ggdraw() +
  draw_label("EL Rank Distribution: Modified vs Background (EL Rank 0-2)",
             fontface = "bold", size = 14, hjust = 0.5)

subtitle_gg <- ggdraw() +
  draw_label("Dashed lines: Strong binder (<0.5%), Weak binder (<2%)",
             size = 10, hjust = 0.5, color = "gray40")

final_plot <- plot_grid(
  title_gg,
  subtitle_gg,
  combined_plot,
  ncol = 1,
  rel_heights = c(0.05, 0.03, 1)
)

# Save
ggsave(file.path(output_dir, "Figure5B_ELRank_Distribution.png"), final_plot,
       width = 4 * n_cols, height = 3.5 * n_rows + 1, units = "in", dpi = 300, bg = "white")

ggsave(file.path(output_dir, "Figure5B_ELRank_Distribution.pdf"), final_plot,
       width = 4 * n_cols, height = 3.5 * n_rows + 1, units = "in", bg = "white")

cat("\nSaved Figure5B_ELRank_Distribution.png/pdf\n")

# =============================================================================
# Summary statistics
# =============================================================================
cat("\n=== Summary ===\n")
cat("Strong binder threshold: EL_Rank < 0.5%\n")
cat("Weak binder threshold: 0.5% <= EL_Rank < 2%\n")
cat("Non-binder: EL_Rank >= 2%\n\n")

# Calculate binder percentages
binder_summary <- binding_data %>%
  filter(!is.na(EL_Rank)) %>%
  group_by(PTM) %>%
  summarise(
    n_total = n(),
    pct_strong = round(100 * sum(EL_Rank < 0.5) / n(), 1),
    pct_weak = round(100 * sum(EL_Rank >= 0.5 & EL_Rank < 2) / n(), 1),
    pct_nonbinder = round(100 * sum(EL_Rank >= 2) / n(), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_strong))

cat("Binder percentages by PTM:\n")
print(binder_summary)
