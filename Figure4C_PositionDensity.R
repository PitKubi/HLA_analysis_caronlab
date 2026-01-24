# =============================================================================
# Figure 4C: PTM Site Distribution vs Background
# =============================================================================
# Shows comparison of where PTMs occur vs where the amino acid occurs in background
# Data source: data_ptm_sites.csv, data_background.csv (from data_loader.R)
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# Load data and config
# -----------------------------------------------------------------------------
config <- readRDS("figure_panels/config.rds")
ptm_sites <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)
background <- read.csv("figure_panels/data_background.csv", stringsAsFactors = FALSE)

output_dir <- "figure_panels"
dir.create(output_dir, showWarnings = FALSE)

# Filter to length 8-14 as per meeting notes
ptm_sites <- ptm_sites %>% filter(Length >= 8 & Length <= 14)
background <- background %>% filter(Length >= 8 & Length <= 14)

# =============================================================================
# Calculate background amino acid position distributions
# =============================================================================

# Function to get AA position distribution from peptide sequences
get_aa_position_dist <- function(peptides, target_aa) {
  positions <- list()
  for (pep in peptides) {
    chars <- strsplit(pep, "")[[1]]
    for (i in seq_along(chars)) {
      if (chars[i] == target_aa) {
        positions[[length(positions) + 1]] <- i
      }
    }
  }
  unlist(positions)
}

# =============================================================================
# Select PTM+Residue combinations to show (most informative ones)
# =============================================================================

# Get PTM+Residue counts at positions 8-14 only
ptm_res_counts <- ptm_sites %>%
  filter(Site >= 8 & Site <= 14) %>%  # Only count modifications at positions 8-14
  group_by(PTM, Residue) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= 50) %>%  # At least 50 observations at positions 8-14
  arrange(desc(n))

# Select top combinations for display
selected <- ptm_res_counts %>%
  head(8)  # Show top 8

cat("Selected PTM+Residue combinations for Figure 4C:\n")
print(selected)

# =============================================================================
# Build comparison data
# =============================================================================

plot_data <- list()

for (i in 1:nrow(selected)) {
  ptm_name <- selected$PTM[i]
  residue <- selected$Residue[i]

  # Get modified positions
  mod_positions <- ptm_sites %>%
    filter(PTM == ptm_name, Residue == residue) %>%
    pull(Site)

  # Get background positions for this AA
  if (residue == "N-term") {
    bg_positions <- rep(1, nrow(background))  # N-term is always position 1
  } else {
    bg_positions <- get_aa_position_dist(background$Peptide, residue)
  }

  # Filter to positions 8-14
  mod_positions <- mod_positions[mod_positions >= 8 & mod_positions <= 14]
  bg_positions <- bg_positions[bg_positions >= 8 & bg_positions <= 14]

  # Create density data for positions 8-14
  mod_counts <- tabulate(mod_positions, nbins = 14)[8:14]
  bg_counts <- tabulate(bg_positions, nbins = 14)[8:14]

  mod_density <- data.frame(
    position = 8:14,
    density = mod_counts / sum(mod_counts),
    type = "Modified",
    panel = paste(residue, ptm_name),
    ptm = ptm_name,
    residue = residue
  )

  bg_density <- data.frame(
    position = 8:14,
    density = bg_counts / sum(bg_counts),
    type = "Background",
    panel = paste(residue, ptm_name),
    ptm = ptm_name,
    residue = residue
  )

  plot_data[[i]] <- bind_rows(mod_density, bg_density)
}

all_plot_data <- bind_rows(plot_data)
all_plot_data$panel <- factor(all_plot_data$panel, levels = unique(all_plot_data$panel))

# =============================================================================
# FIGURE 4C: Position Density Comparison
# =============================================================================

# PTM colors (consistent with circos plot)
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

fig4c <- ggplot(all_plot_data, aes(x = position, y = density)) +
  # Background as gray filled area
  geom_area(data = all_plot_data %>% filter(type == "Background"),
            fill = "gray80", alpha = 0.7) +
  geom_line(data = all_plot_data %>% filter(type == "Background"),
            color = "gray40", linewidth = 0.8) +
  # Modified as colored line
  geom_area(data = all_plot_data %>% filter(type == "Modified"),
            aes(fill = ptm), alpha = 0.5) +
  geom_line(data = all_plot_data %>% filter(type == "Modified"),
            aes(color = ptm), linewidth = 1) +
  scale_fill_manual(values = ptm_colors, guide = "none") +
  scale_color_manual(values = ptm_colors, guide = "none") +
  scale_x_continuous(breaks = 8:14, limits = c(7.5, 14.5)) +
  facet_wrap(~ panel, ncol = 2, scales = "free_y") +
  labs(
    title = "Modified vs Background: PTM Site Position Distribution",
    subtitle = "Gray = background AA distribution, Colored = modification site (Top 8 PTM+Residue with \u226550 occurrences)",
    x = "Peptide Position",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    strip.text = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(output_dir, "Figure4C_Position_Density.png"), fig4c,
       width = 10, height = 12, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure4C_Position_Density.pdf"), fig4c,
       width = 10, height = 12, bg = "white")

cat("Saved Figure4C_Position_Density.png/pdf\n")
