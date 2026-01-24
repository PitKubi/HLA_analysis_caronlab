# =============================================================================
# Figure 2: Peptide Length Distribution by PTM Type (Ridgeline Plot)
# =============================================================================
# Data source: data_ptm_sites.csv, data_background.csv (from data_loader.R)
# =============================================================================

library(ggplot2)
library(ggridges)
library(dplyr)

# -----------------------------------------------------------------------------
# Load data and config
# -----------------------------------------------------------------------------
config <- readRDS("figure_panels/config.rds")
ptm_sites <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)
background <- read.csv("figure_panels/data_background.csv", stringsAsFactors = FALSE)

output_dir <- "figure_panels"
dir.create(output_dir, showWarnings = FALSE)

# =============================================================================
# Prepare length data
# =============================================================================

# Background (unmodified)
bg_lengths <- data.frame(
  PTM = "Unmodified",
  Length = background$Length,
  stringsAsFactors = FALSE
)

# PTM lengths
ptm_lengths <- ptm_sites %>%
  select(PTM, Length)

# Combine
all_lengths <- bind_rows(bg_lengths, ptm_lengths)

# Filter to peptide length <= 14
all_lengths <- all_lengths %>% filter(Length <= 14)

# Calculate median for ordering
median_lengths <- all_lengths %>%
  group_by(PTM) %>%
  summarise(median_len = median(Length, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_len))

# Order PTMs by median length
all_lengths$PTM <- factor(all_lengths$PTM, levels = median_lengths$PTM)

# Get colors - add gray for Unmodified
ptm_colors <- c(config$ptm_colors, "Unmodified" = "#CCCCCC")

# =============================================================================
# FIGURE 2: Ridgeline Plot
# =============================================================================

fig2 <- ggplot(all_lengths, aes(x = Length, y = PTM, fill = PTM)) +
  geom_density_ridges(
    scale = 2.5,
    rel_min_height = 0.01,
    alpha = 0.85,
    color = "white",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = ptm_colors, guide = "none") +
  scale_x_continuous(breaks = seq(8, 14, 1), limits = c(7, 15)) +
  labs(
    title = "Peptide Length Distribution by PTM Type",
    x = "Peptide Length (amino acids)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave(file.path(output_dir, "Figure2_Length_Distribution.png"), fig2,
       width = 10, height = 12, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure2_Length_Distribution.pdf"), fig2,
       width = 10, height = 12, bg = "white")

cat("Saved Figure2_Length_Distribution.png/pdf\n")

# Print summary
cat("\nLength distribution summary:\n")
all_lengths %>%
  group_by(PTM) %>%
  summarise(
    n = n(),
    median = median(Length),
    mean = round(mean(Length), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(median)) %>%
  print()
