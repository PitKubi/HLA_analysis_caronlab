# =============================================================================
# Figure 1: HLA-I Immunopeptidome Composition
# =============================================================================
# Figure 1A: Bar chart showing Modified vs Unmodified
# Figure 1B: Donut chart showing PTM distribution with residue breakdown
#
# Data source: data_summary.csv, data_ptm_sites.csv (from data_loader.R)
# =============================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)

# -----------------------------------------------------------------------------
# Load data and config
# -----------------------------------------------------------------------------
config <- readRDS("figure_panels/config.rds")
summary_data <- read.csv("figure_panels/data_summary.csv", stringsAsFactors = FALSE)
ptm_sites <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)

# Output directory
output_dir <- "figure_panels"
dir.create(output_dir, showWarnings = FALSE)

# Get totals from config
n_total <- config$n_total
n_modified <- config$n_modified
n_background <- config$n_background

cat("Loaded data:\n")
cat("  Total peptides:", n_total, "\n")
cat("  Modified:", n_modified, sprintf("(%.1f%%)\n", n_modified/n_total*100))
cat("  Unmodified:", n_background, sprintf("(%.1f%%)\n", n_background/n_total*100))

# =============================================================================
# COLOR PALETTE (from circos plot)
# =============================================================================
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
  "SUMOylation"     = "#795548",
  "N-Glycosylation" = "#00695C"
)

# =============================================================================
# FIGURE 1A: Horizontal Bar Chart - Overall PTM Composition
# =============================================================================

bar_data <- data.frame(
  category = c("Unmodified", "Modified"),
  count = c(n_background, n_modified),
  percentage = c(n_background/n_total * 100, n_modified/n_total * 100)
)

bar_data$xmin <- c(0, bar_data$percentage[1])
bar_data$xmax <- c(bar_data$percentage[1], 100)
bar_data$label_pos <- (bar_data$xmin + bar_data$xmax) / 2

fig1a <- ggplot(bar_data) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = category),
            color = "white", linewidth = 0.5) +
  geom_text(aes(x = label_pos, y = 0.5,
                label = sprintf("%.1f%%", percentage)),
            fontface = "bold", size = 5,
            color = c("#873600", "white")) +
  scale_fill_manual(values = c("Modified" = "#E74C3C", "Unmodified" = "#FDF2E9"),
                    labels = c(paste0("Modified (", format(n_modified, big.mark=","), ")"),
                              paste0("Unmodified (", format(n_background, big.mark=","), ")"))) +
  scale_x_continuous(breaks = seq(0, 100, 20), labels = seq(0, 100, 20),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = paste0("HLA-I Immunopeptidome Composition (n = ",
                      format(n_total, big.mark=","), " peptides)"),
       x = "Percentage of Total Peptides",
       fill = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(output_dir, "Figure1A_BarChart.png"), fig1a, width = 10, height = 3, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure1A_BarChart.pdf"), fig1a, width = 10, height = 3, bg = "white")
cat("Saved Figure1A_BarChart.png/pdf\n")

# =============================================================================
# FIGURE 1B: Donut Chart - PTM Distribution with Residue Breakdown
# =============================================================================

# Get PTM totals from data
ptm_totals <- ptm_sites %>%
  group_by(PTM) %>%
  summarise(total = n(), .groups = "drop") %>%
  arrange(desc(total))

# Get PTM + Residue breakdown
ptm_residue <- ptm_sites %>%
  group_by(PTM, Residue) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(PTM, desc(count))

# Reorder by abundance
ptm_order <- ptm_totals$PTM
ptm_colors_ordered <- ptm_colors[ptm_order]

# Create inner ring data (PTM types)
inner_df <- ptm_totals %>%
  rename(ptm = PTM, count = total) %>%
  mutate(ptm = factor(ptm, levels = ptm_order))

total_ptms <- sum(inner_df$count)
inner_df$percentage <- inner_df$count / total_ptms * 100

# Flag small PTMs for callout labels (but NO inflation - keep donut closed)
inner_df$is_small <- inner_df$percentage < 1.5

# -----------------------------------------------------------------------------
# Create outer ring data with clustering of small residues (<5% within PTM)
# -----------------------------------------------------------------------------
outer_df <- ptm_residue %>%
  rename(ptm = PTM, residue = Residue) %>%
  mutate(ptm = factor(ptm, levels = ptm_order)) %>%
  group_by(ptm) %>%
  mutate(
    ptm_total = sum(count),
    pct_within_ptm = count / ptm_total * 100
  ) %>%
  ungroup()

# For each PTM, cluster residues - use different thresholds
# Oxidation: <4%
# Others: <6%
outer_clustered <- outer_df %>%
  group_by(ptm) %>%
  mutate(
    threshold = ifelse(ptm == "Oxidation", 4, 6),
    is_small = pct_within_ptm < threshold,
    display_label = ifelse(is_small, "cluster", residue)
  ) %>%
  ungroup() %>%
  select(-threshold)

# Aggregate: keep large residues, combine small ones
# Sort so clustered (small) residues come LAST within each PTM
outer_final <- outer_clustered %>%
  group_by(ptm, display_label) %>%
  summarise(
    count = sum(count),
    residues = paste(sort(unique(residue)), collapse = "/"),
    pct_within_ptm = sum(pct_within_ptm),
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(display_label == "cluster",
                   paste0(residues, " (", sprintf("%.1f", pct_within_ptm), "%)"),
                   residues),
    is_cluster = (display_label == "cluster")
  ) %>%
  arrange(ptm, is_cluster, desc(count)) %>%
  select(-is_cluster)

# Function to lighten colors (keep color saturation, avoid grey)
lighten <- function(color, factor = 0.3) {
  col <- col2rgb(color) / 255
  # Lighter but maintain hue - don't go too close to white/grey
  col <- col + (1 - col) * factor * 0.6  # Reduced lightening
  col <- pmin(col, 0.70)  # Strict cap to keep color visible
  rgb(col[1], col[2], col[3])
}

# Assign colors to outer ring (lighter shades for residues, but keep colored)
outer_final$color <- NA
for (ptm_name in ptm_order) {
  idx <- which(outer_final$ptm == ptm_name)
  if (length(idx) > 0) {
    base_color <- ptm_colors_ordered[ptm_name]
    n_segments <- length(idx)
    for (i in seq_along(idx)) {
      # Scale factor based on position, max out at 0.5 to keep colors visible
      factor <- min(0.5, 0.08 * (i - 1))
      outer_final$color[idx[i]] <- lighten(base_color, factor = factor)
    }
  }
}

# Calculate positions for donut chart (use actual counts - no gaps)
inner_df$ymax <- cumsum(inner_df$count)
inner_df$ymin <- c(0, head(inner_df$ymax, -1))
inner_df$label_pos <- (inner_df$ymin + inner_df$ymax) / 2

outer_final$ymax <- cumsum(outer_final$count)
outer_final$ymin <- c(0, head(outer_final$ymax, -1))
outer_final$label_pos <- (outer_final$ymin + outer_final$ymax) / 2

# Mark clustered labels for special handling
outer_final$is_clustered <- (outer_final$display_label == "cluster")

# Short labels for inner ring
short_labels <- c(
  "Cysteinylation" = "Cys",
  "Deamidation" = "Deam",
  "Oxidation" = "Ox",
  "Phosphorylation" = "Phos",
  "Acetylation" = "Ac",
  "Methylation" = "Me",
  "Ubiquitination" = "Ub",
  "Dimethylation" = "DiMe",
  "Citrullination" = "Cit",
  "SUMOylation" = "Su",
  "N-Glycosylation" = "Glyco"
)
inner_df$short_label <- short_labels[as.character(inner_df$ptm)]

# Separate PTMs by size for labeling
large_ptms <- inner_df[inner_df$percentage >= 5, ]
medium_ptms <- inner_df[inner_df$percentage >= 2.5 & inner_df$percentage < 5, ]
small_ptms <- inner_df[inner_df$is_small, ]

# Pre-compute positions for small PTMs - spread them apart
if (nrow(small_ptms) > 0) {
  small_ptms$row_idx <- 1:nrow(small_ptms)
  # Offset label positions to avoid overlap
  small_ptms$label_y <- small_ptms$label_pos + (small_ptms$row_idx - (nrow(small_ptms)+1)/2) * total_ptms * 0.08
}

# Separate outer labels: regular vs clustered
outer_regular <- outer_final[!outer_final$is_clustered, ]
outer_clustered_labels <- outer_final[outer_final$is_clustered, ]

# Ensure outer ring sums to exactly the same as inner ring
outer_final$ymax[nrow(outer_final)] <- total_ptms

# Plot donut chart
fig1b <- ggplot() +
  # Outer ring (residues) - thinner borders
  geom_rect(data = outer_final,
            aes(xmin = 2.5, xmax = 3.5, ymin = ymin, ymax = ymax, fill = I(color)),
            color = "white", linewidth = 0.3) +
  # Inner ring (PTM types) - thinner borders
  geom_rect(data = inner_df,
            aes(xmin = 1.0, xmax = 2.5, ymin = ymin, ymax = ymax, fill = ptm),
            color = "white", linewidth = 0.5) +
  # Inner ring labels for LARGE PTMs (>= 5%)
  geom_text(data = large_ptms,
            aes(x = 1.75, y = label_pos, label = short_label),
            color = "white", fontface = "bold", size = 4.5) +
  # Medium PTMs (2.5-5%) - smaller font inside ring
  # Adjust Cit label position slightly up
  geom_text(data = medium_ptms,
            aes(x = 1.75,
                y = label_pos + ifelse(short_label == "Cit", total_ptms * 0.012, 0),
                label = short_label),
            color = "white", fontface = "bold", size = 3.2) +
  # Small PTMs - place at small radius, upper portion of circle (staggered angularly)
  annotate("text", x = 0.5, y = total_ptms * 0.92,
           label = paste0(small_ptms$short_label[1], " (", small_ptms$count[1], ")"),
           color = "#333333", fontface = "bold", size = 2.8, hjust = 0.5) +
  annotate("text", x = 0.5, y = total_ptms * 0.85,
           label = if(nrow(small_ptms) > 1) paste0(small_ptms$short_label[2], " (", small_ptms$count[2], ")") else "",
           color = "#333333", fontface = "bold", size = 2.8, hjust = 0.5) +
  # Lines from labels to their segments
  annotate("segment", x = 0.6, xend = 1.0,
           y = total_ptms * 0.92, yend = small_ptms$label_pos[1],
           color = "grey40", linewidth = 0.4) +
  annotate("segment", x = 0.6, xend = 1.0,
           y = total_ptms * 0.85, yend = if(nrow(small_ptms) > 1) small_ptms$label_pos[2] else total_ptms * 0.85,
           color = "grey40", linewidth = 0.4) +
  # Outer ring labels - regular residues (close to ring)
  geom_text_repel(data = outer_regular,
                  aes(x = 3.5, y = label_pos, label = label),
                  color = "#333333", fontface = "bold", size = 3.2,
                  nudge_x = 0.4,
                  direction = "y",
                  segment.size = 0.3,
                  segment.color = "grey50",
                  box.padding = 0.1,
                  point.padding = 0.05,
                  min.segment.length = 0.1,
                  max.overlaps = Inf,
                  force = 1.5) +
  # Outer ring labels - clustered residues (push much further out)
  geom_text_repel(data = outer_clustered_labels,
                  aes(x = 3.5, y = label_pos, label = label),
                  color = "#333333", fontface = "bold", size = 2.5,
                  nudge_x = 1.2,
                  direction = "y",
                  segment.size = 0.3,
                  segment.color = "grey50",
                  box.padding = 0.2,
                  point.padding = 0.1,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  force = 3) +
  # Center text - smaller and lower to avoid overlap
  annotate("text", x = 0, y = total_ptms/2 - total_ptms*0.18,
           label = paste0("\nPTM\n", format(total_ptms, big.mark = ",")),
           fontface = "bold", size = 3.5, color = "#333333", lineheight = 1.1) +
  scale_fill_manual(values = ptm_colors_ordered,
                    labels = paste0(names(ptm_colors_ordered), " (",
                                   format(inner_df$count[match(names(ptm_colors_ordered), inner_df$ptm)], big.mark=","), ", ",
                                   sprintf("%.1f%%", inner_df$percentage[match(names(ptm_colors_ordered), inner_df$ptm)]), ")")) +
  coord_polar(theta = "y", start = 0, direction = 1) +
  xlim(0, 5.5) +
  ylim(0, total_ptms) +
  labs(title = "HLA-I PTM Distribution by Modification Type and Residue",
       fill = "PTM Types") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(file.path(output_dir, "Figure1B_Donut.png"), fig1b, width = 11, height = 9, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure1B_Donut.pdf"), fig1b, width = 11, height = 9, bg = "white")
cat("Saved Figure1B_Donut.png/pdf\n")

# =============================================================================
# Print summary
# =============================================================================
cat("\nPTM breakdown in Figure 1B:\n")
for (i in 1:nrow(inner_df)) {
  cat(sprintf("  %-15s %5d (%5.1f%%)\n",
              as.character(inner_df$ptm[i]),
              inner_df$count[i],
              inner_df$percentage[i]))
}
cat("\nTotal PTM sites:", total_ptms, "\n")

cat("\nOuter ring (residues, clustered <5%):\n")
for (i in 1:nrow(outer_final)) {
  cat(sprintf("  %-15s %-20s %5d (%5.1f%% of PTM)\n",
              as.character(outer_final$ptm[i]),
              outer_final$label[i],
              outer_final$count[i],
              outer_final$pct_within_ptm[i]))
}
