# =============================================================================
# Figure 1: HLA-II Immunopeptidome Composition
# =============================================================================
# Figure 1A: Bar chart showing Modified vs Unmodified
# Figure 1B: Donut chart showing PTM distribution with residue breakdown
#
# Data source: data_summary.csv, data_ptm_sites.csv (from data_loader.R)
# =============================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)
library(geomtextpath)

source("../shared/ptm_config.R")

# -----------------------------------------------------------------------------
# Load data and config
# -----------------------------------------------------------------------------
config <- readRDS("figure_panels/config.rds")
summary_data <- read.csv("figure_panels/data_summary.csv", stringsAsFactors = FALSE)
ptm_sites <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)

# Filter to 10-25 mers (consistent with summary counts in config.rds)
ptm_sites <- ptm_sites[ptm_sites$Length >= 10 & ptm_sites$Length <= 25, ]

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
# COLOR PALETTE (from config)
# =============================================================================
ptm_colors <- config$ptm_colors

# Artifact Oxidation is shown in a separate panel (it dwarfs all other PTMs)
EXCLUDED_PTMS <- c("Artifact Oxidation")

# =============================================================================
# FIGURE 1A: Horizontal Bar Chart - Overall PTM Composition
# =============================================================================

# Exclude artifact oxidation from the composition bar (it is not a biological PTM)
n_artox <- nrow(ptm_sites[ptm_sites$PTM == "Artifact Oxidation", ])
n_modified_bio <- n_modified - n_artox
n_total_bio <- n_background + n_modified_bio

bar_data <- data.frame(
  category = c("Unmodified", "Modified"),
  count = c(n_background, n_modified_bio),
  percentage = c(n_background/n_total_bio * 100, n_modified_bio/n_total_bio * 100)
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
                    labels = c(paste0("Modified (", format(n_modified_bio, big.mark=","), ")"),
                              paste0("Unmodified (", format(n_background, big.mark=","), ")"))) +
  scale_x_continuous(breaks = seq(0, 100, 20), labels = seq(0, 100, 20),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
       x = "Percentage of Total Peptides",
       fill = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
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

# Exclude artifact oxidation from donut (shown separately)
ptm_sites_main <- ptm_sites %>% filter(!PTM %in% EXCLUDED_PTMS)

# Get PTM totals from data
ptm_totals <- ptm_sites_main %>%
  group_by(PTM) %>%
  summarise(total = n(), .groups = "drop") %>%
  arrange(desc(total))

# Get PTM + Residue breakdown
ptm_residue <- ptm_sites_main %>%
  group_by(PTM, Residue) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(PTM, desc(count))

# Canonical shared order so both classes' donuts are directly comparable
ptm_order <- PTM_ORDER[PTM_ORDER %in% ptm_totals$PTM]
ptm_colors_ordered <- ptm_colors[ptm_order]

# Create inner ring data (PTM types)
inner_df <- ptm_totals %>%
  rename(ptm = PTM, count = total) %>%
  mutate(ptm = factor(ptm, levels = ptm_order)) %>%
  arrange(ptm)

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
    n_residues = lengths(strsplit(residues, "/")),
    label = ifelse(display_label == "cluster",
                   ifelse(n_residues > 3,
                          paste0("Other (", sprintf("%.1f", pct_within_ptm), "%)"),
                          paste0(residues, " (", sprintf("%.1f", pct_within_ptm), "%)")),
                   residues),
    is_cluster = (display_label == "cluster")
  ) %>%
  select(-n_residues) %>%
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
  "SUMOylation" = "Sumo",
  "N-Glycosylation" = "Glyco",
  "Carbamidomethylation" = "Carbam",
  "Nitrosylation" = "Nitro"
)
inner_df$short_label <- short_labels[as.character(inner_df$ptm)]

# Separate PTMs by size for labeling
# Large PTMs get labels inside the ring; smaller ones go to center callout
large_ptms <- inner_df[inner_df$percentage >= 4.5, ]
callout_ptms <- inner_df[inner_df$percentage < 4.5, ]

# Prepare individual labels for small PTMs - placed outside the donut
# Get top residue letters for each small PTM
callout_res_list <- sapply(as.character(callout_ptms$ptm), function(p) {
  res <- ptm_residue[ptm_residue$PTM == p, ]
  res <- res[order(-res$count), ]
  top_res <- head(res$Residue, 5)
  paste(top_res, collapse = ", ")
})
# Two-line label: "Name (count)\nresidues"
callout_ptms$label_text <- paste0(
  callout_ptms$short_label, " (",
  format(callout_ptms$count, big.mark = ","), ")\n",
  callout_res_list
)
# Use PTM color for each segment line
callout_ptms$seg_color <- ptm_colors_ordered[as.character(callout_ptms$ptm)]

# Only label outer ring residues for PTMs >= 3.5% to reduce clutter
large_ptm_names <- inner_df$ptm[inner_df$percentage >= 3.5]
outer_labeled <- outer_final[outer_final$ptm %in% large_ptm_names, ]

# Separate outer labels: regular vs clustered
outer_regular <- outer_labeled[!outer_labeled$is_clustered, ]
outer_clustered_labels <- outer_labeled[outer_labeled$is_clustered, ]

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
  # Inner ring labels for LARGE PTMs (>= 4.5%)
  geom_text(data = large_ptms,
            aes(x = 1.75, y = label_pos, label = short_label),
            color = "white", fontface = "bold", size = 4.5) +
  # Small PTMs (Me, DiMe, Glyco, Carbam, Cit, Sumo) are NOT labeled here
  # They are shown in a separate zoomed donut (Figure 1B2) for clarity
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
  # Center text - positioned below callout
  annotate("text", x = 0, y = total_ptms/2 - total_ptms*0.30,
           label = format(total_ptms, big.mark = ","),
           fontface = "bold", size = 7, color = "black") +
  scale_fill_manual(values = ptm_colors_ordered,
                    labels = paste0(names(ptm_colors_ordered), " (",
                                   format(inner_df$count[match(names(ptm_colors_ordered), inner_df$ptm)], big.mark=","), ", ",
                                   sprintf("%.1f%%", inner_df$percentage[match(names(ptm_colors_ordered), inner_df$ptm)]), ")")) +
  coord_polar(theta = "y", start = 0, direction = 1) +
  xlim(0, 5.5) +
  ylim(0, total_ptms) +
  labs(fill = "PTM Types") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(5, 5, 5, 5)
  )

save_panel(fig1b, "Figure1B_Donut", "donut")
cat("Saved Figure1B_Donut.png/pdf/tif\n")

# =============================================================================
# FIGURE 1B2: Zoomed Donut - Small PTMs Only
# =============================================================================
# Shows Me, DiMe, Glyco, Carbam, Cit, Sumo as their own donut for clarity
# User will combine with main donut in Illustrator

# Data-driven: low-abundance PTMs (< 4.5% of modified) get the zoomed donut.
# (HLA-II abundances differ from HLA-I, so this is derived, not hardcoded.)
small_ptm_names <- as.character(inner_df$ptm[inner_df$percentage < 4.5])

# Inner ring: small PTMs only
zoom_inner <- inner_df[as.character(inner_df$ptm) %in% small_ptm_names, ]
zoom_levels <- ptm_order[ptm_order %in% as.character(zoom_inner$ptm)]
zoom_inner <- zoom_inner[match(zoom_levels, as.character(zoom_inner$ptm)), ]
zoom_inner$ptm <- factor(as.character(zoom_inner$ptm), levels = zoom_levels)

zoom_total <- sum(zoom_inner$count)
zoom_inner$percentage_zoom <- zoom_inner$count / zoom_total * 100

# Recalculate positions for closed donut
zoom_inner$ymax <- cumsum(zoom_inner$count)
zoom_inner$ymin <- c(0, head(zoom_inner$ymax, -1))
zoom_inner$label_pos <- (zoom_inner$ymin + zoom_inner$ymax) / 2

# All PTM labels will be curved text outside the donut
zoom_inner$text_color <- ptm_colors[as.character(zoom_inner$ptm)]

# Build arc paths for curved text (geomtextpath) for ALL PTMs
# Each arc traces the angular span of the PTM segment at a fixed radius just outside outer ring
curved_arcs <- do.call(rbind, lapply(1:nrow(zoom_inner), function(i) {
  y_seq <- seq(zoom_inner$ymin[i], zoom_inner$ymax[i], length.out = 50)
  data.frame(
    x = 4.0,
    y = y_seq,
    label = zoom_inner$short_label[i],
    text_color = zoom_inner$text_color[i],
    stringsAsFactors = FALSE
  )
}))

# -----------------------------------------------------------------------------
# Outer ring for zoom: re-derive from raw data with LOWER clustering threshold
# so residues are shown individually; always list residue letters (never "Other")
# -----------------------------------------------------------------------------
zoom_outer_raw <- ptm_residue %>%
  rename(ptm = PTM, residue = Residue) %>%
  filter(ptm %in% small_ptm_names) %>%
  mutate(ptm = factor(ptm, levels = levels(zoom_inner$ptm))) %>%
  group_by(ptm) %>%
  mutate(
    ptm_total = sum(count),
    pct_within_ptm = count / ptm_total * 100
  ) %>%
  ungroup()

# Cluster only residues < 3% within their PTM (much lower than main donut's 6%)
zoom_outer_clust <- zoom_outer_raw %>%
  group_by(ptm) %>%
  mutate(
    is_small = pct_within_ptm < 3,
    display_label = ifelse(is_small, "cluster", residue)
  ) %>%
  ungroup()

zoom_outer <- zoom_outer_clust %>%
  group_by(ptm, display_label) %>%
  summarise(
    count = sum(count),
    residues = paste(sort(unique(residue)), collapse = "/"),
    pct_within_ptm = sum(pct_within_ptm),
    .groups = "drop"
  ) %>%
  mutate(
    # Always list residue letters, never show "Other"
    label = ifelse(display_label == "cluster",
                   paste0(residues, " (", sprintf("%.1f", pct_within_ptm), "%)"),
                   residues),
    is_cluster = (display_label == "cluster")
  ) %>%
  # Sort: within each PTM, non-clustered desc by count, then clustered last
  arrange(ptm, is_cluster, desc(count)) %>%
  select(-is_cluster)

# Assign colors (lighter shades per residue segment)
zoom_outer$color <- NA
for (ptm_name in levels(zoom_inner$ptm)) {
  idx <- which(zoom_outer$ptm == ptm_name)
  if (length(idx) > 0) {
    base_color <- ptm_colors[ptm_name]
    for (i in seq_along(idx)) {
      factor <- min(0.5, 0.08 * (i - 1))
      zoom_outer$color[idx[i]] <- lighten(base_color, factor = factor)
    }
  }
}

# Recalculate positions
zoom_outer$ymax <- cumsum(zoom_outer$count)
zoom_outer$ymin <- c(0, head(zoom_outer$ymax, -1))
zoom_outer$label_pos <- (zoom_outer$ymin + zoom_outer$ymax) / 2

# Separate outer labels: single-residue vs cluster (long text)
zoom_outer_regular <- zoom_outer[zoom_outer$display_label != "cluster", ]
zoom_outer_clusters <- zoom_outer[zoom_outer$display_label == "cluster", ]

# Colors for zoomed inner ring
zoom_colors <- ptm_colors[as.character(levels(zoom_inner$ptm))]

fig1b2 <- ggplot() +
  # Outer ring (residues)
  geom_rect(data = zoom_outer,
            aes(xmin = 2.5, xmax = 3.5, ymin = ymin, ymax = ymax, fill = I(color)),
            color = "white", linewidth = 0.5) +
  # Inner ring (PTM types)
  geom_rect(data = zoom_inner,
            aes(xmin = 1.0, xmax = 2.5, ymin = ymin, ymax = ymax, fill = ptm),
            color = "white", linewidth = 0.5) +
  # All PTM labels as curved text along the donut arc, in their PTM color
  geom_textpath(data = curved_arcs,
                aes(x = x, y = y, label = label, group = label, colour = label),
                fontface = "bold", size = 4.0, text_only = TRUE,
                vjust = -0.3, show.legend = FALSE) +
  scale_colour_manual(values = setNames(zoom_inner$text_color,
                                        zoom_inner$short_label)) +
  # Outer ring labels - single residues (close to ring)
  geom_text_repel(data = zoom_outer_regular,
                  aes(x = 3.5, y = label_pos, label = label),
                  color = "#333333", fontface = "bold", size = 3.2,
                  nudge_x = 0.5,
                  direction = "y",
                  segment.size = 0.3,
                  segment.color = "grey50",
                  box.padding = 0.2,
                  point.padding = 0.05,
                  min.segment.length = 0.1,
                  max.overlaps = Inf,
                  force = 3) +
  # Outer ring labels - clusters (push much further out to avoid overlap)
  geom_text_repel(data = zoom_outer_clusters,
                  aes(x = 3.5, y = label_pos, label = label),
                  color = "#333333", fontface = "bold", size = 2.8,
                  nudge_x = 2.1,
                  direction = "y",
                  segment.size = 0.3,
                  segment.color = "grey50",
                  box.padding = 0.5,
                  point.padding = 0.2,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  force = 8) +
  # Center count
  annotate("text", x = 0, y = zoom_total / 2,
           label = format(zoom_total, big.mark = ","),
           fontface = "bold", size = 7, color = "black") +
  scale_fill_manual(values = zoom_colors,
                    labels = paste0(names(zoom_colors), " (",
                                   format(zoom_inner$count[match(names(zoom_colors), zoom_inner$ptm)], big.mark=","), ", ",
                                   sprintf("%.1f%%", zoom_inner$percentage[match(names(zoom_colors), zoom_inner$ptm)]), ")")) +
  coord_polar(theta = "y", start = 0, direction = 1) +
  xlim(0, 6.8) +
  ylim(0, zoom_total) +
  labs(fill = "PTM Types") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(5, 5, 5, 5)
  )

save_panel(fig1b2, "Figure1B2_Donut_Zoomed", "donut_zoom")
cat("Saved Figure1B2_Donut_Zoomed.png/pdf/tif\n")

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

cat("\nZoomed donut (Figure 1B2) breakdown:\n")
for (i in 1:nrow(zoom_inner)) {
  cat(sprintf("  %-25s %5d (%5.1f%% of total, %5.1f%% of zoomed)\n",
              as.character(zoom_inner$ptm[i]),
              zoom_inner$count[i],
              zoom_inner$percentage[i],
              zoom_inner$percentage_zoom[i]))
}

# =============================================================================
# FIGURE 1C: Artifact Oxidation Breakdown (separate panel)
# =============================================================================

artox_data <- ptm_sites %>%
  filter(PTM == "Artifact Oxidation") %>%
  group_by(Residue) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count)) %>%
  mutate(
    percentage = count / sum(count) * 100,
    label = paste0(Residue, "\n(n=", format(count, big.mark=","), ", ", sprintf("%.1f%%", percentage), ")")
  )

artox_colors <- c("M" = "#FF6F00", "W" = "#FFB300", "F" = "#E65100", "H" = "#BF360C")
n_artox_total <- sum(artox_data$count)

fig1c <- ggplot(artox_data, aes(x = reorder(Residue, -count), y = count, fill = Residue)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(format(count, big.mark=","), "\n(", sprintf("%.1f%%", percentage), ")")),
            vjust = -0.3, fontface = "bold", size = 4) +
  scale_fill_manual(values = artox_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::comma) +
  labs(
    x = "Oxidized Residue",
    y = "Number of Peptides"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(output_dir, "Figure1C_ArtifactOxidation.png"), fig1c, width = 7, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure1C_ArtifactOxidation.pdf"), fig1c, width = 7, height = 5, bg = "white")
cat("Saved Figure1C_ArtifactOxidation.png/pdf\n")
