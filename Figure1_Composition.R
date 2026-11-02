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
library(geomtextpath)

source("shared/ptm_config.R")

# -----------------------------------------------------------------------------
# Load data and config
# -----------------------------------------------------------------------------
config <- readRDS("figure_panels/config.rds")
summary_data <- read.csv("figure_panels/data_summary.csv", stringsAsFactors = FALSE)
ptm_sites <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)

# Filter to 8-14 mers (consistent with summary counts in config.rds)
ptm_sites <- ptm_sites[ptm_sites$Length >= 8 & ptm_sites$Length <= 14, ]

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

# A PTM below this share is too thin to carry a label inside the ring: it gets
# a leader-line name and a protruding outer segment instead
CALLOUT_MAX_PCT <- 4.5
thin_ptm_names <- as.character(inner_df$ptm[inner_df$percentage < CALLOUT_MAX_PCT])

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
    # Every residue label reads the same way: residues on the first line, the
    # share of that PTM on the second. Mixing bare letters with smaller
    # inline-percentage labels made the ring look inconsistent.
    label_text = ifelse(display_label == "cluster" & n_residues > 3,
                        "Other", residues),
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

outer_final$label <- outer_final$label_text

# The thin PTMs are hard to find in the ring, so their outer segment is drawn
# protruding past the others. Only the radius changes; the angular width still
# carries the count, so the proportions stay honest.
POPOUT_RADIUS <- 4.0
outer_final$x_outer <- ifelse(as.character(outer_final$ptm) %in% thin_ptm_names,
                              POPOUT_RADIUS, 3.5)

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
large_ptms <- inner_df[inner_df$percentage >= CALLOUT_MAX_PCT, ]
callout_ptms <- inner_df[inner_df$percentage < CALLOUT_MAX_PCT, ]

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

#' Spread label positions along the angular axis so neighbouring labels keep a
#' minimum gap, then recentre the block on its original span. Thin PTM wedges
#' cluster tightly, and ggrepel's own force cannot separate them inside
#' coord_polar without dragging labels across the donut.
declutter_positions <- function(pos, min_gap, span_total) {
  ord <- order(pos)
  y <- pos[ord]
  # Spread each tight cluster around its own centre, so labels never travel
  # into a neighbouring part of the donut; clamp inside the circle so the
  # last cluster cannot run past the end and be clipped away.
  cluster <- cumsum(c(1, as.integer(diff(y) > 2 * min_gap)))
  for (k in unique(cluster)) {
    idx <- which(cluster == k)
    if (length(idx) == 1) next
    width <- (length(idx) - 1) * min_gap
    start <- mean(y[idx]) - width / 2
    start <- max(0, min(start, span_total - width))
    y[idx] <- start + (seq_along(idx) - 1) * min_gap
  }
  out <- numeric(length(pos))
  out[ord] <- y
  out
}
callout_ptms <- callout_ptms[order(callout_ptms$label_pos), ]
callout_ptms$label_y <- declutter_positions(callout_ptms$label_pos,
                                            total_ptms * 0.033, total_ptms)

# One cutoff decides where a PTM's residue breakdown is shown: at or above it
# the residues are labelled around the main donut, below it they go to the
# zoomed donut. Below this share a wedge sits in the same sector as the
# thin-PTM name callouts and the two label sets cross each other.
MAIN_DONUT_MIN_PCT <- 6
large_ptm_names <- inner_df$ptm[inner_df$percentage >= MAIN_DONUT_MIN_PCT]
outer_labeled <- outer_final[outer_final$ptm %in% large_ptm_names, ]

# One label set for the outer ring; clustered entries carry smaller text
outer_residue_labels <- outer_labeled
outer_residue_labels$label_size <- 3.0

# Ensure outer ring sums to exactly the same as inner ring
outer_final$ymax[nrow(outer_final)] <- total_ptms

# Plot donut chart
fig1b <- ggplot() +
  # Outer ring (residues) - thinner borders
  geom_rect(data = outer_final,
            aes(xmin = 2.5, xmax = x_outer, ymin = ymin, ymax = ymax, fill = I(color)),
            color = "white", linewidth = 0.3) +
  # Inner ring (PTM types) - thinner borders
  geom_rect(data = inner_df,
            aes(xmin = 1.0, xmax = 2.5, ymin = ymin, ymax = ymax, fill = ptm),
            color = "white", linewidth = 0.5) +
  # Inner ring labels for LARGE PTMs (>= 4.5%)
  geom_text(data = large_ptms,
            aes(x = 1.75, y = label_pos, label = short_label),
            color = "white", fontface = "bold", size = 4.5) +
  # Small PTM segments get a leader-line label in their own color: they are
  # too thin for in-ring text, and leaving them unlabeled made the inner ring
  # unreadable there (Moh, 2026-11-02). The zoomed donut (Figure 1B2) still
  # carries their residue breakdown.
  geom_text_repel(data = callout_ptms,
                  aes(x = POPOUT_RADIUS, y = label_pos, label = short_label,
                      colour = as.character(ptm)),
                  fontface = "bold", size = 3.6,
                  nudge_x = 1.7,
                  nudge_y = callout_ptms$label_y - callout_ptms$label_pos,
                  direction = "y",
                  segment.size = 0.4,
                  box.padding = 0.1,
                  point.padding = 0.1,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  force = 0.4,
                  seed = 42,
                  show.legend = FALSE) +
  scale_colour_manual(values = ptm_colors_ordered, guide = "none") +
  # Outer ring residue labels, single repel layer so clustered and single
  # residues are decluttered against each other, and all of them stay in a
  # narrow band next to the ring. Pushing the clustered ones further out put
  # them in the PTM-label band, where the PTM leader lines cut through them.
  geom_text_repel(data = outer_residue_labels,
                  aes(x = 3.5, y = label_pos, label = label, size = label_size),
                  color = "black", fontface = "bold", lineheight = 0.85,
                  nudge_x = 0.45,
                  direction = "y",
                  segment.size = 0.3,
                  segment.color = "grey50",
                  box.padding = 0.12,
                  point.padding = 0.05,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  force = 3,
                  seed = 42) +
  scale_size_identity() +
  # Center text - positioned below callout
  annotate("text", x = 0, y = total_ptms/2 - total_ptms*0.30,
           label = format(total_ptms, big.mark = ","),
           fontface = "bold", size = 7, color = "black") +
  scale_fill_manual(values = ptm_colors_ordered,
                    labels = paste0(names(ptm_colors_ordered), " (",
                                   format(inner_df$count[match(names(ptm_colors_ordered), inner_df$ptm)], big.mark=","), ", ",
                                   sprintf("%.1f%%", inner_df$percentage[match(names(ptm_colors_ordered), inner_df$ptm)]), ")")) +
  coord_polar(theta = "y", start = 0, direction = 1) +
  xlim(0, 6.2) +
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

# The complement of the main donut: every PTM whose residues are not labelled
# there gets its full residue breakdown here
small_ptm_names <- as.character(inner_df$ptm[inner_df$percentage < MAIN_DONUT_MIN_PCT])

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

# Curved text along the arc only where the wedge is wide enough to hold it;
# a sliver (e.g. 5 nitrosylated peptides) would print its label on top of its
# neighbour's. Thin wedges get the same decluttered leader-line callout as the
# main donut.
# Curved arc text is disabled in the zoomed donut: with residue labels now
# running to two lines, curved PTM names collided with them and with the thin
# wedges' callouts. Every PTM here uses the same decluttered leader line.
ZOOM_CURVED_MIN_PCT <- Inf
zoom_curved <- zoom_inner[zoom_inner$percentage_zoom >= ZOOM_CURVED_MIN_PCT, ]
zoom_callout <- zoom_inner[zoom_inner$percentage_zoom < ZOOM_CURVED_MIN_PCT, ]
curved_arcs <- if (nrow(zoom_curved) > 0) {
  do.call(rbind, lapply(seq_len(nrow(zoom_curved)), function(i) {
    y_seq <- seq(zoom_curved$ymin[i], zoom_curved$ymax[i], length.out = 50)
    data.frame(x = 5.8, y = y_seq, label = zoom_curved$short_label[i],
               text_color = zoom_curved$text_color[i], stringsAsFactors = FALSE)
  }))
} else {
  data.frame(x = numeric(0), y = numeric(0), label = character(0),
             text_color = character(0), stringsAsFactors = FALSE)
}
if (nrow(zoom_callout) > 0) {
  zoom_callout <- zoom_callout[order(zoom_callout$label_pos), ]
  zoom_callout$label_y <- declutter_positions(zoom_callout$label_pos,
                                              zoom_total * 0.05, zoom_total)
}

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
    label_text = residues,
    is_cluster = (display_label == "cluster")
  ) %>%
  # Sort: within each PTM, non-clustered desc by count, then clustered last
  arrange(ptm, is_cluster, desc(count)) %>%
  select(-is_cluster)

zoom_outer$label <- zoom_outer$label_text

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
# One label set, so single residues and clustered ones are decluttered against
# each other inside a single narrow band next to the ring
zoom_outer_labels <- zoom_outer
zoom_outer_labels$label_size <- 3.0

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
  geom_text_repel(data = zoom_callout,
                  aes(x = 2.5, y = label_pos, label = short_label,
                      colour = short_label),
                  fontface = "bold", size = 3.6,
                  nudge_x = 4.1,
                  nudge_y = zoom_callout$label_y - zoom_callout$label_pos,
                  direction = "y", segment.size = 0.4,
                  box.padding = 0.1, point.padding = 0.1,
                  min.segment.length = 0, max.overlaps = Inf,
                  force = 0.4, seed = 42, show.legend = FALSE) +
  scale_colour_manual(values = setNames(zoom_inner$text_color,
                                        zoom_inner$short_label)) +
  # Outer ring residue labels, one repel layer kept close to the ring
  geom_text_repel(data = zoom_outer_labels,
                  aes(x = 3.5, y = label_pos, label = label, size = label_size),
                  color = "black", fontface = "bold", lineheight = 0.85,
                  nudge_x = ifelse(zoom_outer_labels$display_label == "cluster", 1.3, 0.5),
                  direction = "y",
                  segment.size = 0.3,
                  segment.color = "grey50",
                  box.padding = 0.18,
                  point.padding = 0.05,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  force = 5,
                  seed = 42) +
  scale_size_identity() +
  # Center count
  annotate("text", x = 0, y = zoom_total / 2,
           label = format(zoom_total, big.mark = ","),
           fontface = "bold", size = 7, color = "black") +
  scale_fill_manual(values = zoom_colors,
                    labels = paste0(names(zoom_colors), " (",
                                   format(zoom_inner$count[match(names(zoom_colors), zoom_inner$ptm)], big.mark=","), ", ",
                                   sprintf("%.1f%%", zoom_inner$percentage[match(names(zoom_colors), zoom_inner$ptm)]), ")")) +
  coord_polar(theta = "y", start = 0, direction = 1) +
  xlim(0, 7.3) +
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
