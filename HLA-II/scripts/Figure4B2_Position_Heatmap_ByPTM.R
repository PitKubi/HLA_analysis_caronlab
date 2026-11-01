# =============================================================================
# FIGURE 4B2: PTM Position Enrichment Heatmap (Aggregated by PTM)
# =============================================================================
# Shows log2 enrichment of PTM positions vs background - all residues clustered per PTM
# Data source: data_ptm_sites.csv, data_background.csv (from data_loader.R)
# =============================================================================

library(dplyr)
library(pheatmap)

# =============================================================================
# CONFIGURATION
# =============================================================================
positions <- 10:25   # Peptide positions to analyze (length 10-25)
min_count <- 50     # Minimum occurrences for a PTM

output_dir <- "figure_panels"
dir.create(output_dir, showWarnings = FALSE)

# =============================================================================
# LOAD DATA
# =============================================================================
cat("=== Figure 4B2: PTM Position Enrichment (Aggregated by PTM) ===\n\n")

# Load PTM data
ptm_data <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)

# Filter to length 10-25
ptm_data <- ptm_data %>% filter(Length >= 10 & Length <= 25)
cat("Loaded", nrow(ptm_data), "PTM records (length 10-25)\n")

# Load background peptides
background <- read.csv("figure_panels/data_background.csv", stringsAsFactors = FALSE)

# Filter to length 10-25
background <- background %>% filter(Length >= 10 & Length <= 25)
bg_peptides <- background$Peptide
cat("Background peptides (length 10-25):", length(bg_peptides), "\n\n")

# =============================================================================
# CALCULATE OVERALL BACKGROUND POSITIONAL DISTRIBUTION
# =============================================================================
# For PTM-level analysis, we compare against overall position distribution
# (not per-AA, since we're aggregating across residues)

# Count all AAs at each position
bg_pos_counts <- numeric(length(positions))
names(bg_pos_counts) <- positions

for (pos in positions) {
  # Count peptides that have this position (i.e., length >= pos)
  n_at_pos <- sum(nchar(bg_peptides) >= pos)
  bg_pos_counts[as.character(pos)] <- n_at_pos
}

# Normalize to get distribution
bg_pos_freq <- bg_pos_counts / sum(bg_pos_counts)

cat("Background positional distribution calculated\n")
cat("Position 10:", sprintf("%.2f%%", 100*bg_pos_freq["10"]), "\n")
cat("Position 11:", sprintf("%.2f%%", 100*bg_pos_freq["11"]), "\n\n")

# =============================================================================
# CALCULATE PTM POSITION FREQUENCIES AND ENRICHMENT
# =============================================================================

# Filter to PTMs with sufficient counts
ptm_counts <- table(ptm_data$PTM)
valid_ptms <- names(ptm_counts[ptm_counts >= min_count])
cat("PTMs with >=", min_count, "occurrences:", length(valid_ptms), "\n")
cat("PTMs:", paste(valid_ptms, collapse = ", "), "\n\n")

# Initialize result matrices
enrichment_matrix <- matrix(NA, nrow = length(valid_ptms), ncol = length(positions))
rownames(enrichment_matrix) <- valid_ptms
colnames(enrichment_matrix) <- positions

correlation_vector <- numeric(length(valid_ptms))
names(correlation_vector) <- valid_ptms

count_vector <- numeric(length(valid_ptms))
names(count_vector) <- valid_ptms

# Calculate enrichment for each PTM
for (ptm in valid_ptms) {
  subset_data <- ptm_data[ptm_data$PTM == ptm, ]
  count_vector[ptm] <- nrow(subset_data)

  # Get position frequencies for this PTM
  pos_freq <- numeric(length(positions))
  names(pos_freq) <- positions

  for (pos in positions) {
    pos_freq[as.character(pos)] <- sum(subset_data$Site == pos, na.rm = TRUE)
  }
  # Normalize over the DISPLAYED window only, matching the background's
  # normalization; dividing by all rows (incl. sites outside the window)
  # imposed a spurious per-row negative offset = log2(window mass)
  n_in_window <- sum(pos_freq)
  if (n_in_window > 0) pos_freq <- pos_freq / n_in_window

  # Calculate log2 enrichment (add small pseudocount to avoid log(0))
  pseudocount <- 0.001
  enrichment <- log2((pos_freq + pseudocount) / (bg_pos_freq + pseudocount))
  enrichment_matrix[ptm, ] <- enrichment

  # Calculate correlation between PTM position distribution and background
  correlation_vector[ptm] <- cor(pos_freq, bg_pos_freq, use = "complete.obs")
}

# Replace infinite values with NA
enrichment_matrix[is.infinite(enrichment_matrix)] <- NA

cat("Calculated enrichment for", nrow(enrichment_matrix), "PTMs\n")

# =============================================================================
# PREPARE DATA FOR HEATMAP
# =============================================================================

# Manuscript rule: canonical PTM order from the shared config (identical in
# HLA-I and HLA-II), never per-figure by local counts.
config <- readRDS(file.path(output_dir, "config.rds"))
stopifnot(all(rownames(enrichment_matrix) %in% config$ptm_order))
sort_order <- order(match(rownames(enrichment_matrix), config$ptm_order))
enrichment_ordered <- enrichment_matrix[sort_order, , drop = FALSE]
correlation_ordered <- correlation_vector[sort_order]
count_ordered <- count_vector[sort_order]

cat("Sorted by canonical PTM order\n")

# Create combined matrix with correlation column
combined_matrix <- cbind(Corr = correlation_ordered, enrichment_ordered)

# Cap extreme values for better visualization
cap_value <- 3
combined_matrix[, -1] <- pmin(pmax(combined_matrix[, -1], -cap_value), cap_value)

# =============================================================================
# CREATE HEATMAP - PNG / PDF / TIF
# =============================================================================
cat("Generating heatmap...\n")

# Create row labels with counts
row_labels <- paste0(rownames(combined_matrix), " (n=", count_ordered, ")")

# Color palettes
enrichment_colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7",
                                         "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100)

# Create row annotation for correlation
row_annotation <- data.frame(
  Corr = combined_matrix[, "Corr"]
)
rownames(row_annotation) <- rownames(combined_matrix)

# Annotation colors
ann_colors <- list(
  Corr = colorRampPalette(c("white", "gray20"))(100)
)

# Calculate figure dimensions based on data
n_rows <- nrow(combined_matrix)
n_cols <- length(positions)
cell_w <- 30
cell_h <- 20
fig_width <- (n_cols * cell_w / 72) + 5
fig_height <- (n_rows * cell_h / 72) + 3

draw_enrichment_heatmap <- function() {
  pheatmap(
    combined_matrix[, -1],  # Exclude Corr column for main heatmap
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = enrichment_colors,
    breaks = seq(-cap_value, cap_value, length.out = 101),
    border_color = "white",
    cellwidth = cell_w,
    cellheight = cell_h,
    fontsize_row = 11,
    fontsize_col = 12,
    main = NA,  # Manuscript panels carry no titles
    labels_row = row_labels,
    labels_col = as.character(positions),
    legend = TRUE,
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3"),
    annotation_row = row_annotation,
    annotation_colors = ann_colors,
    annotation_legend = TRUE,
    annotation_names_row = TRUE
  )
}

open_device <- list(
  png = function(path) png(path, width = fig_width, height = fig_height,
                           units = "in", res = 300),
  pdf = function(path) pdf(path, width = fig_width, height = fig_height),
  tif = function(path) tiff(path, width = fig_width, height = fig_height,
                            units = "in", res = 300, compression = "lzw")
)
for (ext in names(open_device)) {
  out_path <- file.path(output_dir, paste0("Figure4B2_Position_Enrichment_ByPTM.", ext))
  open_device[[ext]](out_path)
  draw_enrichment_heatmap()
  dev.off()
  cat("Saved:", out_path, "\n")
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=== Figure 4B2 Summary ===\n")
cat("PTMs shown:", nrow(combined_matrix), "\n")
cat("Positions analyzed: 10-25\n")
cat("Peptide length filter: 10-25\n")
cat("Minimum count threshold:", min_count, "\n")
cat("Sorted by: canonical PTM order (shared config)\n")
cat("Enrichment = log2(P(position|modified) / P(position|background))\n")
cat("Color scale: capped at +/-", cap_value, "\n")
