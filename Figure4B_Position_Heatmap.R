# =============================================================================
# FIGURE 4B: PTM Position Enrichment Heatmap
# =============================================================================
# Shows log2 enrichment of PTM positions vs background amino acid frequencies
# Similar to mod_position_enrichment.png example

library(readxl)
library(reshape2)
library(ggplot2)
library(pheatmap)

# =============================================================================
# CONFIGURATION
# =============================================================================
xlsx_file <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"
positions <- 1:15  # Peptide positions to analyze (match circos plot)
min_count <- 20    # Minimum occurrences for a PTM+Residue combination

# =============================================================================
# LOAD DATA
# =============================================================================
cat("=== Figure 4B: PTM Position Enrichment Heatmap ===\n\n")

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
# CALCULATE BACKGROUND AA POSITIONAL DISTRIBUTION
# =============================================================================
# For each AA: what fraction of that AA occurs at each position?
# This is the correct comparison: P(position | AA) for both modified and background

aa_list <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

# First, count each AA at each position
bg_counts <- matrix(0, nrow = length(aa_list), ncol = length(positions))
rownames(bg_counts) <- aa_list
colnames(bg_counts) <- positions

for (pos in positions) {
  aa_at_pos <- sapply(bg_peptides, function(pep) {
    if (nchar(pep) >= pos) substr(pep, pos, pos) else NA
  })
  aa_at_pos <- aa_at_pos[!is.na(aa_at_pos)]
  freq_table <- table(aa_at_pos)
  for (aa in names(freq_table)) {
    if (aa %in% aa_list) {
      bg_counts[aa, as.character(pos)] <- freq_table[aa]
    }
  }
}

# Now convert to POSITIONAL DISTRIBUTION for each AA
# bg_freq[aa, pos] = count of AA at position / total count of AA across all positions
bg_freq <- bg_counts
for (aa in aa_list) {
  total_aa <- sum(bg_counts[aa, ])
  if (total_aa > 0) {
    bg_freq[aa, ] <- bg_counts[aa, ] / total_aa
  }
}

# Report some stats
cat("Background AA positional distribution calculated\n")
cat("Example - Cysteine (C) distribution:\n")
cat("  Total C in background:", sum(bg_counts["C", ]), "\n")
cat("  C at pos 1:", sprintf("%.2f%%", 100*bg_freq["C", "1"]), "\n")
cat("  C at pos 9:", sprintf("%.2f%%", 100*bg_freq["C", "9"]), "\n\n")

# =============================================================================
# CALCULATE PTM POSITION FREQUENCIES AND ENRICHMENT
# =============================================================================

# Filter to combinations with sufficient counts
combo_counts <- table(ptm_data$PTM_Residue)
valid_combos <- names(combo_counts[combo_counts >= min_count])
cat("PTM+Residue combinations with >=", min_count, "occurrences:", length(valid_combos), "\n\n")

# Initialize result matrices
enrichment_matrix <- matrix(NA, nrow = length(valid_combos), ncol = length(positions))
rownames(enrichment_matrix) <- valid_combos
colnames(enrichment_matrix) <- positions

correlation_vector <- numeric(length(valid_combos))
names(correlation_vector) <- valid_combos

# Calculate enrichment for each PTM+Residue combination
for (combo in valid_combos) {
  subset_data <- ptm_data[ptm_data$PTM_Residue == combo, ]
  target_residue <- subset_data$Residue[1]

  # Skip non-standard residues (like N-term)
  if (!(target_residue %in% aa_list)) {
    cat("  Skipping", combo, "(non-standard residue)\n")
    next
  }

  # Get position frequencies for this PTM
  pos_freq <- numeric(length(positions))
  names(pos_freq) <- positions

  for (pos in positions) {
    # Count how many modifications occur at this position
    n_at_pos <- sum(subset_data$Site == pos, na.rm = TRUE)
    # Normalize by total for this PTM
    pos_freq[as.character(pos)] <- n_at_pos / nrow(subset_data)
  }

  # Get background frequency for the target residue at each position
  bg_residue_freq <- bg_freq[target_residue, ]

  # Calculate log2 enrichment (add small pseudocount to avoid log(0))
  pseudocount <- 0.001
  enrichment <- log2((pos_freq + pseudocount) / (bg_residue_freq + pseudocount))
  enrichment_matrix[combo, ] <- enrichment

  # Calculate correlation between PTM position distribution and background
  correlation_vector[combo] <- cor(pos_freq, bg_residue_freq, use = "complete.obs")
}

# Replace infinite values with NA
enrichment_matrix[is.infinite(enrichment_matrix)] <- NA

# Remove rows with all NA (skipped combinations)
valid_rows <- apply(enrichment_matrix, 1, function(x) !all(is.na(x)))
enrichment_matrix <- enrichment_matrix[valid_rows, ]
correlation_vector <- correlation_vector[valid_rows]

cat("\nValid combinations after filtering:", nrow(enrichment_matrix), "\n")

# =============================================================================
# PREPARE DATA FOR HEATMAP
# =============================================================================

# Parse PTM and Residue from row names for sorting
row_info <- data.frame(
  combo = rownames(enrichment_matrix),
  stringsAsFactors = FALSE
)
# Extract residue (first character/word) and PTM (rest)
row_info$Residue <- sapply(strsplit(row_info$combo, " "), `[`, 1)
row_info$PTM <- sapply(strsplit(row_info$combo, " "), function(x) paste(x[-1], collapse = " "))

# Sort by PTM first, then by Residue within each PTM
row_info <- row_info[order(row_info$PTM, row_info$Residue), ]
row_order <- match(row_info$combo, rownames(enrichment_matrix))

enrichment_ordered <- enrichment_matrix[row_order, ]
correlation_ordered <- correlation_vector[row_order]

cat("Sorted by PTM type, then by amino acid\n")

# Create combined matrix with correlation column
combined_matrix <- cbind(Corr = correlation_ordered, enrichment_ordered)

# Cap extreme values for better visualization
cap_value <- 3
combined_matrix[, -1] <- pmin(pmax(combined_matrix[, -1], -cap_value), cap_value)

# =============================================================================
# CREATE HEATMAP - PNG
# =============================================================================
cat("Generating heatmap...\n")

# Create annotation for row names (format nicely)
row_labels <- rownames(combined_matrix)
row_labels <- gsub("N-term", "N-term", row_labels)

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
cell_w <- 22
cell_h <- 14
fig_width <- (n_cols * cell_w / 72) + 4  # cells + margins for labels/legend
fig_height <- (n_rows * cell_h / 72) + 2  # cells + margins for title

# Save heatmap
png("figure_panels/Figure4B_Position_Enrichment.png", width = fig_width, height = fig_height, units = "in", res = 300)

# Use pheatmap with correlation annotation
pheatmap(
  combined_matrix[, -1],  # Exclude Corr column for main heatmap
  cluster_rows = FALSE,   # Already ordered by PTM then AA
  cluster_cols = FALSE,
  color = enrichment_colors,
  breaks = seq(-cap_value, cap_value, length.out = 101),
  border_color = "white",
  cellwidth = cell_w,
  cellheight = cell_h,
  fontsize_row = 10,
  fontsize_col = 11,
  main = "PTM Position Enrichment vs Background",
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

dev.off()
cat("Saved: figure_panels/Figure4B_Position_Enrichment.png\n")

# =============================================================================
# CREATE HEATMAP - PDF
# =============================================================================
pdf("figure_panels/Figure4B_Position_Enrichment.pdf", width = fig_width, height = fig_height)

pheatmap(
  combined_matrix[, -1],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = enrichment_colors,
  breaks = seq(-cap_value, cap_value, length.out = 101),
  border_color = "white",
  cellwidth = cell_w,
  cellheight = cell_h,
  fontsize_row = 10,
  fontsize_col = 11,
  main = "PTM Position Enrichment vs Background",
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

dev.off()
cat("Saved: figure_panels/Figure4B_Position_Enrichment.pdf\n")

# =============================================================================
# SAVE DATA
# =============================================================================
write.csv(combined_matrix, "figure_panels/data_figure4B_enrichment.csv")
cat("Saved: figure_panels/data_figure4B_enrichment.csv\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=== Figure 4B Summary ===\n")
cat("PTM+Residue combinations shown:", nrow(combined_matrix), "\n")
cat("Positions analyzed: 1-15\n")
cat("Sorted by: PTM type, then amino acid\n")
cat("Enrichment = log2(P(position|modified) / P(position|AA in background))\n")
cat("Color scale: capped at +/-", cap_value, "\n")
cat("Correlation: Pearson correlation between modified and background positional distributions\n")
