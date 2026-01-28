# =============================================================================
# FIGURE 4X: PTM Position Enrichment at P1-P14
# =============================================================================
# Shows log2 enrichment of PTM positions vs background across ALL peptide
# positions (P1–P14), rather than only the tail positions (8–14).
# Produces two heatmaps:
#   4X1 – PTM+Residue level (like Figure 4B)
#   4X2 – PTM level, all residues aggregated (like Figure 4B2)
#
# Data source: data_ptm_sites.csv, data_background.csv (from data_loader.R)
# =============================================================================

library(dplyr)
library(pheatmap)

# =============================================================================
# CONFIGURATION
# =============================================================================
positions <- 1:14       # Full peptide positions P1-P14
min_count_residue <- 20 # Minimum occurrences for PTM+Residue heatmap
min_count_ptm <- 50     # Minimum occurrences for PTM-level heatmap

output_dir <- "figure_panels"
dir.create(output_dir, showWarnings = FALSE)

# =============================================================================
# LOAD DATA
# =============================================================================
cat("=== Figure 4X: PTM Position Enrichment (P1-P14) ===\n\n")

ptm_data <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)
ptm_data <- ptm_data %>% filter(Length >= 8 & Length <= 14)
ptm_data$PTM_Residue <- paste(ptm_data$Residue, tolower(ptm_data$PTM))
cat("Loaded", nrow(ptm_data), "PTM records (length 8-14)\n")

background <- read.csv("figure_panels/data_background.csv", stringsAsFactors = FALSE)
background <- background %>% filter(Length >= 8 & Length <= 14)
bg_peptides <- background$Peptide
cat("Background peptides (length 8-14):", length(bg_peptides), "\n\n")

# Shared color palette
enrichment_colors <- colorRampPalette(c(
  "#2166AC", "#4393C3", "#92C5DE", "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
))(100)
cap_value <- 3

# Column labels
col_labels <- paste0("P", positions)

# #############################################################################
# FIGURE 4X1 – PTM+Residue enrichment at P1-P14
# #############################################################################

cat("--- Figure 4X1: PTM+Residue enrichment (P1-P14) ---\n")

# ---- Background per-AA positional distribution ----
aa_list <- c("A","C","D","E","F","G","H","I","K","L",
             "M","N","P","Q","R","S","T","V","W","Y")

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

# Normalize: P(position | AA)
bg_freq <- bg_counts
for (aa in aa_list) {
  total_aa <- sum(bg_counts[aa, ])
  if (total_aa > 0) bg_freq[aa, ] <- bg_counts[aa, ] / total_aa
}
cat("Background per-AA positional distribution calculated\n")

# ---- Enrichment per PTM+Residue ----
combo_counts <- table(ptm_data$PTM_Residue)
valid_combos <- names(combo_counts[combo_counts >= min_count_residue])
cat("PTM+Residue combos with >=", min_count_residue, ":", length(valid_combos), "\n")

enrichment_matrix <- matrix(NA, nrow = length(valid_combos), ncol = length(positions))
rownames(enrichment_matrix) <- valid_combos
colnames(enrichment_matrix) <- positions

corr_vec <- setNames(numeric(length(valid_combos)), valid_combos)

for (combo in valid_combos) {
  subset_data <- ptm_data[ptm_data$PTM_Residue == combo, ]
  target_residue <- subset_data$Residue[1]

  if (!(target_residue %in% aa_list)) {
    cat("  Skipping", combo, "(non-standard residue)\n")
    next
  }

  pos_freq <- setNames(numeric(length(positions)), positions)
  for (pos in positions) {
    n_at_pos <- sum(subset_data$Site == pos, na.rm = TRUE)
    pos_freq[as.character(pos)] <- n_at_pos / nrow(subset_data)
  }

  bg_residue_freq <- bg_freq[target_residue, ]
  pseudocount <- 0.001
  enrichment <- log2((pos_freq + pseudocount) / (bg_residue_freq + pseudocount))
  enrichment_matrix[combo, ] <- enrichment
  corr_vec[combo] <- cor(pos_freq, bg_residue_freq, use = "complete.obs")
}

enrichment_matrix[is.infinite(enrichment_matrix)] <- NA
valid_rows <- apply(enrichment_matrix, 1, function(x) !all(is.na(x)))
enrichment_matrix <- enrichment_matrix[valid_rows, , drop = FALSE]
corr_vec <- corr_vec[valid_rows]

cat("Valid combinations:", nrow(enrichment_matrix), "\n")

# Sort by PTM then residue
row_info <- data.frame(combo = rownames(enrichment_matrix), stringsAsFactors = FALSE)
row_info$Residue <- sapply(strsplit(row_info$combo, " "), `[`, 1)
row_info$PTM <- sapply(strsplit(row_info$combo, " "), function(x) paste(x[-1], collapse = " "))
row_info <- row_info[order(row_info$PTM, row_info$Residue), ]
row_order <- match(row_info$combo, rownames(enrichment_matrix))

enrichment_ordered <- enrichment_matrix[row_order, , drop = FALSE]
corr_ordered <- corr_vec[row_order]

combined_res <- cbind(Corr = corr_ordered, enrichment_ordered)
combined_res[, -1] <- pmin(pmax(combined_res[, -1], -cap_value), cap_value)

# Dimensions
n_rows <- nrow(combined_res)
n_cols <- length(positions)
cell_w <- 20
cell_h <- 14
fig_w <- (n_cols * cell_w / 72) + 8
fig_h <- (n_rows * cell_h / 72) + 2.5

row_ann <- data.frame(Corr = combined_res[, "Corr"])
rownames(row_ann) <- rownames(combined_res)
ann_cols <- list(Corr = colorRampPalette(c("white", "gray20"))(100))

# ---- Save PNG ----
png(file.path(output_dir, "Figure4X1_Position_P1P14_ByResidue.png"),
    width = fig_w, height = fig_h, units = "in", res = 300)

pheatmap(
  combined_res[, -1],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = enrichment_colors,
  breaks = seq(-cap_value, cap_value, length.out = 101),
  border_color = "white",
  cellwidth = cell_w,
  cellheight = cell_h,
  fontsize_row = 10,
  fontsize_col = 11,
  main = "PTM Position Enrichment (P1\u2013P14)\nlog2(modified / background), PTM+Residue combinations with \u226520 occurrences",
  labels_row = rownames(combined_res),
  labels_col = col_labels,
  legend = TRUE,
  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
  legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3"),
  annotation_row = row_ann,
  annotation_colors = ann_cols,
  annotation_legend = TRUE,
  annotation_names_row = TRUE
)
dev.off()
cat("Saved: Figure4X1_Position_P1P14_ByResidue.png\n")

# ---- Save PDF ----
pdf(file.path(output_dir, "Figure4X1_Position_P1P14_ByResidue.pdf"),
    width = fig_w, height = fig_h)
pheatmap(
  combined_res[, -1],
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = enrichment_colors,
  breaks = seq(-cap_value, cap_value, length.out = 101),
  border_color = "white",
  cellwidth = cell_w, cellheight = cell_h,
  fontsize_row = 10, fontsize_col = 11,
  main = "PTM Position Enrichment (P1\u2013P14)\nlog2(modified / background), PTM+Residue combinations with \u226520 occurrences",
  labels_row = rownames(combined_res),
  labels_col = col_labels,
  legend = TRUE,
  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
  legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3"),
  annotation_row = row_ann, annotation_colors = ann_cols,
  annotation_legend = TRUE, annotation_names_row = TRUE
)
dev.off()
cat("Saved: Figure4X1_Position_P1P14_ByResidue.pdf\n")

# Save enrichment data
write.csv(combined_res, file.path(output_dir, "data_figure4X1_enrichment.csv"))

# #############################################################################
# FIGURE 4X2 – PTM-level enrichment at P1-P14 (all residues aggregated)
# #############################################################################

cat("\n--- Figure 4X2: PTM-level enrichment (P1-P14) ---\n")

# ---- Background overall positional distribution ----
bg_pos_counts <- setNames(numeric(length(positions)), positions)
for (pos in positions) {
  bg_pos_counts[as.character(pos)] <- sum(nchar(bg_peptides) >= pos)
}
bg_pos_freq <- bg_pos_counts / sum(bg_pos_counts)
cat("Background overall positional distribution calculated\n")

# ---- Enrichment per PTM ----
ptm_counts <- table(ptm_data$PTM)
valid_ptms <- names(ptm_counts[ptm_counts >= min_count_ptm])
cat("PTMs with >=", min_count_ptm, ":", length(valid_ptms), "\n")

enrich_ptm <- matrix(NA, nrow = length(valid_ptms), ncol = length(positions))
rownames(enrich_ptm) <- valid_ptms
colnames(enrich_ptm) <- positions

corr_ptm <- setNames(numeric(length(valid_ptms)), valid_ptms)
count_ptm <- setNames(numeric(length(valid_ptms)), valid_ptms)

for (ptm in valid_ptms) {
  subset_data <- ptm_data[ptm_data$PTM == ptm, ]
  count_ptm[ptm] <- nrow(subset_data)

  pos_freq <- setNames(numeric(length(positions)), positions)
  for (pos in positions) {
    n_at_pos <- sum(subset_data$Site == pos, na.rm = TRUE)
    pos_freq[as.character(pos)] <- n_at_pos / nrow(subset_data)
  }

  pseudocount <- 0.001
  enrich_ptm[ptm, ] <- log2((pos_freq + pseudocount) / (bg_pos_freq + pseudocount))
  corr_ptm[ptm] <- cor(pos_freq, bg_pos_freq, use = "complete.obs")
}

enrich_ptm[is.infinite(enrich_ptm)] <- NA

# Sort by count descending
so <- order(count_ptm[rownames(enrich_ptm)], decreasing = TRUE)
enrich_ptm <- enrich_ptm[so, , drop = FALSE]
corr_ptm <- corr_ptm[so]
count_ptm <- count_ptm[so]

combined_ptm <- cbind(Corr = corr_ptm, enrich_ptm)
combined_ptm[, -1] <- pmin(pmax(combined_ptm[, -1], -cap_value), cap_value)

row_labels_ptm <- paste0(rownames(combined_ptm), " (n=", count_ptm, ")")

n_rows2 <- nrow(combined_ptm)
cell_w2 <- 28
cell_h2 <- 20
fig_w2 <- (n_cols * cell_w2 / 72) + 5
fig_h2 <- (n_rows2 * cell_h2 / 72) + 3

row_ann2 <- data.frame(Corr = combined_ptm[, "Corr"])
rownames(row_ann2) <- rownames(combined_ptm)

# ---- Save PNG ----
png(file.path(output_dir, "Figure4X2_Position_P1P14_ByPTM.png"),
    width = fig_w2, height = fig_h2, units = "in", res = 300)

pheatmap(
  combined_ptm[, -1],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = enrichment_colors,
  breaks = seq(-cap_value, cap_value, length.out = 101),
  border_color = "white",
  cellwidth = cell_w2,
  cellheight = cell_h2,
  fontsize_row = 11,
  fontsize_col = 12,
  main = "PTM Position Enrichment (P1\u2013P14)\nlog2(modified / background), PTMs \u226550 occurrences, all residues aggregated",
  labels_row = row_labels_ptm,
  labels_col = col_labels,
  legend = TRUE,
  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
  legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3"),
  annotation_row = row_ann2,
  annotation_colors = ann_cols,
  annotation_legend = TRUE,
  annotation_names_row = TRUE
)
dev.off()
cat("Saved: Figure4X2_Position_P1P14_ByPTM.png\n")

# ---- Save PDF ----
pdf(file.path(output_dir, "Figure4X2_Position_P1P14_ByPTM.pdf"),
    width = fig_w2, height = fig_h2)
pheatmap(
  combined_ptm[, -1],
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = enrichment_colors,
  breaks = seq(-cap_value, cap_value, length.out = 101),
  border_color = "white",
  cellwidth = cell_w2, cellheight = cell_h2,
  fontsize_row = 11, fontsize_col = 12,
  main = "PTM Position Enrichment (P1\u2013P14)\nlog2(modified / background), PTMs \u226550 occurrences, all residues aggregated",
  labels_row = row_labels_ptm,
  labels_col = col_labels,
  legend = TRUE,
  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
  legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3"),
  annotation_row = row_ann2, annotation_colors = ann_cols,
  annotation_legend = TRUE, annotation_names_row = TRUE
)
dev.off()
cat("Saved: Figure4X2_Position_P1P14_ByPTM.pdf\n")

# Save enrichment data
write.csv(combined_ptm, file.path(output_dir, "data_figure4X2_enrichment.csv"))

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=== Figure 4X Summary ===\n")
cat("4X1 - PTM+Residue combinations:", nrow(combined_res), "\n")
cat("4X2 - PTMs (aggregated):", nrow(combined_ptm), "\n")
cat("Positions analyzed: P1-P14\n")
cat("Peptide length filter: 8-14\n")
cat("Enrichment = log2(P(position|modified) / P(position|background))\n")
cat("Color scale: capped at +/-", cap_value, "\n")
