# =============================================================================
# FIGURE 4X (HLA-II): PTM Position Enrichment at P1-P25
# =============================================================================
# HLA-II port of the HLA-I template Figure4X_Position_P1P14.R: log2 enrichment
# of PTM site positions vs background across full peptide positions P1-P25
# (peptide lengths 10-25). Produces two heatmaps:
#   4X1 - PTM+Residue level -> Figure4X1_Position_P1P25_ByResidue
#   4X2 - PTM level         -> Figure4X2_Position_P1P25_ByPTM
#
# Data source: figure_panels/data_ptm_sites.csv, data_background.csv
# Run from HLA-II/:  Rscript scripts/Figure4X_Position_P1P25.R
# =============================================================================

library(dplyr)
library(pheatmap)

source("../shared/ptm_config.R")

# =============================================================================
# CONFIGURATION
# =============================================================================
positions <- 1:25       # Full peptide positions P1-P25 (HLA-II lengths 10-25)
min_count_residue <- 20 # Minimum occurrences for PTM+Residue heatmap
min_count_ptm <- 50     # Minimum occurrences for PTM-level heatmap

output_dir <- "figure_panels"
dir.create(output_dir, showWarnings = FALSE)

# =============================================================================
# LOAD DATA
# =============================================================================
cat("=== Figure 4X (HLA-II): PTM Position Enrichment (P1-P25) ===\n\n")

ptm_data <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)
ptm_data <- ptm_data %>% filter(Length >= 10 & Length <= 25)
ptm_data$PTM_Residue <- paste(ptm_data$Residue, tolower(ptm_data$PTM))
cat("Loaded", nrow(ptm_data), "PTM records (length 10-25)\n")

background <- read.csv("figure_panels/data_background.csv", stringsAsFactors = FALSE)
background <- background %>% filter(Length >= 10 & Length <= 25)
bg_peptides <- background$Peptide
cat("Background peptides (length 10-25):", length(bg_peptides), "\n\n")

# Shared color palette
enrichment_colors <- colorRampPalette(c(
  "#2166AC", "#4393C3", "#92C5DE", "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
))(100)
cap_value <- 3
pseudocount <- 0.001

col_labels <- paste0("P", positions)
ann_cols <- list(Corr = colorRampPalette(c("white", "gray20"))(100))

# Heatmaps bypass save_panel() (ggplot-only), so the tif is written with an
# explicit tiff() device at the same geometry/dpi as the png.
save_heatmap_png_pdf_tif <- function(name, width_in, height_in, draw_heatmap) {
  png(file.path(output_dir, paste0(name, ".png")),
      width = width_in, height = height_in, units = "in", res = 300)
  draw_heatmap()
  dev.off()
  pdf(file.path(output_dir, paste0(name, ".pdf")),
      width = width_in, height = height_in)
  draw_heatmap()
  dev.off()
  tiff(file.path(output_dir, paste0(name, ".tif")),
       width = width_in, height = height_in, units = "in", res = 300,
       compression = "lzw")
  draw_heatmap()
  dev.off()
  cat("Saved:", name, "(png + pdf + tif)\n")
}

draw_position_enrichment_heatmap <- function(values, row_labels,
                                             cell_w, cell_h,
                                             fontsize_row, fontsize_col) {
  row_annotation <- data.frame(Corr = values[, "Corr"])
  rownames(row_annotation) <- rownames(values)
  pheatmap(
    values[, -1],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = enrichment_colors,
    breaks = seq(-cap_value, cap_value, length.out = 101),
    border_color = "white",
    cellwidth = cell_w,
    cellheight = cell_h,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    main = NA,
    labels_row = row_labels,
    labels_col = col_labels,
    legend = TRUE,
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3"),
    annotation_row = row_annotation,
    annotation_colors = ann_cols,
    annotation_legend = TRUE,
    annotation_names_row = TRUE
  )
}

# #############################################################################
# FIGURE 4X1 - PTM+Residue enrichment at P1-P25
# #############################################################################

cat("--- Figure 4X1: PTM+Residue enrichment (P1-P25) ---\n")

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
  enrichment_matrix[combo, ] <- log2((pos_freq + pseudocount) /
                                       (bg_residue_freq + pseudocount))
  corr_vec[combo] <- cor(pos_freq, bg_residue_freq, use = "complete.obs")
}

enrichment_matrix[is.infinite(enrichment_matrix)] <- NA
valid_rows <- apply(enrichment_matrix, 1, function(x) !all(is.na(x)))
enrichment_matrix <- enrichment_matrix[valid_rows, , drop = FALSE]
corr_vec <- corr_vec[valid_rows]

cat("Valid combinations:", nrow(enrichment_matrix), "\n")

# Sort by canonical PTM order (shared with HLA-I), then residue within PTM
row_info <- data.frame(combo = rownames(enrichment_matrix), stringsAsFactors = FALSE)
row_info$Residue <- sapply(strsplit(row_info$combo, " "), `[`, 1)
row_info$PTM <- sapply(strsplit(row_info$combo, " "), function(x) paste(x[-1], collapse = " "))
row_info <- row_info[order(match(row_info$PTM, tolower(PTM_ORDER)), row_info$Residue), ]
row_order <- match(row_info$combo, rownames(enrichment_matrix))

enrichment_ordered <- enrichment_matrix[row_order, , drop = FALSE]
corr_ordered <- corr_vec[row_order]

combined_res <- cbind(Corr = corr_ordered, enrichment_ordered)
combined_res[, -1] <- pmin(pmax(combined_res[, -1], -cap_value), cap_value)

# Cell geometry and fonts identical to the HLA-I template so the I/II pair
# matches visually; canvas sizes differ with matrix extent (see ptm_config.R).
n_cols <- length(positions)
cell_w <- 20
cell_h <- 14
fig_w <- (n_cols * cell_w / 72) + 8
fig_h <- (nrow(combined_res) * cell_h / 72) + 2.5

save_heatmap_png_pdf_tif(
  "Figure4X1_Position_P1P25_ByResidue", fig_w, fig_h,
  function() draw_position_enrichment_heatmap(
    combined_res, rownames(combined_res),
    cell_w = cell_w, cell_h = cell_h, fontsize_row = 10, fontsize_col = 11
  )
)

write.csv(combined_res, file.path(output_dir, "data_figure4X1_p25_enrichment.csv"))

# #############################################################################
# FIGURE 4X2 - PTM-level enrichment at P1-P25 (all residues aggregated)
# #############################################################################

cat("\n--- Figure 4X2: PTM-level enrichment (P1-P25) ---\n")

# ---- Background overall positional distribution ----
bg_pos_counts <- setNames(numeric(length(positions)), positions)
for (pos in positions) {
  bg_pos_counts[as.character(pos)] <- sum(nchar(bg_peptides) >= pos)
}
bg_pos_freq <- bg_pos_counts / sum(bg_pos_counts)
cat("Background overall positional distribution calculated\n")

# ---- Enrichment per PTM ----
# Collaborator instruction: the PTM-level panel shows biological PTMs only, so
# Artifact Oxidation (sample-prep artifact) is excluded here. The ByResidue
# panel above keeps artifact rows, matching the HLA-I original.
ptm_counts <- table(ptm_data$PTM)
valid_ptms <- names(ptm_counts[ptm_counts >= min_count_ptm])
valid_ptms <- valid_ptms[valid_ptms %in% BIO_PTMS]
cat("Bio PTMs with >=", min_count_ptm, ":", length(valid_ptms), "\n")

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

  enrich_ptm[ptm, ] <- log2((pos_freq + pseudocount) / (bg_pos_freq + pseudocount))
  corr_ptm[ptm] <- cor(pos_freq, bg_pos_freq, use = "complete.obs")
}

enrich_ptm[is.infinite(enrich_ptm)] <- NA

# Sort by canonical PTM order (shared with HLA-I)
so <- order(match(rownames(enrich_ptm), PTM_ORDER))
enrich_ptm <- enrich_ptm[so, , drop = FALSE]
corr_ptm <- corr_ptm[so]
count_ptm <- count_ptm[so]

combined_ptm <- cbind(Corr = corr_ptm, enrich_ptm)
combined_ptm[, -1] <- pmin(pmax(combined_ptm[, -1], -cap_value), cap_value)

row_labels_ptm <- paste0(rownames(combined_ptm), " (n=", count_ptm, ")")

cell_w2 <- 28
cell_h2 <- 20
fig_w2 <- (n_cols * cell_w2 / 72) + 5
fig_h2 <- (nrow(combined_ptm) * cell_h2 / 72) + 3

save_heatmap_png_pdf_tif(
  "Figure4X2_Position_P1P25_ByPTM", fig_w2, fig_h2,
  function() draw_position_enrichment_heatmap(
    combined_ptm, row_labels_ptm,
    cell_w = cell_w2, cell_h = cell_h2, fontsize_row = 11, fontsize_col = 12
  )
)

write.csv(combined_ptm, file.path(output_dir, "data_figure4X2_p25_enrichment.csv"))

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=== Figure 4X (HLA-II) Summary ===\n")
cat("4X1 - PTM+Residue combinations:", nrow(combined_res), "\n")
cat("4X2 - Bio PTMs (aggregated):", nrow(combined_ptm), "\n")
cat("Positions analyzed: P1-P25\n")
cat("Peptide length filter: 10-25\n")
cat("Enrichment = log2(P(position|modified) / P(position|background))\n")
cat("Color scale: capped at +/-", cap_value, "\n")
