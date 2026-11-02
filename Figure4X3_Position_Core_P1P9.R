# =============================================================================
# FIGURE 4X3/4X4 (HLA-I): PTM Position Enrichment within the 9-mer binding core
# =============================================================================
# Core-based counterpart of Figure4X_Position_P1P14.R: positions are NetMHCpan
# core positions P1-P9. Modified counts come from CoreSite (NA = modification
# outside the predicted core, dropped); background counts come from the residue
# composition of the background peptides' predicted 9-mer cores.
#   4X3 - PTM+Residue level -> Figure4X3_Position_CoreP1P9_ByResidue
#   4X4 - PTM level         -> Figure4X4_Position_CoreP1P9_ByPTM
#
# Mirror of HLA-II/scripts/Figure4X3_Position_Core_P1P9.R: statistics, color
# ramp, cell geometry, fonts and thresholds are identical so the HLA-I and
# HLA-II panels can be mounted side by side.
#
# Data source: figure_panels/data_ptm_sites_core.csv, data_hla1_cores.csv,
# data_background.csv (cores built by hla1_netmhcpan_cores.py)
# Run from the repo root:  Rscript Figure4X3_Position_Core_P1P9.R
# =============================================================================

library(dplyr)
library(pheatmap)

source("shared/ptm_config.R")

# =============================================================================
# CONFIGURATION
# =============================================================================
core_positions <- 1:9   # 9-mer binding-core positions P1-P9
min_count_residue <- 20 # Minimum in-core occurrences for PTM+Residue heatmap
min_count_ptm <- 50     # Minimum in-core occurrences for PTM-level heatmap

# NetMHCpan aligns 8-mers to the core by inserting a wildcard '-', so their
# core assignment is not a real residue mapping; the collaborator asked for
# them to be left out. hla1_netmhcpan_cores.py already drops them, and the
# background is filtered to the same 9-14 range here.
shortest_peptide <- 9
longest_peptide <- 14

output_dir <- "figure_panels"
dir.create(output_dir, showWarnings = FALSE)

# =============================================================================
# LOAD DATA
# =============================================================================
cat("=== Figure 4X3/4X4 (HLA-I): PTM Core Position Enrichment (P1-P9) ===\n\n")

ptm_data <- read.csv("figure_panels/data_ptm_sites_core.csv", stringsAsFactors = FALSE)
ptm_data <- ptm_data %>% filter(Length >= shortest_peptide & Length <= longest_peptide)
ptm_data$PTM_Residue <- paste(ptm_data$Residue, tolower(ptm_data$PTM))
cat("Loaded", nrow(ptm_data), "PTM records (length 9-14)\n")

core_table <- read.csv("figure_panels/data_hla1_cores.csv", stringsAsFactors = FALSE)
background <- read.csv("figure_panels/data_background.csv", stringsAsFactors = FALSE)
background <- background %>%
  filter(Length >= shortest_peptide & Length <= longest_peptide) %>%
  left_join(core_table[, c("Peptide", "Core")], by = "Peptide")
bg_cores <- background$Core[!is.na(background$Core) & nchar(background$Core) == 9]
cat("Background cores (9-mer):", length(bg_cores), "of", nrow(background), "peptides\n")

# ---- Per-PTM in-core coverage ----
cat("\nPer-PTM in-core coverage (CoreSite 1-9 vs total sites):\n")
coverage <- ptm_data %>%
  group_by(PTM) %>%
  summarise(in_core = sum(!is.na(CoreSite)), total = n(), .groups = "drop") %>%
  mutate(fraction = in_core / total) %>%
  arrange(match(PTM, PTM_ORDER))
for (i in seq_len(nrow(coverage))) {
  cat(sprintf("  %-22s %5d / %5d  (%.1f%%)\n", coverage$PTM[i],
              coverage$in_core[i], coverage$total[i], 100 * coverage$fraction[i]))
}
cat("\n")

in_core_data <- ptm_data[!is.na(ptm_data$CoreSite), ]

# Shared color palette
enrichment_colors <- colorRampPalette(c(
  "#2166AC", "#4393C3", "#92C5DE", "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
))(100)
cap_value <- 3
pseudocount <- 0.001

col_labels <- paste0("P", core_positions)
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

draw_core_enrichment_heatmap <- function(values, row_labels,
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

# =============================================================================
# BACKGROUND PER-AA POSITIONAL DISTRIBUTION ACROSS CORES
# =============================================================================
aa_list <- c("A","C","D","E","F","G","H","I","K","L",
             "M","N","P","Q","R","S","T","V","W","Y")

bg_counts <- matrix(0, nrow = length(aa_list), ncol = length(core_positions))
rownames(bg_counts) <- aa_list
colnames(bg_counts) <- core_positions

for (pos in core_positions) {
  freq_table <- table(substr(bg_cores, pos, pos))
  for (aa in names(freq_table)) {
    if (aa %in% aa_list) {
      bg_counts[aa, as.character(pos)] <- freq_table[aa]
    }
  }
}

# Normalize: P(core position | AA)
bg_freq <- bg_counts
for (aa in aa_list) {
  total_aa <- sum(bg_counts[aa, ])
  if (total_aa > 0) bg_freq[aa, ] <- bg_counts[aa, ] / total_aa
}
cat("Background per-AA core positional distribution calculated\n")

# #############################################################################
# FIGURE 4X3 - PTM+Residue enrichment at core P1-P9
# #############################################################################

cat("\n--- Figure 4X3: PTM+Residue core enrichment (P1-P9) ---\n")

combo_counts <- table(in_core_data$PTM_Residue)
valid_combos <- names(combo_counts[combo_counts >= min_count_residue])
cat("PTM+Residue combos with >=", min_count_residue, "in-core:", length(valid_combos), "\n")

enrichment_matrix <- matrix(NA, nrow = length(valid_combos), ncol = length(core_positions))
rownames(enrichment_matrix) <- valid_combos
colnames(enrichment_matrix) <- core_positions

corr_vec <- setNames(numeric(length(valid_combos)), valid_combos)

for (combo in valid_combos) {
  subset_data <- in_core_data[in_core_data$PTM_Residue == combo, ]
  target_residue <- subset_data$Residue[1]

  if (!(target_residue %in% aa_list)) {
    cat("  Skipping", combo, "(non-standard residue)\n")
    next
  }

  pos_freq <- setNames(numeric(length(core_positions)), core_positions)
  for (pos in core_positions) {
    n_at_pos <- sum(subset_data$CoreSite == pos, na.rm = TRUE)
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

# Sort by canonical PTM order (shared with HLA-II), then residue within PTM
row_info <- data.frame(combo = rownames(enrichment_matrix), stringsAsFactors = FALSE)
row_info$Residue <- sapply(strsplit(row_info$combo, " "), `[`, 1)
row_info$PTM <- sapply(strsplit(row_info$combo, " "), function(x) paste(x[-1], collapse = " "))
row_info <- row_info[order(match(row_info$PTM, tolower(PTM_ORDER)), row_info$Residue), ]
row_order <- match(row_info$combo, rownames(enrichment_matrix))

enrichment_ordered <- enrichment_matrix[row_order, , drop = FALSE]
corr_ordered <- corr_vec[row_order]

combined_res <- cbind(Corr = corr_ordered, enrichment_ordered)
combined_res[, -1] <- pmin(pmax(combined_res[, -1], -cap_value), cap_value)

# Cell geometry and fonts identical to the HLA-II core panel so the two sit
# side by side; canvas sizes differ with matrix extent (see ptm_config.R).
n_cols <- length(core_positions)
cell_w <- 20
cell_h <- 14
fig_w <- (n_cols * cell_w / 72) + 8
fig_h <- (nrow(combined_res) * cell_h / 72) + 2.5

save_heatmap_png_pdf_tif(
  "Figure4X3_Position_CoreP1P9_ByResidue", fig_w, fig_h,
  function() draw_core_enrichment_heatmap(
    combined_res, rownames(combined_res),
    cell_w = cell_w, cell_h = cell_h, fontsize_row = 10, fontsize_col = 11
  )
)

write.csv(combined_res, file.path(output_dir, "data_figure4X3_core_enrichment.csv"))

# #############################################################################
# FIGURE 4X4 - PTM-level enrichment at core P1-P9 (all residues aggregated)
# #############################################################################

cat("\n--- Figure 4X4: PTM-level core enrichment (P1-P9) ---\n")

# Every 9-mer core covers all of P1-P9, so the full-length template's
# length-availability background is uniform and would leave the Corr strip
# undefined (zero variance). The PTM-level null is therefore the
# residue-composition-weighted average of the per-residue background
# profiles: the direct aggregation of
# log2(P(pos | modified, R) / P(pos | background cores, R)) over the residues
# the PTM actually lands on.
ptm_level_background_profile <- function(subset_data) {
  residue_weights <- table(subset_data$Residue[subset_data$Residue %in% aa_list])
  profile <- setNames(numeric(length(core_positions)), core_positions)
  for (aa in names(residue_weights)) {
    profile <- profile + residue_weights[[aa]] * bg_freq[aa, ]
  }
  profile / sum(residue_weights)
}

# Collaborator instruction: the PTM-level panel shows biological PTMs only, so
# Artifact Oxidation (sample-prep artifact) is excluded here. The ByResidue
# panel above keeps artifact rows.
ptm_counts <- table(in_core_data$PTM)
valid_ptms <- names(ptm_counts[ptm_counts >= min_count_ptm])
valid_ptms <- valid_ptms[valid_ptms %in% BIO_PTMS]
cat("Bio PTMs with >=", min_count_ptm, "in-core:", length(valid_ptms), "\n")

enrich_ptm <- matrix(NA, nrow = length(valid_ptms), ncol = length(core_positions))
rownames(enrich_ptm) <- valid_ptms
colnames(enrich_ptm) <- core_positions

corr_ptm <- setNames(numeric(length(valid_ptms)), valid_ptms)
count_ptm <- setNames(numeric(length(valid_ptms)), valid_ptms)

for (ptm in valid_ptms) {
  subset_data <- in_core_data[in_core_data$PTM == ptm, ]
  count_ptm[ptm] <- nrow(subset_data)

  pos_freq <- setNames(numeric(length(core_positions)), core_positions)
  for (pos in core_positions) {
    n_at_pos <- sum(subset_data$CoreSite == pos, na.rm = TRUE)
    pos_freq[as.character(pos)] <- n_at_pos / nrow(subset_data)
  }

  bg_profile <- ptm_level_background_profile(subset_data)
  enrich_ptm[ptm, ] <- log2((pos_freq + pseudocount) / (bg_profile + pseudocount))
  corr_ptm[ptm] <- cor(pos_freq, bg_profile, use = "complete.obs")
}

enrich_ptm[is.infinite(enrich_ptm)] <- NA

# Sort by canonical PTM order (shared with HLA-II)
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
  "Figure4X4_Position_CoreP1P9_ByPTM", fig_w2, fig_h2,
  function() draw_core_enrichment_heatmap(
    combined_ptm, row_labels_ptm,
    cell_w = cell_w2, cell_h = cell_h2, fontsize_row = 11, fontsize_col = 12
  )
)

write.csv(combined_ptm, file.path(output_dir, "data_figure4X4_core_enrichment.csv"))

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=== Figure 4X3/4X4 (HLA-I core) Summary ===\n")
cat("4X3 - PTM+Residue combinations:", nrow(combined_res), "\n")
cat("4X4 - Bio PTMs (aggregated):", nrow(combined_ptm), "\n")
cat("Positions analyzed: core P1-P9 (modifications outside the core dropped)\n")
cat("Peptide length filter: 9-14 (8-mers excluded: wildcard-padded cores)\n")
cat("Enrichment = log2(P(core position|modified) / P(core position|background cores))\n")
cat("Color scale: capped at +/-", cap_value, "\n")
