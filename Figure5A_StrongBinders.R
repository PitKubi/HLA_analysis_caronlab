# =============================================================================
# FIGURE 5A: % Strong & Weak Binders by PTM(+Residue) and Allele
# =============================================================================
# Heatmaps of % strong / % weak binders per allele, by PTM and by PTM+Residue.
#
# Data source: the loader CSVs (data_binding.csv + data_ptm_sites.csv), NOT the
# xlsx. The earlier direct-xlsx read hardcoded per-sheet column positions for
# the 01062026 layout; the 01262026_3 layout moved those blocks, and grep-based
# column pairing has already silently misaligned EL ranks once upstream. The
# loader CSVs carry correctly paired Allele/EL_Rank/Binder per PTM.
#
# data_binding.csv has no Residue column, so the per-Residue panels take the
# residue of the FIRST PTM site per (PTM, Peptide) from data_ptm_sites.csv.
#
# Artifact Oxidation and Carbamidomethylation stay in these binder heatmaps:
# the collaborator's manuscript Fig 2A/B includes them.
# =============================================================================

library(dplyr)
library(tidyr)
library(pheatmap)

output_dir <- "figure_panels"

cat("=== Figure 5A: % Strong/Weak Binders by PTM(+Residue) and Allele ===\n\n")

binding <- read.csv(file.path(output_dir, "data_binding.csv"), stringsAsFactors = FALSE)
sites <- read.csv(file.path(output_dir, "data_ptm_sites.csv"), stringsAsFactors = FALSE)
stopifnot(all(c("PTM", "Peptide", "Length", "Allele", "Binder") %in% colnames(binding)),
          all(c("PTM", "Peptide", "Residue") %in% colnames(sites)))

# HLA-I manuscript length window (the loader CSV is not pre-filtered)
binding <- binding %>% filter(Length >= 8 & Length <= 14)
cat("Binding records (length 8-14):", nrow(binding), "\n")

first_site_residue <- sites %>%
  group_by(PTM, Peptide) %>%
  slice(1) %>%
  ungroup() %>%
  select(PTM, Peptide, Residue)

merged <- binding %>%
  left_join(first_site_residue, by = c("PTM", "Peptide")) %>%
  filter(!is.na(Binder))

n_residue_join_failures <- sum(is.na(merged$Residue) | merged$Residue == "")
cat("Peptides with no residue joined from data_ptm_sites.csv:",
    n_residue_join_failures, "of", nrow(merged),
    "(kept for the by-PTM panels, dropped from the by-Residue panels)\n")

config <- readRDS(file.path(output_dir, "config.rds"))
ptm_colors <- config$ptm_colors
ptm_order <- config$ptm_order
stopifnot(all(merged$PTM %in% ptm_order))

# SB/WB counts per PTM so they can be eyeballed against the sheet summaries
ptm_binder_counts <- merged %>%
  group_by(PTM) %>%
  summarise(
    n = n(),
    SB = sum(Binder == "Strong"),
    WB = sum(Binder == "Weak"),
    NB = sum(Binder == "Non-binder"),
    n_without_allele = sum(is.na(Allele)),
    .groups = "drop"
  ) %>%
  arrange(match(PTM, ptm_order))
cat("\nSB/WB/NB per PTM (length 8-14, before per-allele split):\n")
print(as.data.frame(ptm_binder_counts))

# The heatmaps are per allele; rows without a Best_Allele in the loader CSV
# (Methylation, Citrullination, Dimethylation blocks carry Binder but no
# allele) cannot be placed in an allele column.
merged <- merged %>% filter(!is.na(Allele))
cat("\nRecords with an allele (used for the heatmaps):", nrow(merged), "\n")

# =============================================================================
# CALCULATE % STRONG AND WEAK BINDERS
# =============================================================================

stats_by_residue <- merged %>%
  filter(!is.na(Residue), Residue != "") %>%
  mutate(PTM_Residue = paste(Residue, tolower(PTM))) %>%
  group_by(PTM, Residue, PTM_Residue, Allele) %>%
  summarise(
    n = n(),
    n_strong = sum(Binder == "Strong"),
    n_weak = sum(Binder == "Weak"),
    n_nonbinder = sum(Binder == "Non-binder"),
    pct_strong = round(100 * n_strong / n, 0),
    pct_weak = round(100 * n_weak / n, 0),
    .groups = "drop"
  ) %>%
  filter(n >= 10)

cat("\nPTM+Residue combinations with n >= 10:",
    length(unique(stats_by_residue$PTM_Residue)), "\n")

stats_by_ptm <- merged %>%
  group_by(PTM, Allele) %>%
  summarise(
    n = n(),
    n_strong = sum(Binder == "Strong"),
    n_weak = sum(Binder == "Weak"),
    pct_strong = round(100 * n_strong / n, 0),
    pct_weak = round(100 * n_weak / n, 0),
    .groups = "drop"
  ) %>%
  filter(n >= 10) %>%
  mutate(Residue = "All", PTM_Residue = PTM)

# =============================================================================
# HEATMAP: rows in canonical PTM order (then residue), png + pdf + tif
# =============================================================================

create_binder_heatmap <- function(stats_df, value_col, color_palette, filename) {
  heatmap_df <- stats_df %>%
    select(PTM, Residue, PTM_Residue, Allele, !!sym(value_col)) %>%
    pivot_wider(names_from = Allele, values_from = !!sym(value_col)) %>%
    as.data.frame() %>%
    arrange(match(PTM, ptm_order), Residue)

  rownames(heatmap_df) <- heatmap_df$PTM_Residue

  allele_cols <- intersect(c("A0201", "B0702", "C0702"), colnames(heatmap_df))
  heatmap_matrix <- heatmap_df[, allele_cols, drop = FALSE]

  # NA = that allele never reached n >= 10 for this row; show as 0%
  heatmap_matrix[is.na(heatmap_matrix)] <- 0

  row_annotation <- data.frame(PTM = heatmap_df$PTM)
  rownames(row_annotation) <- rownames(heatmap_matrix)
  ann_colors <- list(PTM = ptm_colors[unique(heatmap_df$PTM)])

  ptm_rle <- rle(heatmap_df$PTM)
  gap_positions <- cumsum(ptm_rle$lengths)[-length(ptm_rle$lengths)]

  n_rows <- nrow(heatmap_matrix)
  fig_height <- max(8, n_rows * 0.25 + 2)

  draw_heatmap <- function() {
    pheatmap(
      as.matrix(heatmap_matrix),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = color_palette,
      breaks = seq(0, 100, length.out = 101),
      border_color = "white",
      cellwidth = 40,
      cellheight = 14,
      fontsize_row = 8,
      fontsize_col = 10,
      main = NA,  # Manuscript panels carry no titles
      display_numbers = TRUE,
      number_format = "%.0f",
      number_color = "black",
      fontsize_number = 7,
      na_col = "gray90",
      legend = TRUE,
      legend_breaks = c(0, 25, 50, 75, 100),
      legend_labels = c("0%", "25%", "50%", "75%", "100%"),
      annotation_row = row_annotation,
      annotation_colors = ann_colors,
      angle_col = 0,
      gaps_row = gap_positions
    )
  }

  open_device <- list(
    png = function(path) png(path, width = 8, height = fig_height,
                             units = "in", res = 300),
    pdf = function(path) pdf(path, width = 8, height = fig_height),
    tif = function(path) tiff(path, width = 8, height = fig_height,
                              units = "in", res = 300, compression = "lzw")
  )
  for (ext in names(open_device)) {
    open_device[[ext]](file.path(output_dir, paste0(filename, ".", ext)))
    draw_heatmap()
    dev.off()
  }
  cat("Saved:", filename, ".png/.pdf/.tif (", n_rows, "rows )\n")

  heatmap_df
}

strong_colors <- colorRampPalette(c("#FFFFFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C"))(100)
weak_colors <- colorRampPalette(c("#FFFFFF", "#FEE8C8", "#FDBB84", "#E34A33", "#B30000"))(100)

strong_df <- create_binder_heatmap(stats_by_residue, "pct_strong", strong_colors,
                                   "Figure5A_StrongBinders_ByResidue")
weak_df <- create_binder_heatmap(stats_by_residue, "pct_weak", weak_colors,
                                 "Figure5A_WeakBinders_ByResidue")
strong_ptm_df <- create_binder_heatmap(stats_by_ptm, "pct_strong", strong_colors,
                                       "Figure5A_StrongBinders_ByPTM")
weak_ptm_df <- create_binder_heatmap(stats_by_ptm, "pct_weak", weak_colors,
                                     "Figure5A_WeakBinders_ByPTM")

# =============================================================================
# SAVE DATA
# =============================================================================

write.csv(stats_by_residue, file.path(output_dir, "data_figure5A_binders.csv"), row.names = FALSE)
cat("\nSaved: data_figure5A_binders.csv\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== Summary ===\n")
cat("PTM+Residue rows shown:", nrow(strong_df), "\n")
cat("PTM-only rows shown:", nrow(strong_ptm_df), "\n")
cat("Residue join failures:", n_residue_join_failures, "\n")
cat("Sorted by: canonical PTM order (shared config), then amino acid residue\n")

cat("\nTop 10 by % strong binders (mean across alleles):\n")
stats_summary <- stats_by_residue %>%
  group_by(PTM, Residue, PTM_Residue) %>%
  summarise(
    n_total = sum(n),
    pct_strong_mean = round(mean(pct_strong, na.rm = TRUE), 1),
    pct_weak_mean = round(mean(pct_weak, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_strong_mean))

print(head(stats_summary, 10))
