# =============================================================================
# FIGURE 5A: % Strong & Weak Binders by PTM+Residue and Allele
# =============================================================================
# Heatmap showing % strong/weak binders, deconvoluted by amino acid residue
# Sorted by PTM type, then by residue within each PTM
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(pheatmap)

xlsx_file <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"
output_dir <- "figure_panels"

cat("=== Figure 5A: % Strong/Weak Binders by PTM+Residue ===\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

# Load PTM data with residue info
ptm_data <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)
ptm_data$PTM_Residue <- paste(ptm_data$Residue, tolower(ptm_data$PTM))

# Filter to length 8-14
ptm_data <- ptm_data %>% filter(Length >= 8 & Length <= 14)

cat("PTM data loaded:", nrow(ptm_data), "records\n")

# Function to extract binding data from sheets
get_binding_for_sheet <- function(sheet_name) {
  suppressMessages({
    df <- read_excel(xlsx_file, sheet = sheet_name)
  })

  pep_cols <- colnames(df)[grepl("^Peptide", colnames(df))]
  binder_col <- colnames(df)[grepl("^Binder", colnames(df))][1]
  allele_col <- colnames(df)[grepl("^Best_Allele", colnames(df))][1]

  if (length(pep_cols) == 0 || is.na(binder_col)) return(NULL)

  for (pep_col in pep_cols) {
    result <- data.frame(
      Peptide = df[[pep_col]],
      Binder = df[[binder_col]],
      Allele = df[[allele_col]],
      stringsAsFactors = FALSE
    )
    result <- result[!is.na(result$Peptide) & !is.na(result$Binder), ]
    if (nrow(result) > 0) return(result)
  }
  return(NULL)
}

# Get binding data from all sheets
sheets <- c("Phospho.", "Acetyl.", "Cyst.", "Methyl.", "Dimethyl.",
            "Deamid. (NQ)", "bioOxid.", "Citrullination", "G-Ubiq.", "N-glyco")

cat("Loading binding data...\n")
all_binding <- list()
for (sheet in sheets) {
  binding <- tryCatch(get_binding_for_sheet(sheet), error = function(e) NULL)
  if (!is.null(binding)) {
    all_binding[[sheet]] <- binding
    cat("  ", sheet, ":", nrow(binding), "\n")
  }
}

binding_combined <- do.call(rbind, all_binding)

# Merge with PTM data
merged <- ptm_data %>%
  inner_join(binding_combined, by = "Peptide", relationship = "many-to-many")

# Filter out NA alleles
merged <- merged %>% filter(!is.na(Allele))

cat("\nTotal merged records:", nrow(merged), "\n")

# =============================================================================
# CALCULATE % STRONG AND WEAK BINDERS
# =============================================================================

# Calculate stats for each PTM+Residue x Allele
stats <- merged %>%
  group_by(PTM, Residue, PTM_Residue, Allele) %>%
  summarise(
    n = n(),
    n_strong = sum(Binder == "Strong"),
    n_weak = sum(Binder == "Weak"),
    n_nonbinder = sum(Binder == "Non-binder"),
    pct_strong = round(100 * n_strong / n, 0),
    pct_weak = round(100 * n_weak / n, 0),
    .groups = "drop"
  )

# Filter to combinations with at least 10 peptides
stats_filtered <- stats %>%
  filter(n >= 10)

cat("\nPTM+Residue combinations with n >= 10:",
    length(unique(stats_filtered$PTM_Residue)), "\n")

# =============================================================================
# PTM colors (circos palette)
# =============================================================================

ptm_colors <- c(
  "Acetylation"     = "#E65100",
  "Citrullination"  = "#AD1457",
  "Cysteinylation"  = "#2E7D32",
  "Deamidation"     = "#1565C0",
  "Dimethylation"   = "#4527A0",
  "Methylation"     = "#558B2F",
  "N-Glycosylation" = "#00695C",
  "Oxidation"       = "#0277BD",
  "Phosphorylation" = "#6A1B9A",
  "SUMOylation"     = "#4A4A4A",
  "Ubiquitination"  = "#C62828"
)

# =============================================================================
# FUNCTION TO CREATE HEATMAP
# =============================================================================

create_binder_heatmap <- function(stats_df, value_col, title, color_palette, filename) {
  # Create matrix
  heatmap_df <- stats_df %>%
    select(PTM, Residue, PTM_Residue, Allele, !!sym(value_col)) %>%
    pivot_wider(names_from = Allele, values_from = !!sym(value_col)) %>%
    as.data.frame()

  # Sort by PTM first, then by Residue
  heatmap_df <- heatmap_df %>%
    arrange(PTM, Residue)

  # Set row names
  rownames(heatmap_df) <- heatmap_df$PTM_Residue

  # Keep only allele columns that exist
  allele_cols <- intersect(c("A0201", "B0702", "C0702"), colnames(heatmap_df))
  heatmap_matrix <- heatmap_df[, allele_cols, drop = FALSE]

  # Create row annotation showing PTM type
  row_annotation <- data.frame(PTM = heatmap_df$PTM)
  rownames(row_annotation) <- rownames(heatmap_matrix)

  # Get colors for PTMs present
  ptms_present <- unique(heatmap_df$PTM)
  ann_colors <- list(PTM = ptm_colors[ptms_present])

  # Calculate gaps between PTM groups
  ptm_rle <- rle(heatmap_df$PTM)
  gap_positions <- cumsum(ptm_rle$lengths)[-length(ptm_rle$lengths)]

  # Calculate figure height based on rows
  n_rows <- nrow(heatmap_matrix)
  fig_height <- max(8, n_rows * 0.25 + 2)

  # Save PNG
  png(file.path(output_dir, paste0(filename, ".png")),
      width = 8, height = fig_height, units = "in", res = 300)

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
    main = title,
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

  dev.off()

  # Save PDF
  pdf(file.path(output_dir, paste0(filename, ".pdf")),
      width = 8, height = fig_height)

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
    main = title,
    display_numbers = TRUE,
    number_format = "%.0f",
    number_color = "black",
    fontsize_number = 7,
    na_col = "gray90",
    legend = TRUE,
    annotation_row = row_annotation,
    annotation_colors = ann_colors,
    angle_col = 0,
    gaps_row = gap_positions
  )

  dev.off()

  cat("Saved:", filename, ".png/.pdf\n")

  return(heatmap_df)
}

# =============================================================================
# CREATE STRONG BINDERS HEATMAP
# =============================================================================

strong_colors <- colorRampPalette(c("#FFFFFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C"))(100)

strong_df <- create_binder_heatmap(
  stats_filtered,
  "pct_strong",
  "            % Strong Binders (EL_Rank < 0.5) by PTM+Residue and Allele",
  strong_colors,
  "Figure5A_StrongBinders_ByResidue"
)

# =============================================================================
# CREATE WEAK BINDERS HEATMAP
# =============================================================================

weak_colors <- colorRampPalette(c("#FFFFFF", "#FEE8C8", "#FDBB84", "#E34A33", "#B30000"))(100)

weak_df <- create_binder_heatmap(
  stats_filtered,
  "pct_weak",
  "            % Weak Binders (0.5 <= EL_Rank < 2) by PTM+Residue and Allele",
  weak_colors,
  "Figure5A_WeakBinders_ByResidue"
)

# =============================================================================
# CREATE PTM-ONLY VERSIONS (aggregating all residues)
# =============================================================================

ptm_only_stats <- merged %>%
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

# Strong binders by PTM
strong_ptm_df <- create_binder_heatmap(
  ptm_only_stats,
  "pct_strong",
  "            % Strong Binders (EL_Rank < 0.5) by PTM and Allele",
  strong_colors,
  "Figure5A_StrongBinders_ByPTM"
)

# Weak binders by PTM
weak_ptm_df <- create_binder_heatmap(
  ptm_only_stats,
  "pct_weak",
  "            % Weak Binders (0.5 <= EL_Rank < 2) by PTM and Allele",
  weak_colors,
  "Figure5A_WeakBinders_ByPTM"
)

# =============================================================================
# SAVE DATA
# =============================================================================

write.csv(stats_filtered, file.path(output_dir, "data_figure5A_binders.csv"), row.names = FALSE)
cat("\nSaved: data_figure5A_binders.csv\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== Summary ===\n")
cat("PTM+Residue combinations shown:", nrow(strong_df), "\n")
cat("PTM-only combinations shown:", nrow(strong_ptm_df), "\n")
cat("Sorted by: PTM type, then amino acid residue\n")

cat("\nTop 10 by % strong binders (mean across alleles):\n")
stats_summary <- stats_filtered %>%
  group_by(PTM, Residue, PTM_Residue) %>%
  summarise(
    n_total = sum(n),
    pct_strong_mean = round(mean(pct_strong, na.rm = TRUE), 1),
    pct_weak_mean = round(mean(pct_weak, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_strong_mean))

print(head(stats_summary, 10))
