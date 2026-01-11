# =============================================================================
# FIGURE 5A: % Strong Binders by PTM+Residue and Allele
# =============================================================================
# Heatmap showing % strong binders, deconvoluted by amino acid residue
# Sorted by PTM type, then by residue within each PTM

library(readxl)
library(dplyr)
library(tidyr)
library(pheatmap)

xlsx_file <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"

cat("=== Figure 5A: % Strong Binders by PTM+Residue ===\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

# Load our parsed PTM data with residue info
ptm_data <- read.csv("figure_panels/data_figure4A_circos.csv", stringsAsFactors = FALSE)
ptm_data$PTM_Residue <- paste(ptm_data$Residue, tolower(ptm_data$PTM))

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
            "Deamid. (NQ)", "bioOxid.", "Citrullination", "G-Ubiq.")

cat("Loading binding data...\n")
all_binding <- list()
for (sheet in sheets) {
  binding <- get_binding_for_sheet(sheet)
  if (!is.null(binding)) {
    all_binding[[sheet]] <- binding
    cat("  ", sheet, ":", nrow(binding), "\n")
  }
}

binding_combined <- do.call(rbind, all_binding)

# Merge with PTM data
merged <- ptm_data %>%
  inner_join(binding_combined, by = "Peptide", relationship = "many-to-many")

cat("\nTotal merged records:", nrow(merged), "\n")

# =============================================================================
# CALCULATE % STRONG BINDERS
# =============================================================================

# Calculate stats for each PTM+Residue x Allele
stats <- merged %>%
  filter(!is.na(Allele)) %>%
  group_by(PTM, Residue, PTM_Residue, Allele) %>%
  summarise(
    n = n(),
    n_strong = sum(Binder == "Strong"),
    n_weak = sum(Binder == "Weak"),
    n_nonbinder = sum(Binder == "Non-binder"),
    pct_strong = round(100 * n_strong / n, 1),
    .groups = "drop"
  )

# Filter to combinations with at least 10 peptides
stats_filtered <- stats %>%
  filter(n >= 10)

cat("\nPTM+Residue combinations with n >= 10:",
    length(unique(stats_filtered$PTM_Residue)), "\n")

# =============================================================================
# CREATE HEATMAP MATRIX
# =============================================================================

# Create matrix
heatmap_df <- stats_filtered %>%
  select(PTM, Residue, PTM_Residue, Allele, pct_strong) %>%
  pivot_wider(names_from = Allele, values_from = pct_strong) %>%
  as.data.frame()

# Sort by PTM first, then by Residue
heatmap_df <- heatmap_df %>%
  arrange(PTM, Residue)

# Set row names
rownames(heatmap_df) <- heatmap_df$PTM_Residue

# Keep only allele columns
heatmap_matrix <- heatmap_df %>%
  select(A0201, B0702, C0702)

# Handle any NAs
heatmap_matrix[is.na(heatmap_matrix)] <- NA

cat("\nHeatmap dimensions:", nrow(heatmap_matrix), "x", ncol(heatmap_matrix), "\n")

# =============================================================================
# LOAD BACKGROUND FOR COMPARISON
# =============================================================================

suppressMessages({
  bg <- read_excel(xlsx_file, sheet = "Background (wo PTMs)")
})
bg_stats <- bg %>%
  filter(!is.na(Binder) & !is.na(Best_Allele)) %>%
  group_by(Allele = Best_Allele) %>%
  summarise(
    n = n(),
    pct_strong = round(100 * sum(Binder == "Strong") / n, 1),
    .groups = "drop"
  )

cat("\nBackground % Strong Binders:\n")
print(bg_stats)

# =============================================================================
# CREATE HEATMAP
# =============================================================================

# Color palette: white to blue (higher % = more strong binders)
colors <- colorRampPalette(c("#FFFFFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C"))(100)

# Create row annotation showing PTM type
row_annotation <- data.frame(
  PTM = heatmap_df$PTM
)
rownames(row_annotation) <- rownames(heatmap_matrix)

# PTM colors
ptm_types <- unique(heatmap_df$PTM)
ptm_colors <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#FFFF33", "#A65628", "#F781BF", "#999999")[1:length(ptm_types)],
  ptm_types
)
ann_colors <- list(PTM = ptm_colors)

# Calculate figure height based on rows
n_rows <- nrow(heatmap_matrix)
fig_height <- max(8, n_rows * 0.25 + 2)

# Save heatmap
png("figure_panels/Figure5A_StrongBinders_ByResidue.png",
    width = 8, height = fig_height, units = "in", res = 300)

pheatmap(
  as.matrix(heatmap_matrix),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colors,
  breaks = seq(0, 100, length.out = 101),
  border_color = "white",
  cellwidth = 40,
  cellheight = 14,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "% Strong Binders by PTM+Residue and Allele",
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
  gaps_row = cumsum(rle(heatmap_df$PTM)$lengths)[-length(rle(heatmap_df$PTM)$lengths)]
)

dev.off()
cat("\nSaved: figure_panels/Figure5A_StrongBinders_ByResidue.png\n")

# PDF version
pdf("figure_panels/Figure5A_StrongBinders_ByResidue.pdf",
    width = 8, height = fig_height)

pheatmap(
  as.matrix(heatmap_matrix),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colors,
  breaks = seq(0, 100, length.out = 101),
  border_color = "white",
  cellwidth = 40,
  cellheight = 14,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "% Strong Binders by PTM+Residue and Allele",
  display_numbers = TRUE,
  number_format = "%.0f",
  number_color = "black",
  fontsize_number = 7,
  na_col = "gray90",
  legend = TRUE,
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  angle_col = 0,
  gaps_row = cumsum(rle(heatmap_df$PTM)$lengths)[-length(rle(heatmap_df$PTM)$lengths)]
)

dev.off()
cat("Saved: figure_panels/Figure5A_StrongBinders_ByResidue.pdf\n")

# =============================================================================
# SAVE DATA
# =============================================================================

# Add overall stats
stats_summary <- stats_filtered %>%
  group_by(PTM, Residue, PTM_Residue) %>%
  summarise(
    n_total = sum(n),
    pct_strong_mean = round(mean(pct_strong, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(PTM, Residue)

write.csv(stats_filtered, "figure_panels/data_figure5A_strongbinders.csv", row.names = FALSE)
cat("Saved: figure_panels/data_figure5A_strongbinders.csv\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== Summary ===\n")
cat("PTM+Residue combinations shown:", nrow(heatmap_matrix), "\n")
cat("Sorted by: PTM type, then amino acid residue\n")
cat("\nTop 10 by % strong binders (mean across alleles):\n")
print(stats_summary %>% arrange(desc(pct_strong_mean)) %>% head(10))

cat("\nBottom 10 by % strong binders:\n")
print(stats_summary %>% arrange(pct_strong_mean) %>% head(10))
