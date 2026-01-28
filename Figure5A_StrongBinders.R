# =============================================================================
# FIGURE 5A: % Strong & Weak Binders by PTM+Residue and Allele
# =============================================================================
# Heatmap showing % strong/weak binders, deconvoluted by amino acid residue
# Sorted by PTM type, then by residue within each PTM
#
# Method (per collaborator):
#   1. Read binding data directly from each Excel sheet's binding section
#   2. Extract residue info from dedicated columns or parse modifications
#   3. Expand multi-site residues (e.g., "S,T" → both S and T)
#   4. Count by (Binder, Allele, Residue) per PTM
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(pheatmap)

xlsx_file <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"
output_dir <- "figure_panels"

cat("=== Figure 5A: % Strong/Weak Binders by PTM+Residue ===\n\n")

# =============================================================================
# HELPER: Parse residue(s) from Assigned Modifications, filtering by mass
# =============================================================================
#' @param mod_string e.g. "4K(14.0156), 6R(14.0156)"
#' @param target_mass e.g. 14.0156 for methylation
#' @param tol mass tolerance (default 0.02)
#' @return comma-separated residues, or NA
get_residues_by_mass <- function(mod_string, target_mass, tol = 0.02) {
  if (is.na(mod_string) || mod_string == "") return(NA_character_)
  mods <- strsplit(as.character(mod_string), ",\\s*")[[1]]
  residues <- c()
  for (m in mods) {
    m <- trimws(m)
    if (grepl("^N-term\\(", m)) {
      mass <- as.numeric(gsub("N-term\\((.+)\\)", "\\1", m))
      if (!is.na(mass) && abs(mass - target_mass) < tol) residues <- c(residues, "N-term")
      next
    }
    match <- regmatches(m, regexpr("^(\\d+)([A-Z])", m))
    if (length(match) > 0 && nchar(match) >= 2) {
      res <- gsub("[0-9]", "", match)
      mass_str <- regmatches(m, regexpr("\\(([0-9.]+)\\)", m))
      if (length(mass_str) > 0) {
        mass <- as.numeric(gsub("[()]", "", mass_str))
        if (!is.na(mass) && abs(mass - target_mass) < tol) residues <- c(residues, res)
      }
    }
  }
  if (length(residues) == 0) return(NA_character_)
  paste(unique(residues), collapse = ",")
}

# =============================================================================
# LOAD BINDING + RESIDUE DATA DIRECTLY FROM EACH SHEET
# =============================================================================

cat("Loading binding data from Excel sheets...\n")

# Each binding section config: list(pep, allele, binder, and residue source)
# Residue source is one of:
#   fixed_residue = "C"           (constant)
#   residue_col = "Phospho residue"  (read from column)
#   mod_col + mass                (parse from Assigned Modifications)
binding_configs <- list(
  list(sheet = "Phospho.", ptm = "Phosphorylation", sections = list(
    list(pep = "Peptide", allele = "Best_Allele", binder = "Binder",
         residue_col = "Phospho residue")
  )),
  list(sheet = "Cyst.", ptm = "Cysteinylation", sections = list(
    list(pep = "Peptide", allele = "Best_Allele", binder = "Binder",
         fixed_residue = "C")
  )),
  list(sheet = "Deamid. (NQ)", ptm = "Deamidation", sections = list(
    list(pep = "Peptide...7", allele = "Best_Allele...8", binder = "Binder...10",
         fixed_residue = "N"),
    list(pep = "Peptide...16", allele = "Best_Allele...17", binder = "Binder...19",
         fixed_residue = "Q")
  )),
  list(sheet = "Acetyl.", ptm = "Acetylation", sections = list(
    list(pep = "Peptide...12", allele = "Best_Allele...13", binder = "Binder...15",
         fixed_residue = "N-term"),
    list(pep = "Peptide...21", allele = "Best_Allele...22", binder = "Binder...24",
         fixed_residue = "K")
  )),
  list(sheet = "GG-Ubiq.", ptm = "Ubiquitination", sections = list(
    list(pep = "Peptide", allele = "Best_Allele", binder = "Binder",
         residue_col = "Phospho residue")
  )),
  list(sheet = "G-Ubiq.", ptm = "Ubiquitination", sections = list(
    list(pep = "Peptide...13", allele = "Best_Allele", binder = "Binder",
         residue_col = "Phospho residue")
  )),
  list(sheet = "Methyl.", ptm = "Methylation", sections = list(
    list(pep = "Peptide...20", allele = "Best_Allele...21", binder = "Binder...23",
         mod_col = "Assigned Modifications...26", mass = 14.0156),
    list(pep = "Peptide...29", allele = "Best_Allele...30", binder = "Binder...32",
         mod_col = "Assigned Modifications...35", mass = 14.0156)
  )),
  list(sheet = "Dimethyl.", ptm = "Dimethylation", sections = list(
    list(pep = "Peptide...11", allele = "Best_Allele", binder = "Binder",
         mod_col = "Assigned Modifications...17", mass = 28.0313)
  )),
  list(sheet = "Citrullination", ptm = "Citrullination", sections = list(
    list(pep = "Peptide...11", allele = "Best_Allele", binder = "Binder",
         fixed_residue = "R")
  )),
  list(sheet = "bioOxid.", ptm = "Oxidation", sections = list(
    list(pep = "Peptide...142", allele = "Best_Allele", binder = "Binder",
         mod_col = "Assigned Modifications...143", mass = 15.9949)
  )),
  list(sheet = "artOxid.", ptm = "Artifact Oxidation", sections = list(
    list(pep = "Peptide...36", allele = "Best_Allele", binder = "Binder",
         mod_col = "Assigned Modifications...37", mass = 15.9949)
  )),
  list(sheet = "SUMO", ptm = "SUMOylation", sections = list(
    list(pep = "Peptide...30", allele = "Best_Allele", binder = "Binder",
         fixed_residue = "K")
  )),
  list(sheet = "N-glyco", ptm = "N-Glycosylation", sections = list(
    list(pep = "Peptide", allele = "Best_Allele", binder = "Binder",
         fixed_residue = "N")
  )),
  list(sheet = "Carbamid.", ptm = "Carbamidomethylation", sections = list(
    list(pep = "Peptide...11", allele = "Best_Allele", binder = "Binder",
         fixed_residue = "C")
  ))
)

all_binding <- list()
for (cfg in binding_configs) {
  df <- tryCatch(suppressMessages(read_excel(xlsx_file, sheet = cfg$sheet)),
                 error = function(e) NULL)
  if (is.null(df)) { cat("  SKIP:", cfg$sheet, "\n"); next }

  for (sec in cfg$sections) {
    # Check required columns exist
    if (!all(c(sec$pep, sec$allele, sec$binder) %in% colnames(df))) next

    result <- data.frame(
      PTM = cfg$ptm,
      Peptide = df[[sec$pep]],
      Allele = df[[sec$allele]],
      Binder = df[[sec$binder]],
      stringsAsFactors = FALSE
    )

    # Determine residue
    if (!is.null(sec$fixed_residue)) {
      result$Residue <- sec$fixed_residue
    } else if (!is.null(sec$residue_col) && sec$residue_col %in% colnames(df)) {
      result$Residue <- df[[sec$residue_col]]
    } else if (!is.null(sec$mod_col) && sec$mod_col %in% colnames(df)) {
      result$Residue <- sapply(df[[sec$mod_col]], get_residues_by_mass,
                               target_mass = sec$mass)
    }

    result <- result[!is.na(result$Peptide) & !is.na(result$Binder) &
                     !is.na(result$Residue), ]
    if (nrow(result) > 0) {
      all_binding[[length(all_binding) + 1]] <- result
      cat("  ", cfg$ptm, "(", cfg$sheet, sec$pep, "):", nrow(result), "\n")
    }
  }
}

merged <- do.call(rbind, all_binding)

# Expand multi-site residues: "S,T" → one row for S and one row for T
# This follows the collaborator's counting method
expanded <- list()
for (i in 1:nrow(merged)) {
  residues <- strsplit(merged$Residue[i], ",")[[1]]
  for (res in unique(residues)) {
    expanded[[length(expanded) + 1]] <- data.frame(
      PTM = merged$PTM[i],
      Peptide = merged$Peptide[i],
      Allele = merged$Allele[i],
      Binder = merged$Binder[i],
      Residue = res,
      stringsAsFactors = FALSE
    )
  }
}
merged <- do.call(rbind, expanded)
merged$PTM_Residue <- paste(merged$Residue, tolower(merged$PTM))

# Filter out NA alleles
merged <- merged %>% filter(!is.na(Allele))

cat("\nTotal binding records (after expansion):", nrow(merged), "\n")

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

config <- readRDS(file.path(output_dir, "config.rds"))
ptm_colors <- config$ptm_colors

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

  # Replace NA with 0 for consistent display
  heatmap_matrix[is.na(heatmap_matrix)] <- 0

  # Replace NA with 0 for visual consistency (NA = insufficient data for that allele)
  heatmap_matrix[is.na(heatmap_matrix)] <- 0

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
