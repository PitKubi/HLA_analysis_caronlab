# =============================================================================
# FIGURE 6D & 6E: Improved Specific PTM Analysis
# =============================================================================
# 6D: Biological Oxidation - Binding Affinity by Residue (stacked bar)
# 6E: N-Glycosylation - Glycan Types by Complexity with Binding Status
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

xlsx_file <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"
output_dir <- "figure_panels"
dir.create(output_dir, showWarnings = FALSE)

cat("=== Figure 6D & 6E: Improved Analysis ===\n\n")

# Color palette
binder_colors <- c(
  "Strong" = "#2E7D32",
  "Weak" = "#FFA726",
  "Non-binder" = "#90A4AE"
)

# =============================================================================
# FIGURE 6D: Biological Oxidation - Binding Affinity by Residue
# =============================================================================
cat("Creating Figure 6D: Biological Oxidation Binding by Residue...\n")

suppressMessages({
  biooxid <- read_excel(xlsx_file, sheet = "bioOxid.")
})

# Get binding data
binding_data <- biooxid %>%
  select(Peptide = Peptide...142, Length, Best_Allele, EL_Rank, Binder) %>%
  filter(!is.na(Peptide), !is.na(Binder)) %>%
  filter(Length >= 8 & Length <= 14)

# Residue column mapping
residue_cols <- c(
  "P" = "Peptide...2", "I" = "Peptide...11", "L" = "Peptide...19",
  "Q" = "Peptide...27", "S" = "Peptide...36", "T" = "Peptide...45",
  "V" = "Peptide...54", "C" = "Peptide...64", "D" = "Peptide...73",
  "E" = "Peptide...82", "N" = "Peptide...90", "Y" = "Peptide...98",
  "G" = "Peptide...107", "K" = "Peptide...117", "R" = "Peptide...127"
)

# Build peptide-to-residue mapping
peptide_residue_map <- data.frame()
for (res in names(residue_cols)) {
  col <- residue_cols[res]
  peps <- biooxid[[col]]
  peps <- peps[!is.na(peps)]
  if (length(peps) > 0) {
    peptide_residue_map <- rbind(peptide_residue_map,
      data.frame(Peptide = peps, Residue = res, stringsAsFactors = FALSE))
  }
}

# Join with binding data (allow many-to-many since same peptide can appear in multiple contexts)
binding_with_residue <- binding_data %>%
  inner_join(peptide_residue_map, by = "Peptide", relationship = "many-to-many")

cat("  Total peptide-residue-binding records:", nrow(binding_with_residue), "\n")

# Calculate percentages by residue
residue_summary <- binding_with_residue %>%
  group_by(Residue) %>%
  summarise(
    n_total = n(),
    n_strong = sum(Binder == "Strong"),
    n_weak = sum(Binder == "Weak"),
    n_nonbinder = sum(Binder == "Non-binder"),
    pct_strong = round(100 * n_strong / n_total, 1),
    pct_weak = round(100 * n_weak / n_total, 1),
    pct_nonbinder = round(100 * n_nonbinder / n_total, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_strong))

# Add amino acid full names
aa_names <- c(
  "A" = "Ala", "C" = "Cys", "D" = "Asp", "E" = "Glu", "F" = "Phe",
  "G" = "Gly", "H" = "His", "I" = "Ile", "K" = "Lys", "L" = "Leu",
  "M" = "Met", "N" = "Asn", "P" = "Pro", "Q" = "Gln", "R" = "Arg",
  "S" = "Ser", "T" = "Thr", "V" = "Val", "W" = "Trp", "Y" = "Tyr"
)

residue_summary$AA_Name <- aa_names[residue_summary$Residue]
residue_summary$Label <- paste0(residue_summary$Residue, " (", residue_summary$AA_Name, ")\nn=", residue_summary$n_total)

# Create long format for stacked bar
residue_long <- residue_summary %>%
  select(Residue, Label, pct_strong, pct_weak, pct_nonbinder) %>%
  pivot_longer(cols = starts_with("pct_"), names_to = "Binder_Type", values_to = "Percentage") %>%
  mutate(
    Binder_Type = factor(
      gsub("pct_", "", Binder_Type),
      levels = c("nonbinder", "weak", "strong"),
      labels = c("Non-binder", "Weak", "Strong")
    ),
    Residue = factor(Residue, levels = residue_summary$Residue)
  )

# Create stacked bar chart
fig6d <- ggplot(residue_long, aes(x = Residue, y = Percentage, fill = Binder_Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  geom_text(
    data = residue_summary,
    aes(x = Residue, y = 105, label = paste0(pct_strong, "%")),
    inherit.aes = FALSE, size = 3, fontface = "bold", color = "#2E7D32"
  ) +
  scale_fill_manual(values = binder_colors, name = "HLA Binding") +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25)) +
  labs(
    title = "Biological Oxidation: HLA Binding Affinity by Oxidized Residue",
    subtitle = "Peptides length 8-14 | Sorted by % Strong Binders (shown above bars)",
    x = "Oxidized Residue",
    y = "Percentage of Peptides"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
    axis.text.x = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(output_dir, "Figure6D_BioOxid_BindingByResidue.png"), fig6d,
       width = 12, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6D_BioOxid_BindingByResidue.pdf"), fig6d,
       width = 12, height = 6, bg = "white")
cat("  Saved Figure6D_BioOxid_BindingByResidue.png/pdf\n")

# Save data
write.csv(residue_summary, file.path(output_dir, "data_figure6D_biooxid_binding.csv"), row.names = FALSE)

# =============================================================================
# FIGURE 6D-2: Heatmap of Binding by Residue and Allele
# =============================================================================
cat("\nCreating Figure 6D-2: Binding heatmap by Residue and Allele...\n")

# Calculate binding stats by residue AND allele
residue_allele_stats <- binding_with_residue %>%
  filter(!is.na(Best_Allele)) %>%
  group_by(Residue, Allele = Best_Allele) %>%
  summarise(
    n = n(),
    pct_strong = round(100 * sum(Binder == "Strong") / n(), 0),
    .groups = "drop"
  ) %>%
  filter(n >= 5)  # At least 5 observations

# Create heatmap matrix
heatmap_df <- residue_allele_stats %>%
  select(Residue, Allele, pct_strong) %>%
  pivot_wider(names_from = Allele, values_from = pct_strong) %>%
  as.data.frame()

rownames(heatmap_df) <- heatmap_df$Residue
heatmap_df$Residue <- NULL

# Sort rows by mean % strong
row_means <- rowMeans(heatmap_df, na.rm = TRUE)
heatmap_df <- heatmap_df[order(row_means, decreasing = TRUE), ]

# Color palette for heatmap
heatmap_colors <- colorRampPalette(c("#FFFFFF", "#C8E6C9", "#81C784", "#4CAF50", "#2E7D32"))(100)

# Add row labels with AA names
row_labels <- paste0(rownames(heatmap_df), " (", aa_names[rownames(heatmap_df)], ")")

png(file.path(output_dir, "Figure6D2_BioOxid_Heatmap.png"),
    width = 8, height = 8, units = "in", res = 300)

pheatmap(
  as.matrix(heatmap_df),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = heatmap_colors,
  breaks = seq(0, 100, length.out = 101),
  border_color = "white",
  cellwidth = 45,
  cellheight = 25,
  fontsize_row = 11,
  fontsize_col = 11,
  main = "        Biological Oxidation: % Strong Binders by Residue and Allele",
  labels_row = row_labels,
  display_numbers = TRUE,
  number_format = "%.0f",
  number_color = "black",
  fontsize_number = 9,
  na_col = "gray90",
  legend = TRUE,
  legend_breaks = c(0, 25, 50, 75, 100),
  legend_labels = c("0%", "25%", "50%", "75%", "100%")
)

dev.off()

pdf(file.path(output_dir, "Figure6D2_BioOxid_Heatmap.pdf"),
    width = 8, height = 8)

pheatmap(
  as.matrix(heatmap_df),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = heatmap_colors,
  breaks = seq(0, 100, length.out = 101),
  border_color = "white",
  cellwidth = 45,
  cellheight = 25,
  fontsize_row = 11,
  fontsize_col = 11,
  main = "        Biological Oxidation: % Strong Binders by Residue and Allele",
  labels_row = row_labels,
  display_numbers = TRUE,
  number_format = "%.0f",
  number_color = "black",
  fontsize_number = 9,
  na_col = "gray90",
  legend = TRUE,
  legend_breaks = c(0, 25, 50, 75, 100),
  legend_labels = c("0%", "25%", "50%", "75%", "100%")
)

dev.off()

cat("  Saved Figure6D2_BioOxid_Heatmap.png/pdf\n")

# =============================================================================
# FIGURE 6E: N-Glycosylation - Glycan Types with Complexity Analysis
# =============================================================================
cat("\nCreating Figure 6E: Glycan Types by Complexity and Binding...\n")

suppressMessages({
  nglyco <- read_excel(xlsx_file, sheet = "N-glyco")
  glycan_ref <- read_excel(xlsx_file, sheet = "Supp_glycan list ref")
})

# Filter to length 8-14 with binding data
nglyco_filt <- nglyco %>%
  filter(!is.na(StrippedSequences)) %>%
  mutate(Peptide = StrippedSequences, Length = nchar(Peptide)) %>%
  filter(Length >= 8 & Length <= 14) %>%
  filter(!is.na(Binder))

# Extract UniMod number
nglyco_filt$UniMod_Num <- as.numeric(gsub("UniMod:", "", nglyco_filt$UniMod))

# Join with glycan reference for mass shift (complexity indicator)
nglyco_with_ref <- nglyco_filt %>%
  left_join(
    glycan_ref %>% select(`UniMod#`, `Glycan name`, `Mass shift (Da)`),
    by = c("UniMod_Num" = "UniMod#")
  )

# Use Linked glycan column, fall back to Glycan name
nglyco_with_ref$Glycan_Label <- ifelse(
  !is.na(nglyco_with_ref$`Linked glycan`),
  nglyco_with_ref$`Linked glycan`,
  nglyco_with_ref$`Glycan name`
)

# Handle NA glycan labels
nglyco_with_ref$Glycan_Label[is.na(nglyco_with_ref$Glycan_Label)] <- "Unknown"

# Calculate complexity based on monosaccharide count
# Count occurrences of monosaccharide components
count_monosaccharides <- function(glycan_str) {
  if (is.na(glycan_str) || glycan_str == "Unknown") return(0)

  # Extract numbers from patterns like Hex(5), HexNAc(4), etc.
  numbers <- as.numeric(gsub("[^0-9]", "", unlist(regmatches(glycan_str,
    gregexpr("\\([0-9]+\\)", glycan_str)))))

  if (length(numbers) == 0) return(1)  # Simple glycan
  return(sum(numbers))
}

nglyco_with_ref$Complexity <- sapply(nglyco_with_ref$Glycan_Label, count_monosaccharides)

# Categorize complexity
nglyco_with_ref$Complexity_Cat <- cut(
  nglyco_with_ref$Complexity,
  breaks = c(-Inf, 3, 8, 15, Inf),
  labels = c("Simple (1-3)", "Medium (4-8)", "Complex (9-15)", "Very Complex (>15)")
)

# Handle Unknown - convert to character first, then back to factor
nglyco_with_ref$Complexity_Cat <- as.character(nglyco_with_ref$Complexity_Cat)
nglyco_with_ref$Complexity_Cat[nglyco_with_ref$Glycan_Label == "Unknown"] <- "Unknown"
nglyco_with_ref$Complexity_Cat[is.na(nglyco_with_ref$Complexity_Cat)] <- "Unknown"
nglyco_with_ref$Complexity_Cat <- factor(nglyco_with_ref$Complexity_Cat,
  levels = c("Simple (1-3)", "Medium (4-8)", "Complex (9-15)", "Very Complex (>15)", "Unknown"))

cat("  Glycan complexity distribution:\n")
print(table(nglyco_with_ref$Complexity_Cat))

# Summary by complexity category
complexity_summary <- nglyco_with_ref %>%
  group_by(Complexity_Cat) %>%
  summarise(
    n = n(),
    n_strong = sum(Binder == "Strong"),
    n_weak = sum(Binder == "Weak"),
    n_nonbinder = sum(Binder == "Non-binder"),
    pct_binder = round(100 * (n_strong + n_weak) / n, 1),
    pct_strong = round(100 * n_strong / n, 1),
    .groups = "drop"
  )

cat("  Binding by complexity:\n")
print(complexity_summary)

# Create bar chart by complexity with binding status
complexity_long <- nglyco_with_ref %>%
  group_by(Complexity_Cat, Binder) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Binder = factor(Binder, levels = c("Non-binder", "Weak", "Strong")))

# Colors for complexity
complexity_colors <- c(
  "Simple (1-3)" = "#E8F5E9",
  "Medium (4-8)" = "#A5D6A7",
  "Complex (9-15)" = "#66BB6A",
  "Very Complex (>15)" = "#2E7D32",
  "Unknown" = "#BDBDBD"
)

fig6e_1 <- ggplot(complexity_long, aes(x = Complexity_Cat, y = n, fill = Binder)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(
    data = complexity_summary,
    aes(x = Complexity_Cat, y = n + 1, label = paste0("n=", n, "\n(", pct_binder, "% bind)")),
    inherit.aes = FALSE, size = 3, vjust = 0
  ) +
  scale_fill_manual(values = binder_colors, name = "HLA Binding") +
  labs(
    title = "N-Glycosylation: Glycan Complexity vs HLA Binding",
    subtitle = "Complexity = sum of monosaccharide units | Peptides length 8-14",
    x = "Glycan Complexity",
    y = "Number of Peptides"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
    axis.text.x = element_text(size = 9, angle = 15, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  ylim(0, max(complexity_summary$n) * 1.3)

ggsave(file.path(output_dir, "Figure6E_Glycan_Complexity.png"), fig6e_1,
       width = 9, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6E_Glycan_Complexity.pdf"), fig6e_1,
       width = 9, height = 6, bg = "white")
cat("  Saved Figure6E_Glycan_Complexity.png/pdf\n")

# =============================================================================
# FIGURE 6E-2: Detailed Glycan Types Table/Bar
# =============================================================================
cat("\nCreating Figure 6E-2: Detailed Glycan Types...\n")

# Summary by individual glycan type
glycan_type_summary <- nglyco_with_ref %>%
  filter(Glycan_Label != "Unknown") %>%
  group_by(Glycan_Label, Complexity) %>%
  summarise(
    n = n(),
    n_strong = sum(Binder == "Strong"),
    n_weak = sum(Binder == "Weak"),
    n_nonbinder = sum(Binder == "Non-binder"),
    mass_shift = first(`Mass shift (Da)`),
    .groups = "drop"
  ) %>%
  arrange(Complexity, desc(n))

# Create shortened labels for plotting
glycan_type_summary$Short_Label <- gsub("HexNAc", "HN", glycan_type_summary$Glycan_Label)
glycan_type_summary$Short_Label <- gsub("NeuAc", "NA", glycan_type_summary$Short_Label)
glycan_type_summary$Short_Label <- gsub("dHex", "dH", glycan_type_summary$Short_Label)
glycan_type_summary$Short_Label <- gsub("Phos", "P", glycan_type_summary$Short_Label)

cat("  Individual glycan types:\n")
print(glycan_type_summary %>% select(Glycan_Label, Complexity, n, n_strong, n_weak, mass_shift))

# Create horizontal bar chart sorted by complexity
glycan_type_summary$Glycan_Label <- factor(
  glycan_type_summary$Glycan_Label,
  levels = glycan_type_summary$Glycan_Label
)

glycan_long <- glycan_type_summary %>%
  select(Glycan_Label, Short_Label, Complexity, n_strong, n_weak, n_nonbinder) %>%
  pivot_longer(cols = c(n_strong, n_weak, n_nonbinder),
               names_to = "Binder_Type", values_to = "Count") %>%
  mutate(
    Binder_Type = factor(
      gsub("n_", "", Binder_Type),
      levels = c("nonbinder", "weak", "strong"),
      labels = c("Non-binder", "Weak", "Strong")
    )
  )

fig6e_2 <- ggplot(glycan_long, aes(y = Glycan_Label, x = Count, fill = Binder_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    data = glycan_type_summary,
    aes(y = Glycan_Label, x = n + 0.3,
        label = paste0("(", Complexity, " units)")),
    inherit.aes = FALSE, size = 2.5, hjust = 0, color = "gray50"
  ) +
  scale_fill_manual(values = binder_colors, name = "HLA Binding") +
  labs(
    title = "N-Glycosylation: Individual Glycan Types",
    subtitle = "Sorted by glycan complexity (monosaccharide count)",
    x = "Number of Peptides",
    y = "Glycan Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25)))

ggsave(file.path(output_dir, "Figure6E2_Glycan_Types_Detail.png"), fig6e_2,
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6E2_Glycan_Types_Detail.pdf"), fig6e_2,
       width = 10, height = 8, bg = "white")
cat("  Saved Figure6E2_Glycan_Types_Detail.png/pdf\n")

# Save glycan data
write.csv(glycan_type_summary, file.path(output_dir, "data_figure6E_glycan_analysis.csv"), row.names = FALSE)

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=== Figure 6D & 6E Summary ===\n")
cat("Figure 6D: Biological Oxidation binding affinity by residue\n")
cat("  - 15 residue types analyzed\n")
cat("  - Highest % strong binders: I (", residue_summary$pct_strong[1], "%), ",
    residue_summary$Residue[2], " (", residue_summary$pct_strong[2], "%)\n", sep = "")
cat("  - Lowest % strong binders: C (", residue_summary$pct_strong[nrow(residue_summary)], "%)\n", sep = "")
cat("\nFigure 6D-2: Heatmap of % strong binders by residue and allele\n")
cat("\nFigure 6E: Glycan complexity analysis\n")
cat("  - ", nrow(nglyco_with_ref), " glycopeptides analyzed\n", sep = "")
cat("  - Complexity categories: Simple, Medium, Complex, Very Complex\n")
cat("\nFigure 6E-2: Individual glycan types detail\n")
cat("  - ", nrow(glycan_type_summary), " unique glycan types\n", sep = "")
cat("\nAll figures saved to:", output_dir, "\n")
