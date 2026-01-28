# =============================================================================
# FIGURE 6: Specific PTM Analysis Panels
# =============================================================================
# 6A: Doubly phosphorylated peptides - site distribution
# 6B: Ubiquitination G vs GG distribution
# 6C: Oxidation biological vs artifact (pie chart)
# 6D: Biological oxidation residue distribution
# 6E: Glycosylation glycan types distribution
# =============================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

xlsx_file <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"
output_dir <- "figure_panels"
dir.create(output_dir, showWarnings = FALSE)

cat("=== Figure 6: Specific PTM Analysis ===\n\n")

# Color palette
ptm_colors <- c(
  "Cysteinylation"  = "#2E7D32",
  "Deamidation"     = "#1565C0",
  "Phosphorylation" = "#6A1B9A",
  "Acetylation"     = "#E65100",
  "Ubiquitination"  = "#C62828",
  "Citrullination"  = "#AD1457",
  "Dimethylation"   = "#4527A0",
  "Methylation"     = "#558B2F",
  "Oxidation"       = "#0277BD",
  "N-Glycosylation" = "#00695C"
)

# =============================================================================
# FIGURE 6A: Doubly Phosphorylated Peptides - Site Distribution
# =============================================================================
cat("Creating Figure 6A: Doubly phosphorylated peptides...\n")

phospho <- read_excel(xlsx_file, sheet = "Phospho.")

# Get peptides with multiple phospho sites
phospho_multi <- phospho %>%
  filter(!is.na(`All sites`)) %>%
  mutate(
    Peptide = `Stripped Peptide sequences`,
    Length = nchar(Peptide),
    Sites = `All sites`
  ) %>%
  filter(Length >= 8 & Length <= 14) %>%
  filter(grepl(",", Sites))  # Contains comma = multiple sites

cat("  Doubly+ phosphorylated peptides (length 8-14):", nrow(phospho_multi), "\n")

# Parse all sites and count positions
all_sites <- c()
for (sites in phospho_multi$Sites) {
  site_list <- as.numeric(strsplit(as.character(sites), ",")[[1]])
  all_sites <- c(all_sites, site_list)
}

# Filter to positions 1-14
all_sites <- all_sites[all_sites >= 1 & all_sites <= 14]

# Create barplot
site_counts <- as.data.frame(table(all_sites))
colnames(site_counts) <- c("Position", "Count")
site_counts$Position <- as.numeric(as.character(site_counts$Position))

fig6a <- ggplot(site_counts, aes(x = Position, y = Count)) +
  geom_bar(stat = "identity", fill = ptm_colors["Phosphorylation"], width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 3) +
  scale_x_continuous(breaks = 1:14) +
  labs(
    title = "Doubly Phosphorylated Peptides: Site Distribution",
    subtitle = paste0("n = ", nrow(phospho_multi), " peptides (length 8-14)"),
    x = "Phosphorylation Site Position",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "Figure6A_DoublyPhospho_Sites.png"), fig6a,
       width = 8, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6A_DoublyPhospho_Sites.pdf"), fig6a,
       width = 8, height = 5, bg = "white")
cat("  Saved Figure6A_DoublyPhospho_Sites.png/pdf\n")

# =============================================================================
# FIGURE 6B: Ubiquitination G vs GG Distribution
# =============================================================================
cat("\nCreating Figure 6B: Ubiquitination G vs GG...\n")

# Load G-Ubiq
g_ubiq <- read_excel(xlsx_file, sheet = "G-Ubiq.")
g_data <- g_ubiq %>%
  select(Peptide = `Peptide...2`, Length = `Peptide Length`) %>%
  filter(!is.na(Peptide)) %>%
  mutate(Length = as.numeric(Length), Type = "G-Ubiq") %>%
  filter(Length >= 8 & Length <= 14)

# Load GG-Ubiq
gg_ubiq <- read_excel(xlsx_file, sheet = "GG-Ubiq.")
gg_data <- gg_ubiq %>%
  select(Peptide = `stripped peptide`, Length = `Length...7`) %>%
  filter(!is.na(Peptide)) %>%
  mutate(Length = as.numeric(Length), Type = "GG-Ubiq") %>%
  filter(Length >= 8 & Length <= 14)

# Combine
ubiq_combined <- bind_rows(g_data, gg_data)

cat("  G-Ubiq (length 8-14):", nrow(g_data), "\n")
cat("  GG-Ubiq (length 8-14):", nrow(gg_data), "\n")

# Create comparison barplot
ubiq_summary <- ubiq_combined %>%
  group_by(Type) %>%
  summarise(Count = n(), .groups = "drop")

fig6b <- ggplot(ubiq_summary, aes(x = Type, y = Count, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("G-Ubiq" = "#E57373", "GG-Ubiq" = "#C62828")) +
  labs(
    title = "Ubiquitination: G vs GG Distribution",
    subtitle = "Peptides length 8-14",
    x = "Ubiquitination Type",
    y = "Number of Peptides"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  ylim(0, max(ubiq_summary$Count) * 1.15)

ggsave(file.path(output_dir, "Figure6B_Ubiq_G_vs_GG.png"), fig6b,
       width = 6, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6B_Ubiq_G_vs_GG.pdf"), fig6b,
       width = 6, height = 5, bg = "white")
cat("  Saved Figure6B_Ubiq_G_vs_GG.png/pdf\n")

# =============================================================================
# FIGURE 6C: Oxidation - Biological vs Artifact (Pie Chart)
# =============================================================================
cat("\nCreating Figure 6C: Biological vs Artifact Oxidation...\n")

# Count biological oxidation (bioOxid sheet)
bio_ox <- read_excel(xlsx_file, sheet = "bioOxid.")
# Get peptide columns and count unique peptides
bio_pep_cols <- grep("^Peptide", colnames(bio_ox), value = TRUE)
bio_peptides <- c()
for (col in bio_pep_cols) {
  peps <- bio_ox[[col]]
  peps <- peps[!is.na(peps)]
  # Filter by length
  peps <- peps[nchar(peps) >= 8 & nchar(peps) <= 14]
  bio_peptides <- c(bio_peptides, peps)
}
n_bio <- length(unique(bio_peptides))

# Count artifact oxidation (artOxid sheet - M, W, H)
art_ox <- read_excel(xlsx_file, sheet = "artOxid.")
art_pep_cols <- grep("^Peptide", colnames(art_ox), value = TRUE)
art_peptides <- c()
for (col in art_pep_cols) {
  peps <- art_ox[[col]]
  peps <- peps[!is.na(peps)]
  peps <- peps[nchar(peps) >= 8 & nchar(peps) <= 14]
  art_peptides <- c(art_peptides, peps)
}
n_art <- length(unique(art_peptides))

cat("  Biological oxidation peptides:", n_bio, "\n")
cat("  Artifact oxidation peptides:", n_art, "\n")

# Create pie chart data
ox_pie_data <- data.frame(
  Type = c("Biological\nOxidation", "Artifact\nOxidation"),
  Count = c(n_bio, n_art)
) %>%
  mutate(
    Percentage = round(100 * Count / sum(Count), 1),
    Label = paste0(Count, "\n(", Percentage, "%)")
  )

fig6c <- ggplot(ox_pie_data, aes(x = "", y = Count, fill = Type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c("Biological\nOxidation" = "#0277BD", "Artifact\nOxidation" = "#90A4AE")) +
  labs(
    title = "Oxidation: Biological vs Artifact",
    subtitle = "Unique peptides (length 8-14)"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "Figure6C_Oxidation_Bio_vs_Art.png"), fig6c,
       width = 6, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6C_Oxidation_Bio_vs_Art.pdf"), fig6c,
       width = 6, height = 6, bg = "white")
cat("  Saved Figure6C_Oxidation_Bio_vs_Art.png/pdf\n")

# =============================================================================
# FIGURE 6D: Biological Oxidation Residue Distribution
# =============================================================================
cat("\nCreating Figure 6D: Biological oxidation residue distribution...\n")

# Get residue distribution from bioOxid
bio_ox_raw <- read_excel(xlsx_file, sheet = "bioOxid.", col_names = FALSE, n_max = 1)
residue_headers <- as.character(bio_ox_raw[1,])

# Extract peptide counts per residue
bio_ox_data <- read_excel(xlsx_file, sheet = "bioOxid.")
residue_counts <- data.frame(Residue = character(), Count = integer(), stringsAsFactors = FALSE)

# Find columns that start with single letter residues
for (i in seq_along(residue_headers)) {
  res <- residue_headers[i]
  if (!is.na(res) && nchar(res) == 1 && res %in% LETTERS) {
    # Look for the peptide column after this residue header
    next_cols <- colnames(bio_ox_data)[i:(min(i+3, ncol(bio_ox_data)))]
    pep_col <- next_cols[grepl("^Peptide", next_cols)][1]
    if (!is.na(pep_col)) {
      peps <- bio_ox_data[[pep_col]]
      peps <- peps[!is.na(peps)]
      peps <- peps[nchar(peps) >= 8 & nchar(peps) <= 14]
      if (length(peps) > 0) {
        residue_counts <- rbind(residue_counts, data.frame(Residue = res, Count = length(peps)))
      }
    }
  }
}

# Aggregate by residue
residue_agg <- residue_counts %>%
  group_by(Residue) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  mutate(Percentage = round(100 * Count / sum(Count), 1)) %>%
  arrange(desc(Count))

cat("  Residue distribution:\n")
print(residue_agg)

fig6d <- ggplot(residue_agg, aes(x = reorder(Residue, -Count), y = Count, fill = Residue)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = rep(c("#0277BD", "#4FC3F7"), length.out = nrow(residue_agg))) +
  labs(
    title = "Biological Oxidation: Residue Distribution",
    subtitle = "Peptides length 8-14",
    x = "Oxidized Residue",
    y = "Number of Peptides"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold", size = 10)
  ) +
  ylim(0, max(residue_agg$Count) * 1.2)

ggsave(file.path(output_dir, "Figure6D_BioOxid_Residues.png"), fig6d,
       width = 10, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6D_BioOxid_Residues.pdf"), fig6d,
       width = 10, height = 5, bg = "white")
cat("  Saved Figure6D_BioOxid_Residues.png/pdf\n")

# =============================================================================
# FIGURE 6E: Glycosylation - Glycan Types Distribution
# =============================================================================
cat("\nCreating Figure 6E: Glycan types distribution...\n")

nglyco <- read_excel(xlsx_file, sheet = "N-glyco")

# Filter to length 8-14
nglyco_filt <- nglyco %>%
  filter(!is.na(StrippedSequences)) %>%
  mutate(Peptide = StrippedSequences, Length = nchar(Peptide)) %>%
  filter(Length >= 8 & Length <= 14)

# Get glycan types
glycan_counts <- nglyco_filt %>%
  filter(!is.na(`Linked glycan`)) %>%
  group_by(`Linked glycan`) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = round(100 * Count / sum(Count), 1)) %>%
  arrange(desc(Count))

cat("  Glycan types found:", nrow(glycan_counts), "\n")
print(glycan_counts)

# Simplify glycan names for display
glycan_counts$ShortName <- gsub("\\(", "", glycan_counts$`Linked glycan`)
glycan_counts$ShortName <- gsub("\\)", "", glycan_counts$ShortName)
glycan_counts$ShortName <- substr(glycan_counts$ShortName, 1, 20)  # Truncate long names

# Create pie chart
fig6e <- ggplot(glycan_counts, aes(x = "", y = Count, fill = `Linked glycan`)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = ifelse(Percentage >= 3, paste0(Percentage, "%"), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_brewer(palette = "Set3", name = "Glycan Type") +
  labs(
    title = "N-Glycosylation: Glycan Types Distribution",
    subtitle = paste0("n = ", sum(glycan_counts$Count), " peptides (length 8-14)")
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

ggsave(file.path(output_dir, "Figure6E_Glycan_Types.png"), fig6e,
       width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6E_Glycan_Types.pdf"), fig6e,
       width = 10, height = 6, bg = "white")
cat("  Saved Figure6E_Glycan_Types.png/pdf\n")

# Save glycan table
write.csv(glycan_counts, file.path(output_dir, "data_figure6E_glycan_types.csv"), row.names = FALSE)
cat("  Saved data_figure6E_glycan_types.csv\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=== Figure 6 Summary ===\n")
cat("6A: Doubly phosphorylated site distribution -", nrow(phospho_multi), "peptides\n")
cat("6B: Ubiquitination G vs GG -", nrow(g_data), "G,", nrow(gg_data), "GG\n")
cat("6C: Biological vs Artifact oxidation -", n_bio, "bio,", n_art, "artifact\n")
cat("6D: Biological oxidation residues -", nrow(residue_agg), "residue types\n")
cat("6E: Glycan types -", nrow(glycan_counts), "types\n")
cat("\nAll figures saved to:", output_dir, "\n")
