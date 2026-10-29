# =============================================================================
# Figure 4C (HLA-II): Modified vs Background - PTM site position distribution
# One panel per biological PTM (Artifact Oxidation excluded). Colored curve =
# positions of modified sites (10-25mers); grey curve = position distribution
# of that PTM's target residues across background peptides (10-25mers).
# Run from HLA-II/: cd HLA-II && Rscript scripts/Figure4C_PositionDensity.R
# Output: figure_panels/Figure4C_Position_Distribution.*
# =============================================================================

library(ggplot2)
library(dplyr)

source("../shared/ptm_config.R")

ptm_sites <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)
background <- read.csv("figure_panels/data_background.csv", stringsAsFactors = FALSE)

min_peptide_length <- 10
max_peptide_length <- 25
positions <- seq_len(max_peptide_length)

# The HLA-II background export is already length-filtered to 10-25; the site
# table is filtered here for symmetry with the HLA-I script, and NA sites
# (ambiguous localizations) are dropped.
ptm_sites <- ptm_sites %>%
  filter(Length >= min_peptide_length, Length <= max_peptide_length, !is.na(Site))

# Target residue sets per PTM, taken from the data-loader site-assignment
# conventions. N-term-capable PTMs (e.g. Acetylation): modified sites reported
# as "N-term" carry Site = 1 and stay in the modified curve; the background
# uses the side-chain residue set only, since position 1 is where every
# N-terminus sits by definition.
PTM_TARGET_RESIDUES <- list(
  "Phosphorylation"      = c("S", "T", "Y"),
  "Cysteinylation"       = "C",
  "Deamidation"          = c("N", "Q"),
  "Acetylation"          = c("K", "H", "C", "S", "T", "Y"),
  "Ubiquitination"       = c("K", "S", "T", "C"),
  "Methylation"          = c("C", "H", "N", "Q", "I", "L", "D", "E", "S", "T", "K", "R"),
  "Dimethylation"        = c("K", "R", "N", "P"),
  "Citrullination"       = "R",
  "Oxidation"            = c("P", "I", "L", "Q", "S", "T", "V", "C", "D", "E", "N", "Y", "G", "K", "R"),
  "Carbamidomethylation" = "C",
  "SUMOylation"          = "K",
  "N-Glycosylation"      = "N",
  "Nitrosylation"        = c("C", "Y")
)

count_residues_at_each_position <- function(peptides, n_positions) {
  peptide_chars <- strsplit(peptides, "")
  residue <- unlist(peptide_chars, use.names = FALSE)
  position <- unlist(lapply(peptide_chars, seq_along), use.names = FALSE)
  table(factor(residue), factor(position, levels = seq_len(n_positions)))
}

background_residue_position_counts <-
  count_residues_at_each_position(background$Peptide, max_peptide_length)

background_position_distribution <- function(residue_set) {
  observed_residues <- intersect(residue_set, rownames(background_residue_position_counts))
  position_counts <- colSums(background_residue_position_counts[observed_residues, , drop = FALSE])
  position_counts / sum(position_counts)
}

normalize_counts_to_distribution <- function(counts) {
  if (sum(counts) == 0) return(rep(0, length(counts)))
  counts / sum(counts)
}

panel_data <- bind_rows(lapply(BIO_PTMS, function(ptm_name) {
  modified_counts <- tabulate(ptm_sites$Site[ptm_sites$PTM == ptm_name],
                              nbins = max_peptide_length)
  bind_rows(
    data.frame(ptm = ptm_name, position = positions, type = "Modified",
               density = normalize_counts_to_distribution(modified_counts)),
    data.frame(ptm = ptm_name, position = positions, type = "Background",
               density = as.numeric(background_position_distribution(PTM_TARGET_RESIDUES[[ptm_name]])))
  )
}))
panel_data$ptm <- factor(panel_data$ptm, levels = BIO_PTMS)

modified_n_per_ptm <- vapply(BIO_PTMS, function(ptm_name) {
  sum(ptm_sites$PTM == ptm_name)
}, integer(1))
cat("Modified site n per PTM (HLA-II, lengths 10-25):\n")
print(modified_n_per_ptm)
cat("Background peptide n:", nrow(background), "\n")

panel_annotations <- data.frame(
  ptm = factor(BIO_PTMS, levels = BIO_PTMS),
  label = sprintf("%s n = %d", BIO_PTMS, modified_n_per_ptm)
)
# Background n is shown once, in the first panel.
panel_annotations$label[1] <- paste0(
  panel_annotations$label[1],
  "\nbackground n = ", format(nrow(background), big.mark = ",")
)

figure_4c <- ggplot(panel_data, aes(x = position, y = density)) +
  geom_area(data = filter(panel_data, type == "Background"),
            fill = "gray80", alpha = 0.7) +
  geom_line(data = filter(panel_data, type == "Background"),
            color = "gray40", linewidth = 0.8) +
  geom_area(data = filter(panel_data, type == "Modified"),
            aes(fill = ptm), alpha = 0.5) +
  geom_line(data = filter(panel_data, type == "Modified"),
            aes(color = ptm), linewidth = 1) +
  geom_text(data = panel_annotations, aes(label = label),
            x = Inf, y = Inf, hjust = 1.02, vjust = 1.3,
            size = 3.4, color = "black", lineheight = 1.1) +
  scale_fill_manual(values = PTM_COLORS, guide = "none") +
  scale_color_manual(values = PTM_COLORS, guide = "none") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25),
                     limits = c(0.5, max_peptide_length + 0.5)) +
  facet_wrap(~ ptm, nrow = 4, ncol = 4, scales = "free_y", drop = FALSE) +
  labs(x = "Peptide Position", y = "Density") +
  theme_minimal() +
  theme(
    # Panels are labeled by the in-panel "<PTM> n = X" annotation instead.
    strip.text = element_blank(),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  strip_ms_theme()

save_panel(figure_4c, "Figure4C_Position_Distribution", "posdist")
cat("Saved figure_panels/Figure4C_Position_Distribution.png/pdf/tif\n")
