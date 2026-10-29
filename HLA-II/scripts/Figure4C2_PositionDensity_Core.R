# =============================================================================
# Figure 4C2 (HLA-II): Modified vs Background - PTM site distribution across
# the NetMHCIIpan binding core (P1-P9). Colored curve = CoreSite positions of
# modified residues that fall inside the predicted 9-mer core; grey curve =
# per-position frequency of that PTM's target residues across background
# cores, normalized to a distribution.
# Run from HLA-II/: cd HLA-II && Rscript scripts/Figure4C2_PositionDensity_Core.R
# Output: figure_panels/Figure4C2_Position_Distribution_CoreP1P9.*
# =============================================================================

library(ggplot2)
library(dplyr)

source("../shared/ptm_config.R")

ptm_sites <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)
background <- read.csv("figure_panels/data_background.csv", stringsAsFactors = FALSE)

core_length <- 9
core_positions <- seq_len(core_length)

ptm_sites <- ptm_sites %>% filter(!is.na(Site))

# Target residue sets per PTM, taken from the data-loader site-assignment
# conventions (same sets as the full-length Figure 4C scripts).
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

count_residues_at_each_position <- function(sequences, n_positions) {
  sequence_chars <- strsplit(sequences, "")
  residue <- unlist(sequence_chars, use.names = FALSE)
  position <- unlist(lapply(sequence_chars, seq_along), use.names = FALSE)
  table(factor(residue), factor(position, levels = seq_len(n_positions)))
}

# Background cores are all exactly 9-mers, so residue-at-position counts over
# the cores are proportional to residue-at-position frequencies.
background_core_residue_position_counts <-
  count_residues_at_each_position(background$Core, core_length)

background_core_position_distribution <- function(residue_set) {
  observed_residues <- intersect(residue_set, rownames(background_core_residue_position_counts))
  position_counts <- colSums(background_core_residue_position_counts[observed_residues, , drop = FALSE])
  position_counts / sum(position_counts)
}

normalize_counts_to_distribution <- function(counts) {
  if (sum(counts) == 0) return(rep(0, length(counts)))
  counts / sum(counts)
}

panel_data <- bind_rows(lapply(BIO_PTMS, function(ptm_name) {
  in_core_sites <- ptm_sites %>% filter(PTM == ptm_name, !is.na(CoreSite))
  modified_counts <- tabulate(in_core_sites$CoreSite, nbins = core_length)
  bind_rows(
    data.frame(ptm = ptm_name, position = core_positions, type = "Modified",
               density = normalize_counts_to_distribution(modified_counts)),
    data.frame(ptm = ptm_name, position = core_positions, type = "Background",
               density = as.numeric(background_core_position_distribution(PTM_TARGET_RESIDUES[[ptm_name]])))
  )
}))
panel_data$ptm <- factor(panel_data$ptm, levels = BIO_PTMS)

in_core_n_per_ptm <- vapply(BIO_PTMS, function(ptm_name) {
  sum(ptm_sites$PTM == ptm_name & !is.na(ptm_sites$CoreSite))
}, integer(1))
total_n_per_ptm <- vapply(BIO_PTMS, function(ptm_name) {
  sum(ptm_sites$PTM == ptm_name)
}, integer(1))
cat("Modified site n per PTM (HLA-II, in core / total):\n")
print(data.frame(PTM = BIO_PTMS, in_core = in_core_n_per_ptm,
                 total = total_n_per_ptm, row.names = NULL))
cat("Background core n:", nrow(background), "\n")

panel_annotations <- data.frame(
  ptm = factor(BIO_PTMS, levels = BIO_PTMS),
  label = sprintf("%s\nn = %d/%d in core",
                  BIO_PTMS, in_core_n_per_ptm, total_n_per_ptm)
)
# Background n is shown once, in the first panel.
panel_annotations$label[1] <- paste0(
  panel_annotations$label[1],
  "\nbackground n = ", format(nrow(background), big.mark = ",")
)

figure_4c2 <- ggplot(panel_data, aes(x = position, y = density)) +
  geom_area(data = filter(panel_data, type == "Background"),
            fill = "gray80", alpha = 0.7) +
  geom_line(data = filter(panel_data, type == "Background"),
            color = "gray40", linewidth = 0.8) +
  geom_area(data = filter(panel_data, type == "Modified"),
            aes(fill = ptm), alpha = 0.5) +
  geom_line(data = filter(panel_data, type == "Modified"),
            aes(color = ptm), linewidth = 1) +
  geom_text(data = panel_annotations, aes(label = label),
            x = Inf, y = Inf, hjust = 1.02, vjust = 1.2,
            size = 3.4, color = "black", lineheight = 1.1) +
  scale_fill_manual(values = PTM_COLORS, guide = "none") +
  scale_color_manual(values = PTM_COLORS, guide = "none") +
  scale_x_continuous(breaks = core_positions,
                     labels = paste0("P", core_positions),
                     limits = c(0.5, core_length + 0.5)) +
  facet_wrap(~ ptm, nrow = 4, ncol = 4, scales = "free_y", drop = FALSE) +
  labs(x = "Core Position", y = "Density") +
  theme_minimal() +
  theme(
    # Panels are labeled by the in-panel "<PTM> n = X/Y in core" annotation.
    strip.text = element_blank(),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  strip_ms_theme()

save_panel(figure_4c2, "Figure4C2_Position_Distribution_CoreP1P9", "posdist")
cat("Saved figure_panels/Figure4C2_Position_Distribution_CoreP1P9.png/pdf/tif\n")
