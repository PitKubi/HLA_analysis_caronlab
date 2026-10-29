# =============================================================================
# FIGURE 4A: Circos Plot - HLA-I PTM Landscape (manuscript Fig 2C, left panel)
# =============================================================================
# Positional distribution of PTMs across HLA-I immunopeptides (8-14 mers,
# positions 1-15). Sector order and colors are the canonical shared set so the
# plot is directly comparable to the HLA-II counterpart; the peptide-count
# legend and residue annotations live in an external panel (see
# shared/circos_draw.R for the rationale).
#
# Input:  figure_panels/data_ptm_sites.csv (from data_loader.R)
# Output: figure_panels/Figure4A_Circos_PTM_Landscape.{png,pdf,tif}
#         figure_panels/Figure4A2_Circos_ArtifactOxidation.{png,pdf,tif}
# =============================================================================

source("shared/ptm_config.R")
source("shared/circos_draw.R")

ptm_data <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)
ptm_data <- ptm_data[ptm_data$Length >= 8 & ptm_data$Length <= 14, ]

# Artifact oxidation dwarfs all other PTMs and gets its own companion circos
artox_data <- ptm_data[ptm_data$PTM == "Artifact Oxidation", ]
bio_data <- ptm_data[ptm_data$PTM != "Artifact Oxidation", ]

positions <- 1:15

ptm_info <- render_ptm_circos(
  ptm_data = bio_data,
  positions = positions,
  ptm_colors = PTM_COLORS,
  ptm_order_canonical = PTM_ORDER,
  out_base = "Figure4A_Circos_PTM_Landscape",
  dims = PAIR_DIMS$circos
)

render_artifact_circos(
  artox_data = artox_data,
  positions = positions,
  out_base = "Figure4A2_Circos_ArtifactOxidation",
  dims = PAIR_DIMS$circos
)

cat("\n=== Figure 4A: PTM Position Distribution Summary ===\n")
for (p in names(ptm_info)) {
  cat(sprintf("%-22s: %5d peptides | Residues: %s\n",
              p, ptm_info[[p]]$n_peptides,
              paste(ptm_info[[p]]$residues, collapse = "/")))
}
