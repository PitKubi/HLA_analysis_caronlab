# =============================================================================
# FIGURE 4A: Circos Plot - HLA-II PTM Landscape (manuscript Fig 2C, right panel)
# =============================================================================
# Positional distribution of PTMs across HLA-II immunopeptides (10-25 mers,
# positions 1-25). Run from the HLA-II/ directory. Sector order, colors and
# layout come from the shared config so this panel is directly comparable to
# the HLA-I counterpart (see shared/circos_draw.R for the legend rationale).
#
# Input:  figure_panels/data_ptm_sites.csv (from data_loader_hla2.R)
# Output: figure_panels/Figure4A_Circos_PTM_Landscape.{png,pdf,tif}
#         figure_panels/Figure4A2_Circos_ArtifactOxidation.{png,pdf,tif}
# =============================================================================

source("../shared/ptm_config.R")
source("../shared/circos_draw.R")

# Keep rows with NA Site: per-PTM n must equal Moh's canonical peptide counts
# (position segments are built from the localized subset only)
ptm_data <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)

# Artifact oxidation dwarfs most other PTMs and gets its own companion circos
artox_data <- ptm_data[ptm_data$PTM == "Artifact Oxidation", ]
bio_data <- ptm_data[ptm_data$PTM != "Artifact Oxidation", ]

positions <- 1:25

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

cat("\n=== Figure 4A (HLA-II): PTM Position Distribution Summary ===\n")
for (p in names(ptm_info)) {
  cat(sprintf("%-22s: %5d peptides | Residues: %s\n",
              p, ptm_info[[p]]$n_peptides,
              paste(ptm_info[[p]]$residues, collapse = "/")))
}
