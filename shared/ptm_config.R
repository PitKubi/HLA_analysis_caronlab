# =============================================================================
# Shared PTM configuration - single source of truth for BOTH pipelines
# (HLA-I at repo root, HLA-II under HLA-II/). Sourced by each data loader and
# written into the respective config.rds; figure scripts take palette/order
# from config.rds only.
#
# Manuscript requirement (Moh, 2026-10-28, points 4/7/9): HLA-I and HLA-II
# panels are mounted side by side, so PTM order, colors, legend sizes and
# panel dimensions must be identical across classes and across figures.
# =============================================================================

# Canonical PTM order, fixed for both classes (combined HLA-I + HLA-II
# abundance). Do NOT reorder per-figure by local counts.
PTM_ORDER <- c(
  "Deamidation",
  "Ubiquitination",
  "Cysteinylation",
  "Oxidation",
  "Methylation",
  "Citrullination",
  "SUMOylation",
  "Phosphorylation",
  "Acetylation",
  "N-Glycosylation",
  "Dimethylation",
  "Carbamidomethylation",
  "Nitrosylation",
  "Artifact Oxidation"
)

# One PTM = one color, everywhere. These are the hexes both pipelines already
# used (identical in the two loaders before unification); the only addition is
# Nitrosylation, taken from the HLA-II pipeline.
PTM_COLORS <- c(
  "Cysteinylation"       = "#2E7D32",
  "Deamidation"          = "#1565C0",
  "Artifact Oxidation"   = "#FF6F00",
  "Oxidation"            = "#0277BD",
  "Phosphorylation"      = "#6A1B9A",
  "Acetylation"          = "#E65100",
  "Ubiquitination"       = "#C62828",
  "Methylation"          = "#558B2F",
  "Dimethylation"        = "#4527A0",
  "Citrullination"       = "#AD1457",
  "Carbamidomethylation" = "#37474F",
  "SUMOylation"          = "#795548",
  "N-Glycosylation"      = "#00695C",
  "Nitrosylation"        = "#00838F"
)

BIO_PTMS <- setdiff(PTM_ORDER, "Artifact Oxidation")

# Per-pair export dimensions in inches (width, height). Both classes of a
# manuscript pair MUST use the same entry so the panels can sit side by side
# at equal size (Moh point 7). Heatmap pairs with different matrix extents
# (P1-P14 vs P1-P25) instead share cell size + fonts - see figure scripts.
PAIR_DIMS <- list(
  donut      = c(11, 9),
  donut_zoom = c(9, 7),
  ridgeline  = c(10, 12),
  elrank     = c(16, 15),
  posdist    = c(14, 12.5),
  circos     = c(20, 16)
)

# Manuscript panels carry no titles/subtitles/captions (Moh point 8);
# figure titles live in the manuscript legend text.
strip_ms_theme <- function() {
  ggplot2::theme(
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.caption = ggplot2::element_blank()
  )
}

#' Save a ggplot panel as png (300 dpi) + pdf + tif (300 dpi LZW) at the
#' canonical pair dimensions, so every I/II pair is pixel-identical in size.
save_panel <- function(plot, name, pair_id, dir = "figure_panels") {
  dims <- PAIR_DIMS[[pair_id]]
  if (is.null(dims)) stop("Unknown pair_id: ", pair_id)
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), plot,
                  width = dims[1], height = dims[2], dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(dir, paste0(name, ".pdf")), plot,
                  width = dims[1], height = dims[2], bg = "white")
  ggplot2::ggsave(file.path(dir, paste0(name, ".tif")), plot,
                  width = dims[1], height = dims[2], dpi = 300, bg = "white",
                  device = grDevices::tiff, compression = "lzw")
  invisible(plot)
}
