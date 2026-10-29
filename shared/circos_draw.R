# =============================================================================
# Shared circos renderer for the PTM-landscape figure (manuscript Fig 2C pair)
# Used by both pipelines: Figure4A_Circos.R (HLA-I) and
# HLA-II/scripts/Figure4A_Circos.R (HLA-II).
#
# Layout change requested by Moh (2026-10-28, point 2 + annotated circos):
# the center peptide-count legend, the tiny Position-scale numerals and the
# inner residue-ring letters were unreadable at manuscript size. The legend
# now lives in a dedicated panel right of the circle (with per-PTM counts and
# modified residues), position numerals are enlarged, and in-ring residue
# letters are drawn only where the sector is wide enough - the legend always
# carries them.
# =============================================================================

library(circlize)

#' Dominant modified residues (>10% frequency) for one PTM, from site data
dominant_residues <- function(ptm_data, ptm_name) {
  res_table <- sort(table(ptm_data$Residue[ptm_data$PTM == ptm_name]), decreasing = TRUE)
  if (length(res_table) == 0) return(character(0))
  res_freq <- res_table / sum(res_table)
  dominant <- names(res_freq[res_freq > 0.1])
  if (length(dominant) == 0) dominant <- names(res_table)[1]
  dominant
}

#' Per-PTM position counts and residue labels, in canonical order
build_ptm_info <- function(ptm_data, positions, ptm_order_canonical) {
  present <- ptm_order_canonical[ptm_order_canonical %in% unique(ptm_data$PTM)]
  info <- lapply(present, function(ptm) {
    subset <- ptm_data[ptm_data$PTM == ptm, ]
    pos_counts <- vapply(positions, function(p) sum(subset$Site == p, na.rm = TRUE), numeric(1))
    names(pos_counts) <- as.character(positions)
    list(
      n_peptides = nrow(subset),
      pos_counts = pos_counts,
      residues = dominant_residues(ptm_data, ptm)
    )
  })
  names(info) <- present
  info
}

# Short sector labels for narrow sectors (full names live in the legend)
SHORT_LABELS <- c(
  "Methylation" = "Methyl.", "Citrullination" = "Citrul.",
  "Dimethylation" = "DiMe.", "N-Glycosylation" = "N-Glyco",
  "Carbamidomethylation" = "Carbam", "SUMOylation" = "SUMO",
  "Nitrosylation" = "Nitro", "Acetylation" = "Acetyl.",
  "Phosphorylation" = "Phospho."
)

draw_circle <- function(ptm_info, ptm_colors, positions, inner_ring = TRUE) {
  ptm_order_present <- names(ptm_info)
  n_positions <- length(positions)
  all_sectors <- c("Position", ptm_order_present)
  ptm_widths <- vapply(ptm_info, function(x) x$n_peptides, numeric(1))
  # Position wedge scales with the data so its header stays readable for any
  # total count (fixed 400 units shrank to a sliver on the HLA-II totals)
  sector_widths <- c(round(0.055 * sum(ptm_widths)), ptm_widths)
  subring_height <- 0.50 / n_positions
  # With >16 position subrings the numerals collide; label every other one
  pos_label_step <- if (n_positions > 16) 2L else 1L
  pos_label_cex <- if (n_positions > 16) 1.1 else 1.3

  # Stagger the outside labels of consecutive narrow sectors between two radii
  # so they cannot collide (the small PTMs cluster together in canonical order)
  n_gaps <- length(all_sectors)
  approx_deg <- sector_widths / sum(sector_widths) * (360 - 4 - 4 - 2 * (n_gaps - 2))
  names(approx_deg) <- all_sectors
  label_y <- setNames(rep(1.5, length(ptm_order_present)), ptm_order_present)
  toggle_far <- FALSE
  for (s in ptm_order_present) {
    if (approx_deg[[s]] < 8) {
      label_y[[s]] <- if (toggle_far) 2.3 else 1.5
      toggle_far <- !toggle_far
    } else {
      toggle_far <- FALSE
    }
  }

  circos.clear()
  circos.par(
    gap.after = c(4, rep(2, length(ptm_order_present) - 1), 4),
    start.degree = 90,
    track.margin = c(0.001, 0.001),
    cell.padding = c(0, 0, 0, 0),
    canvas.xlim = c(-1.15, 1.15),
    canvas.ylim = c(-1.15, 1.15)
  )
  circos.initialize(
    factors = factor(all_sectors, levels = all_sectors),
    xlim = cbind(rep(0, length(all_sectors)), sector_widths)
  )

  # Outer ring: PTM sectors with counts, names outside
  circos.track(ylim = c(0, 1), track.height = 0.10, bg.border = "white",
    panel.fun = function(x, y) {
      sector.name <- CELL_META$sector.index
      if (sector.name == "Position") {
        circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                    col = "grey60", border = "white")
        circos.text(CELL_META$xcenter, 0.5, "Position",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 1.4, col = "white", font = 2)
      } else {
        circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                    col = ptm_colors[sector.name], border = "white")
        sector_width <- abs(CELL_META$cell.end.degree - CELL_META$cell.start.degree)
        if (sector_width > 8) {
          circos.text(CELL_META$xcenter, 0.5,
                      format(ptm_info[[sector.name]]$n_peptides, big.mark = ","),
                      facing = "bending.inside", niceFacing = TRUE,
                      cex = 0.9, col = "white", font = 2)
        }
        wide <- sector_width > 14
        label <- if (wide || !sector.name %in% names(SHORT_LABELS)) {
          sector.name
        } else {
          SHORT_LABELS[[sector.name]]
        }
        circos.text(CELL_META$xcenter, label_y[[sector.name]], label,
                    facing = "bending.outside", niceFacing = TRUE,
                    cex = if (wide) 1.3 else 1.1,
                    col = ptm_colors[sector.name], font = 2)
      }
    })

  # Middle rings: one subring per position; angular segment = count at position
  for (pos_idx in rev(seq_along(positions))) {
    pos <- positions[pos_idx]
    circos.track(ylim = c(0, 1), track.height = subring_height,
      bg.col = "grey93", bg.border = "grey88",
      panel.fun = function(x, y) {
        sector.name <- CELL_META$sector.index
        if (sector.name == "Position") {
          circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                      col = "grey85", border = "grey80")
          if ((pos - 1) %% pos_label_step == 0) {
            circos.text(CELL_META$xcenter, 0.5, as.character(pos),
                        facing = "bending.inside", niceFacing = TRUE,
                        cex = pos_label_cex, col = "gray20", font = 2)
          }
        } else {
          pos_counts <- ptm_info[[sector.name]]$pos_counts
          cumsum_counts <- cumsum(pos_counts)
          count_at_pos <- pos_counts[as.character(pos)]
          if (count_at_pos > 0) {
            x_start <- if (pos == positions[1]) 0 else cumsum_counts[as.character(pos - 1)]
            circos.rect(x_start, 0, cumsum_counts[as.character(pos)], 1,
                        col = ptm_colors[sector.name], border = NA)
          }
        }
      })
  }

  if (!inner_ring) {
    circos.clear()
    return(invisible(NULL))
  }

  # Inner ring: dominant residues, only where the sector is wide enough to
  # read them; the external legend carries them for every PTM
  circos.track(ylim = c(0, 1), track.height = 0.12,
    bg.col = "white", bg.border = "gray60",
    panel.fun = function(x, y) {
      sector.name <- CELL_META$sector.index
      if (sector.name == "Position") {
        circos.text(CELL_META$xcenter, 0.5, "Residue",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 1.2, col = "gray40", font = 2)
      } else {
        sector_width <- abs(CELL_META$cell.end.degree - CELL_META$cell.start.degree)
        if (sector_width >= 10) {
          residues <- ptm_info[[sector.name]]$residues
          one_line <- paste(residues, collapse = "/")
          lines <- if (length(residues) > 2 || (length(residues) > 1 && nchar(one_line) > 5)) {
            c(residues[1], paste(residues[-1], collapse = "/"))
          } else {
            one_line
          }
          y_pos <- if (length(lines) == 2) c(0.68, 0.28) else 0.5
          for (k in seq_along(lines)) {
            circos.text(CELL_META$xcenter, y_pos[k], lines[k],
                        facing = "bending.inside", niceFacing = TRUE,
                        cex = 1.4, col = ptm_colors[sector.name], font = 2)
          }
        }
      }
    })
  circos.clear()
}

draw_legend <- function(ptm_info, ptm_colors) {
  op <- par(mar = c(0, 0, 0, 0))
  on.exit(par(op))
  plot.new()
  entries <- vapply(names(ptm_info), function(ptm) {
    n_str <- format(ptm_info[[ptm]]$n_peptides, big.mark = ",")
    residues <- ptm_info[[ptm]]$residues
    if (length(residues) > 0) {
      sprintf("%s (%s) @ %s", ptm, n_str, paste(residues, collapse = "/"))
    } else {
      sprintf("%s (%s)", ptm, n_str)
    }
  }, character(1))
  has_residues <- any(vapply(ptm_info, function(x) length(x$residues) > 0, logical(1)))
  legend("left", legend = entries, fill = ptm_colors[names(ptm_info)],
         border = NA, bty = "n", cex = 1.45, y.intersp = 1.35,
         title = if (has_residues) "PTM (n) @ dominant residues" else "Residue (n)",
         title.col = "black", title.cex = 1.5, xpd = TRUE)
}

render_circos_devices <- function(out_base, dims, out_dir, draw_fun) {
  devices <- list(
    png = function(f) png(f, width = dims[1], height = dims[2], units = "in", res = 300),
    pdf = function(f) pdf(f, width = dims[1], height = dims[2]),
    tif = function(f) tiff(f, width = dims[1], height = dims[2], units = "in",
                           res = 300, compression = "lzw")
  )
  for (ext in names(devices)) {
    f <- file.path(out_dir, paste0(out_base, ".", ext))
    devices[[ext]](f)
    layout(matrix(1:2, nrow = 1), widths = c(3, 1.15))
    draw_fun()
    layout(1)
    dev.off()
    cat("Saved", f, "\n")
  }
}

#' Render the PTM circos + external legend to png, pdf and tif
#' dims in inches; legend occupies the right quarter of the canvas
render_ptm_circos <- function(ptm_data, positions, ptm_colors,
                              ptm_order_canonical, out_base, dims,
                              out_dir = "figure_panels") {
  ptm_info <- build_ptm_info(ptm_data, positions, ptm_order_canonical)
  render_circos_devices(out_base, dims, out_dir, function() {
    draw_circle(ptm_info, ptm_colors, positions)
    draw_legend(ptm_info, ptm_colors)
  })
  invisible(ptm_info)
}

#' Artifact-oxidation companion circos: sectors are the oxidized residues
#' (M/W/H/F), no inner residue ring, external legend like the main plot
render_artifact_circos <- function(artox_data, positions, out_base, dims,
                                   out_dir = "figure_panels") {
  artox_colors <- c("M" = "#FF6F00", "W" = "#FFB300", "H" = "#BF360C", "F" = "#E65100")
  residue_order <- names(sort(table(artox_data$Residue), decreasing = TRUE))
  info <- lapply(residue_order, function(res) {
    subset <- artox_data[artox_data$Residue == res, ]
    pos_counts <- vapply(positions, function(p) sum(subset$Site == p, na.rm = TRUE), numeric(1))
    names(pos_counts) <- as.character(positions)
    list(n_peptides = nrow(subset), pos_counts = pos_counts, residues = character(0))
  })
  names(info) <- residue_order
  render_circos_devices(out_base, dims, out_dir, function() {
    draw_circle(info, artox_colors, positions, inner_ring = FALSE)
    draw_legend(info, artox_colors)
  })
  invisible(info)
}
