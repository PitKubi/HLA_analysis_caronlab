# =============================================================================
# FIGURE 4A: Circos Plot - HLA-I PTM Landscape
# =============================================================================
#
# Description:
#   This circos plot visualizes the positional distribution of post-translational
#   modifications (PTMs) across HLA-I immunopeptides. The plot reveals where each
#   PTM type preferentially occurs along the peptide sequence (positions 1-15,
#   N-terminus to C-terminus).
#
# Structure (from outside to inside):
#   - RING 3 (Outer): PTM type names (colored) with peptide counts inside
#   - RING 2 (Middle): 15 concentric subrings representing peptide positions
#                      Angular width of colored segments = number of peptides
#                      modified at that position (creates diagonal clustering)
#   - RING 1 (Inner): Dominant modified amino acid residues for each PTM
#   - CENTER: Color-coded legend showing peptide counts per PTM type
#
# Key Features:
#   - Diagonal pattern shows positional preferences for each PTM
#   - Position 1 (N-terminus) is innermost, Position 15 is outermost
#   - Sector angular width proportional to total peptide count per PTM
#   - Data filtered to 8-14 mer peptides
#
# Input: figure_panels/data_ptm_sites.csv (from data_loader.R)
# Output: figure_panels/Figure4A_Circos_PTM_Landscape.png and .pdf
#
# Author: Generated for HLA-I PTM stability analysis
# =============================================================================

library(circlize)

# =============================================================================
# DATA PREPARATION
# =============================================================================

# Load PTM data (from data_loader.R)
ptm_data <- read.csv("figure_panels/data_ptm_sites.csv", stringsAsFactors = FALSE)

# Filter to 8-14 mers (raw data now includes all lengths from unfiltered columns)
ptm_data <- ptm_data[ptm_data$Length >= 8 & ptm_data$Length <= 14, ]

# Separate artifact oxidation for its own circos (it dwarfs all other PTMs)
artox_data <- ptm_data[ptm_data$PTM == "Artifact Oxidation", ]
ptm_data <- ptm_data[ptm_data$PTM != "Artifact Oxidation", ]

# Order PTMs by abundance (descending)
ptm_counts <- table(ptm_data$PTM)
ptm_order <- names(sort(ptm_counts, decreasing = TRUE))

# Load colors from config
config <- readRDS("figure_panels/config.rds")
ptm_colors <- config$ptm_colors

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Get dominant modified residues for a PTM type
#' @param ptm_name Name of the PTM
#' @return Character vector of residues with >10% frequency
get_dominant_residues <- function(ptm_name) {
  ptm_subset <- ptm_data[ptm_data$PTM == ptm_name, ]
  res_table <- sort(table(ptm_subset$Residue), decreasing = TRUE)
  res_freq <- res_table / sum(res_table)
  dominant <- names(res_freq[res_freq > 0.1])
  if (length(dominant) == 0) dominant <- names(res_table)[1]
  return(dominant)
}

# =============================================================================
# CALCULATE POSITION STATISTICS
# =============================================================================

# Define position range (standard HLA-I peptide positions)
max_pos <- 15
positions <- 1:max_pos
n_positions <- length(positions)

# Calculate statistics for each PTM
ptm_info <- list()
for (ptm in ptm_order) {
  ptm_subset <- ptm_data[ptm_data$PTM == ptm, ]

  # Count modifications at each position
  pos_counts <- sapply(positions, function(p) sum(ptm_subset$Site == p))
  names(pos_counts) <- as.character(positions)

  ptm_info[[ptm]] <- list(
    n_peptides = nrow(ptm_subset),
    pos_counts = pos_counts,
    residues = get_dominant_residues(ptm)
  )
}

# =============================================================================
# N-TERMINAL ACETYLATION: Compute dominant first amino acids
# =============================================================================

nterm_acet <- ptm_data[ptm_data$PTM == "Acetylation" & ptm_data$Residue == "N-term", ]
if (nrow(nterm_acet) > 0) {
  first_aa <- substr(nterm_acet$Peptide, 1, 1)
  aa_table <- sort(table(first_aa), decreasing = TRUE)
  aa_freq <- aa_table / sum(aa_table)
  dominant_nterm_aa <- names(aa_freq[aa_freq > 0.1])
  nterm_label <- paste0("Nt-", paste(dominant_nterm_aa, collapse = "/"))
} else {
  nterm_label <- "N-term"
}

cat("Acetylation N-term dominant residues:", nterm_label, "\n")

# =============================================================================
# CIRCOS PLOT SETUP
# =============================================================================

# Define sectors: Position reference + all PTM types
all_sectors <- c("Position", ptm_order)

# Sector widths: fixed width for Position, proportional to peptide count for PTMs
sector_widths <- c(400, sapply(ptm_order, function(p) ptm_info[[p]]$n_peptides))

# Height of each position subring
subring_height <- 0.50 / n_positions

# Label configuration: stagger small adjacent sectors to avoid overlap
# Large PTMs get full name at standard distance; small ones alternate near/far
label_configs <- list(
  "Cysteinylation"  = list(text = "Cysteinylation",  y = 1.5, cex = 1.3),
  "Deamidation"     = list(text = "Deamidation",     y = 1.5, cex = 1.3),
  "Oxidation"       = list(text = "  Oxidation",     y = 1.5, cex = 1.3),
  "Phosphorylation" = list(text = "Phosphorylation", y = 1.5, cex = 1.3),
  "Ubiquitination"  = list(text = "Ubiquitination",  y = 1.5, cex = 1.3),
  "Acetylation"     = list(text = "Acetylation",     y = 1.5, cex = 1.3),
  "Methylation"     = list(text = "Methyl.",          y = 1.5, cex = 1.1),
  "Citrullination"  = list(text = "Citrul.",          y = 1.5, cex = 1.1),
  "Dimethylation"   = list(text = "DiMe.",            y = 1.5, cex = 1.1),
  "N-Glycosylation"      = list(text = "N-Glyco",    y = 1.5, cex = 1.0),
  "Carbamidomethylation" = list(text = "Carbam",      y = 2.0, cex = 1.0),
  "SUMOylation"          = list(text = "SUMO",        y = 1.5, cex = 1.0)
)

# =============================================================================
# GENERATE PNG OUTPUT
# =============================================================================

png("figure_panels/Figure4A_Circos_PTM_Landscape.png",
    width = 16, height = 16, units = "in", res = 300)

circos.clear()
circos.par(
  gap.after = c(4, rep(2, length(ptm_order) - 1), 4),
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

# -----------------------------------------------------------------------------
# RING 3 (OUTER): PTM names and peptide counts
# -----------------------------------------------------------------------------
circos.track(ylim = c(0, 1), track.height = 0.10, bg.border = "white",
  panel.fun = function(x, y) {
    sector.name <- CELL_META$sector.index

    if (sector.name == "Position") {
      # Position reference sector (grey)
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col = "grey60", border = "white")
      circos.text(CELL_META$xcenter, 0.5, "Position",
                  facing = "bending.inside", niceFacing = TRUE,
                  cex = 1.2, col = "white", font = 2)
    } else {
      # PTM sectors (colored)
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col = ptm_colors[sector.name], border = "white")

      n_count <- ptm_info[[sector.name]]$n_peptides
      sector_width <- CELL_META$cell.end.degree - CELL_META$cell.start.degree

      # Peptide count inside ring (white text)
      if (sector_width > 8) {
        circos.text(CELL_META$xcenter, 0.5, n_count,
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 0.75, col = "white", font = 2)
      }

      # PTM name outside ring
      label_cfg <- label_configs[[sector.name]]
      if (is.null(label_cfg)) label_cfg <- list(text = sector.name, y = 1.5, cex = 1.3)
      circos.text(CELL_META$xcenter, label_cfg$y, label_cfg$text,
                  facing = "bending.outside", niceFacing = TRUE,
                  cex = label_cfg$cex, col = ptm_colors[sector.name], font = 2)
    }
  })

# -----------------------------------------------------------------------------
# RING 2 (MIDDLE): Position subrings showing modification site distribution
# Creates diagonal pattern: angular width = peptide count at each position
# -----------------------------------------------------------------------------
for (pos_idx in rev(seq_along(positions))) {
  pos <- positions[pos_idx]

  circos.track(ylim = c(0, 1), track.height = subring_height,
    bg.col = "grey93", bg.border = "grey88",
    panel.fun = function(x, y) {
      sector.name <- CELL_META$sector.index

      if (sector.name == "Position") {
        # Position number labels
        circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                    col = "grey85", border = "grey80")
        circos.text(CELL_META$xcenter, 0.5, as.character(pos),
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 0.85, col = "gray20", font = 2)
      } else {
        # PTM modification data
        info <- ptm_info[[sector.name]]
        pos_counts <- info$pos_counts
        cumsum_counts <- cumsum(pos_counts)
        count_at_pos <- pos_counts[as.character(pos)]

        # Draw colored segment if modifications exist at this position
        if (count_at_pos > 0) {
          x_start <- if (pos == 1) 0 else cumsum_counts[as.character(pos - 1)]
          x_end <- cumsum_counts[as.character(pos)]
          circos.rect(x_start, 0, x_end, 1,
                      col = ptm_colors[sector.name], border = NA)
        }
      }
    })
}

# -----------------------------------------------------------------------------
# RING 1 (INNER): Dominant modified amino acid residues
# -----------------------------------------------------------------------------
circos.track(ylim = c(0, 1), track.height = 0.12,
  bg.col = "white", bg.border = "gray60",
  panel.fun = function(x, y) {
    sector.name <- CELL_META$sector.index

    if (sector.name == "Position") {
      circos.text(CELL_META$xcenter, 0.5, "Residue",
                  facing = "bending.inside", niceFacing = TRUE,
                  cex = 0.9, col = "gray50", font = 2)
    } else {
      info <- ptm_info[[sector.name]]
      res_text <- paste(info$residues, collapse = "/")
      sector_width <- CELL_META$cell.end.degree - CELL_META$cell.start.degree

      # Handle text placement for small sectors
      if (sector.name == "Acetylation") {
        # Show both N-term and K residues stacked
        circos.text(CELL_META$xcenter, 0.65, "N-term",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 1.0, col = ptm_colors[sector.name], font = 2)
        circos.text(CELL_META$xcenter, 0.25, "K",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 1.0, col = ptm_colors[sector.name], font = 2)
      } else if (sector.name == "Methylation") {
        # Split S/L/D across two lines
        circos.text(CELL_META$xcenter, 0.65, "S/L",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 1.0, col = ptm_colors[sector.name], font = 2)
        circos.text(CELL_META$xcenter, 0.25, "D",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 1.0, col = ptm_colors[sector.name], font = 2)
      } else if (sector.name == "Oxidation") {
        # Show top oxidation residues (M, L, S, P), shift right to avoid overlap
        x_offset <- (CELL_META$xlim[2] - CELL_META$xlim[1]) * 0.3
        circos.text(CELL_META$xcenter + x_offset, 0.65, "M/L",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 1.0, col = ptm_colors[sector.name], font = 2)
        circos.text(CELL_META$xcenter + x_offset, 0.25, "S/P",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 1.0, col = ptm_colors[sector.name], font = 2)
      } else if (sector.name == "Acetylation") {
        # Show K and dominant N-terminal residues
        circos.text(CELL_META$xcenter, 0.7, "K",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 1.0, col = ptm_colors[sector.name], font = 2)
        circos.text(CELL_META$xcenter, 0.3, nterm_label,
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 0.8, col = ptm_colors[sector.name], font = 2)
      } else if (sector.name == "Dimethylation") {
        # Stack R, P, M on separate lines
        circos.text(CELL_META$xcenter, 0.78, "R",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 0.9, col = ptm_colors[sector.name], font = 2)
        circos.text(CELL_META$xcenter, 0.48, "P",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 0.9, col = ptm_colors[sector.name], font = 2)
        circos.text(CELL_META$xcenter, 0.18, "M",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 0.9, col = ptm_colors[sector.name], font = 2)
      } else if (sector_width < 20) {
        circos.text(CELL_META$xcenter, 0.5, res_text,
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 1.0, col = ptm_colors[sector.name], font = 2)
      } else {
        circos.text(CELL_META$xcenter, 0.5, res_text,
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 1.2, col = ptm_colors[sector.name], font = 2)
      }
    }
  })

# -----------------------------------------------------------------------------
# CENTER: Color-coded legend with peptide counts
# -----------------------------------------------------------------------------
n_ptms <- length(ptm_order)
y_positions <- seq(0.18, -0.18, length.out = n_ptms)

for (i in seq_along(ptm_order)) {
  p <- ptm_order[i]
  text(0, y_positions[i], paste0(p, ": ", ptm_info[[p]]$n_peptides),
       cex = 1.0, col = ptm_colors[p], font = 2)
}

dev.off()
circos.clear()

cat("Saved Figure4A_Circos_PTM_Landscape.png\n")

# =============================================================================
# GENERATE PDF OUTPUT (vector format for publication)
# =============================================================================

pdf("figure_panels/Figure4A_Circos_PTM_Landscape.pdf", width = 16, height = 16)

circos.clear()
circos.par(
  gap.after = c(4, rep(2, length(ptm_order) - 1), 4),
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

# Ring 3: PTM names and counts
circos.track(ylim = c(0, 1), track.height = 0.10, bg.border = "white",
  panel.fun = function(x, y) {
    sector.name <- CELL_META$sector.index
    if (sector.name == "Position") {
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = "grey60", border = "white")
      circos.text(CELL_META$xcenter, 0.5, "Position", facing = "bending.inside", niceFacing = TRUE, cex = 1.2, col = "white", font = 2)
    } else {
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = ptm_colors[sector.name], border = "white")
      n_count <- ptm_info[[sector.name]]$n_peptides
      sector_width <- CELL_META$cell.end.degree - CELL_META$cell.start.degree
      if (sector_width > 8) {
        circos.text(CELL_META$xcenter, 0.5, n_count, facing = "bending.inside", niceFacing = TRUE, cex = 0.75, col = "white", font = 2)
      }
      # PTM name outside ring
      label_cfg <- label_configs[[sector.name]]
      if (is.null(label_cfg)) label_cfg <- list(text = sector.name, y = 1.5, cex = 1.3)
      circos.text(CELL_META$xcenter, label_cfg$y, label_cfg$text, facing = "bending.outside", niceFacing = TRUE, cex = label_cfg$cex, col = ptm_colors[sector.name], font = 2)
      # Show residues below PTM name for small sectors
      if (isTRUE(label_cfg$show_res)) {
        res_text <- paste(ptm_info[[sector.name]]$residues, collapse = "/")
        circos.text(CELL_META$xcenter, label_cfg$res_y, res_text, facing = "bending.outside", niceFacing = TRUE, cex = 0.9, col = ptm_colors[sector.name], font = 1)
      }
    }
  })

# Ring 2: Position subrings
for (pos_idx in rev(seq_along(positions))) {
  pos <- positions[pos_idx]
  circos.track(ylim = c(0, 1), track.height = subring_height, bg.col = "grey93", bg.border = "grey88",
    panel.fun = function(x, y) {
      sector.name <- CELL_META$sector.index
      if (sector.name == "Position") {
        circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = "grey85", border = "grey80")
        circos.text(CELL_META$xcenter, 0.5, as.character(pos), facing = "bending.inside", niceFacing = TRUE, cex = 0.85, col = "gray20", font = 2)
      } else {
        info <- ptm_info[[sector.name]]
        pos_counts <- info$pos_counts
        cumsum_counts <- cumsum(pos_counts)
        count_at_pos <- pos_counts[as.character(pos)]
        if (count_at_pos > 0) {
          x_start <- if (pos == 1) 0 else cumsum_counts[as.character(pos - 1)]
          x_end <- cumsum_counts[as.character(pos)]
          circos.rect(x_start, 0, x_end, 1, col = ptm_colors[sector.name], border = NA)
        }
      }
    })
}

# Ring 1: Modified residues
circos.track(ylim = c(0, 1), track.height = 0.12, bg.col = "white", bg.border = "gray60",
  panel.fun = function(x, y) {
    sector.name <- CELL_META$sector.index
    if (sector.name == "Position") {
      circos.text(CELL_META$xcenter, 0.5, "Residue", facing = "bending.inside", niceFacing = TRUE, cex = 0.9, col = "gray50", font = 2)
    } else {
      res_text <- paste(ptm_info[[sector.name]]$residues, collapse = "/")
      sector_width <- CELL_META$cell.end.degree - CELL_META$cell.start.degree
      if (sector.name == "Acetylation") {
        circos.text(CELL_META$xcenter, 0.7, "K", facing = "bending.inside", niceFacing = TRUE, cex = 1.0, col = ptm_colors[sector.name], font = 2)
        circos.text(CELL_META$xcenter, 0.3, nterm_label, facing = "bending.inside", niceFacing = TRUE, cex = 0.8, col = ptm_colors[sector.name], font = 2)
      } else if (sector.name == "Methylation") {
        circos.text(CELL_META$xcenter, 0.65, "S/L", facing = "bending.inside", niceFacing = TRUE, cex = 1.0, col = ptm_colors[sector.name], font = 2)
        circos.text(CELL_META$xcenter, 0.25, "D", facing = "bending.inside", niceFacing = TRUE, cex = 1.0, col = ptm_colors[sector.name], font = 2)
      } else if (sector.name == "Oxidation") {
        x_offset <- (CELL_META$xlim[2] - CELL_META$xlim[1]) * 0.3
        circos.text(CELL_META$xcenter + x_offset, 0.65, "M/L", facing = "bending.inside", niceFacing = TRUE, cex = 1.0, col = ptm_colors[sector.name], font = 2)
        circos.text(CELL_META$xcenter + x_offset, 0.25, "S/P", facing = "bending.inside", niceFacing = TRUE, cex = 1.0, col = ptm_colors[sector.name], font = 2)
      } else if (sector.name == "Dimethylation") {
        circos.text(CELL_META$xcenter, 0.78, "R", facing = "bending.inside", niceFacing = TRUE, cex = 0.9, col = ptm_colors[sector.name], font = 2)
        circos.text(CELL_META$xcenter, 0.48, "P", facing = "bending.inside", niceFacing = TRUE, cex = 0.9, col = ptm_colors[sector.name], font = 2)
        circos.text(CELL_META$xcenter, 0.18, "M", facing = "bending.inside", niceFacing = TRUE, cex = 0.9, col = ptm_colors[sector.name], font = 2)
      } else if (sector_width < 20) {
        circos.text(CELL_META$xcenter, 0.5, res_text, facing = "bending.inside", niceFacing = TRUE, cex = 1.0, col = ptm_colors[sector.name], font = 2)
      } else {
        circos.text(CELL_META$xcenter, 0.5, res_text, facing = "bending.inside", niceFacing = TRUE, cex = 1.2, col = ptm_colors[sector.name], font = 2)
      }
    }
  })

# Center legend
for (i in seq_along(ptm_order)) {
  p <- ptm_order[i]
  text(0, y_positions[i], paste0(p, ": ", ptm_info[[p]]$n_peptides),
       cex = 1.0, col = ptm_colors[p], font = 2)
}

dev.off()
circos.clear()

cat("Saved Figure4A_Circos_PTM_Landscape.pdf\n")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
cat("\n=== Figure 4A: PTM Position Distribution Summary ===\n")
cat("Total PTM types:", length(ptm_order), "\n")
cat("Total modified peptides:", sum(sapply(ptm_info, function(x) x$n_peptides)), "\n\n")

for (p in ptm_order) {
  cat(sprintf("%-15s: %4d peptides | Residues: %s\n",
              p, ptm_info[[p]]$n_peptides,
              paste(ptm_info[[p]]$residues, collapse = "/")))
}

# =============================================================================
# ARTIFACT OXIDATION: Separate Circos Plot
# =============================================================================
cat("\n=== Generating Artifact Oxidation Circos ===\n")

if (nrow(artox_data) > 0) {
  # Compute position stats per residue
  artox_residues <- sort(unique(artox_data$Residue))
  artox_info <- list()
  for (res in artox_residues) {
    res_subset <- artox_data[artox_data$Residue == res, ]
    pos_counts <- sapply(positions, function(p) sum(res_subset$Site == p))
    names(pos_counts) <- as.character(positions)
    artox_info[[res]] <- list(n_peptides = nrow(res_subset), pos_counts = pos_counts)
  }

  # Sort by abundance
  artox_order <- names(sort(sapply(artox_info, function(x) x$n_peptides), decreasing = TRUE))

  artox_colors <- c("M" = "#FF6F00", "W" = "#FFB300", "H" = "#BF360C", "F" = "#E65100")

  # --- PNG ---
  png("figure_panels/Figure4A2_Circos_ArtifactOxidation.png",
      width = 14, height = 14, units = "in", res = 300)

  artox_sectors <- c("Position", artox_order)
  artox_widths <- c(300, sapply(artox_order, function(r) artox_info[[r]]$n_peptides))

  circos.clear()
  circos.par(
    gap.after = c(4, rep(3, length(artox_order) - 1), 4),
    start.degree = 90,
    track.margin = c(0.001, 0.001),
    cell.padding = c(0, 0, 0, 0),
    canvas.xlim = c(-1.0, 1.0),
    canvas.ylim = c(-1.0, 1.0)
  )

  circos.initialize(
    factors = factor(artox_sectors, levels = artox_sectors),
    xlim = cbind(rep(0, length(artox_sectors)), artox_widths)
  )

  # Ring 3: Residue names and counts
  circos.track(ylim = c(0, 1), track.height = 0.10, bg.border = "white",
    panel.fun = function(x, y) {
      sector.name <- CELL_META$sector.index
      if (sector.name == "Position") {
        circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                    col = "grey60", border = "white")
        circos.text(CELL_META$xcenter, 0.5, "Position",
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 1.2, col = "white", font = 2)
      } else {
        circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                    col = artox_colors[sector.name], border = "white")
        n_count <- format(artox_info[[sector.name]]$n_peptides, big.mark = ",")
        circos.text(CELL_META$xcenter, 0.5, n_count,
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 0.9, col = "white", font = 2)
        # Residue label outside
        res_label <- sector.name
        circos.text(CELL_META$xcenter, 1.6, res_label,
                    facing = "bending.outside", niceFacing = TRUE,
                    cex = 1.4, col = artox_colors[sector.name], font = 2)
      }
    })

  # Ring 2: Position subrings
  artox_subring_height <- 0.50 / n_positions
  for (pos_idx in rev(seq_along(positions))) {
    pos <- positions[pos_idx]
    circos.track(ylim = c(0, 1), track.height = artox_subring_height,
      bg.col = "grey93", bg.border = "grey88",
      panel.fun = function(x, y) {
        sector.name <- CELL_META$sector.index
        if (sector.name == "Position") {
          circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                      col = "grey85", border = "grey80")
          circos.text(CELL_META$xcenter, 0.5, as.character(pos),
                      facing = "bending.inside", niceFacing = TRUE,
                      cex = 0.85, col = "gray20", font = 2)
        } else {
          info <- artox_info[[sector.name]]
          pos_counts <- info$pos_counts
          cumsum_counts <- cumsum(pos_counts)
          count_at_pos <- pos_counts[as.character(pos)]
          if (count_at_pos > 0) {
            x_start <- if (pos == 1) 0 else cumsum_counts[as.character(pos - 1)]
            x_end <- cumsum_counts[as.character(pos)]
            circos.rect(x_start, 0, x_end, 1,
                        col = artox_colors[sector.name], border = NA)
          }
        }
      })
  }

  # Center legend
  n_artox_total <- sum(sapply(artox_info, function(x) x$n_peptides))
  text(0, 0.10, paste0("Artifact Oxidation"), cex = 1.3, col = "#FF6F00", font = 2)
  text(0, 0.02, paste0("n = ", format(n_artox_total, big.mark = ",")), cex = 1.1, col = "grey30", font = 1)
  y_pos <- seq(-0.08, -0.08 - 0.06 * (length(artox_order) - 1), by = -0.06)
  for (i in seq_along(artox_order)) {
    r <- artox_order[i]
    text(0, y_pos[i], paste0(r, ": ", format(artox_info[[r]]$n_peptides, big.mark = ",")),
         cex = 1.0, col = artox_colors[r], font = 2)
  }

  dev.off()
  circos.clear()
  cat("Saved Figure4A2_Circos_ArtifactOxidation.png\n")

  # --- PDF ---
  pdf("figure_panels/Figure4A2_Circos_ArtifactOxidation.pdf", width = 14, height = 14)

  circos.clear()
  circos.par(
    gap.after = c(4, rep(3, length(artox_order) - 1), 4),
    start.degree = 90,
    track.margin = c(0.001, 0.001),
    cell.padding = c(0, 0, 0, 0),
    canvas.xlim = c(-1.0, 1.0),
    canvas.ylim = c(-1.0, 1.0)
  )

  circos.initialize(
    factors = factor(artox_sectors, levels = artox_sectors),
    xlim = cbind(rep(0, length(artox_sectors)), artox_widths)
  )

  circos.track(ylim = c(0, 1), track.height = 0.10, bg.border = "white",
    panel.fun = function(x, y) {
      sector.name <- CELL_META$sector.index
      if (sector.name == "Position") {
        circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = "grey60", border = "white")
        circos.text(CELL_META$xcenter, 0.5, "Position", facing = "bending.inside", niceFacing = TRUE, cex = 1.2, col = "white", font = 2)
      } else {
        circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = artox_colors[sector.name], border = "white")
        circos.text(CELL_META$xcenter, 0.5, format(artox_info[[sector.name]]$n_peptides, big.mark=","), facing = "bending.inside", niceFacing = TRUE, cex = 0.9, col = "white", font = 2)
        res_label <- sector.name
        circos.text(CELL_META$xcenter, 1.6, res_label, facing = "bending.outside", niceFacing = TRUE, cex = 1.4, col = artox_colors[sector.name], font = 2)
      }
    })

  for (pos_idx in rev(seq_along(positions))) {
    pos <- positions[pos_idx]
    circos.track(ylim = c(0, 1), track.height = artox_subring_height, bg.col = "grey93", bg.border = "grey88",
      panel.fun = function(x, y) {
        sector.name <- CELL_META$sector.index
        if (sector.name == "Position") {
          circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = "grey85", border = "grey80")
          circos.text(CELL_META$xcenter, 0.5, as.character(pos), facing = "bending.inside", niceFacing = TRUE, cex = 0.85, col = "gray20", font = 2)
        } else {
          info <- artox_info[[sector.name]]
          pos_counts <- info$pos_counts
          cumsum_counts <- cumsum(pos_counts)
          count_at_pos <- pos_counts[as.character(pos)]
          if (count_at_pos > 0) {
            x_start <- if (pos == 1) 0 else cumsum_counts[as.character(pos - 1)]
            x_end <- cumsum_counts[as.character(pos)]
            circos.rect(x_start, 0, x_end, 1, col = artox_colors[sector.name], border = NA)
          }
        }
      })
  }

  text(0, 0.10, paste0("Artifact Oxidation"), cex = 1.3, col = "#FF6F00", font = 2)
  text(0, 0.02, paste0("n = ", format(n_artox_total, big.mark = ",")), cex = 1.1, col = "grey30", font = 1)
  for (i in seq_along(artox_order)) {
    r <- artox_order[i]
    text(0, y_pos[i], paste0(r, ": ", format(artox_info[[r]]$n_peptides, big.mark = ",")),
         cex = 1.0, col = artox_colors[r], font = 2)
  }

  dev.off()
  circos.clear()
  cat("Saved Figure4A2_Circos_ArtifactOxidation.pdf\n")
}
