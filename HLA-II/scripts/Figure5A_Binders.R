# =============================================================================
# FIGURE 5A (HLA-II): % Strong & Weak Binders by PTM(+Residue) and Allele
# =============================================================================
# HLA-II rewrite of the HLA-I Figure5A_StrongBinders.R.
#
# The HLA-I script re-derived binding+residue from the xlsx using HLA-I-
# specific hardcoded sheet names / column positions / 3 alleles. That is
# exactly the kind of mis-mapping Moh iterated on. The HLA-II loader
# (data_loader_hla2.R) already extracted the validated canonical-block
# peptides with their best allele, binder, and left-block residue, and
# reconciled every PTM against Moh's embedded summaries. So here we consume
# the validated loader CSVs instead of re-deriving:
#   data_binding.csv  : PTM, Peptide, Length, Allele(=Best_Allele), EL_Rank, Binder, Subtype
#   data_ptm_sites.csv : PTM, Peptide, Site, Residue, Length, Subtype  (row-aligned)
#
# Method (faithful to HLA-I): group by (PTM[,Residue], Allele); per cell
# %Strong = 100*n_strong/n, %Weak = 100*n_weak/n; keep cells with n >= 10;
# heatmap rows in canonical PTM order (shared config) then residue,
# columns = the 14 HLA-II alleles.
# Length 10-25 already enforced by the loader.
# =============================================================================
suppressMessages({library(dplyr); library(tidyr); library(pheatmap)})

output_dir <- "figure_panels"
cat("=== Figure 5A (HLA-II): % Strong/Weak Binders by PTM+Residue ===\n\n")

bind <- read.csv(file.path(output_dir, "data_binding.csv"),   stringsAsFactors = FALSE)
site <- read.csv(file.path(output_dir, "data_ptm_sites.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(bind) == nrow(site),
          all(paste(bind$PTM, bind$Peptide, bind$Length) ==
              paste(site$PTM, site$Peptide, site$Length)))

merged <- bind %>%
  mutate(Residue = site$Residue) %>%
  filter(!is.na(Allele), !is.na(Binder), !is.na(Residue), Residue != "") %>%
  mutate(PTM_Residue = paste(Residue, tolower(PTM)))

alleles <- sort(unique(merged$Allele))
cat("HLA-II alleles (", length(alleles), "):\n", sep = ""); print(alleles)
cat("\nTotal binding records:", nrow(merged), "\n")

config     <- readRDS(file.path(output_dir, "config.rds"))
ptm_colors <- config$ptm_colors
ptm_order  <- config$ptm_order
stopifnot(all(merged$PTM %in% ptm_order))

strong_colors <- colorRampPalette(c("#FFFFFF","#DEEBF7","#9ECAE1","#4292C6","#08519C"))(100)
weak_colors   <- colorRampPalette(c("#FFFFFF","#FEE8C8","#FDBB84","#E34A33","#B30000"))(100)

# -----------------------------------------------------------------------------
create_binder_heatmap <- function(stats_df, value_col, palette, filename) {
  hm <- stats_df %>%
    select(PTM, Residue, PTM_Residue, Allele, !!sym(value_col)) %>%
    pivot_wider(names_from = Allele, values_from = !!sym(value_col)) %>%
    as.data.frame() %>%
    arrange(match(PTM, ptm_order), Residue)
  rownames(hm) <- hm$PTM_Residue

  allele_cols <- intersect(alleles, colnames(hm))
  mat <- hm[, allele_cols, drop = FALSE]
  mat[is.na(mat)] <- 0

  row_ann <- data.frame(PTM = hm$PTM); rownames(row_ann) <- rownames(mat)
  ann_colors <- list(PTM = ptm_colors[unique(hm$PTM)])
  ptm_rle <- rle(hm$PTM)
  gaps <- cumsum(ptm_rle$lengths)[-length(ptm_rle$lengths)]

  n_rows <- nrow(mat)
  fig_h <- max(7, n_rows * 0.24 + 2.5)
  fig_w <- max(8, length(allele_cols) * 0.85 + 4)

  for (dev_fun in list(
        function() png(file.path(output_dir, paste0(filename, ".png")),
                       width = fig_w, height = fig_h, units = "in", res = 300),
        function() pdf(file.path(output_dir, paste0(filename, ".pdf")),
                       width = fig_w, height = fig_h),
        function() tiff(file.path(output_dir, paste0(filename, ".tif")),
                        width = fig_w, height = fig_h, units = "in", res = 300,
                        compression = "lzw"))) {
    dev_fun()
    pheatmap(as.matrix(mat), cluster_rows = FALSE, cluster_cols = FALSE,
             color = palette, breaks = seq(0, 100, length.out = 101),
             border_color = "white", cellwidth = 22, cellheight = 13,
             fontsize_row = 7, fontsize_col = 8,
             main = NA,  # Manuscript panels carry no titles
             display_numbers = TRUE, number_format = "%.0f",
             number_color = "black", fontsize_number = 6, na_col = "gray90",
             legend = TRUE, legend_breaks = c(0, 25, 50, 75, 100),
             legend_labels = c("0%","25%","50%","75%","100%"),
             annotation_row = row_ann, annotation_colors = ann_colors,
             angle_col = 90, gaps_row = gaps)
    dev.off()
  }
  cat("Saved:", filename, ".png/.pdf/.tif  (", n_rows, "rows )\n")
  hm
}

# --- by PTM + Residue (per-cell n >= 10) -------------------------------------
stats_res <- merged %>%
  group_by(PTM, Residue, PTM_Residue, Allele) %>%
  summarise(n = n(),
            pct_strong = round(100 * sum(Binder == "Strong") / n(), 0),
            pct_weak   = round(100 * sum(Binder == "Weak")   / n(), 0),
            .groups = "drop") %>%
  filter(n >= 10)

# --- by PTM only (residues aggregated; per-cell n >= 10) ---------------------
stats_ptm <- merged %>%
  group_by(PTM, Allele) %>%
  summarise(n = n(),
            pct_strong = round(100 * sum(Binder == "Strong") / n(), 0),
            pct_weak   = round(100 * sum(Binder == "Weak")   / n(), 0),
            .groups = "drop") %>%
  filter(n >= 10) %>%
  mutate(Residue = "All", PTM_Residue = PTM)

create_binder_heatmap(stats_ptm, "pct_strong", strong_colors,
                      "Figure5A_StrongBinders_ByPTM")
create_binder_heatmap(stats_ptm, "pct_weak", weak_colors,
                      "Figure5A_WeakBinders_ByPTM")
create_binder_heatmap(stats_res, "pct_strong", strong_colors,
                      "Figure5A_StrongBinders_ByResidue")
create_binder_heatmap(stats_res, "pct_weak", weak_colors,
                      "Figure5A_WeakBinders_ByResidue")

write.csv(stats_res, file.path(output_dir, "data_figure5A_binders.csv"), row.names = FALSE)
cat("\nSaved data_figure5A_binders.csv\n")
cat("PTM+Residue cells (n>=10):", nrow(stats_res),
    "| PTM cells (n>=10):", nrow(stats_ptm), "\n")
