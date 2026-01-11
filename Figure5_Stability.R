# =============================================================================
# FIGURE 5: PTM Effects on HLA Binding - Matching Example Style
# =============================================================================
# Left: Heatmap (residue-specific PTMs x Alleles) - using mean EL_Rank
# Right: EL Rank histograms (gray filled bg, red outline modified)
# NOTE: Half-life (t1/2) requires NetMHCstab - not available in current data

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(cowplot)
library(grid)

# =============================================================================
# CONFIGURATION
# =============================================================================
xlsx_file <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"

# =============================================================================
# LOAD DATA - Extract residue-specific binding info from each PTM sheet
# =============================================================================
cat("=== Figure 5: PTM Effects on HLA Binding ===\n\n")

# Function to extract residue-specific binding data
extract_residue_binding <- function(sheet_name, ptm_type) {
  suppressMessages({
    df <- read_excel(xlsx_file, sheet = sheet_name)
  })

  # Find columns
  binder_col <- colnames(df)[grepl("^Binder", colnames(df))][1]
  allele_col <- colnames(df)[grepl("^Best_Allele", colnames(df))][1]
  elrank_col <- colnames(df)[grepl("^EL_Rank", colnames(df))][1]

  # Find residue column (varies by sheet)
  residue_col <- NULL
  if (sheet_name == "Phospho.") {
    residue_col <- "1st Phospho residue"
  } else if (sheet_name == "bioOxid.") {
    # Need to parse from modification string
    residue_col <- NULL
  } else if (sheet_name == "Methyl." || sheet_name == "Dimethyl.") {
    residue_col <- colnames(df)[grepl("Residue", colnames(df))][1]
  }

  if (is.na(binder_col) || is.na(allele_col)) return(NULL)

  result <- data.frame(
    PTM_Type = ptm_type,
    Binder = df[[binder_col]],
    Allele = df[[allele_col]],
    EL_Rank = if (!is.na(elrank_col)) df[[elrank_col]] else NA,
    Residue = if (!is.null(residue_col) && residue_col %in% colnames(df)) df[[residue_col]] else NA,
    stringsAsFactors = FALSE
  )

  result <- result[!is.na(result$Binder) & !is.na(result$Allele), ]
  return(result)
}

# Load all PTM data with binding info
cat("Loading PTM binding data...\n")

# Define sheets and their PTM types
sheet_info <- list(
  list(sheet = "Phospho.", ptm = "phosphorylation"),
  list(sheet = "Acetyl.", ptm = "acetylation"),
  list(sheet = "Cyst.", ptm = "cysteinylation"),
  list(sheet = "Methyl.", ptm = "methylation"),
  list(sheet = "Dimethyl.", ptm = "dimethylation"),
  list(sheet = "Deamid. (NQ)", ptm = "deamidation"),
  list(sheet = "bioOxid.", ptm = "oxidation"),
  list(sheet = "Citrullination", ptm = "citrullination"),
  list(sheet = "G-Ubiq.", ptm = "ubiquitination")
)

all_data <- list()
for (info in sheet_info) {
  data <- extract_residue_binding(info$sheet, info$ptm)
  if (!is.null(data) && nrow(data) > 0) {
    # Create PTM label
    if (!all(is.na(data$Residue))) {
      data$PTM <- paste(data$Residue, info$ptm)
    } else {
      data$PTM <- info$ptm
    }
    all_data[[info$ptm]] <- data
    cat("  ", info$sheet, ":", nrow(data), "peptides\n")
  }
}

ptm_combined <- do.call(rbind, all_data)

# Clean up PTM names - use single letter codes
ptm_combined$PTM <- gsub("^NA ", "", ptm_combined$PTM)
ptm_combined$PTM <- trimws(ptm_combined$PTM)

cat("\nPTM types:\n")
print(table(ptm_combined$PTM))

# Load background
cat("\nLoading background...\n")
suppressMessages({
  bg <- read_excel(xlsx_file, sheet = "Background (wo PTMs)")
})
bg_data <- data.frame(
  PTM = "background",
  Binder = bg$Binder,
  Allele = bg$Best_Allele,
  EL_Rank = bg$EL_Rank,
  stringsAsFactors = FALSE
)
bg_data <- bg_data[!is.na(bg_data$Binder) & !is.na(bg_data$Allele), ]
cat("Background:", nrow(bg_data), "peptides\n")

# =============================================================================
# PANEL A: HEATMAP - Mean EL_Rank (lower = better binding)
# =============================================================================
cat("\n=== Creating Heatmap ===\n")

# Calculate mean EL_Rank for each PTM x Allele
ptm_stats <- ptm_combined %>%
  filter(!is.na(EL_Rank)) %>%
  group_by(PTM, Allele) %>%
  summarise(
    n = n(),
    Mean_EL = mean(EL_Rank, na.rm = TRUE),
    Pct_Strong = 100 * sum(Binder == "Strong") / n(),
    .groups = "drop"
  ) %>%
  filter(n >= 10)  # Minimum 10 peptides

bg_stats <- bg_data %>%
  filter(!is.na(EL_Rank)) %>%
  group_by(Allele) %>%
  summarise(
    Bg_EL = mean(EL_Rank, na.rm = TRUE),
    Bg_Strong = 100 * sum(Binder == "Strong") / n(),
    .groups = "drop"
  )

# Calculate log2 ratio of EL_Rank (negative = better than background)
heatmap_data <- ptm_stats %>%
  left_join(bg_stats, by = "Allele") %>%
  mutate(
    # For EL_Rank: lower is better, so ratio < 1 means better binding
    Log2_EL_Ratio = log2(Mean_EL / Bg_EL)
  )

# Create matrix for heatmap
heatmap_matrix <- heatmap_data %>%
  select(PTM, Allele, Log2_EL_Ratio) %>%
  pivot_wider(names_from = Allele, values_from = Log2_EL_Ratio) %>%
  as.data.frame()

rownames(heatmap_matrix) <- heatmap_matrix$PTM
heatmap_matrix$PTM <- NULL

# Remove rows with NA
heatmap_matrix <- heatmap_matrix[complete.cases(heatmap_matrix), ]

# Order by mean effect (most improved binding at top)
row_means <- rowMeans(heatmap_matrix, na.rm = TRUE)
heatmap_matrix <- heatmap_matrix[order(row_means), ]

cat("PTMs in heatmap:", nrow(heatmap_matrix), "\n")

# Color palette: Blue (lower EL = better) - White - Red (higher EL = worse)
heatmap_colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7",
                                      "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100)

# Cap values
cap_val <- 3
heatmap_capped <- pmin(pmax(as.matrix(heatmap_matrix), -cap_val), cap_val)

# Save heatmap
png("figure_panels/Figure5A_Heatmap.png", width = 5, height = 8, units = "in", res = 300)
pheatmap(
  heatmap_capped,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = heatmap_colors,
  breaks = seq(-cap_val, cap_val, length.out = 101),
  border_color = "white",
  cellwidth = 40,
  cellheight = 18,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "",
  legend = TRUE,
  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
  angle_col = 0,
  labels_col = c("A0201", "B0702", "C0702")
)
dev.off()
cat("Saved: figure_panels/Figure5A_Heatmap.png\n")

# =============================================================================
# PANEL B: EL RANK HISTOGRAMS - Per PTM type
# =============================================================================
cat("\n=== Creating EL Rank Histograms ===\n")

# Select PTMs to show (different types)
ptms_to_plot <- c("S phosphorylation", "acetylation", "cysteinylation",
                   "oxidation", "deamidation", "ubiquitination")

# Keep only PTMs that exist in data
ptms_to_plot <- ptms_to_plot[ptms_to_plot %in% unique(ptm_combined$PTM)]

plot_list <- list()
panel_labels <- letters[1:length(ptms_to_plot)]

for (i in seq_along(ptms_to_plot)) {
  ptm <- ptms_to_plot[i]

  # Get PTM EL_Rank (convert to 0-1 scale)
  ptm_el <- ptm_combined$EL_Rank[ptm_combined$PTM == ptm & !is.na(ptm_combined$EL_Rank)]
  ptm_el <- ptm_el / 100
  ptm_el <- ptm_el[ptm_el <= 1]

  # Get background
  bg_el <- bg_data$EL_Rank[!is.na(bg_data$EL_Rank)] / 100
  bg_el <- bg_el[bg_el <= 1]

  if (length(ptm_el) < 10) next

  n_ptm <- length(ptm_el)
  n_bg <- length(bg_el)

  # Create histogram bins
  breaks_seq <- seq(0, 1, by = 0.02)

  bg_hist <- hist(bg_el, breaks = breaks_seq, plot = FALSE)
  ptm_hist <- hist(ptm_el, breaks = breaks_seq, plot = FALSE)

  bg_frac <- bg_hist$counts / sum(bg_hist$counts)
  ptm_frac <- ptm_hist$counts / sum(ptm_hist$counts)

  # Create step coordinates
  n_bins <- length(bg_frac)
  bg_x <- numeric(0)
  bg_y <- numeric(0)
  ptm_x <- numeric(0)
  ptm_y <- numeric(0)

  for (j in 1:n_bins) {
    bg_x <- c(bg_x, breaks_seq[j], breaks_seq[j+1])
    bg_y <- c(bg_y, bg_frac[j], bg_frac[j])
    ptm_x <- c(ptm_x, breaks_seq[j], breaks_seq[j+1])
    ptm_y <- c(ptm_y, ptm_frac[j], ptm_frac[j])
  }

  # Background polygon
  bg_df <- data.frame(x = c(0, bg_x, 1), y = c(0, bg_y, 0))
  ptm_df <- data.frame(x = ptm_x, y = ptm_y)

  y_max <- max(c(bg_frac, ptm_frac), na.rm = TRUE) * 1.15

  # Format title
  title_text <- gsub("^([A-Z]) ", "\\1\n", ptm)

  p <- ggplot() +
    # Background - gray filled
    geom_polygon(data = bg_df, aes(x = x, y = y),
                 fill = "gray75", color = "gray50", linewidth = 0.3) +
    # PTM - RED OUTLINE ONLY
    geom_step(data = ptm_df, aes(x = x, y = y),
              color = "#C62828", linewidth = 0.8, direction = "hv") +
    labs(x = NULL, y = "Fraction") +
    ggtitle(ptm) +
    scale_x_continuous(limits = c(0, 0.8), breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5),
      axis.title.y = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.line = element_line(linewidth = 0.3),
      plot.margin = margin(3, 5, 3, 3)
    )

  plot_list[[ptm]] <- p
  cat("Created histogram for", ptm, "(n =", n_ptm, ")\n")
}

# Combine histograms
n_plots <- length(plot_list)
n_cols <- min(4, n_plots)
n_rows <- ceiling(n_plots / n_cols)

hist_combined <- plot_grid(plotlist = plot_list, nrow = n_rows, ncol = n_cols, align = "hv")

# Add common x-axis label
hist_with_label <- plot_grid(
  hist_combined,
  ggdraw() + draw_label("EL rank", size = 10),
  ncol = 1, rel_heights = c(1, 0.06)
)

ggsave("figure_panels/Figure5B_ELRank_Histograms.png", hist_with_label,
       width = 2.5 * n_cols, height = 2.5 * n_rows + 0.3, units = "in", dpi = 300)
cat("\nSaved: figure_panels/Figure5B_ELRank_Histograms.png\n")

ggsave("figure_panels/Figure5B_ELRank_Histograms.pdf", hist_with_label,
       width = 2.5 * n_cols, height = 2.5 * n_rows + 0.3, units = "in")
cat("Saved: figure_panels/Figure5B_ELRank_Histograms.pdf\n")

# =============================================================================
# SAVE DATA
# =============================================================================
write.csv(heatmap_data, "figure_panels/data_figure5_binding.csv", row.names = FALSE)

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=== Figure 5 Summary ===\n")
cat("Panel A: Heatmap of log2(PTM_EL_Rank / Background_EL_Rank)\n")
cat("  - Blue: better binding than background (lower EL rank)\n")
cat("  - Red: worse binding than background (higher EL rank)\n")
cat("Panel B: EL rank distribution histograms\n")
cat("  - Gray filled: background (unmodified)\n")
cat("  - Red outline: modified peptides\n")
cat("\nNOTE: Half-life (t1/2) stability data requires NetMHCstab predictions\n")
cat("      which are not included in the current Excel dataset.\n")

cat("\nTop PTMs by binding effect:\n")
print(heatmap_data %>%
        group_by(PTM) %>%
        summarise(Mean_Log2_EL = round(mean(Log2_EL_Ratio, na.rm = TRUE), 2),
                  n = sum(n)) %>%
        arrange(Mean_Log2_EL) %>%
        head(15))
