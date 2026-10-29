# =============================================================================
# Figure 5B (HLA-I): EL rank distribution per biological PTM vs background
# =============================================================================
# Manuscript revision (Moh points 5/8/9): every biological PTM gets a panel
# (no minimum-n dropping), Artifact Oxidation is excluded, and the grid is a
# fixed 4 x 4 in PTM_ORDER so the HLA-I and HLA-II exports are identical in
# size and panel positions for side-by-side mounting. Panels carry no
# titles/subtitles; the colored n-annotation identifies each PTM.
# Consumes the loader CSVs (data_binding.csv, data_background.csv).
# Run from the repo root.
# =============================================================================
suppressMessages({library(ggplot2); library(dplyr); library(cowplot)})
source("shared/ptm_config.R")

output_dir <- "figure_panels"
cat("=== Figure 5B (HLA-I): EL Rank Distribution ===\n\n")

binding_data <- read.csv(file.path(output_dir, "data_binding.csv"),    stringsAsFactors = FALSE)
background   <- read.csv(file.path(output_dir, "data_background.csv"), stringsAsFactors = FALSE)

# The HLA-I loader does not length-filter these CSVs; apply the canonical 8-14.
binding_data <- binding_data %>% filter(Length >= 8 & Length <= 14)
background   <- background   %>% filter(Length >= 8 & Length <= 14)

# The panel is a zoom into the binder-relevant EL-rank window; non-binders
# (>= 2) sit off-panel, so the n annotations report the FULL sample size.
EL_WINDOW_MAX <- 2
# Below this n, overlay the individual EL ranks as a rug so tiny samples are
# honest (Moh point 9; e.g. HLA-I Nitrosylation n = 5).
RUG_MAX_N <- 20

background_el_all <- background$EL_Rank[!is.na(background$EL_Rank)]
background_el     <- background_el_all[background_el_all <= EL_WINDOW_MAX]
n_background      <- length(background_el_all)
cat("Background peptides:", n_background,
    "(", length(background_el), "with EL_Rank <=", EL_WINDOW_MAX, ")\n\n")

build_elrank_panel <- function(ptm) {
  ptm_el_all <- binding_data %>%
    filter(PTM == ptm, !is.na(EL_Rank)) %>%
    pull(EL_Rank)
  ptm_el    <- ptm_el_all[ptm_el_all <= EL_WINDOW_MAX]
  n_ptm     <- length(ptm_el_all)
  ptm_color <- PTM_COLORS[[ptm]]
  cat(sprintf("%-22s n = %4d (%4d in window)\n", ptm, n_ptm, length(ptm_el)))

  # binwidth/boundary make the 40 bins tile [0, EL_WINDOW_MAX] exactly;
  # free binning let edge bars cross the scale limits and get silently dropped.
  p <- ggplot() +
    geom_histogram(data = data.frame(EL_Rank = background_el),
                   aes(x = EL_Rank, y = after_stat(density)),
                   binwidth = EL_WINDOW_MAX / 40, boundary = 0, fill = "gray70", color = "gray50",
                   linewidth = 0.2, alpha = 0.7)
  if (length(ptm_el) > 0) {
    p <- p +
      geom_histogram(data = data.frame(EL_Rank = ptm_el),
                     aes(x = EL_Rank, y = after_stat(density)),
                     binwidth = EL_WINDOW_MAX / 40, boundary = 0, fill = NA, color = ptm_color, linewidth = 0.8)
  }
  if (n_ptm < RUG_MAX_N && length(ptm_el) > 0) {
    p <- p +
      geom_rug(data = data.frame(EL_Rank = ptm_el),
               aes(x = EL_Rank), color = ptm_color,
               linewidth = 0.6, length = unit(0.06, "npc"))
  }
  p +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "darkgreen",  linewidth = 0.5) +
    geom_vline(xintercept = 2.0, linetype = "dashed", color = "darkorange", linewidth = 0.5) +
    labs(x = "EL Rank", y = "Density") +
    annotate("text", x = 1.5, y = Inf, vjust = 2,   hjust = 0.5, size = 2.5,
             # black, not grey: Peter's standing rule is no grey text in figures
             label = paste0("Background, n = ", format(n_background, big.mark = ",")),
             color = "black") +
    annotate("text", x = 1.5, y = Inf, vjust = 3.5, hjust = 0.5, size = 2.5,
             label = paste0(ptm, ", n = ", format(n_ptm, big.mark = ",")),
             color = ptm_color) +
    scale_x_continuous(limits = c(0, EL_WINDOW_MAX), breaks = c(0, 0.5, 1, 1.5, 2)) +
    theme_classic() +
    theme(axis.title  = element_text(size = 9),
          axis.text   = element_text(size = 8),
          plot.margin = margin(5, 10, 5, 5)) +
    strip_ms_theme()
}

panels <- lapply(BIO_PTMS, build_elrank_panel)
# 13 panels in a fixed 4 x 4 grid: the 3 trailing cells stay blank so both
# classes export with identical panel positions.
grid_plot <- plot_grid(plotlist = panels, ncol = 4, nrow = 4, align = "hv")

save_panel(grid_plot, "Figure5B_ELRank_Distribution", pair_id = "elrank")
cat("\nSaved Figure5B_ELRank_Distribution.png/pdf/tif\n")

binder_summary <- binding_data %>%
  filter(!is.na(EL_Rank), PTM %in% BIO_PTMS) %>%
  group_by(PTM) %>%
  summarise(n_total = n(),
            pct_strong    = round(100 * sum(EL_Rank < 0.5) / n(), 1),
            pct_weak      = round(100 * sum(EL_Rank >= 0.5 & EL_Rank < 2) / n(), 1),
            pct_nonbinder = round(100 * sum(EL_Rank >= 2) / n(), 1),
            .groups = "drop") %>%
  arrange(desc(pct_strong))
cat("\nBinder percentages by PTM:\n"); print(binder_summary)
