# =============================================================================
# FIGURE 6A-6D (HLA-II): Specific PTM Analysis Panels
# =============================================================================
# HLA-II rewrite of HLA-I Figure6_SpecificPTMs.R + the 6D block of
# Figure6DE_Improved.R. (6E is the separate validated Figure6E_Glycan.R.)
#
# The HLA-I scripts re-derived everything from the xlsx with HLA-I-specific
# hardcoded column positions (Peptide...142 etc.). For HLA-II we instead use
# the validated loader CSVs wherever the loader already captured the quantity
# (6B, 6C, 6D); only 6A needs the raw xlsx, because the loader collapses each
# phospho peptide to its FIRST site and 6A needs the full multi-site list.
#
#  6A Doubly-phosphorylated peptides: site position distribution
#     -> xlsx "Phospho." left block (Peptide...1 stripped + Assigned
#        Modifications), phospho = mass 79.966; doubly+ = >=2 phospho sites.
#        Phospho universe (818) reconciles with data_summary Phosphorylation.
#  6B Ubiquitination G vs GG  -> data_binding Subtype (count-all)
#  6C Oxidation biological vs artifact -> data_summary "All" counts
#  6D Biological oxidation: binding by oxidized residue (stacked %) + 6D2
#     residue x allele %Strong heatmap -> data_binding(PTM==Oxidation) joined
#     row-aligned with data_ptm_sites Residue.
# Length 10-25 throughout (loader-enforced for the CSV-derived panels).
# =============================================================================
suppressMessages({library(readxl); library(dplyr); library(tidyr)
                   library(ggplot2); library(pheatmap)})

xlsx_file  <- "HLA-II_JY_DDA_PTMs_ResultsSummary_05012026_4.xlsx"
output_dir <- "figure_panels"
cat("=== Figure 6A-6D (HLA-II): Specific PTM Analysis ===\n\n")

cfg           <- readRDS(file.path(output_dir, "config.rds"))
ptm_colors    <- cfg$ptm_colors
binder_colors <- c("Strong" = "#2E7D32", "Weak" = "#FFA726", "Non-binder" = "#90A4AE")
bind <- read.csv(file.path(output_dir, "data_binding.csv"),   stringsAsFactors = FALSE)
site <- read.csv(file.path(output_dir, "data_ptm_sites.csv"), stringsAsFactors = FALSE)
summ <- read.csv(file.path(output_dir, "data_summary.csv"),   stringsAsFactors = FALSE)
stopifnot(nrow(bind) == nrow(site),
          all(paste(bind$PTM, bind$Peptide) == paste(site$PTM, site$Peptide)))

# -----------------------------------------------------------------------------
# 6A: Doubly phosphorylated peptides - site position distribution
# -----------------------------------------------------------------------------
cat("6A: Doubly-phosphorylated peptides...\n")
ph <- suppressMessages(read_excel(xlsx_file, sheet = "Phospho."))
pep <- ph[["Peptide...1"]]; am <- ph[["Assigned Modifications"]]
keep <- !is.na(pep) & !is.na(am); pep <- pep[keep]; am <- am[keep]
strip  <- gsub("[^A-Za-z]", "", gsub("\\(UniMod:[0-9]+\\)", "", pep))
Length <- nchar(strip)

ph_sites <- lapply(strsplit(am, ",\\s*"), function(toks) {
  m <- regmatches(toks, regexpr("^[0-9]+[A-Z]\\(79\\.966", toks))
  as.integer(sub("([0-9]+).*", "\\1", m[nchar(m) > 0]))
})
n_site <- vapply(ph_sites, length, integer(1))
inrng  <- Length >= 10 & Length <= 25
multi  <- inrng & n_site >= 2
n_phos_total <- sum(inrng & n_site >= 1)
cat("  phospho peptides (10-25mer):", n_phos_total,
    " | data_summary anchor:", summ$Count[summ$PTM == "Phosphorylation" & summ$Residue == "All"],
    "\n  doubly+ phosphorylated:", sum(multi), "\n")

all_pos <- unlist(ph_sites[multi]); all_pos <- all_pos[all_pos >= 1 & all_pos <= 25]
site_counts <- as.data.frame(table(factor(all_pos, levels = 1:25)))
colnames(site_counts) <- c("Position", "Count")
site_counts$Position <- as.numeric(as.character(site_counts$Position))

fig6a <- ggplot(site_counts, aes(Position, Count)) +
  geom_bar(stat = "identity", fill = ptm_colors["Phosphorylation"], width = 0.7) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")), vjust = -0.5, size = 2.6) +
  scale_x_continuous(breaks = 1:25) +
  labs(title = "HLA-II Doubly Phosphorylated Peptides: Site Position Distribution",
       subtitle = paste0("n = ", sum(multi), " peptides with >=2 phospho sites (length 10-25)"),
       x = "Phosphorylation Site Position", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
        panel.grid.minor = element_blank()) +
  ylim(0, max(site_counts$Count) * 1.15)
ggsave(file.path(output_dir, "Figure6A_DoublyPhospho_Sites.png"), fig6a,
       width = 9, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6A_DoublyPhospho_Sites.pdf"), fig6a,
       width = 9, height = 5, bg = "white")
cat("  Saved Figure6A_DoublyPhospho_Sites.png/pdf\n")

# -----------------------------------------------------------------------------
# 6B: Ubiquitination G vs GG  (count-all, from validated loader Subtype)
# -----------------------------------------------------------------------------
cat("\n6B: Ubiquitination G vs GG...\n")
ubiq <- bind %>% filter(PTM == "Ubiquitination", !is.na(Subtype)) %>%
  mutate(Type = paste0(Subtype, "-Ubiq")) %>%
  group_by(Type) %>% summarise(Count = n(), .groups = "drop")
cat("  ", paste(ubiq$Type, ubiq$Count, collapse = " | "),
    " | data_summary Ubiquitination All =",
    summ$Count[summ$PTM == "Ubiquitination" & summ$Residue == "All"], "\n")

fig6b <- ggplot(ubiq, aes(Type, Count, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("G-Ubiq" = "#E57373", "GG-Ubiq" = "#C62828")) +
  labs(title = "HLA-II Ubiquitination: G vs GG Distribution",
       subtitle = "Peptides length 10-25 (each glyco/mod variant counted distinctly)",
       x = "Ubiquitination Type", y = "Number of Peptides") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
        legend.position = "none", panel.grid.minor = element_blank()) +
  ylim(0, max(ubiq$Count) * 1.15)
ggsave(file.path(output_dir, "Figure6B_Ubiq_G_vs_GG.png"), fig6b,
       width = 6, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6B_Ubiq_G_vs_GG.pdf"), fig6b,
       width = 6, height = 5, bg = "white")
cat("  Saved Figure6B_Ubiq_G_vs_GG.png/pdf\n")

# -----------------------------------------------------------------------------
# 6C: Oxidation biological vs artifact  (validated loader summary counts)
# -----------------------------------------------------------------------------
cat("\n6C: Biological vs Artifact Oxidation...\n")
n_bio <- summ$Count[summ$PTM == "Oxidation"          & summ$Residue == "All"]
n_art <- summ$Count[summ$PTM == "Artifact Oxidation" & summ$Residue == "All"]
cat("  biological:", n_bio, " artifact:", n_art, "\n")
ox <- data.frame(Type = c("Biological\nOxidation", "Artifact\nOxidation"),
                  Count = c(n_bio, n_art)) %>%
  mutate(Percentage = round(100 * Count / sum(Count), 1),
         Label = paste0(Count, "\n(", Percentage, "%)"))
fig6c <- ggplot(ox, aes(x = "", y = Count, fill = Type)) +
  geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c("Biological\nOxidation" = unname(ptm_colors["Oxidation"]),
                               "Artifact\nOxidation"   = unname(ptm_colors["Artifact Oxidation"]))) +
  labs(title = "HLA-II Oxidation: Biological vs Artifact",
       subtitle = "Peptides length 10-25") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
        legend.title = element_blank(), legend.position = "bottom")
ggsave(file.path(output_dir, "Figure6C_Oxidation_Bio_vs_Art.png"), fig6c,
       width = 6, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6C_Oxidation_Bio_vs_Art.pdf"), fig6c,
       width = 6, height = 6, bg = "white")
cat("  Saved Figure6C_Oxidation_Bio_vs_Art.png/pdf\n")

# -----------------------------------------------------------------------------
# 6D + 6D2: Biological oxidation - binding by oxidized residue
# -----------------------------------------------------------------------------
cat("\n6D: Biological oxidation binding by residue...\n")
aa_names <- c(A="Ala",C="Cys",D="Asp",E="Glu",F="Phe",G="Gly",H="His",I="Ile",
              K="Lys",L="Leu",M="Met",N="Asn",P="Pro",Q="Gln",R="Arg",S="Ser",
              T="Thr",V="Val",W="Trp",Y="Tyr")
ox_bind <- bind %>% mutate(Residue = site$Residue) %>%
  filter(PTM == "Oxidation", !is.na(Binder), !is.na(Residue), Residue != "")

res_sum <- ox_bind %>%
  group_by(Residue) %>%
  summarise(n_total = n(),
            pct_strong    = round(100 * sum(Binder == "Strong")     / n(), 1),
            pct_weak      = round(100 * sum(Binder == "Weak")       / n(), 1),
            pct_nonbinder = round(100 * sum(Binder == "Non-binder") / n(), 1),
            .groups = "drop") %>%
  filter(n_total >= 5) %>% arrange(desc(pct_strong))
res_sum$Label <- paste0(res_sum$Residue, " (", aa_names[res_sum$Residue],
                        ")\nn=", res_sum$n_total)
cat("  oxidized-residue summary:\n"); print(res_sum)

res_long <- res_sum %>%
  select(Residue, pct_strong, pct_weak, pct_nonbinder) %>%
  pivot_longer(starts_with("pct_"), names_to = "BT", values_to = "Percentage") %>%
  mutate(BT = factor(gsub("pct_", "", BT),
                      levels = c("nonbinder","weak","strong"),
                      labels = c("Non-binder","Weak","Strong")),
         Residue = factor(Residue, levels = res_sum$Residue))

fig6d <- ggplot(res_long, aes(Residue, Percentage, fill = BT)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  geom_text(data = res_sum, aes(Residue, 105, label = paste0(pct_strong, "%")),
            inherit.aes = FALSE, size = 3, fontface = "bold", color = "#2E7D32") +
  scale_fill_manual(values = binder_colors, name = "HLA-II Binding") +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25)) +
  scale_x_discrete(labels = setNames(res_sum$Label, res_sum$Residue)) +
  labs(title = "HLA-II Biological Oxidation: Binding Affinity by Oxidized Residue",
       subtitle = "Peptides length 10-25 | Sorted by % Strong Binders (shown above bars)",
       x = "Oxidized Residue", y = "Percentage of Peptides") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
        axis.text.x = element_text(face = "bold", size = 9),
        legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
ggsave(file.path(output_dir, "Figure6D_BioOxid_BindingByResidue.png"), fig6d,
       width = 11, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure6D_BioOxid_BindingByResidue.pdf"), fig6d,
       width = 11, height = 6, bg = "white")
cat("  Saved Figure6D_BioOxid_BindingByResidue.png/pdf\n")
write.csv(res_sum, file.path(output_dir, "data_figure6D_biooxid_binding.csv"), row.names = FALSE)

cat("\n6D2: Bio-oxidation %Strong heatmap by residue x allele...\n")
ra <- ox_bind %>%
  group_by(Residue, Allele) %>%
  summarise(n = n(), pct_strong = round(100 * sum(Binder == "Strong") / n(), 0),
            .groups = "drop") %>%
  filter(n >= 5)
hm <- ra %>% select(Residue, Allele, pct_strong) %>%
  pivot_wider(names_from = Allele, values_from = pct_strong) %>% as.data.frame()
rownames(hm) <- hm$Residue; hm$Residue <- NULL
hm <- hm[order(rowMeans(hm, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
hcols <- colorRampPalette(c("#FFFFFF","#C8E6C9","#81C784","#4CAF50","#2E7D32"))(100)
rlab  <- paste0(rownames(hm), " (", aa_names[rownames(hm)], ")")
for (dev_fun in list(
      function() png(file.path(output_dir, "Figure6D2_BioOxid_Heatmap.png"),
                     width = max(8, ncol(hm) * 0.7 + 3), height = max(5, nrow(hm) * 0.45 + 2),
                     units = "in", res = 300),
      function() pdf(file.path(output_dir, "Figure6D2_BioOxid_Heatmap.pdf"),
                     width = max(8, ncol(hm) * 0.7 + 3), height = max(5, nrow(hm) * 0.45 + 2)))) {
  dev_fun()
  pheatmap(as.matrix(hm), cluster_rows = FALSE, cluster_cols = FALSE,
           color = hcols, breaks = seq(0, 100, length.out = 101),
           border_color = "white", cellwidth = 26, cellheight = 20,
           fontsize_row = 10, fontsize_col = 8, angle_col = 90,
           main = "HLA-II Biological Oxidation: % Strong Binders by Residue and Allele",
           labels_row = rlab, display_numbers = TRUE, number_format = "%.0f",
           number_color = "black", fontsize_number = 7, na_col = "gray90",
           legend = TRUE, legend_breaks = c(0, 25, 50, 75, 100),
           legend_labels = c("0%","25%","50%","75%","100%"))
  dev.off()
}
cat("  Saved Figure6D2_BioOxid_Heatmap.png/pdf\n")

cat("\n=== 6A-6D Summary ===\n")
cat("6A doubly+ phospho:", sum(multi), "peptides\n")
cat("6B Ubiq:", paste(ubiq$Type, ubiq$Count, collapse = ", "), "\n")
cat("6C Oxidation bio:", n_bio, " artifact:", n_art, "\n")
cat("6D bio-oxid residues:", nrow(res_sum), " | 6D2 heatmap", nrow(hm), "x", ncol(hm), "\n")
