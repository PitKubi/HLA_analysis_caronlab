# =============================================================================
# Figure 6E: N-Glycosylation glycan types & complexity vs HLA-II binding
# Standalone HLA-II version (the HLA-I Figure6DE_Improved.R reads HLA-I column
# names and bundles a 6D block; this isolates 6E and maps HLA-II columns).
#
# HLA-II specifics:
#  - Sheet "N-Glyco" (left detail block is row-aligned with the canonical
#    binding block; UniMod in left block, Binder in canonical block).
#  - Glycan reference sheet is "Supp" (cols UniMod# / Glycan name /
#    Mass shift (Da)); HLA-II has no "Linked glycan" column.
#  - Length filter 10-25 (per Moh, 2026-05-11). Each glyco-ID counted
#    distinctly (Moh, 2026-01-26: identical sequences w/ different glycans
#    are distinct peptides).
# =============================================================================
suppressMessages({library(readxl); library(dplyr); library(tidyr); library(ggplot2)})

xlsx_file  <- "HLA-II_JY_DDA_PTMs_ResultsSummary_05012026_4.xlsx"
output_dir <- "figure_panels"
dir.create(output_dir, showWarnings = FALSE)
binder_colors <- c("Strong"="#2E7D32","Weak"="#FFA726","Non-binder"="#90A4AE")

ng  <- suppressMessages(read_excel(xlsx_file, sheet = "N-Glyco"))
sup <- suppressMessages(read_excel(xlsx_file, sheet = "Supp"))
cn  <- colnames(ng)
# Glyco identity is embedded in the GlycosylatedPeptideSequence as (UniMod:NN);
# the standalone `UniMod` column is a sparse in-sheet summary (NA for most
# rows) and must NOT be used per-peptide.
glyco_col  <- grep("^GlycosylatedPeptideSequence", cn, value = TRUE)[1]
strip_col  <- grep("^StrippedPeptideSequence", cn, value = TRUE)[1]
binder_col <- grep("^Binder", cn, value = TRUE)[1]
stopifnot(!is.na(glyco_col), !is.na(strip_col), !is.na(binder_col))

nglyco_filt <- ng %>%
  transmute(Glyco   = .data[[glyco_col]],
            Peptide = .data[[strip_col]],
            Binder  = .data[[binder_col]]) %>%
  filter(!is.na(Glyco), !is.na(Peptide), !is.na(Binder)) %>%
  mutate(Length = nchar(gsub("\\(UniMod:[0-9]+\\)", "", Peptide))) %>%
  filter(Length >= 10 & Length <= 25)

# First UniMod tag in the glyco sequence = the glycan on that glyco-ID
nglyco_filt$UniMod_Num <- suppressWarnings(as.integer(
  ifelse(grepl("UniMod:[0-9]+", nglyco_filt$Glyco),
         sub(".*?UniMod:([0-9]+).*", "\\1", nglyco_filt$Glyco), NA)))

glycan_ref <- sup %>% select(`UniMod#`, `Glycan name`, `Mass shift (Da)`) %>%
  filter(!is.na(`UniMod#`))
nglyco_with_ref <- nglyco_filt %>%
  left_join(glycan_ref, by = c("UniMod_Num" = "UniMod#"))
nglyco_with_ref$Glycan_Label <- nglyco_with_ref$`Glycan name`
nglyco_with_ref$Glycan_Label[is.na(nglyco_with_ref$Glycan_Label)] <- "Unknown"

count_monosaccharides <- function(g){
  if (is.na(g) || g == "Unknown") return(0)
  n <- as.numeric(gsub("[^0-9]","", unlist(regmatches(g, gregexpr("\\([0-9]+\\)", g)))))
  if (length(n)==0) return(1)
  sum(n)
}
nglyco_with_ref$Complexity <- sapply(nglyco_with_ref$Glycan_Label, count_monosaccharides)
nglyco_with_ref$Complexity_Cat <- cut(nglyco_with_ref$Complexity,
  breaks=c(-Inf,3,8,15,Inf),
  labels=c("Simple (1-3)","Medium (4-8)","Complex (9-15)","Very Complex (>15)"))
nglyco_with_ref$Complexity_Cat <- as.character(nglyco_with_ref$Complexity_Cat)
nglyco_with_ref$Complexity_Cat[nglyco_with_ref$Glycan_Label=="Unknown"] <- "Unknown"
nglyco_with_ref$Complexity_Cat[is.na(nglyco_with_ref$Complexity_Cat)] <- "Unknown"
nglyco_with_ref$Complexity_Cat <- factor(nglyco_with_ref$Complexity_Cat,
  levels=c("Simple (1-3)","Medium (4-8)","Complex (9-15)","Very Complex (>15)","Unknown"))
nglyco_with_ref$Binder <- factor(nglyco_with_ref$Binder,
  levels=c("Non-binder","Weak","Strong"))

cat("N-Glyco (10-25mer, with binding):", nrow(nglyco_with_ref), "glyco-IDs\n")
print(table(nglyco_with_ref$Complexity_Cat))

complexity_summary <- nglyco_with_ref %>% group_by(Complexity_Cat) %>%
  summarise(n=n(),
            n_strong=sum(Binder=="Strong"), n_weak=sum(Binder=="Weak"),
            n_nonbinder=sum(Binder=="Non-binder"),
            pct_binder=round(100*(n_strong+n_weak)/n,1), .groups="drop")
complexity_long <- nglyco_with_ref %>% group_by(Complexity_Cat, Binder) %>%
  summarise(n=n(), .groups="drop")

fig6e_1 <- ggplot(complexity_long, aes(Complexity_Cat, n, fill=Binder)) +
  geom_bar(stat="identity", position="stack", width=0.7) +
  geom_text(data=complexity_summary,
            aes(Complexity_Cat, n + max(n)*0.05,
                label=paste0("n=",n,"\n(",pct_binder,"% bind)\n")),
            inherit.aes=FALSE, size=3, vjust=0, lineheight=0.9) +
  scale_fill_manual(values=binder_colors, name="HLA-II Binding") +
  labs(title="HLA-II N-Glycosylation: Glycan Complexity vs Binding",
       subtitle="Complexity = sum of monosaccharide units | Peptides length 10-25",
       x="Glycan Complexity", y="Number of Glyco-peptides") +
  theme_minimal() +
  theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
        plot.subtitle=element_text(size=9,hjust=0.5,color="gray40"),
        axis.text.x=element_text(size=9,angle=15,hjust=0.5),
        legend.position="bottom", panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank()) +
  ylim(0, max(complexity_summary$n)*1.3)
ggsave(file.path(output_dir,"Figure6E_Glycan_Complexity.png"), fig6e_1, width=9, height=6, dpi=300, bg="white")
ggsave(file.path(output_dir,"Figure6E_Glycan_Complexity.pdf"), fig6e_1, width=9, height=6, bg="white")
cat("Saved Figure6E_Glycan_Complexity.png/pdf\n")

glycan_type_summary <- nglyco_with_ref %>% filter(Glycan_Label!="Unknown") %>%
  group_by(Glycan_Label, Complexity) %>%
  summarise(n=n(), n_strong=sum(Binder=="Strong"), n_weak=sum(Binder=="Weak"),
            n_nonbinder=sum(Binder=="Non-binder"),
            mass_shift=first(`Mass shift (Da)`), .groups="drop") %>%
  arrange(Complexity, desc(n))
glycan_type_summary$Glycan_Label <- factor(glycan_type_summary$Glycan_Label,
  levels=glycan_type_summary$Glycan_Label)
glycan_long <- glycan_type_summary %>%
  select(Glycan_Label, Complexity, n_strong, n_weak, n_nonbinder) %>%
  pivot_longer(c(n_strong,n_weak,n_nonbinder), names_to="Binder_Type", values_to="Count") %>%
  mutate(Binder_Type=factor(gsub("n_","",Binder_Type),
         levels=c("nonbinder","weak","strong"),
         labels=c("Non-binder","Weak","Strong")))

fig6e_2 <- ggplot(glycan_long, aes(y=Glycan_Label, x=Count, fill=Binder_Type)) +
  geom_bar(stat="identity", position="stack") +
  geom_text(data=glycan_type_summary,
            aes(y=Glycan_Label, x=n+0.3, label=paste0("(",Complexity," units)")),
            inherit.aes=FALSE, size=2.5, hjust=0, color="gray50") +
  scale_fill_manual(values=binder_colors, name="HLA-II Binding") +
  labs(title="HLA-II N-Glycosylation: Individual Glycan Types",
       subtitle="Sorted by glycan complexity (monosaccharide count) | Peptides length 10-25",
       x="Number of Glyco-peptides", y="Glycan Type") +
  theme_minimal() +
  theme(plot.title=element_text(face="bold",size=11,hjust=0.5),
        plot.subtitle=element_text(size=9,hjust=0.5,color="gray40"),
        axis.text.y=element_text(size=7), legend.position="bottom",
        panel.grid.minor=element_blank(), panel.grid.major.y=element_blank()) +
  scale_x_continuous(expand=expansion(mult=c(0,0.25)))
ggsave(file.path(output_dir,"Figure6E2_Glycan_Types_Detail.png"), fig6e_2, width=10, height=8, dpi=300, bg="white")
ggsave(file.path(output_dir,"Figure6E2_Glycan_Types_Detail.pdf"), fig6e_2, width=10, height=8, bg="white")
cat("Saved Figure6E2_Glycan_Types_Detail.png/pdf\n")
write.csv(glycan_type_summary, file.path(output_dir,"data_figure6E_glycan_analysis.csv"), row.names=FALSE)
