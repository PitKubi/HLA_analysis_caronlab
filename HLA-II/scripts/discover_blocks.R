# HLA-II canonical-block discovery + per-PTM count validation
# Read-only probe: confirms we can locate, for every PTM sheet, the canonical
# "Peptide -> Best_Allele/Core/EL_Rank/Binder/Alleles/Tool" block, and reports
# total vs 10-25mer counts + binder split, so we validate BEFORE drawing.
suppressMessages({library(readxl); library(dplyr)})

f <- "HLA-II_JY_DDA_PTMs_ResultsSummary_05012026.xlsx"

# PTM sheet -> display name (14 PTMs incl. Nitrosylation; excludes Background/Supp)
ptm_sheets <- c(
  "Phospho."="Phosphorylation","Cyst."="Cysteinylation","Deamid."="Deamidation",
  "Acetyl."="Acetylation","GG-Ubiq."="Ubiquitination (GG)","G-Ubiq."="Ubiquitination (G)",
  "Methyl."="Methylation","Dimethyl."="Dimethylation","Citrullination"="Citrullination",
  "bioOxid."="Oxidation","artOxid."="Artifact Oxidation","Carbamid."="Carbamidomethylation",
  "SUMO"="SUMOylation","N-Glyco"="N-Glycosylation","Nitrosylation"="Nitrosylation"
)

# Locate canonical block: the Peptide-like column immediately preceding the
# Best_Allele column that is itself followed by Core/EL_Rank/Binder.
find_block <- function(df){
  cn <- colnames(df)
  ba_idx <- grep("^Best_Allele", cn)
  for (bi in ba_idx){
    nxt <- tolower(paste(cn[bi:min(bi+6,length(cn))], collapse="|"))
    if (grepl("el_rank", nxt) && grepl("binder", nxt)){
      pep_i <- rev(grep("[Pp]eptide|Stripped", cn[1:(bi-1)]))
      if (length(pep_i)==0) next
      pcol <- cn[pep_i[1]]
      bcol <- grep("^Binder", cn[bi:length(cn)], value=TRUE)[1]
      ecol <- grep("^EL_Rank", cn[bi:length(cn)], value=TRUE)[1]
      return(list(pep=pcol, allele=cn[bi], binder=bcol, el=ecol))
    }
  }
  NULL
}

cat(sprintf("%-22s %-26s %7s %7s  %5s/%5s/%5s (S/W/N) in 10-25\n",
            "PTM","peptide_col","total","10-25","SB","WB","NB"))
cat(strrep("-",90),"\n")
grand <- 0
for (sh in names(ptm_sheets)){
  df <- suppressMessages(read_excel(f, sheet=sh))
  blk <- find_block(df)
  if (is.null(blk)){ cat(sprintf("%-22s  !! NO CANONICAL BLOCK FOUND\n", ptm_sheets[sh])); next }
  pep <- df[[blk$pep]]
  keep <- !is.na(pep) & nchar(as.character(pep))>0
  pep <- as.character(pep[keep])
  len <- nchar(pep)
  bnd <- if(!is.na(blk$binder)) as.character(df[[blk$binder]][keep]) else rep(NA,length(pep))
  in_rng <- len>=10 & len<=25
  sb <- sum(in_rng & grepl("^Strong", bnd, ignore.case=TRUE))
  wb <- sum(in_rng & grepl("^Weak",   bnd, ignore.case=TRUE))
  nb <- sum(in_rng & grepl("^Non",    bnd, ignore.case=TRUE))
  grand <- grand + sum(in_rng)
  cat(sprintf("%-22s %-26s %7d %7d  %5d/%5d/%5d\n",
              ptm_sheets[sh], substr(blk$pep,1,26), length(pep), sum(in_rng), sb, wb, nb))
}
cat(strrep("-",90),"\n")
cat(sprintf("TOTAL modified peptides (10-25mer, all PTMs): %d\n", grand))

# Background
bg <- suppressMessages(read_excel(f, sheet="Background"))
bcn <- colnames(bg)
bpi <- grep("^Peptide$|^Peptide\\.\\.\\.", bcn)[1]
bp <- as.character(bg[[bpi]]); bp <- bp[!is.na(bp) & nchar(bp)>0]
bl <- nchar(bp)
cat(sprintf("Background peptide col: '%s' | total=%d | 10-25mer=%d\n",
            bcn[bpi], length(bp), sum(bl>=10 & bl<=25)))
