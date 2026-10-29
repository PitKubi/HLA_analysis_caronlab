# =============================================================================
# HLA-II PTM Analysis - Data Loader  (adapted from HLA-I data_loader.R)
# =============================================================================
# Key differences vs HLA-I (do NOT edit the HLA-I repo):
#  - Peptide length filter: 10-25 mer (per Moh, 2026-05-11) instead of 8-14.
#  - 14 PTMs incl. NEW "Nitrosylation".
#  - HLA-II sheets have a different layout: a per-PTM CANONICAL block
#      <Peptide-like col> [Assigned Modifications] Best_Allele Core EL_Rank
#      Binder Alleles Tool
#    plus an embedded summary block (2nd Best_Allele then All/SB/WB/NB/%).
#  - Counting rule (Moh, 2026-01-26): count ALL length-filtered peptides per
#    PTM from the canonical column; keep glyco/phospho variants DISTINCT (no
#    sequence dedup). Binder split only for binding-specific figures.
#  - Output schema identical to HLA-I so figure scripts work unchanged:
#      figure_panels/{data_ptm_sites,data_background,data_binding,data_summary}.csv
#      figure_panels/config.rds
#  - Also writes figure_panels/validation_counts.csv (canonical vs embedded
#    summary) for human sign-off before any figure is trusted.
# =============================================================================
suppressMessages({library(readxl); library(dplyr); library(tidyr)})

XLSX_FILE  <- "HLA-II_JY_DDA_PTMs_ResultsSummary_05012026_4.xlsx"
OUTPUT_DIR <- "figure_panels"
LEN_MIN <- 10L; LEN_MAX <- 25L

# Canonical PTM order + colors shared with the HLA-I pipeline (run from HLA-II/)
source("../shared/ptm_config.R")

# sheet -> display PTM name
PTM_SHEET <- c(
  "Phospho."="Phosphorylation","Cyst."="Cysteinylation","Deamid."="Deamidation",
  "Acetyl."="Acetylation","GG-Ubiq."="Ubiquitination","G-Ubiq."="Ubiquitination",
  "Methyl."="Methylation","Dimethyl."="Dimethylation","Citrullination"="Citrullination",
  "bioOxid."="Oxidation","artOxid."="Artifact Oxidation","Carbamid."="Carbamidomethylation",
  "SUMO"="SUMOylation","N-Glyco"="N-Glycosylation","Nitrosylation"="Nitrosylation"
)
UBIQ_SUBTYPE <- c("GG-Ubiq."="GG","G-Ubiq."="G")

# PTM_COLORS / PTM_ORDER come from ../shared/ptm_config.R (sourced above)

STRONG_BINDER_THRESHOLD <- 0.5; WEAK_BINDER_THRESHOLD <- 2.0
AA <- "ACDEFGHIKLMNPQRSTVWY"

is_pep_like <- function(x){
  x <- as.character(x); x <- x[!is.na(x) & nchar(x)>0]
  if (length(x)==0) return(FALSE)
  x <- gsub("\\(UniMod:[0-9]+\\)","",x)            # strip glyco UniMod tags
  frac <- mean(grepl(sprintf("^[%s]{6,}$",AA), x))
  frac > 0.8
}

# Locate the canonical block: from the Best_Allele that is followed by
# Core/EL_Rank/Binder/Alleles/Tool, walk LEFT past Assigned Modifications /
# Localized / Core / all-NA columns to the first peptide-like column.
find_canonical <- function(df){
  cn <- colnames(df)
  for (bi in grep("^Best_Allele", cn)){
    tail <- tolower(paste(cn[bi:min(bi+6,length(cn))], collapse="|"))
    if (!(grepl("el_rank",tail) && grepl("binder",tail) && grepl("tool",tail))) next
    j <- bi-1
    while (j>=1){
      nm <- cn[j]
      if (!grepl("^Assigned Modifications|^Localized|^Core|^\\.\\.\\.|^Residue", nm) &&
          is_pep_like(df[[j]])) {
        am <- grep("^Assigned Modifications", cn[j:bi], value=TRUE)[1]  # AM inside block, if any
        return(list(pep=cn[j],
                    am = if(length(am)&&!is.na(am)) am else NA_character_,
                    allele=cn[bi],
                    core=if (bi < length(cn) && grepl("^Core", cn[bi+1])) cn[bi+1] else NA_character_,
                    el=grep("^EL_Rank",cn[bi:length(cn)],value=TRUE)[1],
                    binder=grep("^Binder",cn[bi:length(cn)],value=TRUE)[1]))
      }
      j <- j-1
    }
  }
  NULL
}

# Embedded summary: 2nd Best_Allele block then All/SB/WB/NB rows -> sum gives
# Moh's own per-PTM peptide total (each peptide at its best allele).
embedded_total <- function(df){
  cn <- colnames(df)
  ba <- grep("^Best_Allele", cn)
  if (length(ba)<2) return(NA_integer_)
  ai <- grep("^All(\\.\\.\\.|$| )", cn)
  ai <- ai[ai>ba[length(ba)]][1]
  if (is.na(ai)) return(NA_integer_)
  v <- suppressWarnings(as.numeric(df[[ai]]))
  if (all(is.na(v))) return(NA_integer_)
  as.integer(sum(v, na.rm=TRUE))
}

parse_first_site <- function(mod_string, allowed=NULL){
  if (is.na(mod_string) || mod_string=="") return(NULL)
  mods <- trimws(strsplit(as.character(mod_string), ",\\s*")[[1]])
  for (m in mods){
    if (grepl("^N-term\\(", m)) { r<-"N-term"; s<-1 }
    else { mm <- regmatches(m, regexpr("^(\\d+)([A-Z])", m))
           if (length(mm)==0) next
           s <- as.numeric(gsub("[A-Z]","",mm)); r <- gsub("[0-9]","",mm) }
    if (is.null(allowed) || r %in% allowed) return(list(Site=s, Residue=r))
  }
  NULL
}

# Build a stripped-peptide -> "site|residue" map from a sheet's LEFT detail
# block (first ^Peptide column + first ^Assigned Modifications column).
# First relevant-residue site per sequence (mirrors HLA-I "first site").
left_block_site_map <- function(df, allowed){
  cn <- colnames(df)
  pcol <- grep("^Peptide(\\.\\.\\.|$)", cn, value=TRUE)[1]
  if (is.na(pcol)) pcol <- grep("^Peptide", cn, value=TRUE)[1]
  acol <- grep("^Assigned Modifications", cn, value=TRUE)[1]
  if (is.na(pcol) || is.na(acol)) return(NULL)
  p <- as.character(df[[pcol]]); a <- as.character(df[[acol]])
  ok <- !is.na(p) & nchar(p)>0 & !is.na(a)
  p <- p[ok]; a <- a[ok]
  m <- new.env(hash=TRUE, parent=emptyenv())
  for (i in seq_along(p)){
    key <- gsub("\\(UniMod:[0-9]+\\)","",p[i])
    if (!is.null(m[[key]])) next            # keep first occurrence
    ps <- parse_first_site(a[i], allowed)
    if (!is.null(ps)) m[[key]] <- c(ps$Site, ps$Residue)
  }
  m
}
# residue whitelists (mirror HLA-I loader intent; exclude co-occurring artifact oxid M/W/H/F)
RES <- list(
  Phosphorylation=c("S","T","Y"), Cysteinylation=c("C"), Deamidation=c("N","Q"),
  Acetylation=c("K","N-term","H","C","S","T","Y"), Ubiquitination=c("K","S","T","C"),
  Methylation=c("C","H","N","Q","I","L","D","E","S","T","N-term","K","R"),
  Dimethylation=c("K","R","N","P","N-term"), Citrullination=c("R"),
  Oxidation=c("P","I","L","Q","S","T","V","C","D","E","N","Y","G","K","R"),
  "Artifact Oxidation"=c("M","W","H","F"), Carbamidomethylation=c("C"),
  SUMOylation=c("K"), "N-Glycosylation"=c("N"), Nitrosylation=c("C","Y"))

cat("HLA-II data loader |", LEN_MIN, "-", LEN_MAX, "mer | 14 PTMs\n")
dir.create(OUTPUT_DIR, showWarnings=FALSE)

bind_rows_list <- list(); sites_list <- list(); val <- list()

for (sh in names(PTM_SHEET)){
  ptm <- unname(PTM_SHEET[sh])
  df <- suppressMessages(read_excel(XLSX_FILE, sheet=sh))
  blk <- find_canonical(df)
  if (is.null(blk)){ cat(sprintf("  !! %-22s NO BLOCK\n", ptm)); next }

  pep_raw <- as.character(df[[blk$pep]])
  keep <- !is.na(pep_raw) & nchar(pep_raw)>0
  pep_strip <- gsub("\\(UniMod:[0-9]+\\)","", pep_raw[keep])   # for N-glyco length
  L <- nchar(pep_strip)
  binder <- if(!is.na(blk$binder)) as.character(df[[blk$binder]][keep]) else NA
  el     <- if(!is.na(blk$el)) suppressWarnings(as.numeric(df[[blk$el]][keep])) else NA
  allele <- as.character(df[[blk$allele]][keep])
  am     <- if(!is.na(blk$am)) as.character(df[[blk$am]][keep]) else rep(NA, sum(keep))
  inrng  <- L>=LEN_MIN & L<=LEN_MAX

  # NetMHCIIpan 9-mer binding core (column right after Best_Allele in every
  # sheet). CoreStart = 1-based position of the core within the stripped
  # peptide; NA for the ~0.5% deletion/insertion cores that are not a
  # contiguous substring (excluded from core-based positional figures only).
  core <- if(!is.na(blk$core)) as.character(df[[blk$core]][keep]) else rep(NA_character_, sum(keep))
  core_start <- rep(NA_integer_, sum(keep))
  has_core <- !is.na(core) & nchar(core) > 0
  core_start[has_core] <- vapply(which(has_core), function(i) {
    p <- regexpr(core[i], pep_strip[i], fixed = TRUE)
    if (p > 0) as.integer(p) else NA_integer_
  }, integer(1))

  # binding rows (for data_binding + figure 2 background-style + binder figs)
  bind_rows_list[[sh]] <- data.frame(
    PTM=ptm, Peptide=pep_strip[inrng], Length=L[inrng],
    Allele=allele[inrng], EL_Rank=el[inrng], Binder=binder[inrng],
    Core=core[inrng], CoreStart=core_start[inrng],
    Subtype = if (sh %in% names(UBIQ_SUBTYPE)) unname(UBIQ_SUBTYPE[sh]) else NA_character_,
    stringsAsFactors=FALSE)

  # site/residue rows: one row per canonical peptide (keeps Moh's count).
  # Residue/site come from the sheet's LEFT detail block, mapped by stripped
  # sequence (the canonical binding block has no localization for many PTMs).
  allow <- RES[[ptm]]
  pep_s <- pep_strip                       # stripped canonical sequences
  site <- rep(NA_real_, sum(keep)); res <- rep(NA_character_, sum(keep))
  smap <- left_block_site_map(df, allow)
  if (!is.null(smap)){
    for (i in which(inrng)){
      v <- smap[[ pep_s[i] ]]
      if (!is.null(v)){ site[i] <- as.numeric(v[1]); res[i] <- v[2] }
    }
  }
  # Fallback: adjacent Assigned Modifications, if the block had one
  if (!all(is.na(am))){
    for (i in which(inrng & is.na(res))){
      ps <- parse_first_site(am[i], allow)
      if (!is.null(ps)){ site[i]<-ps$Site; res[i]<-ps$Residue }
    }
  }
  # N-Glyco: residue N; site = N position from the GlycosylatedPeptideSequence
  if (ptm=="N-Glycosylation"){
    gcol <- grep("^Glycosylated", colnames(df), value=TRUE)[1]
    if (!is.na(gcol)){
      gseq <- as.character(df[[gcol]][keep])
      for (i in which(inrng)){
        mpos <- regexpr("N\\(UniMod:[0-9]+\\)", gseq[i])
        if (mpos>0){ pre<-substr(gseq[i],1,mpos-1)
          site[i]<-nchar(gsub("\\(UniMod:[0-9]+\\)","",pre))+1; res[i]<-"N" }
      }
    }
  }
  # Core-relative site (P1-P9) where the site falls inside the predicted core
  core_site <- site - core_start + 1
  core_site[is.na(core_start) | core_site < 1 | core_site > 9] <- NA_real_

  sites_list[[sh]] <- data.frame(
    PTM=ptm,
    Peptide=pep_strip[inrng],
    Site=site[inrng],
    Residue=res[inrng],
    Length=L[inrng],
    Core=core[inrng],
    CoreStart=core_start[inrng],
    CoreSite=core_site[inrng],
    Subtype = if (sh %in% names(UBIQ_SUBTYPE)) unname(UBIQ_SUBTYPE[sh]) else NA_character_,
    stringsAsFactors=FALSE)

  emb <- embedded_total(df)
  val[[sh]] <- data.frame(PTM=ptm, sheet=sh, pep_col=blk$pep,
                          canonical_total=sum(keep),
                          canonical_10_25=sum(inrng),
                          embedded_All_sum=emb,
                          SB=sum(inrng & grepl("^Strong",binder,ignore.case=TRUE)),
                          WB=sum(inrng & grepl("^Weak",binder,ignore.case=TRUE)),
                          NB=sum(inrng & grepl("^Non",binder,ignore.case=TRUE)),
                          stringsAsFactors=FALSE)
}

binding <- bind_rows(bind_rows_list)
ptm_sites <- bind_rows(sites_list)
validation <- bind_rows(val)

# ---- Background (unmodified ligandome) ----
bg <- suppressMessages(read_excel(XLSX_FILE, sheet="Background"))
bcn <- colnames(bg)
bpi <- grep("^Peptide$", bcn); if (length(bpi)==0) bpi <- grep("^Peptide", bcn)
bpi <- bpi[1]
bp <- as.character(bg[[bpi]]); ok <- !is.na(bp) & nchar(bp)>0
bg_el <- grep("^EL_Rank", bcn, value=TRUE)[1]
bg_bd <- grep("^Binder", bcn, value=TRUE)[1]
bg_al <- grep("^Best_Allele", bcn, value=TRUE)[1]
bg_co <- grep("^Core", bcn, value=TRUE)[1]
background <- data.frame(
  Peptide=bp[ok], Length=nchar(bp[ok]),
  Allele=if(!is.na(bg_al)) as.character(bg[[bg_al]][ok]) else NA,
  EL_Rank=if(!is.na(bg_el)) suppressWarnings(as.numeric(bg[[bg_el]][ok])) else NA,
  Binder=if(!is.na(bg_bd)) as.character(bg[[bg_bd]][ok]) else NA,
  Core=if(!is.na(bg_co)) as.character(bg[[bg_co]][ok]) else NA,
  stringsAsFactors=FALSE)
background <- background[background$Length>=LEN_MIN & background$Length<=LEN_MAX, ]

# ---- Summary (for Fig 1A/1B) : honour 10-25 + count-all rule ----
n_bg <- nrow(background); n_mod <- nrow(ptm_sites); n_tot <- n_bg + n_mod
ptm_summary <- ptm_sites %>% group_by(PTM) %>% summarise(Count=n(), .groups="drop") %>%
  mutate(Percentage=Count/n_mod*100, Residue="All")
ptm_res <- ptm_sites %>% filter(!is.na(Residue)) %>% group_by(PTM,Residue) %>%
  summarise(Count=n(), .groups="drop") %>%
  mutate(Percentage=Count/n_mod*100, PTM_Residue=paste(Residue,PTM))
overall <- data.frame(PTM=c("Unmodified","Modified"), Residue="All",
  Count=c(n_bg,n_mod), Percentage=c(n_bg/n_tot*100, n_mod/n_tot*100),
  PTM_Residue=c("Unmodified","Modified"), stringsAsFactors=FALSE)
summary_data <- bind_rows(overall, ptm_summary %>% mutate(PTM_Residue=PTM), ptm_res)

write.csv(ptm_sites, file.path(OUTPUT_DIR,"data_ptm_sites.csv"), row.names=FALSE)
write.csv(background, file.path(OUTPUT_DIR,"data_background.csv"), row.names=FALSE)
write.csv(binding,   file.path(OUTPUT_DIR,"data_binding.csv"), row.names=FALSE)
write.csv(summary_data, file.path(OUTPUT_DIR,"data_summary.csv"), row.names=FALSE)
write.csv(validation, file.path(OUTPUT_DIR,"validation_counts.csv"), row.names=FALSE)
saveRDS(list(source_file=XLSX_FILE, len_min=LEN_MIN, len_max=LEN_MAX,
  ptm_colors=PTM_COLORS, ptm_order=PTM_ORDER,
  strong_binder_threshold=STRONG_BINDER_THRESHOLD,
  weak_binder_threshold=WEAK_BINDER_THRESHOLD, n_total=n_tot, n_modified=n_mod,
  n_background=n_bg), file.path(OUTPUT_DIR,"config.rds"))

cat(sprintf("\nBackground %d-%dmer: %d | Modified: %d | Total: %d\n",
            LEN_MIN,LEN_MAX,n_bg,n_mod,n_tot))
cat("\n==== VALIDATION: canonical count vs Moh's embedded summary ====\n")
print(validation[,c("PTM","pep_col","canonical_total","canonical_10_25","embedded_All_sum","SB","WB","NB")], row.names=FALSE)
