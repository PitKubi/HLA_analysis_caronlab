library(readxl)
f <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"

# ============ Critical check: artOxid sub-table residue mismatch ============
cat("========================================================\n")
cat("artOxid. - Cross-contamination between sub-tables\n")
cat("========================================================\n")
cat("The code picks the FIRST residue in (M,W,H,F) from each sub-table's mods.\n")
cat("But sub-tables contain co-occurring modifications on other residues.\n\n")

df <- suppressMessages(read_excel(f, sheet = "artOxid."))

parse_modifications <- function(mod_string) {
  if (is.na(mod_string) || mod_string == "") return(NULL)
  mods <- strsplit(as.character(mod_string), ",\\s*")[[1]]
  results <- lapply(mods, function(m) {
    m <- trimws(m)
    if (grepl("^N-term\\(", m)) {
      mass_val <- as.numeric(gsub("N-term\\(([0-9.]+)\\).*", "\\1", m))
      return(data.frame(Site = 1, Residue = "N-term", Mass = mass_val, stringsAsFactors = FALSE))
    }
    match <- regmatches(m, regexpr("^(\\d+)([A-Z])", m))
    if (length(match) > 0 && nchar(match) >= 2) {
      pos <- as.numeric(gsub("[A-Z]", "", match))
      res <- gsub("[0-9]", "", match)
      mass_val <- as.numeric(gsub(".*\\(([0-9.]+)\\).*", "\\1", m))
      return(data.frame(Site = pos, Residue = res, Mass = mass_val, stringsAsFactors = FALSE))
    }
    return(NULL)
  })
  do.call(rbind, results)
}

# M sub-table: Should be M, but code picks first (M,W,H,F)
cat("M sub-table misclassified residues:\n")
m_wrong <- 0
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...2"]][i]
  mod <- df[["Assigned Modifications...4"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed) && nrow(parsed) > 0) {
      kept <- parsed[parsed$Residue %in% c("M", "W", "H", "F"), , drop=FALSE]
      if (nrow(kept) > 0 && kept$Residue[1] != "M") {
        m_wrong <- m_wrong + 1
        if (m_wrong <= 5) cat(sprintf("  Row %d: mod=%s, first MWHF=%s (should be M)\n", i, mod, kept$Residue[1]))
      }
    }
  }
}
cat(sprintf("  Total misclassified: %d out of 6741\n", m_wrong))

cat("\nW sub-table misclassified residues:\n")
w_wrong <- 0
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...10"]][i]
  mod <- df[["Assigned Modifications...12"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed) && nrow(parsed) > 0) {
      kept <- parsed[parsed$Residue %in% c("M", "W", "H", "F"), , drop=FALSE]
      if (nrow(kept) > 0 && kept$Residue[1] != "W") {
        w_wrong <- w_wrong + 1
        if (w_wrong <= 5) cat(sprintf("  Row %d: mod=%s, first MWHF=%s (should be W)\n", i, mod, kept$Residue[1]))
      }
    }
  }
}
cat(sprintf("  Total misclassified: %d out of 128\n", w_wrong))

cat("\nH sub-table misclassified residues:\n")
h_wrong <- 0
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...18"]][i]
  mod <- df[["Assigned Modifications...20"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed) && nrow(parsed) > 0) {
      kept <- parsed[parsed$Residue %in% c("M", "W", "H", "F"), , drop=FALSE]
      if (nrow(kept) > 0 && kept$Residue[1] != "H") {
        h_wrong <- h_wrong + 1
        if (h_wrong <= 5) cat(sprintf("  Row %d: mod=%s, first MWHF=%s (should be H)\n", i, mod, kept$Residue[1]))
      }
    }
  }
}
cat(sprintf("  Total misclassified: %d out of 42\n", h_wrong))

cat("\nF sub-table misclassified residues:\n")
f_wrong <- 0
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...27"]][i]
  mod <- df[["Assigned Modifications...29"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed) && nrow(parsed) > 0) {
      kept <- parsed[parsed$Residue %in% c("M", "W", "H", "F"), , drop=FALSE]
      if (nrow(kept) > 0 && kept$Residue[1] != "F") {
        f_wrong <- f_wrong + 1
        if (f_wrong <= 5) cat(sprintf("  Row %d: mod=%s, first MWHF=%s (should be F)\n", i, mod, kept$Residue[1]))
      }
    }
  }
}
cat(sprintf("  Total misclassified: %d out of 78\n", f_wrong))

# ============ Check: Methyl KR sub-table also has artifact M(15.9949) ============
cat("\n========================================================\n")
cat("Methyl. KR sub-table - Check for M oxidation co-occurrence\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "Methyl."))
m_cooccur <- 0
for (i in 1:nrow(df)) {
  mod <- df[["Assigned Modifications...10"]][i]
  if (!is.na(mod) && grepl("M\\(15\\.9949\\)", mod)) {
    m_cooccur <- m_cooccur + 1
    if (m_cooccur <= 5) cat(sprintf("  Row %d: %s\n", i, mod))
  }
}
cat(sprintf("KR sub-table rows with M(15.9949) co-occurrence: %d out of %d\n",
  m_cooccur, sum(!is.na(df[["Assigned Modifications...10"]]))))

# Check if any KR mods have dimethyl mass (28.0313) mixed with methyl (14.0156)
cat("\nKR sub-table rows with 28.0313 (dimethyl mass):\n")
dimethyl_in_methyl <- 0
for (i in 1:nrow(df)) {
  mod <- df[["Assigned Modifications...10"]][i]
  if (!is.na(mod) && grepl("28\\.0313", mod)) {
    dimethyl_in_methyl <- dimethyl_in_methyl + 1
    if (dimethyl_in_methyl <= 5) cat(sprintf("  Row %d: %s\n", i, mod))
  }
}
cat(sprintf("Total: %d\n", dimethyl_in_methyl))

# ============ GG-Ubiq residue "n" check ============
cat("\n========================================================\n")
cat("GG-Ubiq. - Lowercase 'n' residue check\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "GG-Ubiq."))
valid <- !is.na(df[["stripped peptide"]]) & !is.na(df[["1st GG site"]])
n_rows <- df[valid & df[["1st Ubiq. residue"]] == "n", ]
if (nrow(n_rows) > 0) {
  cat(sprintf("Found %d rows with lowercase 'n' residue:\n", nrow(n_rows)))
  for (i in 1:min(5, nrow(n_rows))) {
    cat(sprintf("  Peptide: %s, Site: %s, Residue: %s\n",
      n_rows[["stripped peptide"]][i], n_rows[["1st GG site"]][i], n_rows[["1st Ubiq. residue"]][i]))
  }
}

# ============ Binding data extraction check ============
cat("\n========================================================\n")
cat("Binding data - Column identification check\n")
cat("========================================================\n")

# Check what extract_binding_data does for each sheet
sheets_to_check <- list(
  list(sheet="Phospho.", pep_pat="^Peptide$", name="Phosphorylation"),
  list(sheet="Cyst.", pep_pat="^Peptide$", name="Cysteinylation"),
  list(sheet="Deamid. (NQ)", pep_pat="Peptide\\.\\.\\.1", name="Deamidation"),
  list(sheet="Acetyl.", pep_pat="Peptide", name="Acetylation"),
  list(sheet="GG-Ubiq.", pep_pat="^Peptide$", name="Ubiquitination GG"),
  list(sheet="SUMO", pep_pat="Peptide", name="SUMOylation"),
  list(sheet="Carbamid.", pep_pat="Peptide\\.\\.\\.11", name="Carbamidomethylation"),
  list(sheet="artOxid.", pep_pat="Peptide\\.\\.\\.36", name="Artifact Oxidation")
)

for (s in sheets_to_check) {
  df <- suppressMessages(read_excel(f, sheet = s$sheet))
  pep_matches <- grep(s$pep_pat, colnames(df), value = TRUE)
  # Find EL_Rank, Best_Allele, Binder
  el_col <- grep("^EL_Rank|EL_Rank$|^EL.Rank", colnames(df), value = TRUE)
  allele_col <- grep("Best_Allele|Best.Allele", colnames(df), value = TRUE)
  binder_col <- grep("^Binder|Binder$", colnames(df), value = TRUE)

  cat(sprintf("\n%s (%s):\n", s$name, s$sheet))
  cat(sprintf("  Peptide matches for '%s': %s\n", s$pep_pat, paste(pep_matches, collapse=", ")))
  cat(sprintf("  Using (first): %s -> %d non-NA\n", pep_matches[1], sum(!is.na(df[[pep_matches[1]]]))))
  cat(sprintf("  EL_Rank cols: %s\n", paste(el_col, collapse=", ")))
  cat(sprintf("  Best_Allele cols: %s\n", paste(allele_col, collapse=", ")))
  cat(sprintf("  Binder cols: %s\n", paste(binder_col, collapse=", ")))

  # Check: are binding columns aligned with the peptide column?
  if (length(el_col) > 0) {
    pep_idx <- which(colnames(df) == pep_matches[1])
    el_idx <- which(colnames(df) == el_col[1])
    cat(sprintf("  Peptide col index: %d, EL_Rank col index: %d\n", pep_idx, el_idx))
    # Check if they have the same number of non-NA rows
    pep_nonna <- sum(!is.na(df[[pep_matches[1]]]))
    el_nonna <- sum(!is.na(df[[el_col[1]]]))
    cat(sprintf("  Peptide non-NA: %d, EL_Rank non-NA: %d\n", pep_nonna, el_nonna))
  }
}
