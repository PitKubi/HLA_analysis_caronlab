library(readxl)
f <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"

# ============ Binding data: remaining sheets ============
cat("========================================================\n")
cat("Binding data - Remaining sheets column check\n")
cat("========================================================\n")

# G-Ubiq: use_last_pep = TRUE
df <- suppressMessages(read_excel(f, sheet = "G-Ubiq."))
pep_matches <- grep("Peptide", colnames(df), value = TRUE)
cat("G-Ubiq (use_last_pep=TRUE):\n")
cat(sprintf("  All Peptide matches: %s\n", paste(pep_matches, collapse=", ")))
last_pep <- pep_matches[length(pep_matches)]
cat(sprintf("  Using LAST: %s -> %d non-NA\n", last_pep, sum(!is.na(df[[last_pep]]))))
el_col <- grep("^EL_Rank|EL_Rank$|^EL.Rank", colnames(df), value = TRUE)
cat(sprintf("  EL_Rank cols: %s\n", paste(el_col, collapse=", ")))
if (length(el_col) > 0) {
  cat(sprintf("  EL_Rank non-NA: %d\n", sum(!is.na(df[[el_col[1]]]))))
}

# Methyl: pep_pattern = "Peptide"
df <- suppressMessages(read_excel(f, sheet = "Methyl."))
pep_matches <- grep("Peptide", colnames(df), value = TRUE)
cat("\nMethyl (pep_pattern=Peptide):\n")
cat(sprintf("  All Peptide matches: %s\n", paste(pep_matches, collapse=", ")))
cat(sprintf("  Using FIRST: %s -> %d non-NA\n", pep_matches[1], sum(!is.na(df[[pep_matches[1]]]))))
el_col <- grep("^EL_Rank|EL_Rank$|^EL.Rank", colnames(df), value = TRUE)
cat(sprintf("  EL_Rank cols: %s\n", paste(el_col, collapse=", ")))
for (col in el_col) cat(sprintf("    %s: %d non-NA\n", col, sum(!is.na(df[[col]]))))

# Dimethyl: pep_pattern = "Peptide"
df <- suppressMessages(read_excel(f, sheet = "Dimethyl."))
pep_matches <- grep("Peptide", colnames(df), value = TRUE)
cat("\nDimethyl (pep_pattern=Peptide):\n")
cat(sprintf("  All Peptide matches: %s\n", paste(pep_matches, collapse=", ")))
cat(sprintf("  Using FIRST: %s -> %d non-NA\n", pep_matches[1], sum(!is.na(df[[pep_matches[1]]]))))
el_col <- grep("^EL_Rank|EL_Rank$|^EL.Rank", colnames(df), value = TRUE)
cat(sprintf("  EL_Rank cols: %s\n", paste(el_col, collapse=", ")))

# Citrullination: pep_pattern = "Peptide"
df <- suppressMessages(read_excel(f, sheet = "Citrullination"))
pep_matches <- grep("Peptide", colnames(df), value = TRUE)
cat("\nCitrullination (pep_pattern=Peptide):\n")
cat(sprintf("  All Peptide matches: %s\n", paste(pep_matches, collapse=", ")))
cat(sprintf("  Using FIRST: %s -> %d non-NA\n", pep_matches[1], sum(!is.na(df[[pep_matches[1]]]))))
el_col <- grep("^EL_Rank|EL_Rank$|^EL.Rank", colnames(df), value = TRUE)
cat(sprintf("  EL_Rank cols: %s\n", paste(el_col, collapse=", ")))

# Oxidation: pep_pattern = "Peptide"
df <- suppressMessages(read_excel(f, sheet = "bioOxid."))
pep_matches <- grep("Peptide", colnames(df), value = TRUE)
# Remove Peptide Length
pep_matches <- pep_matches[!grepl("Length", pep_matches)]
cat("\nbioOxid (pep_pattern=Peptide):\n")
cat(sprintf("  All Peptide matches (excl Length): %s\n", paste(pep_matches, collapse=", ")))
cat(sprintf("  Using FIRST: %s -> %d non-NA\n", pep_matches[1], sum(!is.na(df[[pep_matches[1]]]))))
el_col <- grep("^EL_Rank|EL_Rank$|^EL.Rank", colnames(df), value = TRUE)
cat(sprintf("  EL_Rank cols: %s\n", paste(el_col, collapse=", ")))

# N-glyco: pep_pattern = "^Peptide$"
df <- suppressMessages(read_excel(f, sheet = "N-glyco"))
pep_matches <- grep("^Peptide$", colnames(df), value = TRUE)
cat("\nN-glyco (pep_pattern=^Peptide$):\n")
cat(sprintf("  Peptide matches: %s\n", paste(pep_matches, collapse=", ")))
if (length(pep_matches) > 0) cat(sprintf("  Using: %s -> %d non-NA\n", pep_matches[1], sum(!is.na(df[[pep_matches[1]]]))))
el_col <- grep("^EL_Rank|EL_Rank$|^EL.Rank", colnames(df), value = TRUE)
cat(sprintf("  EL_Rank cols: %s\n", paste(el_col, collapse=", ")))

# ============ Key cross-check: Are binding peptide counts aligned? ============
cat("\n\n========================================================\n")
cat("PTM site counts vs binding data counts\n")
cat("========================================================\n")
cat("This checks if extract_ptm_sites and extract_binding_data\n")
cat("read from different columns/sub-tables\n\n")

cat("Sheet               | PTM Sites Source     | Count | Binding Source       | Count\n")
cat("--------------------|---------------------|-------|---------------------|------\n")

# Phospho: PTM sites from Stripped Peptide sequences (928), binding from Peptide (702)
cat(sprintf("Phospho.            | Stripped Pep seq     | %5d | Peptide (col 23)    | %5d\n", 928, 702))

# Cyst: PTM sites from stripped peptide (2563), binding from Peptide (2355)
cat(sprintf("Cyst.               | stripped peptide     | %5d | Peptide (col 7)     | %5d\n", 2563, 2355))

# Deamid: PTM sites from Peptide...1 (1816), binding from Peptide...1 (1816) with EL_Rank...9 (958)
cat(sprintf("Deamid.             | Peptide...1          | %5d | Peptide...1 + EL_Rank| %5d*\n", 1816, 958))

# Acetyl: PTM sites from Peptide...2 (455 nterm) + Peptide...7 (54 K), binding from Peptide...2 (455) with EL_Rank...14 (308)
cat(sprintf("Acetyl.             | Pep..2+Pep..7       | %5d | Peptide...2+EL..14  | %5d*\n", 455+54, 308))

# GG-Ubiq: PTM from stripped peptide (532), binding from Peptide (426)
cat(sprintf("GG-Ubiq.            | stripped peptide     | %5d | Peptide (col 19)    | %5d\n", 532, 426))

# G-Ubiq: PTM from Peptide...2 (386), binding from last Peptide match
cat(sprintf("G-Ubiq.             | Peptide...2          | %5d | last Peptide match  | %5d\n", 386, 235))

# SUMO: PTM from Pep..2 (3) + Pep..9 (68) = 71, binding from Peptide...2 (3)
cat(sprintf("SUMO                | Pep..2+Pep..9       | %5d | Peptide...2 (WRONG) | %5d\n", 71, 3))

# Carbamid: PTM from Peptide...2 (59), binding from Peptide...11 (48)
cat(sprintf("Carbamid.           | Peptide...2          | %5d | Peptide...11        | %5d\n", 59, 48))

# N-glyco: PTM from ModifiedPeptideSeq (69), binding from Peptide (42)
cat(sprintf("N-glyco             | ModifiedPepSeq      | %5d | Peptide (col 20)    | %5d\n", 69, 42))

# artOxid: PTM from 4 sub-tables (6989), binding from Peptide...36 (6382)
cat(sprintf("artOxid.            | 4 sub-tables        | %5d | Peptide...36        | %5d\n", 6989, 6382))

# bioOxid: PTM from Peptide...142 (955), binding from Peptide...2 (102)?
cat(sprintf("bioOxid.            | Peptide...142        | %5d | Peptide...2 (WRONG) | %5d\n", 955, 102))

# Methyl: PTM from Pep..2 (168) + Pep..8 (121), binding from Peptide...2 (168)
cat(sprintf("Methyl.             | Pep..2+Pep..8       | %5d | Peptide...2 (WRONG) | %5d\n", 289, 168))

# Dimethyl: PTM from Peptide...2 (192), binding from Peptide...2 (192)
cat(sprintf("Dimethyl.           | Peptide...2          | %5d | Peptide...2         | %5d\n", 192, 192))

# Citrullination: PTM from Peptide...2 (228), binding from Peptide...2 (228)
cat(sprintf("Citrullination      | Peptide...2          | %5d | Peptide...2         | %5d\n", 228, 228))

cat("\n* = Mismatched: binding EL_Rank count differs from peptide count\n")
cat("WRONG = Binding data extraction uses wrong peptide column\n")
