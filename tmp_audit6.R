library(readxl)
f <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"

# ============ Deamidation: Binding data misalignment ============
cat("========================================================\n")
cat("Deamidation - Binding data misalignment\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "Deamid. (NQ)"))
cat("The binding extraction uses Peptide...1 (col 1) with 1816 rows\n")
cat("But EL_Rank...9 (col 9) only has 958 rows (aligned with Peptide...7)\n")
cat("This means binding data extracts 1816 peptides but only 958 have EL_Rank\n")
cat("The remaining 858 get NA for binding info.\n\n")

# Check if Peptide...1 and Peptide...7 overlap
pep1 <- df[["Peptide...1"]][!is.na(df[["Peptide...1"]])]
pep7 <- df[["Peptide...7"]][!is.na(df[["Peptide...7"]])]
cat(sprintf("Peptide...1 unique peptides: %d\n", length(unique(pep1))))
cat(sprintf("Peptide...7 unique peptides: %d\n", length(unique(pep7))))
cat(sprintf("Peptide...7 in Peptide...1: %d\n", sum(pep7 %in% pep1)))

# Check the structure: Peptide...7 is the N-specific sub-table with binding
# Peptide...16 is the Q-specific sub-table with binding
# But Peptide...1 is ALL deamidation peptides (both N and Q)
cat("\nThe structure appears to be:\n")
cat("  Cols 1-3: ALL deamidation peptides (N and Q combined)\n")
cat("  Cols 6-12: Deamid @ N sub-table with binding data\n")
cat("  Cols 14-21: Deamid @ Q sub-table with binding data\n")
cat("The binding extraction should use Peptide...7 for N and Peptide...16 for Q\n")

# ============ Acetylation: Binding data check ============
cat("\n========================================================\n")
cat("Acetylation - Binding data alignment\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "Acetyl."))
cat("The binding extraction uses Peptide...2 (col 2) = 455 N-term peptides\n")
cat("But EL_Rank...14 (col 14) = 308 rows (aligned with Peptide...12)\n")
cat("So 455 peptides get extracted but only 308 have binding info.\n\n")
cat("Structure:\n")
cat("  Cols 1-3: N-term acetylation (raw list, 455 peptides)\n")
cat("  Cols 6-9: K-acetylation (54 peptides)\n")
cat("  Cols 11-17: N-term acetylation BINDING table (308 peptides)\n")
cat("  Cols 19-26: K-acetylation BINDING table (33 peptides)\n")
cat("The binding extraction reads raw Peptide...2 (455) with EL_Rank...14 (308) = misaligned\n")

# ============ SUMO: Binding data critical issue ============
cat("\n========================================================\n")
cat("SUMO - Binding data critical issue\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "SUMO"))
cat("The binding extraction pattern 'Peptide' matches:\n")
pep_matches <- grep("Peptide", colnames(df), value = TRUE)
for (pm in pep_matches) {
  idx <- which(colnames(df) == pm)
  cat(sprintf("  [%d] %s: %d non-NA\n", idx, pm, sum(!is.na(df[[pm]]))))
}
cat("\nUsing FIRST match: Peptide...2 = 3 rows (Other Ubiq. sub-table)\n")
cat("But should be using Peptide...30 = 30 rows (binding table)\n")
cat("EL_Rank column 32 has 30 non-NA (aligned with Peptide...30)\n")
cat("Result: binding data only has 3 peptides instead of 30\n")

# ============ bioOxid: Binding data critical issue ============
cat("\n========================================================\n")
cat("bioOxid - Binding data critical issue\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "bioOxid."))
cat("The binding extraction pattern 'Peptide' matches FIRST:\n")
cat("  Peptide...2 = 102 rows (P sub-table)\n")
cat("But should use Peptide...142 = 955 rows (combined table with binding)\n")
cat("Or the pattern should target the combined/binding table specifically\n")
cat("Result: binding data only has 102 peptides (P sub-table) instead of 955\n")

# ============ Methyl: Binding data issue ============
cat("\n========================================================\n")
cat("Methyl. - Binding data issue\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "Methyl."))
cat("The binding extraction pattern 'Peptide' matches FIRST:\n")
cat("  Peptide...2 = 168 rows (CHNQILDEST sub-table)\n")
cat("But the binding tables are:\n")
cat("  Peptide...20 = 140 rows (CHNQILDEST binding)\n")
cat("  Peptide...29 = 93 rows (KR binding)\n")
cat("EL_Rank columns:\n")
cat("  EL_Rank...22 = 140 non-NA (aligned with Peptide...20)\n")
cat("  EL_Rank...31 = 93 non-NA (aligned with Peptide...29)\n")
cat("Result: binding data uses 168 raw peptides with EL_Rank...22 (140)\n")
cat("        This is partially misaligned, and KR sub-table binding is lost\n")

# ============ Dimethyl and Citrullination binding check ============
cat("\n========================================================\n")
cat("Dimethyl. & Citrullination - Binding data alignment\n")
cat("========================================================\n")

df <- suppressMessages(read_excel(f, sheet = "Dimethyl."))
cat("Dimethyl:\n")
cat(sprintf("  Peptide...2 (raw): %d non-NA\n", sum(!is.na(df[["Peptide...2"]]))))
cat(sprintf("  Peptide...11 (binding): %d non-NA\n", sum(!is.na(df[["Peptide...11"]]))))
cat(sprintf("  EL_Rank: %d non-NA\n", sum(!is.na(df[["EL_Rank"]]))))
cat("  extract_binding_data uses Peptide...2 (192) with EL_Rank (137) = MISALIGNED\n")

df <- suppressMessages(read_excel(f, sheet = "Citrullination"))
cat("\nCitrullination:\n")
cat(sprintf("  Peptide...2 (raw): %d non-NA\n", sum(!is.na(df[["Peptide...2"]]))))
cat(sprintf("  Peptide...11 (binding): %d non-NA\n", sum(!is.na(df[["Peptide...11"]]))))
cat(sprintf("  EL_Rank: %d non-NA\n", sum(!is.na(df[["EL_Rank"]]))))
cat("  extract_binding_data uses Peptide...2 (228) with EL_Rank (138) = MISALIGNED\n")
