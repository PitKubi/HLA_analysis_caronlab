library(readxl)
f <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"

# ============ bioOxid columns 30+ ============
cat("========================================================\n")
cat("bioOxid. - Columns 30-149\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "bioOxid."))
for (i in 30:ncol(df)) {
  nn <- sum(!is.na(df[[i]]))
  cat(sprintf("  [%3d] %s: %d non-NA\n", i, colnames(df)[i], nn))
}

# Sub-table labels
cat("\nSub-table label columns (single uppercase letter):\n")
for (i in seq_along(colnames(df))) {
  cn <- colnames(df)[i]
  if (nchar(cn) == 1 && grepl("^[A-Z]$", cn)) {
    cat(sprintf("  [%3d] %s\n", i, cn))
  }
}

# Total peptides across individual sub-tables vs combined
subtable_pep_cols <- c("Peptide...2","Peptide...11","Peptide...19","Peptide...27",
  "Peptide...36","Peptide...45","Peptide...54","Peptide...64",
  "Peptide...73","Peptide...82","Peptide...90","Peptide...98",
  "Peptide...107","Peptide...117","Peptide...127")
total_st <- 0
for (col in subtable_pep_cols) {
  if (col %in% colnames(df)) {
    n <- sum(!is.na(df[[col]]))
    total_st <- total_st + n
  }
}
combined_n <- sum(!is.na(df[["Peptide...142"]]))
cat(sprintf("\nTotal across individual sub-tables: %d\n", total_st))
cat(sprintf("Combined table (Peptide...142): %d\n", combined_n))
cat(sprintf("Difference: %d\n", combined_n - total_st))

# ============ Deamidation verification ============
cat("\n========================================================\n")
cat("Deamid. (NQ) - Residue parsing verification\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "Deamid. (NQ)"))

parse_modifications <- function(mod_string) {
  if (is.na(mod_string) || mod_string == "") return(NULL)
  mods <- strsplit(as.character(mod_string), ",\\s*")[[1]]
  results <- lapply(mods, function(m) {
    m <- trimws(m)
    if (grepl("^N-term\\(", m)) {
      return(data.frame(Residue = "N-term", stringsAsFactors = FALSE))
    }
    match <- regmatches(m, regexpr("^(\\d+)([A-Z])", m))
    if (length(match) > 0 && nchar(match) >= 2) {
      res <- gsub("[0-9]", "", match)
      return(data.frame(Residue = res, stringsAsFactors = FALSE))
    }
    return(NULL)
  })
  do.call(rbind, results)
}

n_total <- sum(!is.na(df[["Peptide...1"]]) & !is.na(df[["Assigned Modifications"]]))
cat(sprintf("Rows with both Peptide...1 and Assigned Modifications non-NA: %d\n", n_total))

first_residue_counts <- c()
all_residues_in_mods <- c()
n_with_nq <- 0
n_without_nq <- 0
for (i in 1:nrow(df)) {
  if (!is.na(df[["Peptide...1"]][i]) && !is.na(df[["Assigned Modifications"]][i])) {
    parsed <- parse_modifications(df[["Assigned Modifications"]][i])
    if (!is.null(parsed)) {
      all_residues_in_mods <- c(all_residues_in_mods, parsed$Residue)
      nq_only <- parsed[parsed$Residue %in% c("N", "Q"), , drop=FALSE]
      if (nrow(nq_only) > 0) {
        n_with_nq <- n_with_nq + 1
        first_residue_counts <- c(first_residue_counts, nq_only$Residue[1])
      } else {
        n_without_nq <- n_without_nq + 1
      }
    }
  }
}

cat("\nAll residue types in Assigned Modifications:\n")
print(table(all_residues_in_mods))
cat(sprintf("\nPeptides with N or Q (kept by code): %d\n", n_with_nq))
cat(sprintf("Peptides without N or Q (dropped): %d\n", n_without_nq))
cat("\nFirst N/Q residue distribution (what code extracts):\n")
print(table(first_residue_counts))

n_sub <- sum(!is.na(df[["Peptide...7"]]))
q_sub <- sum(!is.na(df[["Peptide...16"]]))
cat(sprintf("\nDeamid. @ N sub-table (Peptide...7): %d\n", n_sub))
cat(sprintf("Deamid. @ Q sub-table (Peptide...16): %d\n", q_sub))
cat(sprintf("Sum of binding sub-tables: %d\n", n_sub + q_sub))
