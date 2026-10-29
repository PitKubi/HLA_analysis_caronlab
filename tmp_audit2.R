library(readxl)
f <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"

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

# ============ Acetylation verification ============
cat("========================================================\n")
cat("Acetyl. - Verification\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "Acetyl."))
cat(sprintf("N-term sub-table (Peptide...2): %d non-NA\n", sum(!is.na(df[["Peptide...2"]]))))
cat(sprintf("K-acetyl sub-table (Peptide...7): %d non-NA\n", sum(!is.na(df[["Peptide...7"]]))))
cat(sprintf("N-term binding table (Peptide...12): %d non-NA\n", sum(!is.na(df[["Peptide...12"]]))))
cat(sprintf("K-acetyl binding table (Peptide...21): %d non-NA\n", sum(!is.na(df[["Peptide...21"]]))))

# Check what happens: code reads Peptide...2 for N-term (455 peptides), Peptide...7 for K (54 peptides)
# but Peptide...2 has 455 rows - these are ALL the raw N-term acetyl peptides
# Let's verify the K-acetyl residue parsing
k_mods <- df[["Assigned Modifications"]][!is.na(df[["Peptide...7"]]) & !is.na(df[["Assigned Modifications"]])]
cat(sprintf("\nK-acetyl modifications to parse: %d\n", length(k_mods)))
k_residues <- c()
for (mod in k_mods) {
  parsed <- parse_modifications(mod)
  if (!is.null(parsed)) {
    k_only <- parsed[parsed$Residue == "K", , drop=FALSE]
    if (nrow(k_only) > 0) k_residues <- c(k_residues, "K")
    else cat(sprintf("  No K found in: %s\n", mod))
  }
}
cat(sprintf("K-acetyl peptides with K residue: %d\n", length(k_residues)))

# ============ Methylation verification ============
cat("\n========================================================\n")
cat("Methyl. - Residue parsing verification\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "Methyl."))

METHYL_RESIDUES_1 <- c("C", "H", "N", "Q", "I", "L", "D", "E", "S", "T", "N-term")
METHYL_RESIDUES_2 <- c("K", "R")

# Sub-table 1: CHNQILDEST,N-term (Peptide...2, Assigned Modifications...4)
cat("Sub-table 1 (CHNQILDEST,N-term):\n")
st1_count <- 0
st1_residues <- c()
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...2"]][i]
  mod <- df[["Assigned Modifications...4"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed) && nrow(parsed) > 0) {
      kept <- parsed[parsed$Residue %in% METHYL_RESIDUES_1, , drop=FALSE]
      if (nrow(kept) > 0) {
        st1_count <- st1_count + 1
        st1_residues <- c(st1_residues, kept$Residue[1])
      } else {
        cat(sprintf("  DROPPED (no matching residue): %s -> residues: %s\n", mod, paste(parsed$Residue, collapse=",")))
      }
    }
  }
}
cat(sprintf("  Peptides extracted: %d (out of %d rows)\n", st1_count, sum(!is.na(df[["Peptide...2"]]))))
cat("  Residue distribution:\n")
print(table(st1_residues))

# Sub-table 2: KR (Peptide...8, Assigned Modifications...10)
cat("\nSub-table 2 (KR):\n")
st2_count <- 0
st2_residues <- c()
st2_dropped <- 0
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...8"]][i]
  mod <- df[["Assigned Modifications...10"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed) && nrow(parsed) > 0) {
      kept <- parsed[parsed$Residue %in% METHYL_RESIDUES_2, , drop=FALSE]
      if (nrow(kept) > 0) {
        st2_count <- st2_count + 1
        st2_residues <- c(st2_residues, kept$Residue[1])
      } else {
        st2_dropped <- st2_dropped + 1
        cat(sprintf("  DROPPED: %s -> residues: %s\n", mod, paste(parsed$Residue, collapse=",")))
      }
    }
  }
}
cat(sprintf("  Peptides extracted: %d (out of %d rows), dropped: %d\n", st2_count, sum(!is.na(df[["Peptide...8"]])), st2_dropped))
cat("  Residue distribution:\n")
print(table(st2_residues))

cat(sprintf("\nTotal methylation: %d\n", st1_count + st2_count))

# Compare with binding tables
cat(sprintf("CHNQILDEST binding table (Peptide...20): %d non-NA\n", sum(!is.na(df[["Peptide...20"]]))))
cat(sprintf("KR binding table (Peptide...29): %d non-NA\n", sum(!is.na(df[["Peptide...29"]]))))

# ============ Dimethylation verification ============
cat("\n========================================================\n")
cat("Dimethyl. - Residue parsing verification\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "Dimethyl."))

DIMETHYL_RESIDUES <- c("K", "R", "N", "P", "N-term")
dm_count <- 0
dm_residues <- c()
dm_dropped <- 0
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...2"]][i]
  mod <- df[["Assigned Modifications...4"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed) && nrow(parsed) > 0) {
      kept <- parsed[parsed$Residue %in% DIMETHYL_RESIDUES, , drop=FALSE]
      if (nrow(kept) > 0) {
        dm_count <- dm_count + 1
        dm_residues <- c(dm_residues, kept$Residue[1])
      } else {
        dm_dropped <- dm_dropped + 1
        cat(sprintf("  DROPPED: %s -> residues: %s\n", mod, paste(parsed$Residue, collapse=",")))
      }
    }
  }
}
cat(sprintf("Peptides extracted: %d (out of %d rows), dropped: %d\n", dm_count, sum(!is.na(df[["Peptide...2"]])), dm_dropped))
cat("Residue distribution:\n")
print(table(dm_residues))
cat(sprintf("Binding table (Peptide...11): %d non-NA\n", sum(!is.na(df[["Peptide...11"]]))))

# ============ Citrullination verification ============
cat("\n========================================================\n")
cat("Citrullination - Residue parsing verification\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "Citrullination"))
cit_count <- 0
cit_dropped <- 0
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...2"]][i]
  mod <- df[["Assigned Modifications"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed)) {
      r_only <- parsed[parsed$Residue == "R", , drop=FALSE]
      if (nrow(r_only) > 0) {
        cit_count <- cit_count + 1
      } else {
        cit_dropped <- cit_dropped + 1
        cat(sprintf("  DROPPED: %s -> residues: %s\n", mod, paste(parsed$Residue, collapse=",")))
      }
    }
  }
}
cat(sprintf("Peptides extracted: %d (out of %d), dropped: %d\n", cit_count, sum(!is.na(df[["Peptide...2"]])), cit_dropped))
cat(sprintf("Binding table (Peptide...11): %d non-NA\n", sum(!is.na(df[["Peptide...11"]]))))
