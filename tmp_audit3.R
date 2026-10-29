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

# ============ bioOxid - Residue parsing from combined table ============
cat("========================================================\n")
cat("bioOxid. - Combined table residue parsing\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "bioOxid."))

ARTIFACT_OXID_RESIDUES <- c("M", "W", "H", "F")
oxid_count <- 0
oxid_residues <- c()
oxid_dropped <- 0
all_residues_in_combined <- c()
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...142"]][i]
  mod <- df[["Assigned Modifications...143"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed) && nrow(parsed) > 0) {
      all_residues_in_combined <- c(all_residues_in_combined, parsed$Residue)
      kept <- parsed[!(parsed$Residue %in% ARTIFACT_OXID_RESIDUES), , drop=FALSE]
      if (nrow(kept) > 0) {
        oxid_count <- oxid_count + 1
        oxid_residues <- c(oxid_residues, kept$Residue[1])
      } else {
        oxid_dropped <- oxid_dropped + 1
      }
    }
  }
}
cat(sprintf("Combined table rows: %d\n", sum(!is.na(df[["Peptide...142"]]))))
cat(sprintf("Peptides with biological oxidation (kept): %d\n", oxid_count))
cat(sprintf("Peptides with only artifact oxidation (dropped): %d\n", oxid_dropped))
cat("\nAll residues in combined modifications:\n")
print(table(all_residues_in_combined))
cat("\nBiological oxidation residue distribution (first non-artifact):\n")
print(table(oxid_residues))

# Check: does combined 955 include peptides that only have artifact oxidation?
cat(sprintf("\nExpected: 955 - %d dropped = %d kept = %d actual\n", oxid_dropped, 955 - oxid_dropped, oxid_count))

# ============ artOxid - Verification ============
cat("\n========================================================\n")
cat("artOxid. - Sub-table verification\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "artOxid."))

# The code reads sub-tables at Peptide...2/...4, ...10/...12, ...18/...20, ...27/...29
subtables <- list(
  list(pep="Peptide...2", mod="Assigned Modifications...4", label="M"),
  list(pep="Peptide...10", mod="Assigned Modifications...12", label="W"),
  list(pep="Peptide...18", mod="Assigned Modifications...20", label="H"),
  list(pep="Peptide...27", mod="Assigned Modifications...29", label="F")
)

art_total <- 0
for (st in subtables) {
  count <- 0
  residues <- c()
  for (i in 1:nrow(df)) {
    pep <- df[[st$pep]][i]
    mod <- df[[st$mod]][i]
    if (!is.na(pep) && !is.na(mod)) {
      parsed <- parse_modifications(mod)
      if (!is.null(parsed) && nrow(parsed) > 0) {
        kept <- parsed[parsed$Residue %in% c("M", "W", "H", "F"), , drop=FALSE]
        if (nrow(kept) > 0) {
          count <- count + 1
          residues <- c(residues, kept$Residue[1])
        }
      }
    }
  }
  art_total <- art_total + count
  cat(sprintf("%s sub-table (%s): %d raw, %d extracted\n", st$label, st$pep,
    sum(!is.na(df[[st$pep]])), count))
  cat("  Residues: "); print(table(residues))
}
cat(sprintf("\nTotal art. oxid. from sub-tables: %d\n", art_total))
cat(sprintf("Combined table (Peptide...36): %d non-NA\n", sum(!is.na(df[["Peptide...36"]]))))

# ============ SUMO - Verification ============
cat("\n========================================================\n")
cat("SUMO - Sub-table verification\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "SUMO"))

# Sub-table 1: Other Ubiq. (Peptide...2, Assigned Modifications...4) - 3 rows
cat("Sub-table 1 (Other Ubiq.):\n")
st1_count <- 0
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...2"]][i]
  mod <- df[["Assigned Modifications...4"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed)) {
      k_only <- parsed[parsed$Residue == "K", , drop=FALSE]
      if (nrow(k_only) > 0) {
        st1_count <- st1_count + 1
        cat(sprintf("  Row %d: pep=%s, mod=%s, K-site=%d\n", i, pep, mod, k_only$Site[1]))
      } else {
        cat(sprintf("  Row %d DROPPED (no K): %s\n", i, mod))
      }
    }
  }
}
cat(sprintf("  Extracted: %d\n", st1_count))

# Sub-table 2: SUMO (Peptide...9, Assigned Modifications...11) - 68 rows
cat("\nSub-table 2 (SUMO):\n")
st2_count <- 0
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...9"]][i]
  mod <- df[["Assigned Modifications...11"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed)) {
      k_only <- parsed[parsed$Residue == "K", , drop=FALSE]
      if (nrow(k_only) > 0) {
        st2_count <- st2_count + 1
      } else {
        cat(sprintf("  Row %d DROPPED (no K): pep=%s, mod=%s, residues=%s\n",
          i, pep, mod, paste(parsed$Residue, collapse=",")))
      }
    }
  }
}
cat(sprintf("  Extracted: %d\n", st2_count))
cat(sprintf("\nTotal SUMO: %d\n", st1_count + st2_count))
cat(sprintf("All SUMO column 27: %d non-NA\n", sum(!is.na(df[["All SUMO-Peptides (except GG@K)*"]]))))
cat(sprintf("Binding table (Peptide...30): %d non-NA\n", sum(!is.na(df[["Peptide...30"]]))))

# ============ G-Ubiq - Verification ============
cat("\n========================================================\n")
cat("G-Ubiq. - Verification\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "G-Ubiq."))
gu_count <- 0
gu_residues <- c()
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...2"]][i]
  mod <- df[["Assigned Modifications"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed) && nrow(parsed) > 0) {
      gu_count <- gu_count + 1
      # The code uses Residue...5 if available, else parsed residue
      res_col_val <- df[["Residue...5"]][i]
      if (!is.na(res_col_val)) {
        # Residue...5 can have comma-separated values like "K,S"
        first_res <- strsplit(res_col_val, ",")[[1]][1]
        gu_residues <- c(gu_residues, first_res)
      } else {
        gu_residues <- c(gu_residues, parsed$Residue[1])
      }
    }
  }
}
cat(sprintf("Peptides extracted: %d (out of %d)\n", gu_count, sum(!is.na(df[["Peptide...2"]]))))
cat("Residue distribution (from Residue...5 col):\n")
print(table(gu_residues))
cat(sprintf("Binding table (Peptide...13): %d non-NA\n", sum(!is.na(df[["Peptide...13"]]))))

# ============ GG-Ubiq - Verification ============
cat("\n========================================================\n")
cat("GG-Ubiq. - Verification\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "GG-Ubiq."))
gg_count <- sum(!is.na(df[["stripped peptide"]]) & !is.na(df[["1st GG site"]]))
cat(sprintf("Peptides with stripped peptide AND 1st GG site: %d\n", gg_count))
cat(sprintf("Binding table (Peptide, col 19): %d non-NA\n", sum(!is.na(df[["Peptide"]]))))
cat("Residue distribution:\n")
valid <- !is.na(df[["stripped peptide"]]) & !is.na(df[["1st GG site"]])
print(table(df[["1st Ubiq. residue"]][valid], useNA = "ifany"))

# ============ Carbamid - Verification ============
cat("\n========================================================\n")
cat("Carbamid. - Verification\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "Carbamid."))
carb_count <- 0
carb_dropped <- 0
for (i in 1:nrow(df)) {
  pep <- df[["Peptide...2"]][i]
  mod <- df[["Assigned Modifications"]][i]
  if (!is.na(pep) && !is.na(mod)) {
    parsed <- parse_modifications(mod)
    if (!is.null(parsed)) {
      c_only <- parsed[parsed$Residue == "C", , drop=FALSE]
      if (nrow(c_only) > 0) {
        carb_count <- carb_count + 1
      } else {
        carb_dropped <- carb_dropped + 1
        cat(sprintf("  DROPPED (no C): %s -> residues: %s\n", mod, paste(parsed$Residue, collapse=",")))
      }
    }
  }
}
cat(sprintf("Extracted: %d, dropped: %d, total: %d\n", carb_count, carb_dropped, sum(!is.na(df[["Peptide...2"]]))))
cat(sprintf("Binding table (Peptide...11): %d non-NA\n", sum(!is.na(df[["Peptide...11"]]))))

# ============ N-glyco - Verification ============
cat("\n========================================================\n")
cat("N-glyco - Verification\n")
cat("========================================================\n")
df <- suppressMessages(read_excel(f, sheet = "N-glyco"))
ng_count <- 0
ng_dropped <- 0
for (i in 1:nrow(df)) {
  mod_seq <- df[["ModifiedPeptideSequence"]][i]
  pep <- df[["StrippedSequences"]][i]
  if (!is.na(mod_seq) && !is.na(pep)) {
    matches <- gregexpr("N\\(UniMod:[0-9]+\\)", mod_seq)[[1]]
    if (matches[1] != -1) {
      ng_count <- ng_count + 1
    } else {
      ng_dropped <- ng_dropped + 1
      cat(sprintf("  DROPPED (no N(UniMod:...)): %s\n", mod_seq))
    }
  }
}
cat(sprintf("Extracted: %d, dropped: %d, total: %d\n", ng_count, ng_dropped, sum(!is.na(df[["ModifiedPeptideSequence"]]))))
cat(sprintf("Binding table (Peptide, col 20): %d non-NA\n", sum(!is.na(df[["Peptide"]]))))
