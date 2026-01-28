# =============================================================================
# HLA-I PTM Analysis - Centralized Data Loader
# =============================================================================
# This script extracts all data from the source Excel file and saves to CSVs.
# Run this ONCE when data changes, then all figure scripts read from CSVs.
#
# Output files:
#   figure_panels/data_ptm_sites.csv   - PTM site-level data (PTM, Peptide, Site, Residue, Length)
#   figure_panels/data_background.csv  - Background peptides (Peptide, Length, Allele, EL_Rank, Binder)
#   figure_panels/data_binding.csv     - All PTM peptides with binding info
#   figure_panels/data_summary.csv     - Aggregated PTM/Residue counts for overview figures
#   figure_panels/config.rds           - Shared settings (colors, PTM order, etc.)
#
# Usage:
#   source("data_loader.R")  # or Rscript data_loader.R
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)

# =============================================================================
# CONFIGURATION - Edit this section when data source changes
# =============================================================================

# Path to source Excel file
XLSX_FILE <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"

# Output directory
OUTPUT_DIR <- "figure_panels"

# Sheet name mapping - update if sheet names change or new PTMs added
SHEET_MAP <- list(
  background    = "Background (wo PTMs)",
  phospho       = "Phospho.",
  cyst          = "Cyst.",
  deamid        = "Deamid. (NQ)",
  acetyl        = "Acetyl.",
  ubiq_gg       = "GG-Ubiq.",
  ubiq_g        = "G-Ubiq.",
  methyl        = "Methyl.",
  dimethyl      = "Dimethyl.",
  citrullination = "Citrullination",
  oxidation     = "bioOxid.",
  art_oxid      = "artOxid.",
  carbamid      = "Carbamid.",
  sumo          = "SUMO",
  nglyco        = "N-glyco"
)

# PTM display names (for figures)
PTM_NAMES <- c(
  "Phosphorylation",
  "Cysteinylation",
  "Deamidation",
  "Acetylation",
  "Ubiquitination",
  "Methylation",
  "Dimethylation",
  "Citrullination",
  "Oxidation",
  "Artifact Oxidation",
  "Carbamidomethylation",
  "SUMOylation",
  "N-Glycosylation"
)

# Color palette - consistent across all figures
# Ordered by typical abundance (can reorder as needed)
PTM_COLORS <- c(
  "Cysteinylation"       = "#2E7D32",
  "Deamidation"          = "#1565C0",
  "Artifact Oxidation"   = "#FF6F00",
  "Oxidation"            = "#0277BD",
  "Phosphorylation"      = "#6A1B9A",
  "Acetylation"          = "#E65100",
  "Ubiquitination"       = "#C62828",
  "Methylation"          = "#558B2F",
  "Dimethylation"        = "#4527A0",
  "Citrullination"       = "#AD1457",
  "Carbamidomethylation" = "#37474F",
  "SUMOylation"          = "#795548",
  "N-Glycosylation"      = "#00695C"
)

# Binding thresholds (NetMHCpan EL_Rank)
STRONG_BINDER_THRESHOLD <- 0.5
WEAK_BINDER_THRESHOLD <- 2.0

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Parse "Assigned Modifications" column
#' Format examples: "5N(0.9840)", "3M(15.9949), 4N(0.9840)"
#' Returns data.frame with Site and Residue columns, or NULL
parse_modifications <- function(mod_string, peptide = NULL) {
  if (is.null(mod_string) || length(mod_string) == 0) return(NULL)
  mod_string <- as.character(mod_string)[1]
  if (is.na(mod_string) || mod_string == "") return(NULL)

  # Split by comma if multiple modifications

mods <- strsplit(mod_string, ",\\s*")[[1]]

  results <- lapply(mods, function(m) {
    m <- trimws(m)
    # Handle N-terminal modifications: "N-term(mass)"
    if (grepl("^N-term\\(", m)) {
      return(data.frame(Site = 1, Residue = "N-term", stringsAsFactors = FALSE))
    }
    match <- regmatches(m, regexpr("^(\\d+)([A-Z])", m))
    if (length(match) > 0 && nchar(match) >= 2) {
      pos <- as.numeric(gsub("[A-Z]", "", match))
      res <- gsub("[0-9]", "", match)
      return(data.frame(Site = pos, Residue = res, stringsAsFactors = FALSE))
    }
    return(NULL)
  })

  do.call(rbind, results)
}

#' Safely read an Excel sheet, return NULL on error
safe_read_excel <- function(file, sheet) {
  tryCatch({
    suppressMessages(read_excel(file, sheet = sheet))
  }, error = function(e) {
    warning(paste("Could not read sheet:", sheet, "-", e$message))
    return(NULL)
  })
}

#' Extract binding info from a data frame
#' Looks for common column patterns for EL_Rank, Allele, Binder classification
extract_binding_info <- function(df) {
  # Find EL_Rank column (try various patterns)
  el_col <- grep("^EL_Rank|EL_Rank$|^EL.Rank", colnames(df), value = TRUE)[1]

  # Find Best_Allele column
  allele_col <- grep("Best_Allele|Best.Allele", colnames(df), value = TRUE)[1]

  # Find Binder column
  binder_col <- grep("^Binder|Binder$", colnames(df), value = TRUE)[1]

  list(
    el_rank = if (!is.na(el_col)) el_col else NULL,
    allele = if (!is.na(allele_col)) allele_col else NULL,
    binder = if (!is.na(binder_col)) binder_col else NULL
  )
}

# =============================================================================
# DATA EXTRACTION FUNCTIONS
# =============================================================================

#' Extract PTM site data from all sheets
#' Returns data.frame: PTM, Peptide, Site, Residue, Length
extract_ptm_sites <- function(xlsx_file, sheet_map) {
  cat("Extracting PTM site data...\n")
  all_data <- list()

  # --- Phosphorylation ---
  cat("  - Phosphorylation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$phospho)
  if (!is.null(df)) {
    pep_col <- df$`Stripped Peptide sequences`
    site_col <- df$`1st Phosphosite`
    res_col <- df$`1st Phospho residue`
    valid <- !is.na(pep_col) & !is.na(site_col)
    if (sum(valid) > 0) {
      all_data$phospho <- data.frame(
        PTM = "Phosphorylation",
        Peptide = pep_col[valid],
        Site = as.numeric(site_col[valid]),
        Residue = res_col[valid],
        stringsAsFactors = FALSE
      )
    }
  }

  # --- Cysteinylation ---
  cat("  - Cysteinylation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$cyst)
  if (!is.null(df)) {
    valid <- !is.na(df$`stripped peptide`) & !is.na(df$`PTM-site`)
    if (sum(valid) > 0) {
      all_data$cyst <- data.frame(
        PTM = "Cysteinylation",
        Peptide = df$`stripped peptide`[valid],
        Site = as.numeric(df$`PTM-site`[valid]),
        Residue = "C",
        stringsAsFactors = FALSE
      )
    }
  }

  # --- Deamidation ---
  cat("  - Deamidation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$deamid)
  if (!is.null(df)) {
    deamid_list <- list()
    for (i in 1:nrow(df)) {
      if (!is.na(df$`Peptide...1`[i]) && !is.na(df$`Assigned Modifications`[i])) {
        parsed <- parse_modifications(df$`Assigned Modifications`[i])
        if (!is.null(parsed)) {
          parsed <- parsed[parsed$Residue %in% c("N", "Q"), ]
          if (nrow(parsed) > 0) {
            # One entry per peptide (first site only)
            first_site <- parsed[1, , drop = FALSE]
            first_site$PTM <- "Deamidation"
            first_site$Peptide <- df$`Peptide...1`[i]
            deamid_list[[length(deamid_list) + 1]] <- first_site[, c("PTM", "Peptide", "Site", "Residue")]
          }
        }
      }
    }
    if (length(deamid_list) > 0) all_data$deamid <- do.call(rbind, deamid_list)
  }

  # --- Acetylation ---
  cat("  - Acetylation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$acetyl)
  if (!is.null(df)) {
    acetyl_list <- list()
    # N-term acetylation
    if ("Peptide...2" %in% colnames(df)) {
      nterm <- df[!is.na(df$`Peptide...2`), ]
      if (nrow(nterm) > 0) {
        acetyl_list$nterm <- data.frame(
          PTM = "Acetylation",
          Peptide = nterm$`Peptide...2`,
          Site = 1,
          Residue = "N-term",
          stringsAsFactors = FALSE
        )
      }
    }
    # K acetylation
    if ("Peptide...7" %in% colnames(df) && "Assigned Modifications" %in% colnames(df)) {
      k_acetyl <- df[!is.na(df$`Peptide...7`) & !is.na(df$`Assigned Modifications`), ]
      for (i in 1:nrow(k_acetyl)) {
        parsed <- parse_modifications(k_acetyl$`Assigned Modifications`[i])
        if (!is.null(parsed)) {
          parsed <- parsed[parsed$Residue == "K", ]
          if (nrow(parsed) > 0) {
            # One entry per peptide (first K site only)
            first_site <- parsed[1, , drop = FALSE]
            first_site$PTM <- "Acetylation"
            first_site$Peptide <- k_acetyl$`Peptide...7`[i]
            acetyl_list[[length(acetyl_list) + 1]] <- first_site[, c("PTM", "Peptide", "Site", "Residue")]
          }
        }
      }
    }
    if (length(acetyl_list) > 0) all_data$acetyl <- do.call(rbind, acetyl_list)
  }

  # --- Ubiquitination (GG) ---
  cat("  - Ubiquitination (GG)\n")
  df <- safe_read_excel(xlsx_file, sheet_map$ubiq_gg)
  if (!is.null(df)) {
    pep_col <- df$`stripped peptide`
    site_col <- df$`1st GG site`
    res_col <- df$`1st Ubiq. residue`
    valid <- !is.na(pep_col) & !is.na(site_col)
    if (sum(valid) > 0) {
      all_data$ubiq_gg <- data.frame(
        PTM = "Ubiquitination",
        Peptide = pep_col[valid],
        Site = as.numeric(site_col[valid]),
        Residue = res_col[valid],
        Subtype = "GG",
        stringsAsFactors = FALSE
      )
    }
  }

  # --- Ubiquitination (G) ---
  cat("  - Ubiquitination (G)\n")
  df <- safe_read_excel(xlsx_file, sheet_map$ubiq_g)
  if (!is.null(df)) {
    gubiq_list <- list()
    pep_col_name <- "Peptide...2"
    mod_col_name <- grep("^Assigned Modifications", colnames(df), value = TRUE)[1]
    res_col_name <- grep("^Residue", colnames(df), value = TRUE)[1]
    if (pep_col_name %in% colnames(df) && !is.na(mod_col_name)) {
      for (i in 1:nrow(df)) {
        pep_val <- df[[pep_col_name]][i]
        mod_val <- df[[mod_col_name]][i]
        if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
          parsed <- parse_modifications(mod_val)
          if (!is.null(parsed) && nrow(parsed) > 0) {
            # One entry per peptide (first site only)
            first_site <- parsed[1, , drop = FALSE]
            # Use Residue column if available, otherwise use parsed residue
            if (!is.na(res_col_name) && !is.na(df[[res_col_name]][i])) {
              first_site$Residue <- df[[res_col_name]][i]
            }
            first_site$PTM <- "Ubiquitination"
            first_site$Peptide <- pep_val
            first_site$Subtype <- "G"
            gubiq_list[[length(gubiq_list) + 1]] <- first_site[, c("PTM", "Peptide", "Site", "Residue", "Subtype")]
          }
        }
      }
    }
    if (length(gubiq_list) > 0) all_data$ubiq_g <- do.call(rbind, gubiq_list)
  }

  # --- Methylation ---
  cat("  - Methylation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$methyl)
  if (!is.null(df)) {
    methyl_list <- list()
    cols_to_check <- c("Peptide...2", "Peptide...8")
    mod_cols <- c("Assigned Modifications...4", "Assigned Modifications...10")
    for (j in seq_along(cols_to_check)) {
      if (cols_to_check[j] %in% colnames(df) && mod_cols[j] %in% colnames(df)) {
        for (i in 1:nrow(df)) {
          pep_val <- df[[cols_to_check[j]]][i]
          mod_val <- df[[mod_cols[j]]][i]
          if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
            parsed <- parse_modifications(mod_val)
            if (!is.null(parsed) && nrow(parsed) > 0) {
              # One entry per peptide (first site only)
              first_site <- parsed[1, , drop = FALSE]
              first_site$PTM <- "Methylation"
              first_site$Peptide <- pep_val
              methyl_list[[length(methyl_list) + 1]] <- first_site[, c("PTM", "Peptide", "Site", "Residue")]
            }
          }
        }
      }
    }
    if (length(methyl_list) > 0) all_data$methyl <- do.call(rbind, methyl_list)
  }

  # --- Dimethylation ---
  cat("  - Dimethylation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$dimethyl)
  if (!is.null(df)) {
    dimethyl_list <- list()
    if ("Peptide...2" %in% colnames(df) && "Assigned Modifications...4" %in% colnames(df)) {
      for (i in 1:nrow(df)) {
        pep_val <- df$`Peptide...2`[i]
        mod_val <- df$`Assigned Modifications...4`[i]
        if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
          parsed <- parse_modifications(mod_val)
          if (!is.null(parsed) && nrow(parsed) > 0) {
            # One entry per peptide (first site only)
            first_site <- parsed[1, , drop = FALSE]
            first_site$PTM <- "Dimethylation"
            first_site$Peptide <- pep_val
            dimethyl_list[[length(dimethyl_list) + 1]] <- first_site[, c("PTM", "Peptide", "Site", "Residue")]
          }
        }
      }
    }
    if (length(dimethyl_list) > 0) all_data$dimethyl <- do.call(rbind, dimethyl_list)
  }

  # --- Citrullination ---
  cat("  - Citrullination\n")
  df <- safe_read_excel(xlsx_file, sheet_map$citrullination)
  if (!is.null(df)) {
    citrul_list <- list()
    if ("Peptide...2" %in% colnames(df) && "Assigned Modifications" %in% colnames(df)) {
      for (i in 1:nrow(df)) {
        pep_val <- df$`Peptide...2`[i]
        mod_val <- df$`Assigned Modifications`[i]
        if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
          parsed <- parse_modifications(mod_val)
          if (!is.null(parsed)) {
            parsed <- parsed[parsed$Residue == "R", ]
            if (nrow(parsed) > 0) {
              # One entry per peptide (first R site only)
              first_site <- parsed[1, , drop = FALSE]
              first_site$PTM <- "Citrullination"
              first_site$Peptide <- pep_val
              citrul_list[[length(citrul_list) + 1]] <- first_site[, c("PTM", "Peptide", "Site", "Residue")]
            }
          }
        }
      }
    }
    if (length(citrul_list) > 0) all_data$citrul <- do.call(rbind, citrul_list)
  }

  # --- Oxidation ---
  cat("  - Oxidation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$oxidation)
  if (!is.null(df)) {
    oxid_list <- list()
    if ("Peptide...142" %in% colnames(df) && "Assigned Modifications...143" %in% colnames(df)) {
      for (i in 1:nrow(df)) {
        pep_val <- df$`Peptide...142`[i]
        mod_val <- df$`Assigned Modifications...143`[i]
        if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
          parsed <- parse_modifications(mod_val)
          if (!is.null(parsed) && nrow(parsed) > 0) {
            # One entry per peptide (first site only)
            first_site <- parsed[1, , drop = FALSE]
            first_site$PTM <- "Oxidation"
            first_site$Peptide <- pep_val
            oxid_list[[length(oxid_list) + 1]] <- first_site[, c("PTM", "Peptide", "Site", "Residue")]
          }
        }
      }
    }
    if (length(oxid_list) > 0) all_data$oxid <- do.call(rbind, oxid_list)
  }

  # --- SUMOylation ---
  cat("  - SUMOylation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$sumo)
  if (!is.null(df)) {
    sumo_list <- list()
    # Sub-table pairs: (peptide_col, modifications_col)
    sumo_subtables <- list(
      c("Peptide...2", "Assigned Modifications...4"),    # Other Ubiq. sub-table
      c("Peptide...9", "Assigned Modifications...11")    # SUMO sub-table
    )
    for (subtable in sumo_subtables) {
      pep_col_name <- subtable[1]
      mod_col_name <- subtable[2]
      if (pep_col_name %in% colnames(df) && mod_col_name %in% colnames(df)) {
        for (i in 1:nrow(df)) {
          pep_val <- df[[pep_col_name]][i]
          mod_val <- df[[mod_col_name]][i]
          if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
            parsed <- parse_modifications(mod_val)
            if (!is.null(parsed)) {
              parsed <- parsed[parsed$Residue == "K", ]
              if (nrow(parsed) > 0) {
                # One entry per peptide (first K site only)
                first_k <- parsed[1, , drop = FALSE]
                first_k$PTM <- "SUMOylation"
                first_k$Peptide <- pep_val
                sumo_list[[length(sumo_list) + 1]] <- first_k[, c("PTM", "Peptide", "Site", "Residue")]
              }
            }
          }
        }
      }
    }
    if (length(sumo_list) > 0) all_data$sumo <- do.call(rbind, sumo_list)
  }

  # --- N-Glycosylation ---
  cat("  - N-Glycosylation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$nglyco)
  if (!is.null(df)) {
    nglyco_list <- list()
    if ("ModifiedPeptideSequence" %in% colnames(df) && "StrippedSequences" %in% colnames(df)) {
      for (i in 1:nrow(df)) {
        mod_seq <- df$ModifiedPeptideSequence[i]
        pep_val <- df$StrippedSequences[i]
        if (!is.null(mod_seq) && !is.na(mod_seq) && !is.null(pep_val) && !is.na(pep_val)) {
          # Find N(UniMod:xxx) patterns and their positions
          # The position in the modified sequence needs to account for removed UniMod text
          clean_seq <- gsub("\\(UniMod:[0-9]+\\)", "", mod_seq)
          # Find positions of N that had UniMod annotations
          matches <- gregexpr("N\\(UniMod:[0-9]+\\)", mod_seq)[[1]]
          if (matches[1] != -1) {
            # One entry per peptide (first glycosylation site only)
            m <- matches[1]
            prefix <- substr(mod_seq, 1, m - 1)
            prefix_clean <- gsub("\\(UniMod:[0-9]+\\)", "", prefix)
            site <- nchar(prefix_clean) + 1
            nglyco_list[[length(nglyco_list) + 1]] <- data.frame(
              PTM = "N-Glycosylation",
              Peptide = pep_val,
              Site = site,
              Residue = "N",
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    if (length(nglyco_list) > 0) all_data$nglyco <- do.call(rbind, nglyco_list)
  }

  # --- Carbamidomethylation ---
  cat("  - Carbamidomethylation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$carbamid)
  if (!is.null(df)) {
    carbamid_list <- list()
    if ("Peptide...2" %in% colnames(df) && "Assigned Modifications" %in% colnames(df)) {
      for (i in 1:nrow(df)) {
        pep_val <- df$`Peptide...2`[i]
        mod_val <- df$`Assigned Modifications`[i]
        if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
          parsed <- parse_modifications(mod_val)
          if (!is.null(parsed)) {
            # Keep only C residues (carbamidomethylation targets cysteine)
            parsed_c <- parsed[parsed$Residue == "C", ]
            if (nrow(parsed_c) > 0) {
              first_site <- parsed_c[1, , drop = FALSE]
              first_site$PTM <- "Carbamidomethylation"
              first_site$Peptide <- pep_val
              carbamid_list[[length(carbamid_list) + 1]] <- first_site[, c("PTM", "Peptide", "Site", "Residue")]
            }
          }
        }
      }
    }
    if (length(carbamid_list) > 0) all_data$carbamid <- do.call(rbind, carbamid_list)
  }

  # --- Artifact Oxidation ---
  cat("  - Artifact Oxidation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$art_oxid)
  if (!is.null(df)) {
    art_oxid_list <- list()
    # 4 sub-tables by residue: M, W, H, F
    art_oxid_subtables <- list(
      c("Peptide...2",  "Assigned Modifications...4"),   # M
      c("Peptide...10", "Assigned Modifications...12"),  # W
      c("Peptide...18", "Assigned Modifications...20"),  # H
      c("Peptide...27", "Assigned Modifications...29")   # F
    )
    for (subtable in art_oxid_subtables) {
      pep_col_name <- subtable[1]
      mod_col_name <- subtable[2]
      if (pep_col_name %in% colnames(df) && mod_col_name %in% colnames(df)) {
        for (i in 1:nrow(df)) {
          pep_val <- df[[pep_col_name]][i]
          mod_val <- df[[mod_col_name]][i]
          if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
            parsed <- parse_modifications(mod_val)
            if (!is.null(parsed) && nrow(parsed) > 0) {
              # Keep only oxidation-relevant residues (M, W, H, F)
              parsed_oxid <- parsed[parsed$Residue %in% c("M", "W", "H", "F"), ]
              if (nrow(parsed_oxid) > 0) {
                first_site <- parsed_oxid[1, , drop = FALSE]
                first_site$PTM <- "Artifact Oxidation"
                first_site$Peptide <- pep_val
                art_oxid_list[[length(art_oxid_list) + 1]] <- first_site[, c("PTM", "Peptide", "Site", "Residue")]
              }
            }
          }
        }
      }
    }
    if (length(art_oxid_list) > 0) all_data$art_oxid <- do.call(rbind, art_oxid_list)
  }

  # Combine all and add Length
  result <- bind_rows(all_data)
  result <- result[!is.na(result$Site) & !is.na(result$Residue), ]
  result$Length <- nchar(result$Peptide)

  # Ensure Subtype column exists (NA for non-ubiquitination PTMs)
  if (!"Subtype" %in% colnames(result)) result$Subtype <- NA_character_

  return(result)
}

#' Extract background peptides
#' Returns data.frame: Peptide, Length, Allele, EL_Rank, Binder
extract_background <- function(xlsx_file, sheet_map) {
  cat("Extracting background peptides...\n")

  df <- safe_read_excel(xlsx_file, sheet_map$background)
  if (is.null(df)) return(NULL)

  # Find the peptide column
  pep_col <- grep("^Peptide$|Peptide\\.\\.\\.1", colnames(df), value = TRUE)[1]
  if (is.na(pep_col)) {
    warning("Could not find Peptide column in background sheet")
    return(NULL)
  }

  # Get binding info columns
  binding_cols <- extract_binding_info(df)

  result <- data.frame(
    Peptide = df[[pep_col]],
    stringsAsFactors = FALSE
  )
  result$Length <- nchar(result$Peptide)

  # Add binding info if available
  if (!is.null(binding_cols$allele)) {
    result$Allele <- df[[binding_cols$allele]]
  }
  if (!is.null(binding_cols$el_rank)) {
    result$EL_Rank <- as.numeric(df[[binding_cols$el_rank]])
  }
  if (!is.null(binding_cols$binder)) {
    result$Binder <- df[[binding_cols$binder]]
  }

  # Remove rows with NA peptides
  result <- result[!is.na(result$Peptide), ]

  return(result)
}

#' Extract binding data for all PTM peptides
#' Returns data.frame: PTM, Peptide, Length, Allele, EL_Rank, Binder
extract_binding_data <- function(xlsx_file, sheet_map) {
  cat("Extracting binding data for PTM peptides...\n")
  all_data <- list()

  # Configuration for each sheet: sheet_key, PTM_name, peptide_col pattern, binding col patterns
  # use_last_pep: if TRUE, use the last matching Peptide column (for sheets with multiple Peptide columns)
  sheet_configs <- list(
    list(key = "phospho", ptm = "Phosphorylation", pep_pattern = "^Peptide$"),
    list(key = "cyst", ptm = "Cysteinylation", pep_pattern = "^Peptide$"),
    list(key = "deamid", ptm = "Deamidation", pep_pattern = "Peptide\\.\\.\\.1"),
    list(key = "acetyl", ptm = "Acetylation", pep_pattern = "Peptide"),
    list(key = "ubiq_gg", ptm = "Ubiquitination", pep_pattern = "^Peptide$"),
    list(key = "ubiq_g", ptm = "Ubiquitination", pep_pattern = "Peptide", use_last_pep = TRUE),
    list(key = "methyl", ptm = "Methylation", pep_pattern = "Peptide"),
    list(key = "dimethyl", ptm = "Dimethylation", pep_pattern = "Peptide"),
    list(key = "citrullination", ptm = "Citrullination", pep_pattern = "Peptide"),
    list(key = "oxidation", ptm = "Oxidation", pep_pattern = "Peptide"),
    list(key = "sumo", ptm = "SUMOylation", pep_pattern = "Peptide"),
    list(key = "nglyco", ptm = "N-Glycosylation", pep_pattern = "^Peptide$"),
    list(key = "carbamid", ptm = "Carbamidomethylation", pep_pattern = "Peptide\\.\\.\\.11"),
    list(key = "art_oxid", ptm = "Artifact Oxidation", pep_pattern = "Peptide\\.\\.\\.36")
  )

  for (cfg in sheet_configs) {
    cat("  -", cfg$ptm, "\n")
    df <- safe_read_excel(xlsx_file, sheet_map[[cfg$key]])
    if (is.null(df)) next

    # Find peptide column (use last match for sheets with multiple Peptide columns)
    pep_matches <- grep(cfg$pep_pattern, colnames(df), value = TRUE)
    if (length(pep_matches) == 0) next
    pep_col <- if (!is.null(cfg$use_last_pep) && cfg$use_last_pep) {
      pep_matches[length(pep_matches)]
    } else {
      pep_matches[1]
    }
    if (is.na(pep_col)) next

    # Get binding columns
    binding_cols <- extract_binding_info(df)

    # Build result for this PTM
    valid <- !is.na(df[[pep_col]])
    if (sum(valid) == 0) next

    ptm_df <- data.frame(
      PTM = cfg$ptm,
      Peptide = df[[pep_col]][valid],
      stringsAsFactors = FALSE
    )
    ptm_df$Length <- nchar(ptm_df$Peptide)

    if (!is.null(binding_cols$allele)) {
      ptm_df$Allele <- df[[binding_cols$allele]][valid]
    }
    if (!is.null(binding_cols$el_rank)) {
      ptm_df$EL_Rank <- as.numeric(df[[binding_cols$el_rank]][valid])
    }
    if (!is.null(binding_cols$binder)) {
      ptm_df$Binder <- df[[binding_cols$binder]][valid]
    }

    all_data[[cfg$key]] <- ptm_df
  }

  result <- do.call(rbind, all_data)
  result <- result[!is.na(result$Peptide), ]

  # Ensure Binder column exists and classify if not
  if (!"Binder" %in% colnames(result) && "EL_Rank" %in% colnames(result)) {
    result$Binder <- case_when(
      result$EL_Rank < STRONG_BINDER_THRESHOLD ~ "Strong",
      result$EL_Rank < WEAK_BINDER_THRESHOLD ~ "Weak",
      TRUE ~ "Non-binder"
    )
  }

  return(result)
}

#' Generate summary statistics for overview figures
#' Returns data.frame: PTM, Residue, Count, Percentage, PTM_Residue
generate_summary <- function(ptm_sites, background) {
  cat("Generating summary statistics...\n")

  # Filter to 8-14 mers before counting (per collaborator: count all 8-14 mer peptides)
  ptm_sites <- ptm_sites[ptm_sites$Length >= 8 & ptm_sites$Length <= 14, ]
  background <- background[background$Length >= 8 & background$Length <= 14, ]

  # Total counts
  n_background <- nrow(background)
  n_modified <- nrow(ptm_sites)
  n_total <- n_background + n_modified

  # PTM-level summary
  ptm_summary <- ptm_sites %>%
    group_by(PTM) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(
      Percentage = Count / n_modified * 100,
      Residue = "All"
    )

  # PTM + Residue summary
  ptm_residue_summary <- ptm_sites %>%
    group_by(PTM, Residue) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(
      Percentage = Count / n_modified * 100,
      PTM_Residue = paste(Residue, PTM)
    )

  # Overall summary (modified vs unmodified)
  overall <- data.frame(
    PTM = c("Unmodified", "Modified"),
    Residue = c("All", "All"),
    Count = c(n_background, n_modified),
    Percentage = c(n_background / n_total * 100, n_modified / n_total * 100),
    PTM_Residue = c("Unmodified", "Modified"),
    stringsAsFactors = FALSE
  )

  # Combine
  result <- bind_rows(
    overall,
    ptm_summary %>% mutate(PTM_Residue = PTM),
    ptm_residue_summary
  )

  # Add total counts as attributes for reference
  attr(result, "n_total") <- n_total

attr(result, "n_modified") <- n_modified
  attr(result, "n_background") <- n_background

  return(result)
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

run_data_loader <- function() {
  cat("\n")
  cat("=============================================================================\n")
  cat("HLA-I PTM Analysis - Data Loader\n")
  cat("=============================================================================\n")
  cat("Source file:", XLSX_FILE, "\n")
  cat("Output directory:", OUTPUT_DIR, "\n\n")

  # Check source file exists
  if (!file.exists(XLSX_FILE)) {
    stop("Source Excel file not found: ", XLSX_FILE)
  }

  # Create output directory
  dir.create(OUTPUT_DIR, showWarnings = FALSE)

  # Extract all data
  ptm_sites <- extract_ptm_sites(XLSX_FILE, SHEET_MAP)
  background <- extract_background(XLSX_FILE, SHEET_MAP)
  binding_data <- extract_binding_data(XLSX_FILE, SHEET_MAP)
  summary_data <- generate_summary(ptm_sites, background)

  # Save CSVs
  cat("\nSaving output files...\n")

  write.csv(ptm_sites, file.path(OUTPUT_DIR, "data_ptm_sites.csv"), row.names = FALSE)
  cat("  - data_ptm_sites.csv:", nrow(ptm_sites), "records\n")

  write.csv(background, file.path(OUTPUT_DIR, "data_background.csv"), row.names = FALSE)
  cat("  - data_background.csv:", nrow(background), "records\n")

  write.csv(binding_data, file.path(OUTPUT_DIR, "data_binding.csv"), row.names = FALSE)
  cat("  - data_binding.csv:", nrow(binding_data), "records\n")

  write.csv(summary_data, file.path(OUTPUT_DIR, "data_summary.csv"), row.names = FALSE)
  cat("  - data_summary.csv:", nrow(summary_data), "records\n")

  # Save config as RDS (R native format for complex objects)
  config <- list(
    source_file = XLSX_FILE,
    extraction_date = Sys.time(),
    sheet_map = SHEET_MAP,
    ptm_names = PTM_NAMES,
    ptm_colors = PTM_COLORS,
    strong_binder_threshold = STRONG_BINDER_THRESHOLD,
    weak_binder_threshold = WEAK_BINDER_THRESHOLD,
    n_total = attr(summary_data, "n_total"),
    n_modified = attr(summary_data, "n_modified"),
    n_background = attr(summary_data, "n_background")
  )
  saveRDS(config, file.path(OUTPUT_DIR, "config.rds"))
  cat("  - config.rds: shared settings and color palettes\n")

  # Print summary
  cat("\n")
  cat("=============================================================================\n")
  cat("Summary\n")
  cat("=============================================================================\n")
  cat("Total peptides:", config$n_total, "\n")
  cat("  - Unmodified:", config$n_background, sprintf("(%.1f%%)\n", config$n_background/config$n_total*100))
  cat("  - Modified:", config$n_modified, sprintf("(%.1f%%)\n", config$n_modified/config$n_total*100))
  cat("\nPTM breakdown (8-14 mers):\n")

  ptm_counts <- ptm_sites %>%
    filter(Length >= 8 & Length <= 14) %>%
    group_by(PTM) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n))

  for (i in 1:nrow(ptm_counts)) {
    cat(sprintf("  %-20s %5d (%.1f%%)\n",
                ptm_counts$PTM[i],
                ptm_counts$n[i],
                ptm_counts$n[i] / config$n_modified * 100))
  }

  cat("\n=============================================================================\n")
  cat("Data extraction complete!\n")
  cat("=============================================================================\n\n")

  invisible(config)
}

# Run if executed directly (not sourced)
if (!interactive() || identical(environment(), globalenv())) {
  run_data_loader()
}
