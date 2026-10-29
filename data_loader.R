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
XLSX_FILE <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01262026_3.xlsx"

# Output directory
OUTPUT_DIR <- "figure_panels"

# Canonical PTM order + colors shared with the HLA-II pipeline
source("shared/ptm_config.R")

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
  nglyco        = "N-glyco",
  nitro         = "Nitrosylation"
)

# PTM display names (for figures) - canonical shared order
PTM_NAMES <- PTM_ORDER

# Binding thresholds (NetMHCpan EL_Rank)
STRONG_BINDER_THRESHOLD <- 0.5
WEAK_BINDER_THRESHOLD <- 2.0

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Parse "Assigned Modifications" column
#' Format examples: "5N(0.9840)", "3M(15.9949), 4N(0.9840)"
#' Returns data.frame with Site, Residue, and Mass columns, or NULL
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

#' Parse a "Localized <PTM>" column value
#' Format examples: "N-term", "7S", "4T, N-term", "4K, 6R"
#' Returns data.frame with Site and Residue for the FIRST localized site
#' (same first-site-per-peptide convention as parse_modifications usage).
parse_localized <- function(loc_string) {
  if (is.null(loc_string) || length(loc_string) == 0) return(NULL)
  loc_string <- as.character(loc_string)[1]
  if (is.na(loc_string) || trimws(loc_string) == "") return(NULL)
  tok <- trimws(strsplit(loc_string, ",\\s*")[[1]][1])
  if (tok == "N-term") {
    return(data.frame(Site = 1, Residue = "N-term", stringsAsFactors = FALSE))
  }
  m <- regmatches(tok, regexpr("^(\\d+)([A-Z])$", tok))
  if (length(m) == 1) {
    return(data.frame(
      Site = as.numeric(gsub("[A-Z]", "", m)),
      Residue = gsub("[0-9]", "", m),
      stringsAsFactors = FALSE
    ))
  }
  NULL
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
  # The sheet's own PTM-site column is computed on the MASS-ANNOTATED peptide
  # string (e.g. "AILSYIM[15.9949]C[119.0041]V" -> site 17 on a 9-mer), so for
  # peptides with co-occurring oxidation it points past the stripped sequence.
  # Every affected peptide carries exactly one C: recompute Site as the C
  # position whenever the sheet site does not land on a C.
  cat("  - Cysteinylation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$cyst)
  if (!is.null(df)) {
    valid <- !is.na(df$`stripped peptide`) & !is.na(df$`PTM-site`)
    if (sum(valid) > 0) {
      pep <- df$`stripped peptide`[valid]
      site <- as.numeric(df$`PTM-site`[valid])
      bad_site <- is.na(site) | site < 1 | site > nchar(pep) |
        substr(pep, site, site) != "C"
      site[bad_site] <- regexpr("C", pep[bad_site], fixed = TRUE)
      keep <- site > 0
      cat(sprintf("    (repaired %d annotated-string site offsets)\n", sum(bad_site & keep)))
      all_data$cyst <- data.frame(
        PTM = "Cysteinylation",
        Peptide = pep[keep],
        Site = site[keep],
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
  # 01262026_3 layout: seven per-residue sub-tables (N-term, K, H, C, S, T, Y)
  # plus a combined binding block (Peptide...43) carrying Moh's own
  # length-filtered union with a "Localized acetytilation (residue)" column
  # (sic, typo is in the source file). We read the combined block: it is the
  # collaborator's canonical acetyl set (381 unique 8-14mer peptides, vs the
  # old N-term+K block sum of 341), avoids double-counting peptides shared
  # between residue sub-tables, and includes the new H/C/S/T/Y acetyl sites.
  ACETYL_RESIDUES <- c("K", "H", "C", "S", "T", "Y", "N-term")
  cat("  - Acetylation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$acetyl)
  if (!is.null(df)) {
    stopifnot(all(c("Peptide...43", "Localized acetytilation (residue)",
                    "Assigned Modifications...45") %in% colnames(df)))
    acetyl_list <- list()
    for (i in 1:nrow(df)) {
      pep_val <- df$`Peptide...43`[i]
      if (is.null(pep_val) || is.na(pep_val)) next
      parsed <- parse_localized(df$`Localized acetytilation (residue)`[i])
      if (is.null(parsed)) {
        # Fallback: parse Assigned Modifications, restricted to acetyl residues
        am <- parse_modifications(df$`Assigned Modifications...45`[i])
        if (!is.null(am)) {
          am <- am[am$Residue %in% ACETYL_RESIDUES, ]
          if (nrow(am) > 0) parsed <- am[1, c("Site", "Residue"), drop = FALSE]
        }
      }
      if (!is.null(parsed) && parsed$Residue %in% ACETYL_RESIDUES) {
        parsed$PTM <- "Acetylation"
        parsed$Peptide <- pep_val
        acetyl_list[[length(acetyl_list) + 1]] <- parsed[, c("PTM", "Peptide", "Site", "Residue")]
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
  # Site/residue come from parsing Assigned Modifications restricted to the
  # ubiquitination target residues: the strings also carry co-occurring
  # artifact M oxidation, and the sheet's Residue column holds comma-joined
  # multi-site strings ("T,K,K") that are not per-site residues.
  UBIQ_RESIDUES <- c("K", "S", "T", "C")
  cat("  - Ubiquitination (G)\n")
  df <- safe_read_excel(xlsx_file, sheet_map$ubiq_g)
  if (!is.null(df)) {
    gubiq_list <- list()
    pep_col_name <- "Peptide...2"
    mod_col_name <- grep("^Assigned Modifications", colnames(df), value = TRUE)[1]
    if (pep_col_name %in% colnames(df) && !is.na(mod_col_name)) {
      for (i in 1:nrow(df)) {
        pep_val <- df[[pep_col_name]][i]
        mod_val <- df[[mod_col_name]][i]
        if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
          parsed <- parse_modifications(mod_val)
          if (!is.null(parsed) && nrow(parsed) > 0) {
            parsed <- parsed[parsed$Residue %in% UBIQ_RESIDUES, ]
            if (nrow(parsed) > 0) {
              first_site <- parsed[1, , drop = FALSE]
              first_site$PTM <- "Ubiquitination"
              first_site$Peptide <- pep_val
              first_site$Subtype <- "G"
              gubiq_list[[length(gubiq_list) + 1]] <- first_site[, c("PTM", "Peptide", "Site", "Residue", "Subtype")]
            }
          }
        }
      }
    }
    if (length(gubiq_list) > 0) all_data$ubiq_g <- do.call(rbind, gubiq_list)
  }

  # --- Methylation ---
  # Sub-tables (01262026_3): "Methyl. on CHNQILDEST" and "Methyl. on KR, N-term"
  # (N-term moved from sub-table 1 to sub-table 2 vs the 01062026 file).
  # The Assigned Modifications column can contain co-occurring artifact M
  # oxidation (mass 15.9949) which must be excluded. Methylation mass is
  # 14.0156 — filter by expected residues per sub-table.
  METHYL_RESIDUES <- list(
    c("C", "H", "N", "Q", "I", "L", "D", "E", "S", "T"),  # sub-table 1
    c("K", "R", "N-term")                                  # sub-table 2
  )
  cat("  - Methylation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$methyl)
  if (!is.null(df)) {
    methyl_list <- list()
    cols_to_check <- c("Peptide...2", "Peptide...8")
    mod_cols <- c("Assigned Modifications...4", "Assigned Modifications...10")
    loc_cols <- c("Localized Methyl....5", "Localized Methyl....11")
    for (j in seq_along(cols_to_check)) {
      if (cols_to_check[j] %in% colnames(df) && mod_cols[j] %in% colnames(df)) {
        allowed_res <- METHYL_RESIDUES[[j]]
        has_loc <- loc_cols[j] %in% colnames(df)
        for (i in 1:nrow(df)) {
          pep_val <- df[[cols_to_check[j]]][i]
          mod_val <- df[[mod_cols[j]]][i]
          if (!is.null(pep_val) && !is.na(pep_val)) {
            # Prefer the curated Localized column (01262026_3); fall back to
            # parsing Assigned Modifications with the per-sub-table whitelist.
            parsed <- if (has_loc) parse_localized(df[[loc_cols[j]]][i]) else NULL
            if (!is.null(parsed) && !parsed$Residue %in% allowed_res) parsed <- NULL
            if (is.null(parsed) && !is.null(mod_val) && !is.na(mod_val)) {
              am <- parse_modifications(mod_val)
              if (!is.null(am) && nrow(am) > 0) {
                am <- am[am$Residue %in% allowed_res, ]
                if (nrow(am) > 0) parsed <- am[1, c("Site", "Residue"), drop = FALSE]
              }
            }
            if (!is.null(parsed)) {
              parsed$PTM <- "Methylation"
              parsed$Peptide <- pep_val
              methyl_list[[length(methyl_list) + 1]] <- parsed[, c("PTM", "Peptide", "Site", "Residue")]
            }
          }
        }
      }
    }
    if (length(methyl_list) > 0) all_data$methyl <- do.call(rbind, methyl_list)
  }

  # --- Dimethylation ---
  # Sub-table: "KRNP, N-term". Exclude co-occurring artifact M oxidation.
  DIMETHYL_RESIDUES <- c("K", "R", "N", "P", "N-term")
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
            # Keep only residues expected for dimethylation
            parsed <- parsed[parsed$Residue %in% DIMETHYL_RESIDUES, ]
            if (nrow(parsed) > 0) {
              first_site <- parsed[1, , drop = FALSE]
              first_site$PTM <- "Dimethylation"
              first_site$Peptide <- pep_val
              dimethyl_list[[length(dimethyl_list) + 1]] <- first_site[, c("PTM", "Peptide", "Site", "Residue")]
            }
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
  # The bioOxid sheet lists biological oxidation on P,I,L,Q,S,T,V,C,D,E,N,Y,G,K,R.
  # The combined table (cols 142-143) includes ALL modifications per peptide,
  # including co-occurring artifact M/W/H/F oxidation. We must exclude artifact
  # oxidation residues and pick the first biological oxidation site.
  ARTIFACT_OXID_RESIDUES <- c("M", "W", "H", "F")
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
            # Exclude artifact oxidation residues (M, W, H, F)
            parsed <- parsed[!parsed$Residue %in% ARTIFACT_OXID_RESIDUES, ]
            if (nrow(parsed) > 0) {
              first_site <- parsed[1, , drop = FALSE]
              first_site$PTM <- "Oxidation"
              first_site$Peptide <- pep_val
              oxid_list[[length(oxid_list) + 1]] <- first_site[, c("PTM", "Peptide", "Site", "Residue")]
            }
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
    # 01262026_3: the binding block duplicates the detail-block headers, so
    # readxl disambiguates them to ...1/...2 (detail) and ...17/...18 (binding)
    nglyco_list <- list()
    stopifnot(all(c("ModifiedPeptideSequence...1", "StrippedSequences...2") %in% colnames(df)))
    if (TRUE) {
      for (i in 1:nrow(df)) {
        mod_seq <- df$`ModifiedPeptideSequence...1`[i]
        pep_val <- df$`StrippedSequences...2`[i]
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

  # --- Nitrosylation ---
  # New sheet in 01262026_3 (n = 5: 4x C, 1x Y). Union block = Peptide...1 +
  # Localized Nitro....4. The Assigned Modifications strings also carry
  # co-occurring artifact oxidation (15.9949), so use the Localized column.
  cat("  - Nitrosylation\n")
  df <- safe_read_excel(xlsx_file, sheet_map$nitro)
  if (!is.null(df)) {
    stopifnot(all(c("Peptide...1", "Localized Nitro....4") %in% colnames(df)))
    nitro_list <- list()
    for (i in 1:nrow(df)) {
      pep_val <- df$`Peptide...1`[i]
      if (is.null(pep_val) || is.na(pep_val)) next
      parsed <- parse_localized(df$`Localized Nitro....4`[i])
      if (!is.null(parsed) && parsed$Residue %in% c("C", "Y")) {
        parsed$PTM <- "Nitrosylation"
        parsed$Peptide <- pep_val
        nitro_list[[length(nitro_list) + 1]] <- parsed[, c("PTM", "Peptide", "Site", "Residue")]
      }
    }
    if (length(nitro_list) > 0) all_data$nitro <- do.call(rbind, nitro_list)
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

  # Invariant: every single-letter residue must sit at its claimed position.
  # Guards against source columns computed on mass-annotated strings
  # (the Cyst. sheet shipped such offsets in 01262026_3).
  single <- nchar(result$Residue) == 1
  mismatch <- single &
    (result$Site > result$Length |
       substr(result$Peptide, result$Site, result$Site) != result$Residue)
  if (any(mismatch)) {
    print(head(result[mismatch, ], 20))
    stop(sum(mismatch), " PTM sites do not match their residue position - fix the extraction")
  }

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
#'
#' Each sheet's binding block is addressed by EXPLICIT column names. The old
#' grep-pattern approach paired the first "Peptide"-matching column (often the
#' left detail sub-table) with the sheet-wide first EL_Rank column (belonging
#' to a different, differently-sorted block), silently misaligning EL ranks
#' for Deamid./Methyl./SUMO/bioOxid./Dimethyl./Citrullination/Acetyl.
#' Explicit specs + stopifnot make any future layout drift fail loudly.
extract_binding_data <- function(xlsx_file, sheet_map) {
  cat("Extracting binding data for PTM peptides...\n")
  all_data <- list()

  # Per sheet: one or more binding blocks, each with exact column names.
  # allele/binder may be NA when the block does not carry them.
  blk <- function(pep, el, allele = NA, binder = NA) {
    list(pep = pep, el = el, allele = allele, binder = binder)
  }
  sheet_configs <- list(
    list(key = "phospho", ptm = "Phosphorylation",
         blocks = list(blk("Peptide...23", "EL_Rank", "Best_Allele...24", "Binder"))),
    list(key = "cyst", ptm = "Cysteinylation",
         blocks = list(blk("Peptide", "EL_Rank", "Best_Allele...8", "Binder"))),
    list(key = "deamid", ptm = "Deamidation",
         blocks = list(blk("Peptide...7", "EL_Rank...9", "Best_Allele...8", "Binder...10"),
                       blk("Peptide...16", "EL_Rank...18", "Best_Allele...17", "Binder...19"))),
    list(key = "acetyl", ptm = "Acetylation",
         blocks = list(blk("Peptide...43", "EL_Rank", "Best_Allele...47", "Binder"))),
    list(key = "ubiq_gg", ptm = "Ubiquitination",
         blocks = list(blk("Peptide", "EL_Rank", "Best_Allele...20", "Binder"))),
    list(key = "ubiq_g", ptm = "Ubiquitination",
         blocks = list(blk("Peptide...13", "EL_Rank", "Best_Allele...14", "Binder"))),
    list(key = "methyl", ptm = "Methylation",
         blocks = list(blk("Peptide...20", "EL_Rank...22", "Best_Allele...21", "Binder...23"),
                       blk("Peptide...29", "EL_Rank...31", "Best_Allele...30", "Binder...32"))),
    list(key = "dimethyl", ptm = "Dimethylation",
         blocks = list(blk("Peptide...11", "EL_Rank", "Best_Allele...12", "Binder"))),
    list(key = "citrullination", ptm = "Citrullination",
         blocks = list(blk("Peptide...11", "EL_Rank", "Best_Allele...12", "Binder"))),
    list(key = "oxidation", ptm = "Oxidation",
         blocks = list(blk("Peptide...142", "EL_Rank", "Best_Allele...145", "Binder"))),
    list(key = "sumo", ptm = "SUMOylation",
         blocks = list(blk("Peptide...31", "EL_Rank", "Best_Allele...32", "Binder"))),
    list(key = "nglyco", ptm = "N-Glycosylation",
         blocks = list(blk("StrippedSequences...18", "EL_Rank", "Best_Allele...20", "Binder"))),
    list(key = "carbamid", ptm = "Carbamidomethylation",
         blocks = list(blk("Peptide...11", "EL_Rank", "Best_Allele...12", "Binder"))),
    list(key = "art_oxid", ptm = "Artifact Oxidation",
         blocks = list(blk("Peptide...36", "EL_Rank", "Best_Allele...39", "Binder"))),
    list(key = "nitro", ptm = "Nitrosylation",
         blocks = list(blk("Peptide...22", "EL_Rank", "Best_Allele...23", "Binder")))
  )

  for (cfg in sheet_configs) {
    cat("  -", cfg$ptm, "(", sheet_map[[cfg$key]], ")\n")
    df <- safe_read_excel(xlsx_file, sheet_map[[cfg$key]])
    if (is.null(df)) next

    for (b in cfg$blocks) {
      required <- c(b$pep, b$el)
      if (!all(required %in% colnames(df))) {
        stop(sprintf("Sheet '%s': expected binding columns not found: %s (layout drift?)",
                     sheet_map[[cfg$key]],
                     paste(setdiff(required, colnames(df)), collapse = ", ")))
      }
      valid <- !is.na(df[[b$pep]])
      if (sum(valid) == 0) next

      ptm_df <- data.frame(
        PTM = cfg$ptm,
        Peptide = df[[b$pep]][valid],
        stringsAsFactors = FALSE
      )
      ptm_df$Length <- nchar(ptm_df$Peptide)
      ptm_df$Allele <- if (!is.na(b$allele) && b$allele %in% colnames(df)) {
        df[[b$allele]][valid]
      } else NA_character_
      ptm_df$EL_Rank <- as.numeric(df[[b$el]][valid])
      ptm_df$Binder <- if (!is.na(b$binder) && b$binder %in% colnames(df)) {
        df[[b$binder]][valid]
      } else NA_character_

      all_data[[length(all_data) + 1]] <- ptm_df
    }
  }

  result <- do.call(rbind, all_data)
  result <- result[!is.na(result$Peptide), ]

  # Classify binder from EL_Rank where the block carried no Binder column
  missing_binder <- is.na(result$Binder) & !is.na(result$EL_Rank)
  result$Binder[missing_binder] <- dplyr::case_when(
    result$EL_Rank[missing_binder] < STRONG_BINDER_THRESHOLD ~ "Strong",
    result$EL_Rank[missing_binder] < WEAK_BINDER_THRESHOLD ~ "Weak",
    TRUE ~ "Non-binder"
  )

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
    ptm_order = PTM_ORDER,
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
