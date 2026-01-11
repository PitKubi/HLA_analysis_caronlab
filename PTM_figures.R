# =============================================================================
# HLA-I PTM Analysis Figures
# =============================================================================
# Required packages
# install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "ggridges", "scales", "patchwork"))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(scales)

# Set file path
xlsx_file <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"

# Create output directory
dir.create("figure_panels", showWarnings = FALSE)

# =============================================================================
# DATA EXTRACTION: Generate all required CSV files from source Excel
# =============================================================================
# This section extracts PTM data from the Excel file and saves to figure_panels/
# Run once to generate data files, then figures can be regenerated anytime

cat("=== Extracting PTM data from Excel file ===\n")

# -----------------------------------------------------------------------------
# Helper function: Parse "Assigned Modifications" column
# Format examples: "5N(0.9840)", "3M(15.9949), 4N(0.9840)"
# -----------------------------------------------------------------------------
parse_modifications <- function(mod_string, peptide) {
  # Handle NULL, NA, empty, or non-character input
  if (is.null(mod_string) || length(mod_string) == 0) return(NULL)
  mod_string <- as.character(mod_string)[1]
  if (is.na(mod_string) || mod_string == "") return(NULL)

  # Split by comma if multiple modifications
  mods <- strsplit(mod_string, ",\\s*")[[1]]

  results <- lapply(mods, function(m) {
    # Extract position and residue: pattern like "5N(" or "3M("
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

# -----------------------------------------------------------------------------
# Extract PTM data for Figure 4A Circos plot
# -----------------------------------------------------------------------------
circos_data <- list()

# 1. Phosphorylation
cat("  Extracting Phosphorylation...\n")
suppressMessages({
  df <- read_excel(xlsx_file, sheet = "Phospho.")
})
phospho <- df[!is.na(df$Peptide) & !is.na(df$`1st Phosphosite`), ]
circos_data$Phosphorylation <- data.frame(
  PTM = "Phosphorylation",
  Peptide = phospho$Peptide,
  Site = as.numeric(phospho$`1st Phosphosite`),
  Residue = phospho$`1st Phospho residue`,
  stringsAsFactors = FALSE
)

# 2. Cysteinylation
cat("  Extracting Cysteinylation...\n")
suppressMessages({
  df <- read_excel(xlsx_file, sheet = "Cyst.")
})
cyst <- df[!is.na(df$Peptide) & !is.na(df$`PTM-site`), ]
circos_data$Cysteinylation <- data.frame(
  PTM = "Cysteinylation",
  Peptide = cyst$Peptide,
  Site = as.numeric(cyst$`PTM-site`),
  Residue = "C",
  stringsAsFactors = FALSE
)

# 3. Deamidation (parse from Assigned Modifications)
cat("  Extracting Deamidation...\n")
suppressMessages({
  df <- read_excel(xlsx_file, sheet = "Deamid. (NQ)")
})
deamid_list <- list()
for (i in 1:nrow(df)) {
  if (!is.na(df$`Peptide...1`[i]) && !is.na(df$`Assigned Modifications`[i])) {
    parsed <- parse_modifications(df$`Assigned Modifications`[i], df$`Peptide...1`[i])
    if (!is.null(parsed)) {
      # Only keep N and Q (deamidation residues)
      parsed <- parsed[parsed$Residue %in% c("N", "Q"), ]
      if (nrow(parsed) > 0) {
        parsed$PTM <- "Deamidation"
        parsed$Peptide <- df$`Peptide...1`[i]
        deamid_list[[length(deamid_list) + 1]] <- parsed[, c("PTM", "Peptide", "Site", "Residue")]
      }
    }
  }
}
circos_data$Deamidation <- do.call(rbind, deamid_list)

# 4. Acetylation
cat("  Extracting Acetylation...\n")
suppressMessages({
  df <- read_excel(xlsx_file, sheet = "Acetyl.")
})
# N-term acetylation: data is in Peptide...2 (site always 1)
nterm <- df[!is.na(df$`Peptide...2`), ]
if (nrow(nterm) > 0) {
  nterm_data <- data.frame(
    PTM = "Acetylation",
    Peptide = nterm$`Peptide...2`,
    Site = 1,
    Residue = "N-term",
    stringsAsFactors = FALSE
  )
} else {
  nterm_data <- NULL
}
# K acetylation: data in Peptide...7 with Assigned Modifications
k_acetyl <- df[!is.na(df$`Peptide...7`) & !is.na(df$`Assigned Modifications`), ]
k_list <- list()
for (i in 1:nrow(k_acetyl)) {
  pep_val <- k_acetyl$`Peptide...7`[i]
  mod_val <- k_acetyl$`Assigned Modifications`[i]
  parsed <- parse_modifications(mod_val, pep_val)
  if (!is.null(parsed)) {
    parsed <- parsed[parsed$Residue == "K", ]
    if (nrow(parsed) > 0) {
      parsed$PTM <- "Acetylation"
      parsed$Peptide <- pep_val
      k_list[[length(k_list) + 1]] <- parsed[, c("PTM", "Peptide", "Site", "Residue")]
    }
  }
}
k_data <- if (length(k_list) > 0) do.call(rbind, k_list) else NULL
circos_data$Acetylation <- rbind(nterm_data, k_data)

# 5. Ubiquitination (GG)
cat("  Extracting Ubiquitination...\n")
suppressMessages({
  df <- read_excel(xlsx_file, sheet = "GG-Ubiq.")
})
ubiq <- df[!is.na(df$Peptide) & !is.na(df$`1st GG site`), ]
circos_data$Ubiquitination <- data.frame(
  PTM = "Ubiquitination",
  Peptide = ubiq$Peptide,
  Site = as.numeric(ubiq$`1st GG site`),
  Residue = ubiq$`1st Ubiq. residue`,
  stringsAsFactors = FALSE
)

# 6. Methylation
cat("  Extracting Methylation...\n")
suppressMessages({
  df <- read_excel(xlsx_file, sheet = "Methyl.")
})
methyl_list <- list()
# Correct column names: Peptide...2/Peptide...8 and Assigned Modifications...4/...10
cols_to_check <- c("Peptide...2", "Peptide...8")
mod_cols <- c("Assigned Modifications...4", "Assigned Modifications...10")
for (j in seq_along(cols_to_check)) {
  if (cols_to_check[j] %in% colnames(df) && mod_cols[j] %in% colnames(df)) {
    for (i in 1:nrow(df)) {
      pep_val <- df[[cols_to_check[j]]][i]
      mod_val <- df[[mod_cols[j]]][i]
      if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
        parsed <- parse_modifications(mod_val, pep_val)
        if (!is.null(parsed) && nrow(parsed) > 0) {
          parsed$PTM <- "Methylation"
          parsed$Peptide <- pep_val
          methyl_list[[length(methyl_list) + 1]] <- parsed[, c("PTM", "Peptide", "Site", "Residue")]
        }
      }
    }
  }
}
circos_data$Methylation <- if (length(methyl_list) > 0) do.call(rbind, methyl_list) else NULL

# 7. Dimethylation
cat("  Extracting Dimethylation...\n")
suppressMessages({
  df <- read_excel(xlsx_file, sheet = "Dimethyl.")
})
dimethyl_list <- list()
# Correct column names: Peptide...2 and Assigned Modifications...4
if ("Peptide...2" %in% colnames(df) && "Assigned Modifications...4" %in% colnames(df)) {
  for (i in 1:nrow(df)) {
    pep_val <- df$`Peptide...2`[i]
    mod_val <- df$`Assigned Modifications...4`[i]
    if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
      parsed <- parse_modifications(mod_val, pep_val)
      if (!is.null(parsed) && nrow(parsed) > 0) {
        parsed$PTM <- "Dimethylation"
        parsed$Peptide <- pep_val
        dimethyl_list[[length(dimethyl_list) + 1]] <- parsed[, c("PTM", "Peptide", "Site", "Residue")]
      }
    }
  }
}
circos_data$Dimethylation <- if (length(dimethyl_list) > 0) do.call(rbind, dimethyl_list) else NULL

# 8. Citrullination
cat("  Extracting Citrullination...\n")
suppressMessages({
  df <- read_excel(xlsx_file, sheet = "Citrullination")
})
citrul_list <- list()
# Correct column names: Peptide...2 and Assigned Modifications
if ("Peptide...2" %in% colnames(df) && "Assigned Modifications" %in% colnames(df)) {
  for (i in 1:nrow(df)) {
    pep_val <- df$`Peptide...2`[i]
    mod_val <- df$`Assigned Modifications`[i]
    if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
      parsed <- parse_modifications(mod_val, pep_val)
      if (!is.null(parsed)) {
        parsed <- parsed[parsed$Residue == "R", ]
        if (nrow(parsed) > 0) {
          parsed$PTM <- "Citrullination"
          parsed$Peptide <- pep_val
          citrul_list[[length(citrul_list) + 1]] <- parsed[, c("PTM", "Peptide", "Site", "Residue")]
        }
      }
    }
  }
}
circos_data$Citrullination <- if (length(citrul_list) > 0) do.call(rbind, citrul_list) else NULL

# 9. Oxidation (biological)
# bioOxid. sheet: 955 rows, 624 unique peptides across multiple oxidized residue types
# Use consolidated column (Peptide...142, Assigned Modifications...143) for all biological oxidation
cat("  Extracting Oxidation...\n")
suppressMessages({
  df <- read_excel(xlsx_file, sheet = "bioOxid.")
})
oxid_list <- list()
# Use consolidated columns for all biological oxidation data
if ("Peptide...142" %in% colnames(df) && "Assigned Modifications...143" %in% colnames(df)) {
  seen_peptide_site <- list()  # Track unique peptide+site combinations
  for (i in 1:nrow(df)) {
    pep_val <- df$`Peptide...142`[i]
    mod_val <- df$`Assigned Modifications...143`[i]
    if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
      parsed <- parse_modifications(mod_val, pep_val)
      if (!is.null(parsed) && nrow(parsed) > 0) {
        parsed$PTM <- "Oxidation"
        parsed$Peptide <- pep_val
        # Deduplicate by peptide+site combination
        for (r in 1:nrow(parsed)) {
          key <- paste(pep_val, parsed$Site[r], parsed$Residue[r], sep = "_")
          if (is.null(seen_peptide_site[[key]])) {
            seen_peptide_site[[key]] <- TRUE
            oxid_list[[length(oxid_list) + 1]] <- parsed[r, c("PTM", "Peptide", "Site", "Residue")]
          }
        }
      }
    }
  }
}
circos_data$Oxidation <- if (length(oxid_list) > 0) do.call(rbind, oxid_list) else NULL

# 10. SUMOylation
cat("  Extracting SUMOylation...\n")
suppressMessages({
  df <- read_excel(xlsx_file, sheet = "SUMO")
})
sumo_list <- list()
# Correct column names: Peptide...2 and Assigned Modifications...4
if ("Peptide...2" %in% colnames(df) && "Assigned Modifications...4" %in% colnames(df)) {
  for (i in 1:nrow(df)) {
    pep_val <- df$`Peptide...2`[i]
    mod_val <- df$`Assigned Modifications...4`[i]
    if (!is.null(pep_val) && !is.na(pep_val) && !is.null(mod_val) && !is.na(mod_val)) {
      parsed <- parse_modifications(mod_val, pep_val)
      if (!is.null(parsed)) {
        parsed <- parsed[parsed$Residue == "K", ]
        if (nrow(parsed) > 0) {
          parsed$PTM <- "SUMOylation"
          parsed$Peptide <- pep_val
          sumo_list[[length(sumo_list) + 1]] <- parsed[, c("PTM", "Peptide", "Site", "Residue")]
        }
      }
    }
  }
}
circos_data$SUMOylation <- if (length(sumo_list) > 0) do.call(rbind, sumo_list) else NULL

# Combine all PTM data and save
all_circos_data <- do.call(rbind, circos_data)
all_circos_data <- all_circos_data[!is.na(all_circos_data$Site) & !is.na(all_circos_data$Residue), ]
write.csv(all_circos_data, "figure_panels/data_figure4A_circos.csv", row.names = FALSE)
cat("  Saved: figure_panels/data_figure4A_circos.csv (", nrow(all_circos_data), " records)\n")

# -----------------------------------------------------------------------------
# Extract PTM + Residue + Allele data for Figure 3B Heatmap
# -----------------------------------------------------------------------------
cat("  Extracting PTM/Residue/Allele data for heatmap...\n")

heatmap_data <- list()

# Function to extract allele data from a sheet
extract_allele_data <- function(df, ptm_name, peptide_col, allele_col, residue_val = NULL) {
  if (!peptide_col %in% colnames(df) || !allele_col %in% colnames(df)) return(NULL)

  valid <- !is.na(df[[peptide_col]]) & !is.na(df[[allele_col]])
  if (sum(valid) == 0) return(NULL)

  alleles <- df[[allele_col]][valid]

  # Count by allele
  allele_counts <- as.data.frame(table(alleles), stringsAsFactors = FALSE)
  colnames(allele_counts) <- c("Allele", "n")
  allele_counts$PTM <- ptm_name
  allele_counts$Residue <- if (!is.null(residue_val)) residue_val else "All"
  allele_counts$PTM_Residue <- paste(allele_counts$Residue, ptm_name)

  return(allele_counts[, c("PTM", "Residue", "Allele", "n", "PTM_Residue")])
}

# Extract from each sheet (simplified - using Best_Allele where available)
sheets_config <- list(
  list(sheet = "Acetyl.", ptm = "acetylation", peptide = "N-term Acetyl....1", allele = "Best_Allele...4", residue = "N-term"),
  list(sheet = "Acetyl.", ptm = "acetylation", peptide = "Acetyl. @ K...6", allele = "Best_Allele...8", residue = "K"),
  list(sheet = "Cyst.", ptm = "cysteinylation", peptide = "Peptide", allele = "Best_Allele", residue = "C"),
  list(sheet = "Deamid. (NQ)", ptm = "deamidation", peptide = "Peptide...1", allele = "Best_Allele...4", residue = "N/Q"),
  list(sheet = "Phospho.", ptm = "phosphorylation", peptide = "Peptide", allele = "Best_Allele", residue = "S/T/Y"),
  list(sheet = "GG-Ubiq.", ptm = "ubiquitination", peptide = "Peptide", allele = "Best_Allele", residue = "K/S/T"),
  list(sheet = "Methyl.", ptm = "methylation", peptide = "Peptide...1", allele = "Best_Allele...5", residue = "K/R"),
  list(sheet = "Dimethyl.", ptm = "dimethylation", peptide = "Peptide...1", allele = "Best_Allele", residue = "R/K"),
  list(sheet = "Citrullination", ptm = "citrullination", peptide = "Peptide...1", allele = "Best_Allele", residue = "R"),
  list(sheet = "bioOxid.", ptm = "oxidation", peptide = "Peptide...1", allele = "Best_Allele", residue = "M/P")
)

for (cfg in sheets_config) {
  suppressMessages({
    df <- tryCatch(read_excel(xlsx_file, sheet = cfg$sheet), error = function(e) NULL)
  })
  if (!is.null(df)) {
    result <- extract_allele_data(df, cfg$ptm, cfg$peptide, cfg$allele, cfg$residue)
    if (!is.null(result)) {
      heatmap_data[[length(heatmap_data) + 1]] <- result
    }
  }
}

all_heatmap_data <- do.call(rbind, heatmap_data)
if (!is.null(all_heatmap_data) && nrow(all_heatmap_data) > 0) {
  write.csv(all_heatmap_data, "figure_panels/data_figure3B_heatmap.csv", row.names = FALSE)
  cat("  Saved: figure_panels/data_figure3B_heatmap.csv (", nrow(all_heatmap_data), " records)\n")
}

cat("=== Data extraction complete ===\n\n")

# =============================================================================
# FIGURE 1A: Horizontal Bar Chart - Overall PTM Composition
# =============================================================================

# Read data and count peptides per sheet
sheet_names <- excel_sheets(xlsx_file)
counts <- sapply(sheet_names, function(s) {
  df <- read_excel(xlsx_file, sheet = s)
  nrow(df)
})
names(counts) <- sheet_names

# Separate background and modified
background <- counts["Background (wo PTMs)"]
modified_counts <- counts[names(counts) != "Background (wo PTMs)"]
modified_counts <- modified_counts[modified_counts > 0]  # Remove empty sheets

total_modified <- sum(modified_counts)
total_all <- background + total_modified

# Create data for bar plot
bar_data <- data.frame(
  category = c("Unmodified", "Modified"),
  count = c(background, total_modified),
  percentage = c(background/total_all * 100, total_modified/total_all * 100)
)

# Create stacked bar data
bar_data$xmin <- c(0, bar_data$percentage[1])
bar_data$xmax <- c(bar_data$percentage[1], 100)
bar_data$label_pos <- (bar_data$xmin + bar_data$xmax) / 2

# Plot Figure 1A
fig1a <- ggplot(bar_data) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = category),
            color = "white", size = 0.5) +
  geom_text(aes(x = label_pos, y = 0.5, 
                label = sprintf("%.1f%%", percentage)),
            fontface = "bold", size = 5,
            color = c("#873600", "white")) +
  scale_fill_manual(values = c("Modified" = "#E74C3C", "Unmodified" = "#FDF2E9"),
                    labels = c(paste0("Modified (", format(total_modified, big.mark=","), ")"),
                              paste0("Unmodified (", format(background, big.mark=","), ")"))) +
  scale_x_continuous(breaks = seq(0, 100, 20), labels = seq(0, 100, 20),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = paste0("HLA-I Immunopeptidome Composition (n = ", 
                      format(total_all, big.mark=","), " peptides)"),
       x = "Percentage of Total Peptides",
       fill = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("Figure1A_BarChart.png", fig1a, width = 10, height = 3, dpi = 300, bg = "white")
ggsave("Figure1A_BarChart.pdf", fig1a, width = 10, height = 3, bg = "white")

cat("Saved Figure1A_BarChart.png and .pdf\n")

# =============================================================================
# FIGURE 1B: Donut Chart - PTM Distribution with Residue Breakdown
# =============================================================================

# Define PTM data with residue breakdown
ptm_data <- list(
  Cysteinylation = list(total = 2563, residues = c(C = 2563)),
  Deamidation = list(total = 1816, residues = c(N = 1178, Q = 638)),
  Phosphorylation = list(total = 960, residues = c(S = 715, T = 198, Y = 47)),
  `Ubiquitination (GG)` = list(total = 541, residues = c(K = 197, S = 197, T = 131, C = 16)),
  Acetylation = list(total = 509, residues = c(`N-term` = 455, K = 54)),
  `Ubiquitination (G)` = list(total = 386, residues = c(S = 179, K = 156, T = 51)),
  Dimethylation = list(total = 192, residues = c(R = 140, P = 35, K = 14, `N-term` = 3)),
  Methylation = list(total = 168, residues = c(K = 85, R = 50, Other = 33)),
  SUMOylation = list(total = 68, residues = c(K = 68))
)

# Warm color palette
ptm_colors <- c(
  "Cysteinylation" = "#8B1A1A",
  "Deamidation" = "#CD3333",
  "Phosphorylation" = "#D35400",
  "Ubiquitination (GG)" = "#E67E22",
  "Acetylation" = "#F39C12",
  "Ubiquitination (G)" = "#F4D03F",
  "Dimethylation" = "#D4AC0D",
  "Methylation" = "#B7950B",
  "SUMOylation" = "#6E4B3A"
)

# Create inner ring data (PTM types)
inner_df <- data.frame(
  ptm = names(ptm_data),
  count = sapply(ptm_data, function(x) sum(x$residues)),
  stringsAsFactors = FALSE
)
inner_df$ptm <- factor(inner_df$ptm, levels = names(ptm_data))
total_ptms <- sum(inner_df$count)
inner_df$percentage <- inner_df$count / total_ptms * 100

# Create outer ring data (residues)
outer_list <- list()
for (ptm_name in names(ptm_data)) {
  residues <- ptm_data[[ptm_name]]$residues
  for (i in seq_along(residues)) {
    outer_list[[length(outer_list) + 1]] <- data.frame(
      ptm = ptm_name,
      residue = names(residues)[i],
      count = residues[i],
      stringsAsFactors = FALSE
    )
  }
}
outer_df <- do.call(rbind, outer_list)
outer_df$ptm <- factor(outer_df$ptm, levels = names(ptm_data))

# Function to lighten colors
lighten <- function(color, factor = 0.3) {
  col <- col2rgb(color) / 255
  col <- col + (1 - col) * factor
  rgb(col[1], col[2], col[3])
}

# Assign colors to outer ring
outer_df$color <- NA
for (ptm_name in names(ptm_data)) {
  idx <- which(outer_df$ptm == ptm_name)
  base_color <- ptm_colors[ptm_name]
  n_res <- length(idx)
  for (i in seq_along(idx)) {
    outer_df$color[idx[i]] <- lighten(base_color, factor = 0.15 * (i - 1))
  }
}

# Calculate positions for donut chart
inner_df$ymax <- cumsum(inner_df$count)
inner_df$ymin <- c(0, head(inner_df$ymax, -1))
inner_df$label_pos <- (inner_df$ymin + inner_df$ymax) / 2

outer_df$ymax <- cumsum(outer_df$count)
outer_df$ymin <- c(0, head(outer_df$ymax, -1))
outer_df$label_pos <- (outer_df$ymin + outer_df$ymax) / 2

# Short labels for inner ring
short_labels <- c(
  "Cysteinylation" = "Cys",
  "Deamidation" = "Deam",
  "Phosphorylation" = "Ph",
  "Ubiquitination (GG)" = "Ub",
  "Acetylation" = "Ac",
  "Ubiquitination (G)" = "gUb",
  "Dimethylation" = "Dim",
  "Methylation" = "Me",
  "SUMOylation" = "Su"
)
inner_df$short_label <- short_labels[as.character(inner_df$ptm)]

# Plot donut chart
fig1b <- ggplot() +
  # Outer ring (residues)
  geom_rect(data = outer_df, 
            aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax, fill = I(color)),
            color = "white", size = 0.5) +
  # Inner ring (PTM types)
  geom_rect(data = inner_df,
            aes(xmin = 1.5, xmax = 3, ymin = ymin, ymax = ymax, fill = ptm),
            color = "white", size = 0.8) +
  # Inner ring labels
  geom_text(data = inner_df[inner_df$percentage > 2.5, ],
            aes(x = 2.25, y = label_pos, label = short_label),
            color = "white", fontface = "bold", size = 4) +
  # Outer ring labels
  geom_text(data = outer_df[outer_df$count / total_ptms * 100 > 1.5, ],
            aes(x = 3.5, y = label_pos, label = residue),
            color = "#333333", fontface = "bold", size = 3) +
  # Center text
  annotate("text", x = 0, y = total_ptms/2, label = "PTM sites", 
           fontface = "bold", size = 5, color = "#333333") +
  annotate("text", x = 0, y = total_ptms/2 - total_ptms*0.08, 
           label = paste0("n = ", format(total_ptms, big.mark = ",")),
           size = 4, color = "#555555", fontface = "italic") +
  scale_fill_manual(values = ptm_colors,
                    labels = paste0(names(ptm_colors), " (", 
                                   format(inner_df$count, big.mark=","), ", ",
                                   sprintf("%.1f%%", inner_df$percentage), ")")) +
  coord_polar(theta = "y") +
  xlim(0, 4.5) +
  labs(title = "HLA-I PTM Distribution by Modification Type and Residue",
       fill = "PTM Types") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("Figure1B_Donut.png", fig1b, width = 12, height = 10, dpi = 300, bg = "white")
ggsave("Figure1B_Donut.pdf", fig1b, width = 12, height = 10, bg = "white")

cat("Saved Figure1B_Donut.png and .pdf\n")

# =============================================================================
# FIGURE 2: Ridgeline Plot - Peptide Length Distribution by PTM Type
# =============================================================================

# Extract length data from each sheet
length_list <- list()

# Unmodified
df <- read_excel(xlsx_file, sheet = "Background (wo PTMs)")
length_list[["Unmodified"]] <- df$Length

# Cysteinylation
df <- read_excel(xlsx_file, sheet = "Cyst.")
length_list[["Cysteinylation"]] <- df$Length

# Deamidation
df <- read_excel(xlsx_file, sheet = "Deamid. (NQ)")
length_list[["Deamidation"]] <- df$`Peptide Length`

# Phosphorylation
df <- read_excel(xlsx_file, sheet = "Phospho.")
length_list[["Phosphorylation"]] <- df$Length

# Ubiquitination (GG)
df <- read_excel(xlsx_file, sheet = "GG-Ubiq.")
length_list[["Ubiquitination (GG)"]] <- df$Length

# Ubiquitination (G)
df <- read_excel(xlsx_file, sheet = "G-Ubiq.")
length_list[["Ubiquitination (G)"]] <- df$`Peptide Length`

# Acetylation
df <- read_excel(xlsx_file, sheet = "Acetyl.")
length_list[["Acetylation"]] <- c(df$Length, df$`Peptide Length`)

# Methylation
df <- read_excel(xlsx_file, sheet = "Methyl.")
length_list[["Methylation"]] <- c(df$`Peptide Length`, df$`Peptide Length.1`)

# Dimethylation
df <- read_excel(xlsx_file, sheet = "Dimethyl.")
length_list[["Dimethylation"]] <- df$`Peptide Length`

# SUMOylation
df <- read_excel(xlsx_file, sheet = "SUMO")
length_list[["SUMOylation"]] <- c(df$`Peptide Length`, df$`Peptide Length.1`)

# Combine into single dataframe
length_df <- do.call(rbind, lapply(names(length_list), function(ptm) {
  data.frame(
    PTM = ptm,
    Length = as.numeric(na.omit(length_list[[ptm]])),
    stringsAsFactors = FALSE
  )
}))

# Order by median length
median_order <- length_df %>%
  group_by(PTM) %>%
  summarise(median_len = median(Length, na.rm = TRUE)) %>%
  arrange(median_len) %>%
  pull(PTM)

length_df$PTM <- factor(length_df$PTM, levels = median_order)

# Color palette including Unmodified
ridge_colors <- c(
  "Unmodified" = "#BEBEBE",
  "Cysteinylation" = "#8B1A1A",
  "Deamidation" = "#CD3333",
  "Phosphorylation" = "#D35400",
  "Ubiquitination (GG)" = "#E67E22",
  "Acetylation" = "#F39C12",
  "Ubiquitination (G)" = "#F4D03F",
  "Dimethylation" = "#D4AC0D",
  "Methylation" = "#B7950B",
  "SUMOylation" = "#6E4B3A"
)

# Plot ridgeline
fig2 <- ggplot(length_df, aes(x = Length, y = PTM, fill = PTM)) +
  geom_density_ridges(scale = 2.5, rel_min_height = 0.01, alpha = 0.85,
                      color = "black", size = 0.3) +
  scale_fill_manual(values = ridge_colors) +
  scale_x_continuous(breaks = seq(8, 24, 2), limits = c(7, 25),
                     expand = c(0, 0)) +
  labs(title = "HLA-I Peptide Length Distribution by PTM Type",
       x = "Peptide length",
       y = "Modification type") +
  theme_ridges() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.text.x = element_text(size = 10),
    legend.position = "none",
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave("Figure2_Length_Distribution.png", fig2, width = 10, height = 12, dpi = 300, bg = "white")
ggsave("Figure2_Length_Distribution.pdf", fig2, width = 10, height = 12, bg = "white")

cat("Saved Figure2_Length_Distribution.png and .pdf\n")

# =============================================================================
# FIGURE 3: Phosphopeptide Binding Affinity (EL_Rank) by HLA Allele
# =============================================================================

# Read phosphopeptide data
phospho_df <- read_excel(xlsx_file, sheet = "Phospho.")

# Extract relevant columns and clean data
phospho_binding <- phospho_df %>%
  select(Peptide, Best_Allele, EL_Rank, Binder) %>%
  filter(!is.na(Best_Allele) & !is.na(EL_Rank)) %>%
  mutate(
    EL_Rank = as.numeric(EL_Rank),
    # Map allele codes to HLA nomenclature
    HLA_Allele = case_when(
      Best_Allele == "A0201" ~ "HLA-A*02:01",
      Best_Allele == "B0702" ~ "HLA-B*07:02",
      Best_Allele == "C0702" ~ "HLA-C*07:02",
      TRUE ~ Best_Allele
    ),
    # Create allele group for ordering
    Allele_Group = case_when(
      grepl("^A", Best_Allele) ~ "HLA-A",
      grepl("^B", Best_Allele) ~ "HLA-B",
      grepl("^C", Best_Allele) ~ "HLA-C",
      TRUE ~ "Other"
    )
  )

# Set factor order for plotting
phospho_binding$HLA_Allele <- factor(
  phospho_binding$HLA_Allele,
  levels = c("HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02")
)

# Calculate summary statistics for labels
allele_summary <- phospho_binding %>%
  group_by(HLA_Allele) %>%
  summarise(
    n = n(),
    median_rank = median(EL_Rank, na.rm = TRUE),
    .groups = "drop"
  )

# Perform pairwise Wilcoxon tests for significance
allele_levels <- levels(phospho_binding$HLA_Allele)
pairwise_tests <- list()
comparisons <- list(
  c("HLA-A*02:01", "HLA-B*07:02"),
  c("HLA-A*02:01", "HLA-C*07:02"),
  c("HLA-B*07:02", "HLA-C*07:02")
)

for (comp in comparisons) {
  group1 <- phospho_binding$EL_Rank[phospho_binding$HLA_Allele == comp[1]]
  group2 <- phospho_binding$EL_Rank[phospho_binding$HLA_Allele == comp[2]]
  test_result <- wilcox.test(group1, group2)
  pairwise_tests[[paste(comp, collapse = " vs ")]] <- test_result$p.value
}

# Create significance labels
get_sig_label <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("n.s.")
}

sig_labels <- sapply(pairwise_tests, get_sig_label)

# Color palette consistent with other figures
allele_colors <- c(
  "HLA-A*02:01" = "#D35400",
  "HLA-B*07:02" = "#8B1A1A",
  "HLA-C*07:02" = "#F39C12"
)

# Define y positions for significance bars (log scale) - staggered to avoid overlap
y_log_max <- max(phospho_binding$EL_Rank, na.rm = TRUE)
y_log_min <- min(phospho_binding$EL_Rank[phospho_binding$EL_Rank > 0], na.rm = TRUE)

# Create significance annotation dataframe - staggered heights to avoid overlap
# Bar 1: A vs B (positions 1-2) - lowest
# Bar 2: B vs C (positions 2-3) - middle
# Bar 3: A vs C (positions 1-3) - highest (widest span)
sig_df <- data.frame(
  xmin = c(1, 2, 1),
  xmax = c(2, 3, 3),
  y = c(y_log_max * 1.8, y_log_max * 4, y_log_max * 10),
  label = sig_labels
)

# Significance legend text
sig_legend <- "Wilcoxon test: *** p<0.001, ** p<0.01, * p<0.05, n.s. p>=0.05"

# Plot Figure 3
fig3 <- ggplot(phospho_binding, aes(x = HLA_Allele, y = EL_Rank, fill = HLA_Allele)) +
  # Box plot
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8, color = "black", linewidth = 0.4) +
  # Individual points with jitter
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.4, color = "black") +
  # Significance bars
  geom_segment(data = sig_df, aes(x = xmin, xend = xmax, y = y, yend = y),
               inherit.aes = FALSE, color = "black", linewidth = 0.4) +
  geom_segment(data = sig_df, aes(x = xmin, xend = xmin, y = y, yend = y * 0.75),
               inherit.aes = FALSE, color = "black", linewidth = 0.4) +
  geom_segment(data = sig_df, aes(x = xmax, xend = xmax, y = y, yend = y * 0.75),
               inherit.aes = FALSE, color = "black", linewidth = 0.4) +
  geom_text(data = sig_df, aes(x = (xmin + xmax) / 2, y = y * 1.3, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  # Scales and colors
  scale_fill_manual(values = allele_colors) +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
    labels = c("0.001", "0.01", "0.1", "1", "10", "100"),
    expand = expansion(mult = c(0.05, 0.3))
  ) +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
  # Labels
  labs(
    title = "Phosphopeptide Binding Affinity by HLA Allele",
    subtitle = "JY Cell Line (HLA-I Immunopeptidome)",
    x = NULL,
    y = "EL Rank (%, log scale)",
    caption = paste0("Wilcoxon test: *** p<0.001, ** p<0.01, * p<0.05, n.s. p>=0.05\n",
                     "n: ", paste(paste0(allele_summary$HLA_Allele, "=", allele_summary$n), collapse = ", "))
  ) +
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#555555"),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 11, margin = margin(t = 5)),
    axis.text.y = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 9, color = "#555555", hjust = 0),
    plot.margin = margin(20, 20, 10, 20)
  )

# Save Figure 3
ggsave("figure_panels/Figure3_Phospho_ELRank.png", fig3,
       width = 8, height = 7, dpi = 300, bg = "white")
ggsave("figure_panels/Figure3_Phospho_ELRank.pdf", fig3,
       width = 8, height = 7, bg = "white")

cat("Saved Figure3_Phospho_ELRank.png and .pdf\n")

# Print statistical summary
cat("\n--- Figure 3 Statistics ---\n")
cat("Sample sizes per allele:\n")
print(allele_summary)
cat("\nPairwise Wilcoxon test p-values:\n")
for (name in names(pairwise_tests)) {
  cat(sprintf("  %s: p = %.4g (%s)\n", name, pairwise_tests[[name]],
              get_sig_label(pairwise_tests[[name]])))
}

# =============================================================================
# FIGURE 3B: PTM Fraction Heatmap by HLA Allele (Residue-Level)
# =============================================================================
# Note: This section uses data extracted at the beginning of this script
# Data file: figure_panels/data_figure3B_heatmap.csv

# Check if CSV exists
if (!file.exists("figure_panels/data_figure3B_heatmap.csv")) {
  cat("WARNING: figure_panels/data_figure3B_heatmap.csv not found.\n")
  cat("Run the data extraction section at the beginning of this script first.\n")
} else {
  # Read the comprehensive PTM data
  ptm_data <- read.csv("figure_panels/data_figure3B_heatmap.csv", stringsAsFactors = FALSE)

  # Remove NA rows
  ptm_data <- ptm_data[!is.na(ptm_data$PTM_Residue) & ptm_data$PTM_Residue != "NA NA", ]
  ptm_data <- ptm_data[ptm_data$Residue != "" & !is.na(ptm_data$Residue), ]

  # Map alleles to HLA nomenclature
  ptm_data$HLA_Allele <- ifelse(ptm_data$Allele == "A0201", "HLA-A*02:01",
                         ifelse(ptm_data$Allele == "B0702", "HLA-B*07:02", "HLA-C*07:02"))

  # Calculate fractions (rows sum to 1)
  totals <- aggregate(n ~ PTM_Residue, data = ptm_data, sum)
  names(totals)[2] <- "total"
  ptm_data <- merge(ptm_data, totals, by = "PTM_Residue")
  ptm_data$fraction <- ptm_data$n / ptm_data$total

  # Filter out rows with very low counts (< 5 total) for statistical reliability
  ptm_data <- ptm_data[ptm_data$total >= 5, ]

  # Clean PTM names for better display
  ptm_data$PTM_clean <- ptm_data$PTM
  ptm_data$PTM_clean <- gsub("bio\\. oxidation", "bio-oxidation", ptm_data$PTM_clean)
  ptm_data$PTM_clean <- gsub("gg-ubiq\\.", "ubiq-GG", ptm_data$PTM_clean)
  ptm_data$PTM_clean <- gsub("g-ubiq\\.", "ubiq-G", ptm_data$PTM_clean)
  ptm_data$PTM_clean <- gsub("n-glycosylation", "N-glycosylation", ptm_data$PTM_clean)

  # Create clean label: Residue (PTM)
  ptm_data$Label <- paste0(ptm_data$Residue, " (", ptm_data$PTM_clean, ")")

  # Order by PTM alphabetically, then by Residue within each PTM
  order_df <- unique(ptm_data[, c("Label", "PTM_clean", "Residue")])
  order_df <- order_df[order(order_df$PTM_clean, order_df$Residue), ]
  ptm_data$Label <- factor(ptm_data$Label, levels = rev(order_df$Label))

  ptm_data$HLA_Allele <- factor(ptm_data$HLA_Allele,
                                 levels = c("HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"))

  # Plot Figure 3B - Compact Heatmap
  fig3b <- ggplot(ptm_data, aes(x = HLA_Allele, y = Label, fill = fraction)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.2f", fraction)), size = 2.2, color = "black") +
    scale_fill_gradient(low = "white", high = "#8B1A1A",
                        limits = c(0, 1), name = "Fraction") +
    scale_x_discrete(position = "top") +
    labs(
      title = "PTM Distribution Across HLA-I Alleles",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text.x = element_text(face = "bold", size = 9),
      axis.text.y = element_text(size = 7, hjust = 1),
      legend.position = "right",
      legend.key.height = unit(0.8, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )

  # Save Figure 3B
  ggsave("figure_panels/Figure3B_PTM_Heatmap.png", fig3b,
         width = 5, height = 8, dpi = 300, bg = "white")
  ggsave("figure_panels/Figure3B_PTM_Heatmap.pdf", fig3b,
         width = 5, height = 8, bg = "white")

  cat("Saved Figure3B_PTM_Heatmap.png and .pdf\n")

  # Print summary
  cat("\n--- Figure 3B Summary ---\n")
  cat("PTM+Residue combinations included:", length(unique(ptm_data$Label)), "\n")
  cat("HLA alleles: A*02:01, B*07:02, C*07:02 (JY cell line)\n")
}

# =============================================================================
# FIGURE 4: Circos Plot - PTM Landscape with Position-Relative AA Composition
# =============================================================================
# Note: Main Figure 4A circos plot is in separate file: Figure4A_Circos.R
# This section generates supplementary Figure 4 panels
# Data file: figure_panels/data_figure4A_circos.csv (generated above)

if (file.exists("figure_panels/data_figure4A_circos.csv")) {

  if (!require(circlize)) {
    cat("Installing circlize package...\n")
    install.packages("circlize", repos = "https://cloud.r-project.org")
  }
  library(circlize)

  # Read data
  ptm_circos <- read.csv("figure_panels/data_figure4A_circos.csv", stringsAsFactors = FALSE)

  # Define PTM order (by count, descending)
  ptm_counts <- table(ptm_circos$PTM)
  ptm_order <- names(sort(ptm_counts, decreasing = TRUE))

  # Color palette matching publication style
  ptm_colors_circos <- c(
    "Oxidation" = "#0D47A1",
    "Deamidation" = "#1E88E5",
    "Cysteinylation" = "#43A047",
    "Methylation" = "#8BC34A",
    "Acetylation" = "#FF9800",
    "Ubiquitination" = "#E53935",
    "Citrullination" = "#EC407A",
    "Dimethylation" = "#8E24AA",
    "Phosphorylation" = "#673AB7",
    "SUMOylation" = "#455A64"
  )

  # AA property colors
  get_aa_color <- function(aa) {
    hydrophobic <- c("A", "V", "I", "L", "M", "F", "W", "P", "G")
    polar <- c("S", "T", "C", "Y", "N", "Q")
    positive <- c("K", "R", "H")
    negative <- c("D", "E")
    if (aa %in% hydrophobic) return("#FFA726")
    if (aa %in% polar) return("#42A5F5")
    if (aa %in% positive) return("#EF5350")
    if (aa %in% negative) return("#66BB6A")
    return("#9E9E9E")
  }

  # Sample representative peptides
  set.seed(42)
  max_peptides_per_ptm <- 30
  sampled_data <- do.call(rbind, lapply(ptm_order, function(ptm) {
    subset_df <- ptm_circos[ptm_circos$PTM == ptm, ]
    if (nrow(subset_df) > max_peptides_per_ptm) {
      subset_df <- subset_df[sample(nrow(subset_df), max_peptides_per_ptm), ]
    }
    subset_df$PTM_factor <- ptm
    return(subset_df)
  }))

  sector_sizes <- table(sampled_data$PTM_factor)[ptm_order]

  # Create circos plot
  png("figure_panels/Figure4_Circos_PTM.png", width = 14, height = 14, units = "in", res = 300)

  circos.clear()
  circos.par(
    gap.after = c(rep(4, length(ptm_order) - 1), 15),
    start.degree = 90,
    track.margin = c(0.005, 0.005),
    cell.padding = c(0.005, 0, 0.005, 0)
  )

  circos.initialize(
    factors = factor(ptm_order, levels = ptm_order),
    xlim = cbind(rep(0, length(ptm_order)), as.numeric(sector_sizes))
  )

  # Track 1: PTM labels
  circos.track(
    ylim = c(0, 1), track.height = 0.06,
    bg.col = ptm_colors_circos[ptm_order], bg.border = NA,
    panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, 0.5, CELL_META$sector.index,
                  facing = "bending.inside", niceFacing = TRUE,
                  cex = 0.8, col = "white", font = 2)
    }
  )

  # Track 2: Peptide sequences as colored tiles
  circos.track(
    ylim = c(0, 15), track.height = 0.35,
    bg.col = "white", bg.border = "gray90",
    panel.fun = function(x, y) {
      sector.name <- CELL_META$sector.index
      sector_peptides <- sampled_data[sampled_data$PTM_factor == sector.name, ]
      n_pep <- nrow(sector_peptides)
      if (n_pep > 0) {
        x_positions <- seq(0.5, n_pep - 0.5, length.out = n_pep)
        for (i in 1:n_pep) {
          pep <- sector_peptides$Peptide[i]
          site <- sector_peptides$Site[i]
          y_start <- 7.5 - (site - 1)
          for (j in 1:nchar(pep)) {
            aa <- substr(pep, j, j)
            y_pos <- y_start + (j - 1)
            if (y_pos >= 0 && y_pos <= 15) {
              col <- if (j == site) ptm_colors_circos[sector.name] else get_aa_color(aa)
              circos.rect(x_positions[i] - 0.4, y_pos,
                          x_positions[i] + 0.4, y_pos + 0.8,
                          col = col, border = NA)
            }
          }
        }
      }
    }
  )

  # Track 3: Position numbers
  circos.track(
    ylim = c(-8, 8), track.height = 0.08,
    bg.col = "#F5F5F5", bg.border = "gray80",
    panel.fun = function(x, y) {
      for (pos in seq(-6, 6, by = 2)) {
        circos.text(CELL_META$xcenter, pos, as.character(pos),
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 0.45, col = "gray30")
      }
      circos.segments(CELL_META$xlim[1], 0, CELL_META$xlim[2], 0, col = "red", lwd = 1)
    }
  )

  # Track 4: AA frequency bars
  circos.track(
    ylim = c(-7.5, 7.5), track.height = 0.2,
    bg.col = "white", bg.border = "gray80",
    panel.fun = function(x, y) {
      sector.name <- CELL_META$sector.index
      ptm_pos <- position_freq[position_freq$PTM == sector.name, ]
      for (pos in -5:5) {
        pos_data <- ptm_pos[ptm_pos$Position == pos, ]
        if (nrow(pos_data) > 0) {
          pos_data <- pos_data[order(-pos_data$Frequency), ]
          x_start <- CELL_META$xlim[1]
          bar_width <- CELL_META$xlim[2] - CELL_META$xlim[1]
          cum_width <- 0
          for (j in 1:min(3, nrow(pos_data))) {
            freq <- pos_data$Frequency[j]
            aa <- pos_data$AA[j]
            if (freq > 0.05) {
              width <- freq * bar_width
              circos.rect(x_start + cum_width, pos - 0.4,
                          x_start + cum_width + width, pos + 0.4,
                          col = get_aa_color(aa), border = NA)
              if (freq > 0.15) {
                circos.text(x_start + cum_width + width/2, pos, aa,
                            facing = "bending.inside", niceFacing = TRUE,
                            cex = 0.35, col = "white", font = 2)
              }
              cum_width <- cum_width + width
            }
          }
        }
      }
    }
  )

  # Track 5: Residue distribution
  circos.track(
    ylim = c(0, 1), track.height = 0.08,
    bg.col = "gray95", bg.border = "gray80",
    panel.fun = function(x, y) {
      sector.name <- CELL_META$sector.index
      ptm_peps <- ptm_circos[ptm_circos$PTM == sector.name, ]
      res_counts <- sort(table(ptm_peps$Residue), decreasing = TRUE)
      total <- sum(res_counts)
      if (length(res_counts) > 0) {
        x_start <- 0
        xlim_range <- CELL_META$xlim[2] - CELL_META$xlim[1]
        for (i in seq_along(res_counts)) {
          width <- as.numeric(res_counts[i]) / total * xlim_range
          res <- names(res_counts)[i]
          circos.rect(x_start, 0.1, x_start + width, 0.9,
                      col = get_aa_color(res), border = "white")
          if (width > xlim_range * 0.08) {
            circos.text(x_start + width/2, 0.5, res,
                        facing = "bending.inside", niceFacing = TRUE,
                        cex = 0.4, col = "white", font = 2)
          }
          x_start <- x_start + width
        }
      }
    }
  )

  title("HLA-I PTM Landscape\nPosition-Relative Amino Acid Composition",
        line = -2, cex.main = 1.3, font.main = 2)

  legend("bottomright",
         legend = c(names(ptm_colors_circos), "", "AA Properties:",
                    "Hydrophobic", "Polar", "Positive", "Negative"),
         fill = c(ptm_colors_circos, NA, NA, "#FFA726", "#42A5F5", "#EF5350", "#66BB6A"),
         border = c(rep("white", length(ptm_colors_circos)), NA, NA, "white", "white", "white", "white"),
         title = "Legend", cex = 0.7, bty = "n", ncol = 2)

  mtext(paste0("Total PTM peptides: ", format(nrow(ptm_circos), big.mark = ",")),
        side = 1, line = -2, cex = 0.9)

  dev.off()
  circos.clear()

  cat("Saved Figure4_Circos_PTM.png\n")

  # Figure 4D: S phosphorylation position frequency
  s_phospho <- ptm_circos[ptm_circos$PTM == "Phosphorylation" & ptm_circos$Residue == "S", ]

  positions <- -5:5
  p_freq <- sapply(positions, function(pos) {
    aa_at_pos <- sapply(1:nrow(s_phospho), function(i) {
      pep <- s_phospho$Peptide[i]
      site <- s_phospho$Site[i]
      target <- site + pos
      if (target >= 1 && target <= nchar(pep)) return(substr(pep, target, target))
      return(NA)
    })
    aa_at_pos <- aa_at_pos[!is.na(aa_at_pos)]
    freq_table <- table(aa_at_pos)
    if ("P" %in% names(freq_table)) return(freq_table["P"] / length(aa_at_pos))
    return(0)
  })

  png("figure_panels/Figure4D_S_phospho_position.png", width = 4, height = 4, units = "in", res = 300)
  barplot(p_freq, names.arg = positions, col = "#673AB7", border = NA,
          main = "S phosphorylation\n(Proline frequency)",
          xlab = "Position relative to pS", ylab = "Frequency",
          ylim = c(0, max(p_freq) * 1.2))
  dev.off()

  cat("Saved Figure4D_S_phospho_position.png\n")

} else {
  cat("WARNING: figure_panels/data_figure4A_circos.csv not found.\n")
  cat("Run the data extraction section at the top of this script first.\n")
}

# =============================================================================
# FIGURE 5A: % Strong Binders by PTM+Residue and Allele
# =============================================================================
# Heatmap showing % strong binders, deconvoluted by amino acid residue
# Sorted by PTM type, then by residue within each PTM

cat("\n=== Figure 5A: % Strong Binders by PTM+Residue ===\n")

if (file.exists("figure_panels/data_figure4A_circos.csv")) {

  # Load parsed PTM data with residue info
  ptm_data_fig5 <- read.csv("figure_panels/data_figure4A_circos.csv", stringsAsFactors = FALSE)
  ptm_data_fig5$PTM_Residue <- paste(ptm_data_fig5$Residue, tolower(ptm_data_fig5$PTM))

  # Function to extract binding data from sheets
  get_binding_for_sheet <- function(sheet_name) {
    suppressMessages({
      df <- read_excel(xlsx_file, sheet = sheet_name)
    })

    pep_cols <- colnames(df)[grepl("^Peptide", colnames(df))]
    binder_col <- colnames(df)[grepl("^Binder", colnames(df))][1]
    allele_col <- colnames(df)[grepl("^Best_Allele", colnames(df))][1]

    if (length(pep_cols) == 0 || is.na(binder_col)) return(NULL)

    for (pep_col in pep_cols) {
      result <- data.frame(
        Peptide = df[[pep_col]],
        Binder = df[[binder_col]],
        Allele = df[[allele_col]],
        stringsAsFactors = FALSE
      )
      result <- result[!is.na(result$Peptide) & !is.na(result$Binder), ]
      if (nrow(result) > 0) return(result)
    }
    return(NULL)
  }

  # Get binding data from all sheets
  sheets_fig5 <- c("Phospho.", "Acetyl.", "Cyst.", "Methyl.", "Dimethyl.",
                   "Deamid. (NQ)", "bioOxid.", "Citrullination", "G-Ubiq.")

  cat("Loading binding data for Figure 5...\n")
  all_binding <- list()
  for (sheet in sheets_fig5) {
    binding <- get_binding_for_sheet(sheet)
    if (!is.null(binding)) {
      all_binding[[sheet]] <- binding
    }
  }

  binding_combined <- do.call(rbind, all_binding)

  # Merge with PTM data
  suppressWarnings({
    merged_fig5 <- ptm_data_fig5 %>%
      inner_join(binding_combined, by = "Peptide", relationship = "many-to-many")
  })

  cat("Merged records:", nrow(merged_fig5), "\n")

  # Calculate stats for each PTM+Residue x Allele
  stats_fig5 <- merged_fig5 %>%
    filter(!is.na(Allele)) %>%
    group_by(PTM, Residue, PTM_Residue, Allele) %>%
    summarise(
      n = n(),
      n_strong = sum(Binder == "Strong"),
      pct_strong = round(100 * n_strong / n, 1),
      .groups = "drop"
    ) %>%
    filter(n >= 10)

  # Create heatmap matrix
  heatmap_df_fig5 <- stats_fig5 %>%
    select(PTM, Residue, PTM_Residue, Allele, pct_strong) %>%
    pivot_wider(names_from = Allele, values_from = pct_strong) %>%
    arrange(PTM, Residue) %>%
    as.data.frame()

  rownames(heatmap_df_fig5) <- heatmap_df_fig5$PTM_Residue
  heatmap_matrix_fig5 <- heatmap_df_fig5 %>% select(any_of(c("A0201", "B0702", "C0702")))

  cat("PTM+Residue combinations:", nrow(heatmap_matrix_fig5), "\n")

  # Color palette
  colors_fig5 <- colorRampPalette(c("#FFFFFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C"))(100)

  # Row annotation
  row_annotation_fig5 <- data.frame(PTM = heatmap_df_fig5$PTM)
  rownames(row_annotation_fig5) <- rownames(heatmap_matrix_fig5)

  ptm_types_fig5 <- unique(heatmap_df_fig5$PTM)
  ptm_colors_fig5 <- setNames(
    c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
      "#FFFF33", "#A65628", "#F781BF", "#999999")[1:length(ptm_types_fig5)],
    ptm_types_fig5
  )
  ann_colors_fig5 <- list(PTM = ptm_colors_fig5)

  # Calculate gaps between PTM groups
  gaps_fig5 <- cumsum(rle(heatmap_df_fig5$PTM)$lengths)
  gaps_fig5 <- gaps_fig5[-length(gaps_fig5)]

  # Figure height
  fig_height_5 <- max(8, nrow(heatmap_matrix_fig5) * 0.25 + 2)

  # Save PNG
  png("figure_panels/Figure5A_StrongBinders_ByResidue.png",
      width = 8, height = fig_height_5, units = "in", res = 300)

  pheatmap(
    as.matrix(heatmap_matrix_fig5),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colors_fig5,
    breaks = seq(0, 100, length.out = 101),
    border_color = "white",
    cellwidth = 40,
    cellheight = 14,
    fontsize_row = 8,
    fontsize_col = 10,
    main = "% Strong Binders by PTM+Residue and Allele",
    display_numbers = TRUE,
    number_format = "%.0f",
    number_color = "black",
    fontsize_number = 7,
    na_col = "gray90",
    legend = TRUE,
    legend_breaks = c(0, 25, 50, 75, 100),
    legend_labels = c("0%", "25%", "50%", "75%", "100%"),
    annotation_row = row_annotation_fig5,
    annotation_colors = ann_colors_fig5,
    angle_col = 0,
    gaps_row = gaps_fig5
  )

  dev.off()
  cat("Saved: figure_panels/Figure5A_StrongBinders_ByResidue.png\n")

  # Save PDF
  pdf("figure_panels/Figure5A_StrongBinders_ByResidue.pdf",
      width = 8, height = fig_height_5)

  pheatmap(
    as.matrix(heatmap_matrix_fig5),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colors_fig5,
    breaks = seq(0, 100, length.out = 101),
    border_color = "white",
    cellwidth = 40,
    cellheight = 14,
    fontsize_row = 8,
    fontsize_col = 10,
    main = "% Strong Binders by PTM+Residue and Allele",
    display_numbers = TRUE,
    number_format = "%.0f",
    number_color = "black",
    fontsize_number = 7,
    na_col = "gray90",
    annotation_row = row_annotation_fig5,
    annotation_colors = ann_colors_fig5,
    angle_col = 0,
    gaps_row = gaps_fig5
  )

  dev.off()
  cat("Saved: figure_panels/Figure5A_StrongBinders_ByResidue.pdf\n")

  # Save data
  write.csv(stats_fig5, "figure_panels/data_figure5A_strongbinders.csv", row.names = FALSE)
  cat("Saved: figure_panels/data_figure5A_strongbinders.csv\n")

} else {
  cat("WARNING: figure_panels/data_figure4A_circos.csv not found for Figure 5.\n")
}

# =============================================================================
# SESSION INFO
# =============================================================================
cat("\n=== All figures generated successfully! ===\n")
cat("Files created:\n")
cat("  - Figure1A_BarChart.png/.pdf\n")
cat("  - Figure1B_Donut.png/.pdf\n")
cat("  - Figure2_Length_Distribution.png/.pdf\n")
cat("  - Figure3_Phospho_ELRank.png/.pdf\n")
cat("  - Figure3B_PTM_Heatmap.png/.pdf\n")
cat("  - Figure4A_Circos_PTM_Landscape.png/.pdf (run Figure4A_Circos.R)\n")
cat("  - Figure4D_S_phospho_position.png\n")
cat("  - Figure5A_StrongBinders_ByResidue.png/.pdf\n")
cat("\nNote: Figure 4A circos plot is in separate file: Figure4A_Circos.R\n")
cat("Note: Figure 4B/4C position analysis in: Figure4B_Position_Heatmap.R, Figure4C_Position_Density.R\n")

