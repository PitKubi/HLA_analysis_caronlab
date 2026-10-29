library(readxl)
f <- "HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx"

# ============ Check actual output CSVs ============
cat("========================================================\n")
cat("Actual output CSV verification\n")
cat("========================================================\n")

if (file.exists("figure_panels/data_ptm_sites.csv")) {
  ptm <- read.csv("figure_panels/data_ptm_sites.csv")
  cat(sprintf("data_ptm_sites.csv: %d rows\n", nrow(ptm)))
  cat("\nPTM breakdown:\n")
  print(table(ptm$PTM))
  cat("\nResidue distribution by PTM:\n")
  for (ptm_name in unique(ptm$PTM)) {
    cat(sprintf("\n  %s:\n", ptm_name))
    print(table(ptm$Residue[ptm$PTM == ptm_name]))
  }
  cat("\nLength distribution:\n")
  print(table(ptm$Length))
} else {
  cat("data_ptm_sites.csv not found\n")
}

if (file.exists("figure_panels/data_background.csv")) {
  bg <- read.csv("figure_panels/data_background.csv")
  cat(sprintf("\ndata_background.csv: %d rows\n", nrow(bg)))
  cat(sprintf("Length range: %d - %d\n", min(bg$Length), max(bg$Length)))
} else {
  cat("data_background.csv not found\n")
}

if (file.exists("figure_panels/data_binding.csv")) {
  bd <- read.csv("figure_panels/data_binding.csv")
  cat(sprintf("\ndata_binding.csv: %d rows\n", nrow(bd)))
  cat("\nPTM breakdown:\n")
  print(table(bd$PTM))
  cat("\nRows with NA EL_Rank:\n")
  print(table(bd$PTM[is.na(bd$EL_Rank)]))
} else {
  cat("data_binding.csv not found\n")
}
