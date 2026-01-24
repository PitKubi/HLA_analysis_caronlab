# HLA-I PTM Analysis Pipeline

Analysis pipeline for post-translational modifications (PTMs) on HLA class I presented peptides from the JY cell line.

## Overview

This repository contains R scripts for analyzing immunopeptidomics data focusing on PTM characterization, including:

- PTM composition and distribution analysis
- Peptide length distribution by modification type
- Positional enrichment of modifications
- HLA binding affinity analysis (NetMHCpan EL Rank)
- Specific analyses for phosphorylation, ubiquitination, oxidation, and glycosylation

## Quick Start

```bash
# Clone repository
git clone https://github.com/PitKubi/HLA_analysis_caronlab.git
cd HLA_analysis_caronlab

# Run complete pipeline
Rscript run_all_figures.R

# Generate HTML report
Rscript render_report.R
```

## Requirements

### Software
- R ≥ 4.0

### R Packages
```r
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2",
                   "circlize", "pheatmap", "cowplot",
                   "ggridges", "scales", "base64enc"))
```

### Data
- `HLA-I_JY_DDA_PTMs_ResultsSummary_01062026.xlsx` (not tracked in git)
- Contact lab for data access

## Repository Structure

```
├── run_all_figures.R          # Master script - runs entire pipeline
├── render_report.R            # Generate HTML report with figures
├── data_loader.R              # Data extraction from Excel
├── Figure1_Composition.R      # PTM composition (bar + donut)
├── Figure2_LengthDistribution.R  # Length distribution ridgeline
├── Figure4A_Circos.R          # PTM landscape circos plot
├── Figure4B_Position_Heatmap.R   # Position enrichment by residue
├── Figure4B2_Position_Heatmap_ByPTM.R  # Position enrichment by PTM
├── Figure4C_PositionDensity.R # PTM position density curves
├── Figure5A_StrongBinders.R   # Binding affinity heatmaps
├── Figure5B_ELRank_Distribution.R  # EL Rank distributions
├── Figure6_SpecificPTMs.R     # Specific PTM analyses (6A-6C)
├── Figure6DE_Improved.R       # Enhanced 6D and 6E analyses
├── HLA-I_PTM_Analysis_Report.md   # Full documentation (Markdown)
├── HLA-I_PTM_Analysis_Report.html # Self-contained HTML report
└── figure_panels/             # Output directory
    ├── *.png                  # Figure files
    ├── *.pdf                  # Vector graphics
    └── *.csv                  # Data tables
```

## Figures Generated

| Figure | Description |
|--------|-------------|
| 1A | PTM distribution bar chart |
| 1B | PTM composition donut chart |
| 2 | Peptide length distribution ridgeline |
| 4A | Circos plot: PTM-residue-allele relationships |
| 4B | Position enrichment heatmap (by residue) |
| 4B2 | Position enrichment heatmap (by PTM) |
| 4C | PTM position density curves |
| 5A | Binding affinity heatmaps (strong/weak) |
| 5B | EL Rank distribution histograms |
| 6A | Doubly phosphorylated site distribution |
| 6B | Ubiquitination G vs GG comparison |
| 6C | Biological vs artifact oxidation |
| 6D | Biological oxidation binding by residue |
| 6D2 | Biological oxidation by residue × allele |
| 6E | Glycan complexity analysis |
| 6E2 | Individual glycan types |

## Analysis Parameters

- **Peptide length filter:** 8-14 amino acids
- **Binding thresholds (EL Rank):**
  - Strong binder: < 0.5%
  - Weak binder: 0.5% - 2%
  - Non-binder: ≥ 2%
- **Minimum counts:**
  - Position enrichment by residue: ≥ 20
  - Position enrichment by PTM: ≥ 50
  - Binding heatmaps: ≥ 10

## PTM Types Analyzed

| PTM | Target Residues |
|-----|-----------------|
| Phosphorylation | S, T, Y |
| Cysteinylation | C |
| Deamidation | N, Q |
| Acetylation | K, N-term |
| Methylation | K, R |
| Dimethylation | K, R |
| Citrullination | R |
| Biological Oxidation | P, I, L, Q, S, T, V, C, D, E, N, Y, G, K, R |
| Artifact Oxidation | M, W, H |
| Ubiquitination (G/GG) | K |
| N-Glycosylation | N |
| SUMOylation | K |

## Data Updates

The pipeline is designed for reproducibility with updated data:

1. Replace the Excel file with new data (same format)
2. Run `Rscript run_all_figures.R`
3. All figures regenerate automatically

## Documentation

- **HLA-I_PTM_Analysis_Report.md** - Full scientific documentation
- **HLA-I_PTM_Analysis_Report.html** - Self-contained HTML with embedded figures

## Citation

If using this pipeline, please cite:
- Caron Lab immunopeptidomics methodology
- NetMHCpan for binding predictions

## License

[Specify license]

## Contact

Caron Laboratory
CHU Sainte-Justine Research Center
University of Montreal
