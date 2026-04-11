# BGI DNBSeq Amplicon Analysis Pipeline

An R-based pipeline for reproducing BGI's 16S/ITS amplicon analysis workflow from raw OTU tables to publication-ready figures. Implements all 18 statistical modules described in BGI's standard amplicon report, with automated group-wise comparison execution across 11 sample groupings.

## Overview

This pipeline takes BGI's intermediary deliverables (OTU tables, taxonomy assignments, phylogenetic trees, PICRUSt2 predictions) and performs a comprehensive downstream analysis including:

- **Alpha Diversity** — Sobs, Chao1, ACE, Shannon, Simpson, Good's coverage with rarefaction
- **Beta Diversity** — Bray-Curtis, weighted/unweighted UniFrac, Pearson dissimilarity with bootstrapped PCoA
- **Taxa Composition** — Phylum-through-species level barplots and heatmaps
- **Differential Analysis** — Wilcoxon/Kruskal-Wallis with FDR correction at OTU and multi-taxonomy levels
- **Functional Prediction** — PICRUSt2-based KEGG pathway barplots
- **Advanced Analysis** — DCA model selection → RDA/CCA ordination, Random Forest with 10×10 CV + ROC/AUC
- **Rarefaction Curves** — Species richness and Shannon rarefaction
- **PCA/PCoA** — At OTU and L2–L7 taxonomy levels
- **Similarity Tests** — ANOSIM and MRPP on Bray-Curtis distances
- **Venn/Flower Diagrams** — Shared and unique OTU visualization
- **Rank Abundance** — Rank abundance and cumulative proportion curves
- **PLS-DA** — Partial least squares discriminant analysis
- **Network Analysis** — FDR-corrected Spearman correlation co-occurrence networks
- **Enterotypes** — PAM clustering on Jensen-Shannon divergence with silhouette width optimization
- **Multi-level Taxa** — Stacked barplots at all taxonomic ranks (L2–L7)
- **Functional Expansion** — Differential testing across KO, COG, EC, and MetaCyc pathways
- **UniFrac Beta** — Phylogenetic beta diversity with tree auto-detection per comparison
- **NMDS** — Non-metric multidimensional scaling ordination

## Repository Structure

```
BGI_Amplicon_Workflow/
├── README.md
├── metadata.tsv                        # Sample metadata (51 samples, 17 groups)
├── .gitignore
│
├── analysis/
│   ├── install_packages.R              # One-time R package installer
│   ├── 00_run_all_groups.R             # Master wrapper — runs all scripts × all comparisons
│   ├── 01_alpha_diversity.R            # Alpha diversity indices + boxplots
│   ├── 02_beta_diversity.R             # Beta diversity (Bray-Curtis, Pearson, PCoA)
│   ├── 03_taxa_composition.R           # Taxonomic barplots + heatmaps
│   ├── 04_differential_analysis.R      # OTU + multi-level differential tests
│   ├── 05_function_prediction.R        # PICRUSt2 KEGG visualization
│   ├── 06_advanced_analysis.R          # DCA/RDA/CCA + Random Forest + ROC
│   ├── 07_rarefaction_curves.R         # Species + Shannon rarefaction
│   ├── 08_pca_analysis.R              # PCA at OTU + taxonomy levels
│   ├── 09_similarity_tests.R           # ANOSIM + MRPP
│   ├── 10_venn_flower.R               # Venn diagrams + flower plots
│   ├── 11_rank_abundance.R            # Rank abundance + cumulative curves
│   ├── 12_plsda.R                     # PLS-DA ordination
│   ├── 13_network.R                   # Spearman network + FDR heatmap
│   ├── 14_enterotypes.R               # PAM enterotyping on JSD
│   ├── 15_multilevel_taxa.R           # Stacked barplots (L2–L7)
│   ├── 16_function_expansion.R         # COG/EC/MetaCyc differential
│   ├── 17_unifrac_beta.R             # UniFrac w/ phylogenetic tree
│   └── 18_nmds.R                     # NMDS ordination
│
├── BGI_Result/                         # BGI's original intermediary data (gitignored)
│   ├── OTU/                            #   OTU tables, taxonomy, L2-L7 tables
│   ├── Beta/                           #   Phylogenetic trees per comparison
│   ├── Genus_Tree/                     #   Genus-level trees
│   └── Picrust/                        #   PICRUSt2 predictions (KO/COG/EC/MetaCyc)
│
└── BGI_Reproduced/                     # Pipeline output directory (gitignored)
```

## Prerequisites

### R ≥ 4.5

The pipeline requires **R 4.5+** with the following packages:

**CRAN:**
`vegan`, `ggplot2`, `reshape2`, `ggpubr`, `pheatmap`, `randomForest`, `caret`, `pROC`, `ape`, `igraph`, `cluster`, `ade4`, `VennDiagram`, `UpSetR`, `scales`, `futile.logger`, `psych`, `MLmetrics`

**Bioconductor:**
`phyloseq`, `mixOmics`

### Installation

```r
# From the analysis/ directory:
source("install_packages.R")
```

This automatically installs all required packages and reports any failures.

## Usage

### Quick Start

```bash
cd analysis

# 1. Install dependencies (first time only)
Rscript install_packages.R

# 2. Run the full pipeline (18 scripts × 11 comparisons)
Rscript 00_run_all_groups.R
```

### What `00_run_all_groups.R` Does

1. Reads `metadata.tsv` and iterates through 11 predefined group comparisons (A-B, A-C, ..., all 17 groups)
2. For each comparison, creates a subset metadata file
3. Sources each of the 18 analysis scripts in a sandboxed R environment
4. Injects output directory overrides to route results to `BGI_Reproduced/`
5. Reports OK/ERROR status for each script-comparison pair

### Running Individual Scripts

Each script can also run standalone:

```r
setwd("analysis")

# Override defaults before sourcing:
meta_file <- "../metadata.tsv"
output_dir <- "../BGI_Reproduced/Alpha_Box"

source("01_alpha_diversity.R")
```

All configuration variables are guarded with `if (!exists(...) || is.null(...))`, allowing the wrapper or user to inject overrides without modifying script source.

### Input Data

The pipeline reads from `BGI_Result/` (gitignored, not included in this repository):

| File | Format | Description |
|------|--------|-------------|
| `OTU/OTU_table_for_biom.txt` | TSV, 2-line header | 5,358 OTUs × 51 samples + taxonomy |
| `OTU/OTU_taxonomy.xls` | TSV | OTU → taxonomic lineage mapping |
| `OTU/OTU_table_L{2-7}.txt` | TSV | Taxonomy-aggregated abundance tables |
| `Beta/{comparison}/*.tree.txt` | Newick | Phylogenetic trees per comparison |
| `Picrust/Function_Prdeict/` | TSV (.xls) | PICRUSt2 KO/COG/EC/MetaCyc predictions |

## Key Design Decisions

- **Path injection:** All scripts use `if (!exists("var") || is.null(var))` guards, enabling the master wrapper to redirect outputs to `BGI_Reproduced/` while keeping `BGI_Result/` read-only.
- **Input/output separation:** The wrapper explicitly excludes input-only directories (`otu_dir`, `base_dir`, `tree_dir_*`) from output redirection to avoid breaking file discovery.
- **Explicit namespacing:** `vegan::diversity()` and `igraph::degree()` are explicitly namespaced to prevent function masking when `caret`, `randomForest`, or `psych` load competing generics in the same session.
- **Defensive loading:** `make.unique()` is used for row names in PICRUSt tables where duplicate pathway IDs exist in BGI deliverables.
- **DCA fallback:** `decorana()` is wrapped in `tryCatch` — comparisons where DCA fails to converge (high sparsity) gracefully fall back to RDA.

## Group Comparisons

The pipeline runs the following 11 group comparisons as defined in the BGI report:

| Comparison | Groups | Samples |
|-----------|--------|---------|
| A-B | A, B | 6 |
| A-C | A, C | 6 |
| B-C | B, C | 6 |
| A-B-C | A, B, C | 9 |
| A-B-C-D-E-P | A–E, P | 18 |
| F-G-H-I-J-P | F–J, P | 18 |
| K-L-M-N-O-P-Q | K–Q | 21 |
| A-B-C-D-E-F-G-H-I-J-P | A–J, P | 33 |
| A-B-C-D-E-K-L-M-N-O-P-Q | A–E, K–Q | 36 |
| P-Q | P, Q | 6 |
| All | A–Q (17 groups) | 51 |

## License

This project is for academic and research purposes.
