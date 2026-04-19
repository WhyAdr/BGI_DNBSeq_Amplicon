# BGI DNBSeq Amplicon Analysis Pipeline

An R-based pipeline for reproducing BGI's 16S amplicon analysis workflow from OTU tables to publication-ready figures. Implements 18 statistical modules described in BGI's standard amplicon report, with 17 running automatically in a group-wise execution loop across 11 sample comparisons.

## Overview

This pipeline takes BGI's intermediary deliverables (OTU tables, taxonomy assignments, phylogenetic trees, PICRUSt2 predictions) and performs comprehensive downstream analysis including:

- **Alpha Diversity** — Sobs, Chao1, ACE, Shannon, Simpson, Good's coverage with rarefaction
- **Beta Diversity** — Bray-Curtis, weighted/unweighted UniFrac, Pearson dissimilarity with bootstrapped PCoA
- **Taxa Composition** — Phylum-through-species level barplots and heatmaps
- **Differential Analysis** — Wilcoxon/Kruskal-Wallis with FDR correction at OTU and multi-taxonomy levels; LEfSe input preparation
- **Functional Prediction** — PICRUSt2-based KEGG pathway barplots
- **Advanced Analysis** — DCA model selection → RDA/CCA ordination, Random Forest with 10×10 CV + ROC/AUC
- **Rarefaction Curves** — Species richness and Shannon rarefaction
- **PCA** — At OTU and L2–L7 taxonomy levels
- **Similarity Tests** — ANOSIM and MRPP on Bray-Curtis distances
- **Venn/Flower Diagrams** — Shared and unique OTU visualization
- **Rank Abundance** — Rank abundance and cumulative proportion curves
- **PLS-DA** — Partial least squares discriminant analysis
- **Network Analysis** — Species-aggregated Spearman correlation co-occurrence networks with Cytoscape-style visualization
- **Enterotypes** — PAM clustering on Jensen-Shannon divergence with Calinski-Harabasz optimization
- **Multi-level Taxa** — Stacked barplots at all taxonomic ranks (L2–L7)
- **Functional Expansion** — Differential testing across KO, COG, EC, and MetaCyc pathways
- **UniFrac Beta** — Phylogenetic beta diversity with per-comparison tree auto-detection
- **NMDS** — Non-metric multidimensional scaling ordination with ANOSIM overlay

## Repository Structure

```
BGI_Amplicon_Workflow/
├── README.md
├── metadata.tsv                        # Sample metadata (51 samples, 17 groups A–Q)
├── .gitignore
│
├── analysis/
│   ├── install_packages.R              # One-time R package installer
│   ├── 00_run_all_groups.R             # Master wrapper — runs scripts × comparisons
│   ├── 01_alpha_diversity.R            # Alpha diversity indices + boxplots
│   ├── 02_beta_diversity.R             # Beta diversity (Bray-Curtis, Pearson, PCoA)
│   ├── 03_taxa_composition.R           # Taxonomic barplots + heatmaps
│   ├── 04_differential_analysis.R      # OTU + multi-level differential tests + LEfSe prep
│   ├── 05_function_prediction.R        # PICRUSt2 KEGG visualization
│   ├── 06_advanced_analysis.R          # DCA/RDA/CCA + Random Forest (manual only)
│   ├── 07_rarefaction_curves.R         # Species + Shannon rarefaction
│   ├── 08_pca_analysis.R              # PCA at OTU + taxonomy levels (L2–L7)
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
├── BGI_Result/                         # BGI's original data + reference outputs (gitignored)
│   ├── OTU/                            #   ⬤ INPUT: OTU tables, taxonomy, L2–L7 tables
│   ├── Beta/                           #   ⬤ INPUT: Per-comparison phylogenetic trees
│   ├── Genus_Tree/                     #   ⬤ INPUT: Genus-level phylogenetic trees
│   ├── Picrust/Function_Prdeict/       #   ⬤ INPUT: PICRUSt2 predictions (KO/COG/EC/MetaCyc)
│   ├── Alpha_Box/                      #   ○ REFERENCE: BGI alpha diversity outputs
│   ├── Alpha_Rarefaction/              #   ○ REFERENCE: BGI rarefaction outputs
│   ├── Barplot/                        #   ○ REFERENCE: BGI taxonomic barplots
│   ├── Heatmap/                        #   ○ REFERENCE: BGI heatmap outputs
│   ├── Network/                        #   ○ REFERENCE: BGI network + correlation outputs
│   ├── Diff/                           #   ○ REFERENCE: BGI differential analysis outputs
│   ├── Lefse/                          #   ○ REFERENCE: BGI LEfSe results
│   ├── Graphlan/                       #   ○ REFERENCE: BGI GraPhlAn circular trees
│   ├── Enterotypes/                    #   ○ REFERENCE: BGI enterotype results
│   ├── NMDS/, PCA/, PLSDA/            #   ○ REFERENCE: BGI ordination outputs
│   ├── SimilarityAnalysis/             #   ○ REFERENCE: BGI ANOSIM/MRPP outputs
│   ├── Venn/, Flower/                  #   ○ REFERENCE: BGI Venn/Flower outputs
│   ├── Cumulative_Curve/, OTU_Rank/    #   ○ REFERENCE: BGI rank/curve outputs
│   ├── Cleandata/, Rawdata/, Tag/      #   ○ UPSTREAM: Sequencing reads (not used by pipeline)
│   └── readme.en.txt                   #   BGI's own file manifest
│
└── BGI_Reproduced/                     # Pipeline output directory (gitignored)
```

> **⬤ INPUT** = read by the R pipeline as data dependencies.
> **○ REFERENCE** = BGI's own analysis results, used only for visual parity comparison.

## Prerequisites

### R ≥ 4.1

The pipeline requires **R 4.1+** with the following packages:

**CRAN:**
`vegan`, `ggplot2`, `reshape2`, `ggpubr`, `pheatmap`, `randomForest`, `caret`, `pROC`, `ape`, `igraph`, `cluster`, `clusterSim`, `ade4`, `VennDiagram`, `UpSetR`, `scales`, `futile.logger`, `psych`, `MLmetrics`

**Bioconductor:**
`phyloseq`, `mixOmics`

### Installation

```r
# From the analysis/ directory:
source("install_packages.R")
```

This automatically installs all required CRAN and Bioconductor packages and reports any failures.

## Usage

### Quick Start

```bash
cd analysis

# 1. Install dependencies (first time only)
Rscript install_packages.R

# 2. Run the full pipeline (17 scripts × 11 comparisons)
Rscript 00_run_all_groups.R
```

> **Note:** `06_advanced_analysis.R` (DCA/RDA/CCA + Random Forest) is excluded from the automated loop by design. Run it manually for specific comparisons as needed.

### What `00_run_all_groups.R` Does

1. Reads `metadata.tsv` and iterates through 11 predefined group comparisons
2. For each comparison, writes a temporary subset metadata file (`metadata_{comparison}.tsv`)
3. Sources 17 analysis scripts in a sandboxed `new.env()` R environment
4. Parses each script for `*_dir` variable assignments and redirects output paths from `BGI_Result/` to `BGI_Reproduced/`, while preserving input-only directories (`otu_dir`, `tree_dir_beta`, `tree_dir_genus`, `base_dir`)
5. Applies flat-directory routing for modules that BGI outputs without group subdirectories
6. Reports OK/ERROR status for each script × comparison pair

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

The pipeline reads from 4 directories inside `BGI_Result/` (gitignored, not included in this repository):

| File | Format | Description | Consumers |
|------|--------|-------------|-----------|
| `OTU/OTU_table_for_biom.txt` | TSV, 2-line header | 5,358 OTUs × 51 samples (raw counts) | 13 scripts |
| `OTU/OTU_taxonomy.xls` | TSV | OTU → full taxonomic lineage mapping | 03, 04, 13 |
| `OTU/OTU_table_L{2-7}.txt` | TSV | Pre-aggregated abundance at Phylum–Species | 04, 08, 14, 15 |
| `Beta/{comparison}/*.tree.txt` | Newick | Per-comparison pruned phylogenetic trees | 17 |
| `Genus_Tree/{comparison}.genus.phylogeny.tree` | Newick | Per-comparison genus-level trees (fallback) | 17 |
| `Picrust/Function_Prdeict/KO/ko_Level2_Function.xls` | TSV | PICRUSt2 KEGG Level 2 pathway abundances | 05 |
| `Picrust/Function_Prdeict/{KO,COG,EC,METACYC}/*` | TSV (.xls) | Multi-level functional predictions | 16 |

## Key Design Decisions

- **Path injection:** All scripts use `if (!exists("var") || is.null(var))` guards, enabling the master wrapper to redirect outputs to `BGI_Reproduced/` while keeping `BGI_Result/` read-only.
- **Input/output separation:** The wrapper explicitly excludes input-only directories (`otu_dir`, `base_dir`, `tree_dir_beta`, `tree_dir_genus`) from output redirection to avoid breaking file discovery.
- **Flat-directory whitelist:** Modules whose BGI reference outputs are flat (no group subdirectories) — `Alpha_Rarefaction`, `Alpha_Box`, `Cumulative_Curve`, `OTU_Rank`, `PCA`, `PLSDA`, `NMDS`, `Picrust`, `Function_Diff` — are routed flat via the `flat_dirs` whitelist in `00_run_all_groups.R`.
- **Explicit namespacing:** Function calls like `vegan::diversity()` and `vegan::estimateR()` are explicitly namespaced to prevent function masking when `igraph`, `caret`, `randomForest`, or `psych` load competing generics in the same session.
- **Defensive loading:** `make.unique()` is used for row names in PICRUSt tables where duplicate pathway IDs exist in BGI deliverables.
- **DCA fallback:** `decorana()` is wrapped in `tryCatch` — comparisons where DCA fails to converge (high sparsity) gracefully fall back to RDA.

## Known Issues

- **Lefse write leak:** `04_differential_analysis.R` (line 240) writes LEfSe input files directly to `../BGI_Result/Lefse/` instead of routing through `output_dir`. This violates the read-only contract on `BGI_Result/`.
- **Hardcoded L6 path:** `14_enterotypes.R` (line 36) hardcodes `../BGI_Result/OTU/OTU_table_L6.txt` instead of using the `otu_dir` variable, bypassing the orchestrator's input-directory protection.
- **CRAN mirror error:** `14_enterotypes.R` may report a CRAN mirror warning when `clusterSim` is not pre-installed.
- **"Too few points" warnings:** Expected behavior in beta diversity and NMDS modules for group subsets with very small sample sizes (e.g., 2-group comparisons with 3 samples each).

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
| All 17 groups | A–Q | 51 |

## Changelog

### 2026-04-13: Network Module Parity
- Implemented strict group-based sample subsetting (fixed "all samples" bug)
- Added mandatory species-level aggregation to prevent duplicate-label crashes
- Generated Species Abundance matrices with "Other" tier
- Applied Cytoscape-style visualization (phylum-based node colors, curved edges)
- Added dual PDF+PNG output for all plots

### 2026-04-13: Alpha Rarefaction Stability
- Flattened directory structure to match BGI output layout
- Added `vegan::` namespace protection for `diversity()` and `estimateR()` to prevent `igraph` masking crashes

### 2026-04-12: Beta Diversity & Ordination Fixes
- Fixed out-of-bounds PCoA matrix dimension mapping for low-N subgroups
- Integrated ANOSIM test into `09_similarity_tests.R`
- Matched BGI aesthetic frames in `18_nmds.R`

### 2026-04-11: Pipeline Orchestration
- Removed `06_advanced_analysis.R` from automated loop (manual-only)
- Fixed PICRUSt path sandboxing and `diff_dir` routing
- Expanded `flat_dirs` whitelist for correct directory structure parity

## License

This project is for academic and research purposes.
