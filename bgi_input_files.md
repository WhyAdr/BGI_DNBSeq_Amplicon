# BGI_Result Input File Mapping

This document maps every file in the `BGI_Result/` directory that serves as a **read-only input** to the R pipeline. Files that exist in `BGI_Result/` purely as BGI's own analysis outputs (which we reproduce into `BGI_Reproduced/`) are excluded — only upstream data dependencies are catalogued here.

> [!IMPORTANT]
> The `BGI_Result/` directory contains **24 subdirectories**, but only **4 of them** supply input data to the pipeline. The remaining 20 are BGI's own analysis outputs that serve as visual references for parity comparison. This distinction is critical for understanding the pipeline's actual dependency surface.

---

## 1. Core Abundance & Taxonomy — `OTU/`

The `OTU/` directory is the single most important input source. Nearly every analysis module depends on it.

### `OTU/OTU_table_for_biom.txt`

| Property | Value |
|:--|:--|
| **Format** | TSV with 2-line header (line 1: `# Constructed from biom file`, line 2: column headers). Trailing `taxonomy` column stripped at load time. |
| **Dimensions** | 5,358 OTUs × 51 samples |
| **Loading pattern** | `read.table(..., skip = 1, row.names = 1, comment.char = "")` followed by `if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL` |

**Consuming scripts** (13 of 17 automated):

| Script | Usage |
|:--|:--|
| [01_alpha_diversity.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/01_alpha_diversity.R) | Sobs, Chao1, ACE, Shannon, Simpson, Good's coverage |
| [02_beta_diversity.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/02_beta_diversity.R) | Bray-Curtis dissimilarity, Pearson PCoA |
| [04_differential_analysis.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/04_differential_analysis.R) | OTU-level Wilcoxon/Kruskal-Wallis, LEfSe prep |
| [06_advanced_analysis.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/06_advanced_analysis.R) | DCA → RDA/CCA, Random Forest (manual only) |
| [07_rarefaction_curves.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/07_rarefaction_curves.R) | Rarefaction subsampling for richness/Shannon curves |
| [08_pca_analysis.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/08_pca_analysis.R) | OTU-level PCA |
| [09_similarity_tests.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/09_similarity_tests.R) | ANOSIM + MRPP on Bray-Curtis |
| [10_venn_flower.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/10_venn_flower.R) | Shared/unique OTU sets |
| [11_rank_abundance.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/11_rank_abundance.R) | Rank abundance + cumulative curves |
| [12_plsda.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/12_plsda.R) | PLS-DA ordination |
| [13_network.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/13_network.R) | Species-aggregated Spearman correlation network |
| [14_enterotypes.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/14_enterotypes.R) | Fallback if L6 file is missing |
| [17_unifrac_beta.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/17_unifrac_beta.R) | UniFrac distance matrix construction |
| [18_nmds.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/18_nmds.R) | NMDS ordination |

---

### `OTU/OTU_taxonomy.xls`

| Property | Value |
|:--|:--|
| **Format** | TSV. Row names = OTU IDs, column `Taxonomy` = semicolon-delimited lineage (Domain through Species) |
| **Size** | ~425 KB, covers all 5,358 OTUs |

**Consuming scripts** (3 scripts):

| Script | Usage |
|:--|:--|
| [03_taxa_composition.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/03_taxa_composition.R) | Aggregate OTUs into taxonomic tiers for barplots/heatmaps |
| [04_differential_analysis.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/04_differential_analysis.R) | Multi-level differential testing + LEfSe input prep |
| [13_network.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/13_network.R) | Resolve OTU IDs → species names for heatmap labeling; extract Phylum for node coloring |

---

### `OTU/OTU_table_L{2-7}.txt`

Pre-aggregated abundance tables at six taxonomic depths:

| File | Level | Approx. Size |
|:--|:--|:--|
| `OTU_table_L2.txt` | Phylum | 11 KB |
| `OTU_table_L3.txt` | Class | 28 KB |
| `OTU_table_L4.txt` | Order | 58 KB |
| `OTU_table_L5.txt` | Family | 121 KB |
| `OTU_table_L6.txt` | Genus | 292 KB |
| `OTU_table_L7.txt` | Species | 442 KB |

**Consuming scripts** (4 scripts):

| Script | Files Used | Access Pattern |
|:--|:--|:--|
| [04_differential_analysis.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/04_differential_analysis.R) | L2–L7 | Via `otu_dir` variable, loops through levels |
| [08_pca_analysis.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/08_pca_analysis.R) | L2–L7 | Via `otu_dir` variable, loops through levels |
| [14_enterotypes.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/14_enterotypes.R) | L6 only | **Hardcoded** path `../BGI_Result/OTU/OTU_table_L6.txt` (not via `otu_dir`) |
| [15_multilevel_taxa.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/15_multilevel_taxa.R) | L2–L7 | Via `otu_dir` variable, loops through levels |

> [!WARNING]
> `14_enterotypes.R` hardcodes the L6 path instead of using the `otu_dir` variable. This means the orchestrator's input-dir protection doesn't apply to it — a potential issue if the `OTU/` directory is ever relocated.

---

## 2. Functional Predictions — `Picrust/Function_Prdeict/`

> [!NOTE]
> The directory name `Function_Prdeict` is a BGI typo (should be `Function_Predict`). All pipeline code preserves this misspelling to maintain path compatibility.

### `Function_Prdeict/KO/ko_Level2_Function.xls`

| Property | Value |
|:--|:--|
| **Format** | TSV. Rows = KEGG pathway IDs at Level 2, Columns = samples |
| **Consuming script** | [05_function_prediction.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/05_function_prediction.R) — generates KEGG pathway overview barplots |

### `Function_Prdeict/{DB}/*_Function.xls` (Full inventory)

`16_function_expansion.R` iterates over four databases, each with multiple hierarchy levels:

| Database | Files | Total Size |
|:--|:--|:--|
| **KO** | `ko_Level1_Function.xls`, `ko_Level2_Function.xls`, `ko_Level3_Function.xls` | ~196 KB |
| **COG** | `cog_Level1_Function.xls`, `cog_Level2_Function.xls`, `cog_Level3_Function.xls` | ~2.7 MB |
| **EC** | `EC.Function.descrip.xls` (auto-detected via `*.xls` glob) | ~1.3 MB |
| **METACYC** | `metacyc_Level1_Function.xls`, `metacyc_Level2_Function.xls`, `metacyc_Level3_Function.xls` | ~457 KB |

**Access pattern**: KO files are explicitly named in the script; COG, EC, and METACYC use `list.files(..., pattern = "\\.xls$")` auto-detection.

---

## 3. Phylogenetic Trees — `Beta/` and `Genus_Tree/`

### `Beta/{comparison}/*.tree.txt`

| Property | Value |
|:--|:--|
| **Format** | Newick phylogenetic trees, one per comparison group |
| **Example** | `Beta/A-B/A-B.OTU_final_phylogeny_tree.txt` (~130 KB) |
| **Consuming script** | [17_unifrac_beta.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/17_unifrac_beta.R) — weighted + unweighted UniFrac, PCoA, UPGMA trees |
| **Discovery** | Dynamic: the script searches `tree_dir_beta` using `comp_suffix` to find the matching tree for the current group permutation |

> [!NOTE]
> BGI provides **per-comparison pruned trees** rather than a single global tree. This means each tree contains only the OTUs observed in that subset of samples, which improves UniFrac precision but requires the dynamic discovery mechanism.

### `Beta/{comparison}.Mapping.txt` and `.Mapping.Box.txt`

These are BGI's sample-to-group mapping files for each comparison. They are **not read** by the pipeline (which generates its own mappings from `metadata.tsv`), but exist as part of BGI's original output structure.

### `Genus_Tree/{comparison}.genus.phylogeny.tree`

| Property | Value |
|:--|:--|
| **Format** | Newick trees at genus level, one per comparison |
| **Consuming script** | [17_unifrac_beta.R](file:///d:/W/BGI%20Amplicon%20Workflow/analysis/17_unifrac_beta.R) — accessed via `tree_dir_genus` variable as a **fallback** when the OTU-level tree in `Beta/` fails to load or match |

---

## 4. Sample Metadata — `metadata.tsv` (root directory)

While not inside `BGI_Result/`, this is a critical pipeline-wide input:

| Property | Value |
|:--|:--|
| **Format** | TSV, 2 columns: `SampleID`, `Group` |
| **Dimensions** | 51 samples across 17 groups (A–Q) |
| **Consuming scripts** | All 17 automated scripts + the orchestrator `00_run_all_groups.R` |
| **Orchestrator behavior** | The wrapper subsets this file per comparison group and writes temporary `metadata_{comparison}.tsv` files that individual scripts read via the injected `meta_file` variable |

---

## Summary: Input vs. Reference Directories

| Directory | Role | Read by Pipeline? |
|:--|:--|:--|
| `OTU/` | **Input** — abundance tables, taxonomy | ✅ Yes (13+ scripts) |
| `Picrust/Function_Prdeict/` | **Input** — functional predictions | ✅ Yes (2 scripts) |
| `Beta/` (trees only) | **Input** — phylogenetic trees | ✅ Yes (1 script) |
| `Genus_Tree/` | **Input** — genus-level trees (fallback) | ✅ Yes (1 script) |
| `Alpha_Box/`, `Alpha_Rarefaction/`, `Barplot/`, `Heatmap/`, `Cumulative_Curve/`, `OTU_Rank/`, `PCA/`, `PLSDA/`, `NMDS/`, `SimilarityAnalysis/`, `Venn/`, `Flower/`, `Network/`, `Enterotypes/`, `Diff/`, `Lefse/`, `Graphlan/` | **Reference** — BGI outputs for parity comparison | ❌ Not read (reproduced into `BGI_Reproduced/`) |
| `Picrust/Function_Diff/` | **Reference** — BGI's differential functional results | ❌ Not read (reproduced) |
| `Cleandata/`, `Rawdata/`, `Tag/` | **Upstream QC** — raw/filtered sequencing reads | ❌ Not read by R pipeline |

> [!CAUTION]
> **Known leak**: `04_differential_analysis.R` line 240 writes LEfSe input files directly to `../BGI_Result/Lefse/` instead of routing through `output_dir`. This violates the read-only contract on `BGI_Result/` and should be addressed in a future refactor.
