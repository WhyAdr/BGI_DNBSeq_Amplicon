# BGI Amplicon Workflow — HPC Setup Guide (Ubuntu)

Complete environment setup for all 19 R scripts + local LEfSe, GraPhlAn, and FastTree2.

---

## 1. Conda Base (Recommended)

If you don't have `sudo`, Conda is the cleanest way to manage everything:

```bash
# Install Miniconda (if not already present)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

---

## 2. R Environment

```bash
conda create -n amplicon-r -c conda-forge r-base=4.2 r-essentials
conda activate amplicon-r
```

### System libraries (may need admin)

Some R packages require C/C++ libraries. If installs fail, ask your admin:

```bash
sudo apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev \
    libglpk-dev libgmp3-dev libfontconfig1-dev libfreetype6-dev \
    libpng-dev libtiff5-dev libjpeg-dev
```

### CRAN packages (14 total)

```r
install.packages(c(
    "ggplot2",       # Core visualization (all scripts)
    "reshape2",      # Data reshaping (01, 03, 05, 15, 16)
    "vegan",         # Diversity indices, ANOSIM, MRPP, NMDS (01, 02, 07, 09, 11, 17, 18)
    "ggpubr",        # Stat comparison annotations (01)
    "pheatmap",      # Heatmaps (03, 13, 16)
    "randomForest",  # Random Forest classifier (06)
    "VennDiagram",   # Venn diagrams (10)
    "UpSetR",        # UpSet intersection plots (10)
    "cluster",       # PAM clustering for enterotypes (14)
    "ade4",          # BCA/PCoA ordination (14)
    "igraph",        # Network analysis (13)
    "ape",           # Phylogenetic tree I/O (17)
    "scales"         # Rescale for network node sizing (13)
), repos = "https://cloud.r-project.org")
```

### Bioconductor packages (2 total)

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "phyloseq",      # UniFrac distance calculations (17)
    "mixOmics"       # PLS-DA ordination (12)
))
```

> [!TIP]
> **Verification**: After installing, run this in R to confirm all packages load:
> ```r
> pkgs <- c("ggplot2","reshape2","vegan","ggpubr","pheatmap","randomForest",
>           "VennDiagram","UpSetR","cluster","ade4","igraph","ape","scales",
>           "phyloseq","mixOmics")
> sapply(pkgs, require, character.only = TRUE)
> ```
> All should return `TRUE`.

---

## 3. FastTree2 (Local Installation)

FastTree builds phylogenetic trees from aligned sequences. You only need this if you:
- Re-run OTU clustering from scratch (new `OTU_final.fasta`)
- Want to regenerate per-comparison genus-level trees (like BGI's `Genus_Tree/` directory)

> [!NOTE]
> For the current pipeline using BGI's pre-computed data, FastTree is **optional** — the existing tree files are sufficient to compute UniFrac distances.

### Install via Conda (easiest)

```bash
conda activate amplicon-r   # or create a separate env
conda install -c bioconda fasttree
```

### Install from source (alternative)

```bash
wget http://www.microbesonline.org/fasttree/FastTree.c
gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
mv FastTree ~/bin/  # or any directory in your $PATH
```

### Usage example

```bash
# Build a tree from aligned OTU representative sequences
# Input must be a multiple sequence alignment (FASTA format)
# Use MUSCLE or MAFFT to align first if needed:
mafft --auto OTU_final.fasta > OTU_final_aligned.fasta

# Then build the tree:
FastTree -gtr -nt OTU_final_aligned.fasta > OTU_phylogeny_tree.nwk
```

> [!IMPORTANT]
> FastTree requires a **multiple sequence alignment** as input, not raw sequences. BGI uses MUSCLE for alignment before FastTree. You'll need MUSCLE or MAFFT:
> ```bash
> conda install -c bioconda muscle   # or mafft
> ```

---

## 4. LEfSe (Local Installation)

LEfSe performs biomarker discovery via LDA effect sizes.

```bash
# Create a dedicated environment (Python 2/3 compatibility issues)
conda create -n lefse -c bioconda -c conda-forge lefse
conda activate lefse
```

### Usage with our pipeline

Script `04_differential_analysis.R` generates the input file. Then:

```bash
conda activate lefse

# 1. Format the input
lefse_format_input.py lefse_input.in lefse_formatted.in \
    -c 1 -s 2 -u 3 -o 1000000

# 2. Run the statistical analysis
run_lefse.py lefse_formatted.in lefse_result.res

# 3. Plot the bar chart (LDA scores)
plot_res.py lefse_result.res lefse_barplot.png \
    --format png --dpi 300

# 4. Plot the cladogram
plot_cladogram.py lefse_result.res lefse_cladogram.png \
    --format png --dpi 300
```

---

## 5. GraPhlAn (Local Installation)

GraPhlAn renders circular taxonomic cladograms.

```bash
# Install in the same env as LEfSe, or a new one
conda create -n graphlan -c bioconda -c conda-forge graphlan
conda activate graphlan
```

### Usage

```bash
conda activate graphlan

# 1. Create annotation file (from your taxonomy data)
# 2. Annotate the tree
graphlan_annotate.py input_tree.txt annotated_tree.xml \
    --annot annotation.txt

# 3. Render the cladogram
graphlan.py annotated_tree.xml cladogram.png \
    --dpi 300 --size 7
```

> [!NOTE]
> GraPhlAn input preparation is manual — you need to create annotation files from the taxonomy tables. This is not automated in our R pipeline yet.

---

## 6. Running the Full Pipeline

```bash
# 1. Activate the R environment
conda activate amplicon-r

# 2. Navigate to the analysis directory
cd /path/to/BGI_Amplicon_Workflow/analysis/

# 3. (Optional) Edit metadata.tsv to remove outliers
#    Just delete the rows for samples you want to exclude.

# 4. Run individual scripts
Rscript 01_alpha_diversity.R
Rscript 17_unifrac_beta.R
# ... etc.

# 5. Or run ALL analyses across ALL 11 group comparisons
Rscript 00_run_all_groups.R

# 6. Run LEfSe on the output from script 04
conda activate lefse
lefse_format_input.py ../BGI_Result/Differential/lefse_input.in formatted.in -c 1 -s 2 -u 3
run_lefse.py formatted.in result.res
plot_res.py result.res lefse_barplot.png
```

---

## Summary: What's handled where

| Analysis | Tool | Where it runs |
|----------|------|---------------|
| Alpha diversity (Shannon, Chao1, etc.) | R (`vegan`) | `01_alpha_diversity.R` |
| Beta diversity (Bray-Curtis PCoA) | R (`vegan`) | `02_beta_diversity.R` |
| UniFrac (weighted + unweighted) | R (`phyloseq`) | `17_unifrac_beta.R` |
| UPGMA clustering trees | R (`hclust`) | `17_unifrac_beta.R` |
| PERMANOVA | R (`vegan::adonis2`) | `17_unifrac_beta.R` |
| NMDS | R (`vegan::metaMDS`) | `18_nmds.R` |
| PLS-DA | R (`mixOmics`) | `12_plsda.R` |
| Rarefaction curves | R (`vegan`) | `07_rarefaction_curves.R` |
| ANOSIM + MRPP | R (`vegan`) | `09_similarity_tests.R` |
| Co-occurrence network | R (`igraph`) | `13_network.R` |
| Enterotypes (JSD + PAM) | R (`cluster`) | `14_enterotypes.R` |
| PICRUSt2 visualization | R | `05_function_prediction.R` + `16_function_expansion.R` |
| LEfSe biomarker discovery | Python (`lefse`) | Command line (local) |
| GraPhlAn cladograms | Python (`graphlan`) | Command line (local) |
| Phylogenetic tree building | FastTree2 | Command line (only if re-clustering) |
| OTU clustering | USEARCH/UPARSE | **Not needed** (using BGI's OTU table) |
| Taxonomy classification | RDP Classifier | **Not needed** (using BGI's taxonomy) |
| Read merging | FLASH | **Not needed** (using BGI's OTU table) |

> [!IMPORTANT]
> **USEARCH, RDP Classifier, FLASH, and QIIME** are upstream tools that BGI already ran. You do **not** need to install these unless you plan to re-process from raw FASTQ reads.
