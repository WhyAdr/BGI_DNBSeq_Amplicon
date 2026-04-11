# ==============================================================================
# 01_alpha_diversity.R
# BGI Amplicon Workflow - Alpha Diversity Analysis (Optimized)
# ==============================================================================
# Calculates and visualizes intra-sample diversity (Sections 8 and M7)
# Evaluates species richness, evenness, and sequencing depth coverage.
# ==============================================================================

# Load necessary libraries
# Note: Install missing packages using install.packages(c("vegan", "ggplot2", "reshape2", "ggpubr"))
library(vegan)
library(ggplot2)
library(reshape2)
library(ggpubr)

# --- Configuration ---
if (!exists("otu_file") || is.null(otu_file)) otu_file <- "../BGI_Result/OTU/OTU_table_for_biom.txt"
if (!exists("meta_file") || is.null(meta_file)) meta_file <- "../metadata.tsv"
if (!exists("output_dir") || is.null(output_dir)) output_dir <- "../BGI_Result/Alpha_Box"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading & Preprocessing ---
# Read OTU table (samples in columns, OTUs in rows)
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t",
                  comment.char = "", skip = 1)  # Skip "# Construct" junk line; preserve #OTU ID header
# Remove the 'taxonomy' column (present in OTU_table_for_biom.txt)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
# Read metadata
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

# Align samples between OTU table and metadata
common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- Alpha Diversity Calculation (Section 8) ---
# NOTE: BGI reports Simpson as dominance D = Σp_i² (not diversity 1−D).
#       vegan::diversity(index="simpson") returns 1−Σp², so we invert.
alpha_metrics <- data.frame(
    Sample_Name = common_samples,
    Sobs = specnumber(otu, MARGIN = 2), # Observed species richness
    Chao1 = estimateR(t(otu))["S.chao1", ], # Estimated richness (Chao1)
    ACE = estimateR(t(otu))["S.ACE", ], # Estimated richness (ACE)
    Shannon = diversity(t(otu), index = "shannon"), # Shannon diversity (accounts for richness & evenness)
    Simpson = 1 - diversity(t(otu), index = "simpson"), # Dominance = Σp² (BGI Table 7 convention)
    Coverage = NA  # Placeholder — computed on rarefied data below
)

# --- Good's Coverage on Rarefied Data (per BGI Section 9) ---
# BGI: "Normalisation is first applied by randomly sub-sampling all samples
#       to the minimum sequencing depth (rarefaction)."
min_depth <- min(colSums(otu))
set.seed(42)
otu_rare <- as.data.frame(t(rrarefy(t(otu), sample = min_depth)))
# Good's coverage: C = 1 − (n_singletons / total_reads_after_rarefaction)
alpha_metrics$Coverage <- 1 - (colSums(otu_rare == 1) / colSums(otu_rare))
cat(sprintf("Good's Coverage computed on rarefied data (depth = %d).\n", min_depth))

# Merge calculations with metadata
alpha_data <- merge(alpha_metrics, metadata, by.x = "Sample_Name", by.y = 1)
write.table(alpha_data, file = file.path(output_dir, "alpha_diversity_summary.xls"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Group-level Boxplots (Figure 20 / Section 8) ---
# Assuming metadata has a column named "Group"
group_col <- "Group"

if (!group_col %in% colnames(alpha_data)) {
    stop("Error: 'Group' column not found in metadata.")
}

# Dynamic statistical test based on number of groups in the metadata
num_groups <- length(unique(alpha_data[[group_col]]))
if (num_groups < 2) {
    print("Only 1 group found. Skipping statistical testing.")
    test_method <- NULL
} else if (num_groups == 2) {
    # Wilcoxon Rank-Sum for exactly 2 groups
    test_method <- "wilcox.test"
} else {
    # Kruskal-Wallis Test for >= 3 groups
    test_method <- "kruskal.test"
}

metrics_to_plot <- c("Sobs", "Chao1", "ACE", "Shannon", "Simpson")

for (metric in metrics_to_plot) {
    p <- ggplot(alpha_data, aes_string(x = group_col, y = metric, fill = group_col)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.8) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        theme_bw() +
        labs(title = paste(metric, "Index by Group"), x = "Group", y = metric)
    
    # Overlay statistical significance if groups are valid
    if (!is.null(test_method)) {
        p <- p + stat_compare_means(method = test_method)
    }
    
    ggsave(file.path(output_dir, paste0(metric, "_boxplot.png")), p, width = 8, height = 6)
    ggsave(file.path(output_dir, paste0(metric, "_boxplot.pdf")), p, width = 8, height = 6)
}

print("Optimized Alpha diversity analysis complete.")
