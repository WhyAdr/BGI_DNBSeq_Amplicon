# ==============================================================================
# 02_beta_diversity.R
# BGI Amplicon Workflow - Beta Diversity Analysis (Optimized)
# ==============================================================================
# Quantifies compositional dissimilarity between samples (Sections 9 and M8)
# Performs normalization, distance calculation, PCoA, NMDS, and PERMANOVA.
# ==============================================================================

library(vegan)
library(ggplot2)

# --- Configuration ---
if (!exists("otu_file") || is.null(otu_file)) otu_file <- "../BGI_Result/OTU/OTU_table_for_biom.txt"
if (!exists("meta_file") || is.null(meta_file)) meta_file <- "../metadata.tsv"
if (!exists("output_dir") || is.null(output_dir)) output_dir <- "../BGI_Result/Beta"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t",
                  comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

# Align samples
common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- Normalization: Rarefaction (Section 9) ---
# Sub-sampling all samples to the minimum sequencing depth as per BGI report
min_depth <- min(colSums(otu))
set.seed(42)  # For reproducibility
otu_rare <- as.data.frame(t(rrarefy(t(otu), sample = min_depth)))

# --- Beta Diversity Metrics (Bray-Curtis) ---
dist_bray <- vegdist(t(otu_rare), method = "bray")

# --- PCoA Analysis (Section 9) ---
# BGI: "PCoA is computed using QIIME's iterative algorithm: 100 iterations
#       of 75% sub-sampling at the minimum sequencing depth."
# We replicate with Procrustes-aligned consensus across 100 iterations.

n_iter <- 100
sub_frac <- 0.75
sub_depth <- floor(min_depth * sub_frac)

cat(sprintf("Bootstrapped PCoA: %d iterations, sub-depth = %d (%.0f%% of min %d)\n",
            n_iter, sub_depth, sub_frac * 100, min_depth))

set.seed(42)
pcoa_list <- vector("list", n_iter)

for (i in seq_len(n_iter)) {
    r <- rrarefy(t(otu), sub_depth)
    d <- vegdist(r, method = "bray")
    pcoa_list[[i]] <- cmdscale(d, k = 2, eig = FALSE)
}

# Use first iteration as reference, Procrustes-align all others
ref <- pcoa_list[[1]]
aligned <- lapply(pcoa_list, function(coords) {
    proc <- procrustes(ref, coords, symmetric = TRUE)
    proc$Yrot  # Rotated coordinates
})

# Consensus: element-wise mean across all aligned iterations
consensus_coords <- Reduce("+", aligned) / n_iter

pcoa_data <- as.data.frame(consensus_coords)
colnames(pcoa_data) <- c("PCoA1", "PCoA2")
pcoa_data$Group <- metadata$Group
pcoa_data$Sample <- rownames(pcoa_data)

# Variance explained from the full-depth single PCoA (for axis labels)
# Note: Bray-Curtis is non-Euclidean and can produce negative eigenvalues;
# divide by sum of positive eigenvalues only (standard PCoA convention).
pcoa_full <- cmdscale(dist_bray, k = 3, eig = TRUE)
pos_eig <- pcoa_full$eig[pcoa_full$eig > 0]
var_exp <- round(100 * pcoa_full$eig[1:2] / sum(pos_eig), 2)

p_pcoa <- ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, color = Group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95, linetype = 2) +
    theme_bw() +
    labs(title = sprintf("PCoA (Bray-Curtis, %d-iter bootstrap)", n_iter),
         x = paste0("PCoA 1 (", var_exp[1], "%)"),
         y = paste0("PCoA 2 (", var_exp[2], "%)"))

ggsave(file.path(output_dir, "PCoA_BrayCurtis.png"), p_pcoa, width = 8, height = 6)
ggsave(file.path(output_dir, "PCoA_BrayCurtis.pdf"), p_pcoa, width = 8, height = 6)

# --- PERMANOVA (ADONIS) Analysis ---
# Evaluates whether group identity significantly explains total beta diversity variance
if ("Group" %in% colnames(metadata) && length(unique(metadata$Group)) > 1) {
    adonis_res <- adonis2(dist_bray ~ Group, data = metadata, permutations = 999)
    write.table(as.data.frame(adonis_res), file = file.path(output_dir, "PERMANOVA_result.txt"), sep = "\t", quote = FALSE)
}

# --- Pearson Dissimilarity (BGI Section 9, Table 5) ---
# BGI specifies four beta diversity metrics: Bray-Curtis, Weighted UniFrac,
# Unweighted UniFrac, and Pearson. Pearson distance = 1 - Pearson correlation.
# Uses a 10-iteration Procrustes-aligned consensus PCoA for stability.
cat("Computing Pearson dissimilarity...\n")
cor_pearson <- cor(otu_rare, method = "pearson")
dist_pearson <- as.dist(1 - cor_pearson)

# Export Pearson distance matrix
write.table(as.matrix(dist_pearson),
            file.path(output_dir, "Pearson.Beta_diversity.txt"),
            sep = "\t", quote = FALSE)

# 10-iteration bootstrapped PCoA (Pearson)
n_iter_p <- 10
sub_depth_p <- floor(min_depth * 0.75)
cat(sprintf("Bootstrapped Pearson PCoA: %d iterations, sub-depth = %d\n",
            n_iter_p, sub_depth_p))

set.seed(42)
pcoa_p_list <- vector("list", n_iter_p)
for (i in seq_len(n_iter_p)) {
    r <- rrarefy(t(otu), sub_depth_p)
    cor_r <- cor(t(r), method = "pearson")  # r is 51×5358; t(r) is 5358×51; cor on columns → 51×51
    d_r <- as.dist(1 - cor_r)
    pcoa_p_list[[i]] <- cmdscale(d_r, k = 2, eig = FALSE)
}

# Procrustes-align and average
ref_p <- pcoa_p_list[[1]]
aligned_p <- lapply(pcoa_p_list, function(coords) {
    proc <- procrustes(ref_p, coords, symmetric = TRUE)
    proc$Yrot
})
consensus_p <- Reduce("+", aligned_p) / n_iter_p

# Variance explained from single full-depth Pearson PCoA
pcoa_p_full <- cmdscale(dist_pearson, k = 2, eig = TRUE)
pos_eig_p <- pcoa_p_full$eig[pcoa_p_full$eig > 0]
var_exp_p <- round(100 * pcoa_p_full$eig[1:2] / sum(pos_eig_p), 2)

pcoa_p_data <- data.frame(
    PCoA1 = consensus_p[,1],
    PCoA2 = consensus_p[,2],
    Group = metadata$Group,
    Sample = rownames(consensus_p)
)

p_pcoa_p <- ggplot(pcoa_p_data, aes(x = PCoA1, y = PCoA2, color = Group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95, linetype = 2) +
    theme_bw() +
    labs(title = sprintf("PCoA (Pearson, %d-iter bootstrap)", n_iter_p),
         x = paste0("PCoA 1 (", var_exp_p[1], "%)"),
         y = paste0("PCoA 2 (", var_exp_p[2], "%)"))

ggsave(file.path(output_dir, "PCoA_Pearson.png"), p_pcoa_p, width = 8, height = 6)
ggsave(file.path(output_dir, "PCoA_Pearson.pdf"), p_pcoa_p, width = 8, height = 6)

# Pearson PERMANOVA
if ("Group" %in% colnames(metadata) && length(unique(metadata$Group)) >= 2) {
    adonis_pearson <- adonis2(dist_pearson ~ Group, data = metadata, permutations = 999)
    write.table(as.data.frame(adonis_pearson),
                file.path(output_dir, "Pearson_PERMANOVA.txt"),
                sep = "\t", quote = FALSE)
}

# Export Pearson PCoA coordinates
write.table(pcoa_p_data, file.path(output_dir, "Pearson.coordinate.xls"),
            sep = "\t", row.names = FALSE, quote = FALSE)

print("Optimized Beta diversity analysis complete.")
