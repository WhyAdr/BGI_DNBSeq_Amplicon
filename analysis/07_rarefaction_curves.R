# ==============================================================================
# 07_rarefaction_curves.R
# BGI Amplicon Workflow - Rarefaction / Species Accumulation Curves
# ==============================================================================
# Replicates BGI Alpha_Rarefaction directory.
# Software reference: mothur v1.31.2 (BGI M7); R vegan::rarefy()
# ==============================================================================

library(vegan)
library(ggplot2)

# --- Configuration ---
if (!exists("otu_file") || is.null(otu_file)) otu_file <- "../BGI_Result/OTU/OTU_table_for_biom.txt"
if (!exists("meta_file") || is.null(meta_file)) meta_file <- "../metadata.tsv"
if (!exists("output_dir") || is.null(output_dir)) output_dir <- "../BGI_Result/Alpha_Rarefaction"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- Rarefaction Curve Calculation ---
# Transpose: samples as rows for vegan
otu_t <- t(otu)
max_depth <- max(rowSums(otu_t))

# Generate rarefaction data at stepped intervals
step_size <- max(1, floor(max_depth / 50))
steps <- seq(1, max_depth, by = step_size)

rare_data <- do.call(rbind, lapply(rownames(otu_t), function(samp) {
    row <- as.numeric(otu_t[samp, ])
    total <- sum(row)
    valid_steps <- steps[steps <= total]
    sobs <- sapply(valid_steps, function(d) rarefy(row, d))
    data.frame(Sample = samp, Depth = valid_steps, Sobs = sobs,
               Group = metadata[samp, "Group"], stringsAsFactors = FALSE)
}))

# --- Rarefaction Curve Plot ---
p_rare <- ggplot(rare_data, aes(x = Depth, y = Sobs, color = Group, group = Sample)) +
    geom_line(alpha = 0.7, linewidth = 0.5) +
    theme_bw() +
    labs(title = "Rarefaction Curves (Observed Species)",
         x = "Sequencing Depth", y = "Observed OTUs") +
    theme(legend.position = "right")

ggsave(file.path(output_dir, "Rarefaction_Sobs.png"), p_rare, width = 12, height = 7)
ggsave(file.path(output_dir, "Rarefaction_Sobs.pdf"), p_rare, width = 12, height = 7)

# --- Shannon rarefaction ---
rare_shannon <- do.call(rbind, lapply(rownames(otu_t), function(samp) {
    row <- as.numeric(otu_t[samp, ])
    total <- sum(row)
    valid_steps <- steps[steps <= total]
    set.seed(42)  # seed once per sample for reproducible rarefaction
    shannon_vals <- sapply(valid_steps, function(d) {
        sub <- rrarefy(row, d)
        diversity(sub, index = "shannon")
    })
    data.frame(Sample = samp, Depth = valid_steps, Shannon = shannon_vals,
               Group = metadata[samp, "Group"], stringsAsFactors = FALSE)
}))

p_shannon <- ggplot(rare_shannon, aes(x = Depth, y = Shannon, color = Group, group = Sample)) +
    geom_line(alpha = 0.7, linewidth = 0.5) +
    theme_bw() +
    labs(title = "Rarefaction Curves (Shannon Index)",
         x = "Sequencing Depth", y = "Shannon Diversity") +
    theme(legend.position = "right")

ggsave(file.path(output_dir, "Rarefaction_Shannon.png"), p_shannon, width = 12, height = 7)
ggsave(file.path(output_dir, "Rarefaction_Shannon.pdf"), p_shannon, width = 12, height = 7)

print("Rarefaction curve analysis complete.")
