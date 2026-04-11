# ==============================================================================
# 18_nmds.R
# BGI Amplicon Workflow - Non-Metric Multidimensional Scaling (NMDS)
# ==============================================================================
# Replicates BGI NMDS directory.
# Software reference: R v3.1.1 (vegan package) per BGI report.
# NMDS stress < 0.2 indicates acceptable representation.
# ==============================================================================

library(vegan)
library(ggplot2)

# --- Configuration ---
otu_file <- "../BGI_Result/OTU/OTU_table_for_biom.txt"
meta_file <- "../metadata.tsv"
output_dir <- "../BGI_Result/NMDS"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "#")
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- Rarefaction ---
min_depth <- min(colSums(otu))
set.seed(42)
otu_rare <- as.data.frame(t(rrarefy(t(otu), sample = min_depth)))

# --- NMDS ---
dist_bray <- vegdist(t(otu_rare), method = "bray")
nmds_res <- metaMDS(dist_bray, k = 2, trymax = 100, trace = FALSE)

cat(sprintf("NMDS Stress: %.4f\n", nmds_res$stress))
if (nmds_res$stress >= 0.2) {
    warning("NMDS stress >= 0.2 — ordination may not adequately represent distances.")
}

# --- Plot ---
nmds_df <- data.frame(NMDS1 = nmds_res$points[, 1], NMDS2 = nmds_res$points[, 2],
                       Sample = rownames(nmds_res$points),
                       Group = metadata[rownames(nmds_res$points), "Group"])

p <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95, linetype = 2) +
    annotate("text", x = Inf, y = -Inf,
             label = paste0("Stress = ", round(nmds_res$stress, 4)),
             hjust = 1.1, vjust = -0.5, size = 4) +
    theme_bw() +
    labs(title = "NMDS (Bray-Curtis)")

ggsave(file.path(output_dir, "NMDS_BrayCurtis.png"), p, width = 8, height = 6)
ggsave(file.path(output_dir, "NMDS_BrayCurtis.pdf"), p, width = 8, height = 6)

# Export coordinates
write.table(nmds_df, file.path(output_dir, "NMDS_coordinates.xls"),
            sep = "\t", row.names = FALSE, quote = FALSE)

print("NMDS analysis complete.")
