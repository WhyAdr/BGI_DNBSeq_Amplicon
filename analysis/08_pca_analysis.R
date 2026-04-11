# ==============================================================================
# 08_pca_analysis.R
# BGI Amplicon Workflow - Principal Component Analysis
# ==============================================================================
# Replicates BGI PCA directory.
# Software reference: R v3.1.1 (ade4 package) per BGI report.
# We use prcomp() (base R, equivalent) for OTU-level PCA,
# plus taxon-level PCA across L2-L7.
# ==============================================================================

library(ggplot2)

# --- Configuration ---
if (!exists("otu_file") || is.null(otu_file)) otu_file <- "../BGI_Result/OTU/OTU_table_for_biom.txt"
if (!exists("meta_file") || is.null(meta_file)) meta_file <- "../metadata.tsv"
if (!exists("otu_dir") || is.null(otu_dir)) otu_dir <- "../BGI_Result/OTU"
if (!exists("output_dir") || is.null(output_dir)) output_dir <- "../BGI_Result/PCA"
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

# --- Helper: PCA plotting function ---
run_pca_plot <- function(data_mat, samples, meta, title_suffix, out_prefix) {
    # Relative abundance
    data_rel <- sweep(data_mat, 2, colSums(data_mat), "/")
    # Remove zero-variance features
    data_rel <- data_rel[apply(data_rel, 1, var) > 0, , drop = FALSE]

    pca_res <- prcomp(t(data_rel), center = TRUE, scale. = TRUE)
    var_exp <- round(100 * summary(pca_res)$importance[2, 1:2], 2)

    pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
                         Sample = rownames(pca_res$x),
                         Group = meta[rownames(pca_res$x), "Group"])

    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(size = 3) +
        stat_ellipse(level = 0.95, linetype = 2) +
        theme_bw() +
        labs(title = paste0("PCA (", title_suffix, ")"),
             x = paste0("PC1 (", var_exp[1], "%)"),
             y = paste0("PC2 (", var_exp[2], "%)"))

    ggsave(file.path(output_dir, paste0(out_prefix, "_PCA.png")), p, width = 8, height = 6)
    ggsave(file.path(output_dir, paste0(out_prefix, "_PCA.pdf")), p, width = 8, height = 6)
}

# --- OTU-level PCA ---
run_pca_plot(otu, common_samples, metadata, "OTU-level", "OTU")

# --- Taxon-level PCA (L2-L7) ---
level_names <- c("L2" = "Phylum", "L3" = "Class", "L4" = "Order",
                 "L5" = "Family", "L6" = "Genus", "L7" = "Species")

for (lvl in names(level_names)) {
    lvl_file <- file.path(otu_dir, paste0("OTU_table_", lvl, ".txt"))
    if (!file.exists(lvl_file)) next

    taxa <- read.table(lvl_file, header = TRUE, row.names = 1, check.names = FALSE,
                       sep = "\t", comment.char = "")
    if ("taxonomy" %in% colnames(taxa)) taxa$taxonomy <- NULL

    common <- intersect(colnames(taxa), rownames(metadata))
    if (length(common) < 3) next
    taxa <- taxa[, common, drop = FALSE]

    run_pca_plot(taxa, common, metadata, paste0(level_names[lvl], " Level"),
                 paste0(lvl, "_", level_names[lvl]))
}

print("PCA analysis complete.")
