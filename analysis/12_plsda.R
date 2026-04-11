# ==============================================================================
# 12_plsda.R
# BGI Amplicon Workflow - Partial Least Squares Discriminant Analysis
# ==============================================================================
# Replicates BGI PLSDA directory.
# Software reference: R v3.2.1 (mixOmics package) per BGI report.
# PLS-DA is a supervised linear classification model that maximises group
# separation (Barker & Rayens 2003).
# ==============================================================================

if (!requireNamespace("mixOmics", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("mixOmics")
}
library(mixOmics)
library(ggplot2)

# --- Configuration ---
if (!exists("otu_file") || is.null(otu_file)) otu_file <- "../BGI_Result/OTU/OTU_table_for_biom.txt"
if (!exists("meta_file") || is.null(meta_file)) meta_file <- "../metadata.tsv"
if (!exists("output_dir") || is.null(output_dir)) output_dir <- "../BGI_Result/PLSDA"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
X <- t(otu[, common_samples])
Y <- as.factor(metadata[common_samples, "Group"])

# --- PLS-DA ---
plsda_res <- plsda(X, Y, ncomp = 2)
var_exp <- round(plsda_res$prop_expl_var$X * 100, 2)

plsda_df <- data.frame(Comp1 = plsda_res$variates$X[, 1],
                        Comp2 = plsda_res$variates$X[, 2],
                        Group = Y)

p <- ggplot(plsda_df, aes(x = Comp1, y = Comp2, color = Group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95, linetype = 2) +
    theme_bw() +
    labs(title = "PLS-DA (OTU-level)",
         x = paste0("Component 1 (", var_exp[1], "%)"),
         y = paste0("Component 2 (", var_exp[2], "%)"))

ggsave(file.path(output_dir, "PLSDA_plot.png"), p, width = 8, height = 6)
ggsave(file.path(output_dir, "PLSDA_plot.pdf"), p, width = 8, height = 6)

print("PLS-DA analysis complete.")
