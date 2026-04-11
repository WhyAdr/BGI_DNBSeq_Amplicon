# ==============================================================================
# 17_unifrac_beta.R
# BGI Amplicon Workflow - UniFrac Beta Diversity + UPGMA Trees
# ==============================================================================
# Replicates BGI Beta/ UniFrac analyses.
# Software reference: QIIME v1.80 per BGI report; FastTree v2.1.3 for
# phylogenetic tree construction.
# Uses the BGI-generated phylogenetic tree to compute weighted & unweighted
# UniFrac distances, PCoA, UPGMA trees, and beta boxplots.
# ==============================================================================

# Required: phyloseq or GUniFrac for UniFrac calculation
if (!requireNamespace("phyloseq", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("phyloseq")
}
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)

# --- Configuration ---
if (!exists("otu_file") || is.null(otu_file)) otu_file <- "../BGI_Result/OTU/OTU_table_for_biom.txt"
if (!exists("meta_file") || is.null(meta_file)) meta_file <- "../metadata.tsv"
# BGI provides FastTree-generated phylogenetic trees in two locations:
#   Beta/<comparison>/ — OTU-level trees
#   Genus_Tree/        — genus-level trees (per comparison)
if (!exists("tree_dir_beta") || is.null(tree_dir_beta)) tree_dir_beta <- "../BGI_Result/Beta"
if (!exists("tree_dir_genus") || is.null(tree_dir_genus)) tree_dir_genus <- "../BGI_Result/Genus_Tree"
if (!exists("output_dir") || is.null(output_dir)) output_dir <- "../BGI_Result/Beta"

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- Find the phylogenetic tree ---
# Auto-select comparison-specific tree when run via 00_run_all_groups.R
# (which injects 'comp_suffix' into the environment)
find_tree <- function(comp_name = NULL) {
    # 1) Try comparison-specific OTU tree in Beta/
    if (!is.null(comp_name) && comp_name != "") {
        specific_beta <- file.path(tree_dir_beta, comp_name,
                                   paste0(comp_name, ".OTU_final_phylogeny_tree.txt"))
        if (file.exists(specific_beta)) {
            cat(sprintf("Using comparison-specific OTU tree: %s\n", specific_beta))
            return(specific_beta)
        }
        # 2) Try comparison-specific genus tree in Genus_Tree/
        specific_genus <- file.path(tree_dir_genus,
                                    paste0(comp_name, ".genus.phylogeny.tree"))
        if (file.exists(specific_genus)) {
            cat(sprintf("Using comparison-specific genus tree: %s\n", specific_genus))
            return(specific_genus)
        }
    }

    # 3) Fall back to all-sample OTU tree
    all_beta <- file.path(tree_dir_beta,
                          "A-B-C-D-E-F-G-H-I-J-K-L-M-N-O-P-Q",
                          "A-B-C-D-E-F-G-H-I-J-K-L-M-N-O-P-Q.OTU_final_phylogeny_tree.txt")
    if (file.exists(all_beta)) return(all_beta)

    # 4) Fall back to all-sample genus tree
    all_genus <- file.path(tree_dir_genus,
                           "A-B-C-D-E-F-G-H-I-J-K-L-M-N-O-P-Q.genus.phylogeny.tree")
    if (file.exists(all_genus)) return(all_genus)

    # 5) Search both directories for any available tree
    candidates <- c(
        list.files(tree_dir_beta, pattern = "phylogeny_tree\\.txt$",
                   recursive = TRUE, full.names = TRUE),
        list.files(tree_dir_genus, pattern = "\\.phylogeny\\.tree$",
                   recursive = TRUE, full.names = TRUE)
    )
    if (length(candidates) > 0) {
        cat(sprintf("Using fallback tree: %s\n", candidates[1]))
        return(candidates[1])
    }

    stop("No phylogenetic tree file found in Beta/ or Genus_Tree/. Cannot compute UniFrac.")
}

# comp_suffix is injected by 00_run_all_groups.R; defaults to NULL for standalone runs
comp_name_env <- if (exists("comp_suffix")) comp_suffix else NULL
tree_file <- find_tree(comp_name_env)
cat(sprintf("Final tree selection: %s\n", tree_file))

# --- Read tree ---
tree <- read.tree(tree_file)

# --- Prune OTU table to match tree tips ---
shared_otus <- intersect(rownames(otu), tree$tip.label)
cat(sprintf("OTU table: %d OTUs. Tree tips: %d. Shared: %d.\n",
            nrow(otu), length(tree$tip.label), length(shared_otus)))

if (length(shared_otus) < 10) {
    stop("Too few shared OTUs between table and tree. Check OTU ID consistency.")
}

otu_pruned <- otu[shared_otus, ]
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, shared_otus))

# --- Rarefaction (per BGI Section 9: subsample to minimum depth) ---
min_depth <- min(colSums(otu_pruned))
set.seed(42)
otu_rare <- as.data.frame(t(rrarefy(t(otu_pruned), sample = min_depth)))

# --- Build phyloseq object ---
ps <- phyloseq(otu_table(as.matrix(otu_rare), taxa_are_rows = TRUE),
               sample_data(metadata),
               phy_tree(tree_pruned))

# --- Compute UniFrac distances ---
cat("Computing Unweighted UniFrac...\n")
dist_uw <- UniFrac(ps, weighted = FALSE)

cat("Computing Weighted UniFrac...\n")
dist_w <- UniFrac(ps, weighted = TRUE)

# --- Also compute Bray-Curtis for comparison ---
dist_bc <- vegdist(t(otu_rare), method = "bray")

# --- Helper: PCoA + PERMANOVA per metric ---
plot_pcoa_unifrac <- function(dist_mat, metric_name, meta, out_dir) {
    pcoa <- cmdscale(dist_mat, k = min(3, nrow(as.matrix(dist_mat)) - 1), eig = TRUE)
    var_exp <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 2)

    pcoa_df <- data.frame(PCoA1 = pcoa$points[, 1], PCoA2 = pcoa$points[, 2],
                           Sample = rownames(pcoa$points),
                           Group = meta[rownames(pcoa$points), "Group"])

    p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
        geom_point(size = 3) +
        stat_ellipse(level = 0.95, linetype = 2) +
        theme_bw() +
        labs(title = paste0("PCoA (", metric_name, ")"),
             x = paste0("PCoA1 (", var_exp[1], "%)"),
             y = paste0("PCoA2 (", var_exp[2], "%)"))

    ggsave(file.path(out_dir, paste0(metric_name, ".PCoA.png")), p, width = 8, height = 6)
    ggsave(file.path(out_dir, paste0(metric_name, ".PCoA.pdf")), p, width = 8, height = 6)

    # PERMANOVA
    if (length(unique(meta$Group)) >= 2) {
        adonis_res <- adonis2(dist_mat ~ Group, data = meta, permutations = 999)
        write.table(as.data.frame(adonis_res),
                    file.path(out_dir, paste0(metric_name, ".permanova.test.xls")),
                    sep = "\t", quote = FALSE)
    }

    # Export distance matrix
    write.table(as.matrix(dist_mat),
                file.path(out_dir, paste0(metric_name, ".Beta_diversity.txt")),
                sep = "\t", quote = FALSE)

    # Export coordinates
    write.table(pcoa_df, file.path(out_dir, paste0(metric_name, ".coordinate.xls")),
                sep = "\t", row.names = FALSE, quote = FALSE)
}

# --- UPGMA Tree ---
plot_upgma <- function(dist_mat, metric_name, meta, out_dir) {
    hc <- hclust(dist_mat, method = "average")  # UPGMA = average linkage
    dend <- as.dendrogram(hc)

    # Color by group
    sample_groups <- meta[hc$labels[hc$order], "Group"]
    n_groups <- length(unique(meta$Group))
    group_cols <- rainbow(n_groups)
    names(group_cols) <- sort(unique(meta$Group))

    png(file.path(out_dir, paste0(metric_name, ".UPGMA_tree.png")),
        width = 1200, height = 800, res = 120)
    plot(hc, main = paste0("UPGMA Clustering Tree (", metric_name, ")"),
         xlab = "", sub = "", cex = 0.7)
    dev.off()

    pdf(file.path(out_dir, paste0(metric_name, ".UPGMA_tree.pdf")),
        width = 12, height = 8)
    plot(hc, main = paste0("UPGMA Clustering Tree (", metric_name, ")"),
         xlab = "", sub = "", cex = 0.7)
    dev.off()
}

# --- Beta Boxplot ---
plot_beta_box <- function(dist_mat, metric_name, meta, out_dir) {
    dist_m <- as.matrix(dist_mat)
    groups <- sort(unique(meta$Group))

    box_data <- do.call(rbind, lapply(groups, function(g) {
        samps <- rownames(meta)[meta$Group == g]
        if (length(samps) < 2) return(NULL)
        pairs <- combn(samps, 2)
        dists <- apply(pairs, 2, function(x) dist_m[x[1], x[2]])
        data.frame(Group = g, Distance = dists, stringsAsFactors = FALSE)
    }))

    if (is.null(box_data) || nrow(box_data) == 0) return()

    p <- ggplot(box_data, aes(x = Group, y = Distance, fill = Group)) +
        geom_boxplot(alpha = 0.8) +
        theme_bw() +
        labs(title = paste0("Intra-group Beta Diversity (", metric_name, ")"),
             y = paste0(metric_name, " Distance"))

    ggsave(file.path(out_dir, paste0(metric_name, ".Beta.Box.png")), p, width = 10, height = 6)

    write.table(box_data, file.path(out_dir, paste0(metric_name, ".Beta_Box.xls")),
                sep = "\t", row.names = FALSE, quote = FALSE)
}

# --- Run all metrics ---
all_comp_dir <- file.path(output_dir, "R_Analysis")
dir.create(all_comp_dir, showWarnings = FALSE, recursive = TRUE)

metrics <- list(
    "unweighted_unifrac" = dist_uw,
    "weighted_unifrac" = dist_w,
    "bray_curtis" = dist_bc
)

for (metric_name in names(metrics)) {
    cat(sprintf("\n=== %s ===\n", metric_name))
    plot_pcoa_unifrac(metrics[[metric_name]], metric_name, metadata, all_comp_dir)
    plot_upgma(metrics[[metric_name]], metric_name, metadata, all_comp_dir)
    plot_beta_box(metrics[[metric_name]], metric_name, metadata, all_comp_dir)
}

print("UniFrac beta diversity analysis (weighted + unweighted + Bray-Curtis) complete.")
