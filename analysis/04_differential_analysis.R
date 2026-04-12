# ==============================================================================
# 04_differential_analysis.R
# BGI Amplicon Workflow - Differential Species Analysis (Optimized)
# ==============================================================================
# Identifies biomarkers and significant features (Sections 10 and M9)
# Performs Wilcoxon or Kruskal-Wallis non-parametric tests on relative abundances.
# ==============================================================================

library(ggplot2)

# --- Configuration ---
if (!exists("otu_file") || is.null(otu_file)) otu_file <- "../BGI_Result/OTU/OTU_table_for_biom.txt"
if (!exists("tax_file") || is.null(tax_file)) tax_file <- "../BGI_Result/OTU/OTU_taxonomy.xls"
if (!exists("meta_file") || is.null(meta_file)) meta_file <- "../metadata.tsv"
if (!exists("output_dir") || is.null(output_dir)) output_dir <- "../BGI_Result/Diff"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t",
                  comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
tax <- read.table(tax_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t",
                  comment.char = "")
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

groups <- unique(metadata$Group)
if (length(groups) < 2) {
    stop("Error: Differential analysis requires at least 2 distinct groups in metadata.")
}

# --- Pre-filtering for Statistical Power (Optimization) ---
# Removing extremely rare OTUs drastically reduces the multiple testing penalty 
# (FDR correction via Benjamini-Hochberg).
# Keep OTUs present in at least 10% of samples or with mean relative abundance > 0.01%
rel_abund <- sweep(otu, 2, colSums(otu), "/")
prevalence <- rowSums(otu > 0) / ncol(otu)
keep_otus <- rownames(otu)[rowMeans(rel_abund) > 0.0001 | prevalence >= 0.1]

cat(sprintf("Filtering: Retaining %d out of %d OTUs for testing.\n", length(keep_otus), nrow(otu)))
otu_filtered <- otu[keep_otus, , drop = FALSE]

# --- Statistical Testing (Section 10) ---
# The Wilcoxon Rank-Sum test (Mann-Whitney U) is used to identify significantly 
# different species between two independent groups.
# The Kruskal-Wallis test is used across three or more groups.
results <- data.frame(OTU = keep_otus)

if (length(groups) == 2) {
    # Wilcoxon test (2 groups)
    p_values <- apply(otu_filtered, 1, function(x) wilcox.test(x ~ metadata$Group)$p.value)
} else {
    # Kruskal-Wallis test (3+ groups)
    p_values <- apply(otu_filtered, 1, function(x) kruskal.test(x ~ metadata$Group)$p.value)
}

results$p_value <- p_values
# FDR Correction only on the filtered high-quality OTUs to preserve actual biomarkers
results$FDR <- p.adjust(p_values, method = "BH")
results$Taxonomy <- tax[keep_otus, "Taxonomy"]

# --- Log₂ Fold-Change (BGI Section 10 requirement) ---
group_lvls <- sort(unique(metadata$Group))
group_means <- sapply(group_lvls, function(g)
    rowMeans(rel_abund[keep_otus, metadata$Group == g, drop = FALSE]))
pseudo <- 1e-9  # Pseudocount to avoid log(0)

if (length(group_lvls) == 2) {
    # Two-group: log2(Group2 / Group1)
    results$log2FC <- log2((group_means[, 2] + pseudo) / (group_means[, 1] + pseudo))
    results$Comparison <- paste0(group_lvls[2], "_vs_", group_lvls[1])
} else {
    # Multi-group: report max |log2FC| across all pairwise comparisons
    n_grp <- length(group_lvls)
    max_lfc <- rep(0, length(keep_otus))
    max_pair <- rep("", length(keep_otus))
    for (a in 1:(n_grp - 1)) {
        for (b in (a + 1):n_grp) {
            lfc <- log2((group_means[, b] + pseudo) / (group_means[, a] + pseudo))
            update <- abs(lfc) > abs(max_lfc)
            max_lfc[update] <- lfc[update]
            max_pair[update] <- paste0(group_lvls[b], "_vs_", group_lvls[a])
        }
    }
    results$log2FC <- max_lfc
    results$MaxPair <- max_pair
}

results <- results[order(results$p_value), ]
write.table(results, file = file.path(output_dir, "differential_species_stats.xls"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Histogram of Key Species Difference Comparison (Section 10) ---
# Top significant abundance comparisons are boxed/bar graphed contextually.
sig_results <- results[!is.na(results$FDR) & results$FDR < 0.05, ]
if (nrow(sig_results) > 0) {
    # Extract top 10 for visualization
    top_diff <- head(sig_results$OTU, 10)
    
    for (td_otu in top_diff) {
        dt <- data.frame(Abundance = rel_abund[td_otu, ], Group = metadata$Group)
        p_diff <- ggplot(dt, aes(x = Group, y = Abundance, fill = Group)) +
            geom_boxplot() + theme_bw() +
            labs(title = paste("Differential Abundance:", td_otu), y = "Relative Abundance")
        ggsave(file.path(output_dir, paste0(td_otu, "_diff_boxplot.png")), p_diff, width = 6, height = 5)
    }
}

# --- Structural preparation for LEfSe Analysis (Strict TXT generation) ---
lefse_dir <- "../BGI_Result/Lefse"
dir.create(lefse_dir, showWarnings = FALSE, recursive = TRUE)

# Compute multi-level hierarchical aggregation
# Create full taxonomy strings
tax_strings <- tax[rownames(rel_abund), "Taxonomy"]
# Replace ; with | and remove whitespace
tax_strings <- gsub("; *", "|", tax_strings)

agg_abund_list <- list()

for (i in seq_along(tax_strings)) {
    full_tax <- tax_strings[i]
    if (is.na(full_tax) || full_tax == "") next
    
    parts <- unlist(strsplit(full_tax, "\\|"))
    for (j in seq_along(parts)) {
        sub_tax <- paste(parts[1:j], collapse = "|")
        if (is.null(agg_abund_list[[sub_tax]])) {
            agg_abund_list[[sub_tax]] <- rel_abund[i, ]
        } else {
            agg_abund_list[[sub_tax]] <- agg_abund_list[[sub_tax]] + rel_abund[i, ]
        }
    }
}

agg_abund <- do.call(rbind, agg_abund_list)
agg_abund <- agg_abund * 100 # Multiply relative abundance * 100 for LEfSe

# BGI schema format
out_lefse <- rbind(
    c("Group", as.character(metadata$Group)),
    c("Sample", rownames(metadata)),
    cbind(rownames(agg_abund), format(agg_abund, scientific=FALSE))
)
comp_suffix <- paste(sort(unique(metadata$Group)), collapse = "-")
write.table(out_lefse, file.path(lefse_dir, paste0(comp_suffix, ".OTU_tax_assignments.txt")), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# ==============================================================================
# --- Multi-level Differential Analysis (BGI Section 10) ---
# ==============================================================================
# BGI output path: Diff/wilcoxon.test/each_level/{Phylum,Class,...,Species}/
# Tests are run on taxonomically aggregated abundance tables at each level.
if (!exists("otu_dir") || is.null(otu_dir)) otu_dir <- "../BGI_Result/OTU"
level_map <- c("L2" = "Phylum", "L3" = "Class", "L4" = "Order",
               "L5" = "Family", "L6" = "Genus", "L7" = "Species")

for (lvl in names(level_map)) {
    lvl_name <- level_map[lvl]
    lvl_file <- file.path(otu_dir, paste0("OTU_table_", lvl, ".txt"))

    if (!file.exists(lvl_file)) {
        cat(sprintf("  Skipping differential at %s: file not found.\n", lvl_name))
        next
    }

    lvl_otu <- read.table(lvl_file, header = TRUE, row.names = 1,
                          check.names = FALSE, sep = "\t", comment.char = "")
    if ("taxonomy" %in% colnames(lvl_otu)) lvl_otu$taxonomy <- NULL

    lvl_common <- intersect(colnames(lvl_otu), common_samples)
    if (length(lvl_common) < 3) next
    lvl_otu <- lvl_otu[, lvl_common, drop = FALSE]
    lvl_rel <- sweep(lvl_otu, 2, colSums(lvl_otu), "/")

    # Filter: mean relative abundance > 0.01% or prevalence >= 10%
    prev <- rowSums(lvl_otu > 0) / ncol(lvl_otu)
    keep <- rownames(lvl_otu)[rowMeans(lvl_rel) > 0.0001 | prev >= 0.1]
    if (length(keep) < 2) next

    lvl_filtered <- lvl_otu[keep, , drop = FALSE]

    lvl_meta <- metadata[lvl_common, , drop = FALSE]
    lvl_groups <- unique(lvl_meta$Group)

    if (length(lvl_groups) == 2) {
        pvals <- apply(lvl_filtered, 1, function(x)
            tryCatch(wilcox.test(x ~ lvl_meta$Group)$p.value,
                     error = function(e) NA))
        test_label <- "wilcoxon"
    } else {
        pvals <- apply(lvl_filtered, 1, function(x)
            tryCatch(kruskal.test(x ~ lvl_meta$Group)$p.value,
                     error = function(e) NA))
        test_label <- "kruskal"
    }

    lvl_results <- data.frame(
        Feature = keep,
        p_value = pvals,
        FDR = p.adjust(pvals, method = "BH")
    )

    # Log2FC for this level
    lvl_grp_means <- sapply(lvl_groups, function(g)
        rowMeans(lvl_rel[keep, lvl_meta$Group == g, drop = FALSE]))
    if (length(lvl_groups) == 2) {
        lvl_results$log2FC <- log2((lvl_grp_means[, 2] + 1e-9) /
                                    (lvl_grp_means[, 1] + 1e-9))
    } else {
        # Multi-group: max |LFC| pairwise
        n_g <- length(lvl_groups)
        mlfc <- rep(0, length(keep))
        mpair <- rep("", length(keep))
        for (a in 1:(n_g - 1)) {
            for (b in (a + 1):n_g) {
                lfc_lvl <- log2((lvl_grp_means[, b] + 1e-9) / (lvl_grp_means[, a] + 1e-9))
                upd <- abs(lfc_lvl) > abs(mlfc)
                mlfc[upd] <- lfc_lvl[upd]
                mpair[upd] <- paste0(lvl_groups[b], "_vs_", lvl_groups[a])
            }
        }
        lvl_results$log2FC <- mlfc
        lvl_results$MaxPair <- mpair
    }

    lvl_results <- lvl_results[order(lvl_results$p_value), ]

    lvl_out <- file.path(output_dir, paste0(test_label, ".test"), "each_level", lvl_name)
    dir.create(lvl_out, showWarnings = FALSE, recursive = TRUE)
    write.table(lvl_results, file.path(lvl_out, paste0(lvl_name, "_diff.xls")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    cat(sprintf("  Differential at %s: %d features tested, %d significant (FDR<0.05).\n",
                lvl_name, nrow(lvl_results), sum(lvl_results$FDR < 0.05, na.rm = TRUE)))
}

print("Differential analysis (OTU-level + multi-taxonomy) complete.")
