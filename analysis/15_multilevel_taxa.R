# ==============================================================================
# 15_multilevel_taxa.R
# BGI Amplicon Workflow - Multi-level Taxonomic Barplots (L2-L7)
# ==============================================================================
# Replicates full BGI Barplot directory across all classification levels.
# Software reference: R v3.4.1 per BGI report.
# Produces sample-level AND group-level barplots at every level.
# Species < 0.5% relative abundance consolidated to "Others" (per BGI).
# ==============================================================================

library(ggplot2)
library(reshape2)

# --- Configuration ---
meta_file <- "../metadata.tsv"
otu_dir <- "../BGI_Result/OTU"
output_dir <- "../BGI_Result/Barplot"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

level_names <- c("L2" = "Phylum", "L3" = "Class", "L4" = "Order",
                 "L5" = "Family", "L6" = "Genus", "L7" = "Species")

for (lvl in names(level_names)) {
    lvl_file <- file.path(otu_dir, paste0("OTU_table_", lvl, ".txt"))
    if (!file.exists(lvl_file)) {
        cat(sprintf("  Skipping %s: file not found.\n", lvl))
        next
    }

    taxa <- read.table(lvl_file, header = TRUE, row.names = 1, check.names = FALSE,
                       sep = "\t", comment.char = "#")
    if ("taxonomy" %in% colnames(taxa)) taxa$taxonomy <- NULL

    common <- intersect(colnames(taxa), rownames(metadata))
    if (length(common) == 0) next
    taxa <- taxa[, common, drop = FALSE]
    meta_sub <- metadata[common, , drop = FALSE]

    # Relative abundance
    taxa_rel <- sweep(taxa, 2, colSums(taxa), "/")

    # Top N taxa, others < 0.5% consolidated (per BGI Section 7)
    taxa_rel$Mean <- rowMeans(taxa_rel[, common])
    taxa_rel <- taxa_rel[order(-taxa_rel$Mean), ]

    # Determine top taxa (those with mean > 0.5%, up to 15 max)
    above_thresh <- rownames(taxa_rel)[taxa_rel$Mean >= 0.005]
    top_taxa <- head(above_thresh, 15)
    if (length(top_taxa) < 5) top_taxa <- rownames(taxa_rel)[1:min(10, nrow(taxa_rel))]

    others <- setdiff(rownames(taxa_rel), top_taxa)
    plot_mat <- taxa_rel[top_taxa, common, drop = FALSE]
    if (length(others) > 0) {
        plot_mat["Others", ] <- colSums(taxa_rel[others, common, drop = FALSE])
    }

    # --- Sample-level Barplot ---
    plot_df <- melt(as.matrix(plot_mat))
    colnames(plot_df) <- c("Taxa", "Sample", "Abundance")
    plot_df <- merge(plot_df, meta_sub, by.x = "Sample", by.y = 1)

    p_sample <- ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Taxa)) +
        geom_bar(stat = "identity", position = "stack") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) +
        labs(title = paste0("Taxonomic Composition (", level_names[lvl], " Level)"),
             y = "Relative Abundance")

    if ("Group" %in% colnames(plot_df)) {
        p_sample <- p_sample + facet_grid(~ Group, scales = "free_x", space = "free_x")
    }

    ggsave(file.path(output_dir, paste0("Barplot.Sample.", lvl, ".", level_names[lvl], ".png")),
           p_sample, width = 16, height = 8)
    ggsave(file.path(output_dir, paste0("Barplot.Sample.", lvl, ".", level_names[lvl], ".pdf")),
           p_sample, width = 16, height = 8)

    # --- Group-level Barplot ---
    if ("Group" %in% colnames(meta_sub)) {
        group_means <- aggregate(t(plot_mat[, common]), 
                                  by = list(Group = meta_sub$Group), FUN = mean)
        group_melt <- melt(group_means, id.vars = "Group")
        colnames(group_melt) <- c("Group", "Taxa", "Abundance")

        p_group <- ggplot(group_melt, aes(x = Group, y = Abundance, fill = Taxa)) +
            geom_bar(stat = "identity", position = "stack") +
            theme_bw() +
            labs(title = paste0("Taxonomic Composition - Group Mean (", level_names[lvl], ")"),
                 y = "Mean Relative Abundance")

        ggsave(file.path(output_dir, paste0("Barplot.Group.", lvl, ".", level_names[lvl], ".png")),
               p_group, width = 10, height = 7)
        ggsave(file.path(output_dir, paste0("Barplot.Group.", lvl, ".", level_names[lvl], ".pdf")),
               p_group, width = 10, height = 7)
    }

    # --- Export data tables ---
    write.table(data.frame(Taxa = rownames(plot_mat), plot_mat),
                file.path(output_dir, paste0("Sample.", lvl, ".", level_names[lvl], ".Barplot.xls")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    cat(sprintf("  Completed %s (%s) - %d taxa plotted.\n", lvl, level_names[lvl], nrow(plot_mat)))
}

print("Multi-level taxonomic barplot analysis complete.")
