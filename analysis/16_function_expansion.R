# ==============================================================================
# 16_function_expansion.R
# BGI Amplicon Workflow - COG, MetaCyc, EC Pathway Visualization + Diff
# ==============================================================================
# Extends 05_function_prediction.R to cover all PICRUSt2 databases.
# Replicates BGI Picrust/Function_Prdeict/{COG,EC,METACYC} and Function_Diff.
# Software reference: PICRUSt2 v2.3.0-b; R v3.4.1 per BGI report.
# ==============================================================================

library(ggplot2)
library(reshape2)
library(pheatmap)

# --- Configuration ---
meta_file <- "../metadata.tsv"
base_dir <- "../BGI_Result/Picrust/Function_Prdeict"
diff_dir <- "../BGI_Result/Picrust/Function_Diff"
output_dir <- "../BGI_Result/Picrust"
dir.create(diff_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

# --- Pathway databases to process ---
pathways <- list(
    KO = list(dir = file.path(base_dir, "KO"),
              files = c("ko_Level1_Function.xls", "ko_Level2_Function.xls", "ko_Level3_Function.xls")),
    COG = list(dir = file.path(base_dir, "COG"),
               files = NULL),  # will auto-detect
    EC = list(dir = file.path(base_dir, "EC"),
              files = NULL),
    METACYC = list(dir = file.path(base_dir, "METACYC"),
                   files = NULL)
)

for (pw_name in names(pathways)) {
    pw_info <- pathways[[pw_name]]
    if (!dir.exists(pw_info$dir)) {
        cat(sprintf("  Skipping %s: directory not found.\n", pw_name))
        next
    }

    # Find available files
    if (is.null(pw_info$files)) {
        pw_files <- list.files(pw_info$dir, pattern = "\\.xls$", full.names = TRUE)
    } else {
        pw_files <- file.path(pw_info$dir, pw_info$files)
        pw_files <- pw_files[file.exists(pw_files)]
    }

    if (length(pw_files) == 0) {
        cat(sprintf("  Skipping %s: no .xls files found.\n", pw_name))
        next
    }

    for (pw_file in pw_files) {
        level_name <- gsub("\\.xls$", "", basename(pw_file))
        cat(sprintf("  Processing %s / %s...\n", pw_name, level_name))

        pw_data <- tryCatch(
            read.table(pw_file, header = TRUE, row.names = 1,
                       check.names = FALSE, sep = "\t", comment.char = "#"),
            error = function(e) { cat(sprintf("    Error reading: %s\n", e$message)); NULL }
        )
        if (is.null(pw_data)) next

        common <- intersect(colnames(pw_data), rownames(metadata))
        if (length(common) == 0) next
        pw_data <- pw_data[, common, drop = FALSE]
        meta_sub <- metadata[common, , drop = FALSE]

        # Relative abundance
        pw_rel <- sweep(pw_data, 2, colSums(pw_data), "/")
        pw_rel$Mean <- rowMeans(pw_rel[, common])
        pw_rel <- pw_rel[order(-pw_rel$Mean), ]
        top_n <- min(20, nrow(pw_rel))
        top_names <- rownames(pw_rel)[1:top_n]

        # --- Barplot (sample level) ---
        barplot_dir <- file.path(pw_info$dir, "barplot")
        dir.create(barplot_dir, showWarnings = FALSE, recursive = TRUE)

        plot_df <- melt(as.matrix(pw_rel[top_names, common, drop = FALSE]))
        colnames(plot_df) <- c("Pathway", "Sample", "Abundance")
        plot_df <- merge(plot_df, meta_sub, by.x = "Sample", by.y = 1)

        p <- ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Pathway)) +
            geom_bar(stat = "identity", position = "stack") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
                  legend.text = element_text(size = 6)) +
            labs(title = paste0("Top ", top_n, " ", pw_name, " (", level_name, ")"),
                 y = "Relative Abundance")

        ggsave(file.path(barplot_dir, paste0("Sample.", pw_name, ".", level_name, ".png")),
               p, width = 14, height = 8)

        # --- Group-level barplot ---
        if ("Group" %in% colnames(meta_sub)) {
            group_means <- aggregate(t(pw_rel[top_names, common]),
                                      by = list(Group = meta_sub$Group), FUN = mean)
            group_melt <- melt(group_means, id.vars = "Group")
            colnames(group_melt) <- c("Group", "Pathway", "Abundance")

            p_g <- ggplot(group_melt, aes(x = Group, y = Abundance, fill = Pathway)) +
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                labs(title = paste0("Group Mean - ", pw_name, " (", level_name, ")"),
                     y = "Mean Relative Abundance")
            ggsave(file.path(barplot_dir, paste0("Group.", pw_name, ".", level_name, ".png")),
                   p_g, width = 10, height = 7)
        }

        # --- Heatmap (log10) ---
        heatmap_dir <- file.path(pw_info$dir, "heatmap")
        dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

        heatmap_mat <- pw_rel[top_names, common, drop = FALSE]
        pseudo <- min(heatmap_mat[heatmap_mat > 0]) * 0.5
        if (is.na(pseudo) || pseudo == 0) pseudo <- 1e-6

        tryCatch({
            pheatmap(log10(heatmap_mat + pseudo),
                     cluster_cols = TRUE, cluster_rows = TRUE,
                     clustering_method = "complete",
                     main = paste0(pw_name, " ", level_name, " (log10)"),
                     filename = file.path(heatmap_dir,
                                          paste0("Heatmap.", pw_name, ".", level_name, ".png")),
                     width = 12, height = 8)
        }, error = function(e) cat(sprintf("    Heatmap error: %s\n", e$message)))

        # --- Differential testing (BGI Section 11: Function Diff) ---
        if ("Group" %in% colnames(meta_sub) && length(unique(meta_sub$Group)) >= 2) {
            diff_subdir <- file.path(diff_dir, pw_name, level_name)
            dir.create(diff_subdir, showWarnings = FALSE, recursive = TRUE)

            groups <- unique(meta_sub$Group)
            if (length(groups) == 2) {
                p_vals <- apply(pw_data[, common], 1, function(x)
                    tryCatch(wilcox.test(x ~ meta_sub$Group)$p.value, error = function(e) NA))
                test_name <- "wilcox.test"
            } else {
                p_vals <- apply(pw_data[, common], 1, function(x)
                    tryCatch(kruskal.test(x ~ meta_sub$Group)$p.value, error = function(e) NA))
                test_name <- "kruskal.test"
            }

            diff_results <- data.frame(
                Pathway = rownames(pw_data),
                P_value = p_vals,
                FDR = p.adjust(p_vals, method = "BH")
            )
            diff_results <- diff_results[order(diff_results$P_value), ]
            write.table(diff_results,
                        file.path(diff_subdir, paste0(test_name, "_result.xls")),
                        sep = "\t", row.names = FALSE, quote = FALSE)
        }
    }
}

print("Functional pathway expansion (COG/EC/MetaCyc + differential) complete.")
