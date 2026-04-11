# ==============================================================================
# 13_network.R
# BGI Amplicon Workflow - Co-occurrence Network & Correlation Heatmap
# ==============================================================================
# Replicates BGI Network directory.
# Software reference: R v3.4.1; Cytoscape per BGI report.
# Spearman rank correlation, |rho| > 0.2 threshold, network + heatmap.
# ==============================================================================

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(igraph)
library(pheatmap)

# --- Configuration ---
if (!exists("otu_file") || is.null(otu_file)) otu_file <- "../BGI_Result/OTU/OTU_table_for_biom.txt"
if (!exists("meta_file") || is.null(meta_file)) meta_file <- "../metadata.tsv"
if (!exists("output_dir") || is.null(output_dir)) output_dir <- "../BGI_Result/Network"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL

# --- Filter: keep OTUs with relative abundance > 0.5% (per BGI Section 12) ---
otu_rel <- sweep(otu, 2, colSums(otu), "/")
keep <- rowMeans(otu_rel) > 0.005
otu_filt <- otu[keep, ]

if (nrow(otu_filt) < 5) {
    print("Too few OTUs passed filtering. Lowering threshold to 0.1%.")
    keep <- rowMeans(otu_rel) > 0.001
    otu_filt <- otu[keep, ]
}

cat(sprintf("Network analysis using %d OTUs (rel. abund. > 0.5%%).\n", nrow(otu_filt)))

# --- Spearman Correlation (vectorized) ---
if (!requireNamespace("psych", quietly = TRUE)) install.packages("psych")
library(psych)
cor_result <- corr.test(t(otu_filt), method = "spearman", adjust = "none")
cor_mat <- cor_result$r
p_mat <- cor_result$p
n_otus <- nrow(otu_filt)

# --- FDR Correction (Benjamini-Hochberg) ---
# With many pairwise tests, FDR correction prevents false-positive correlations.
p_vec <- p_mat[upper.tri(p_mat)]
p_adj_vec <- p.adjust(p_vec, method = "BH")
p_adj_mat <- matrix(1, n_otus, n_otus)
p_adj_mat[upper.tri(p_adj_mat)] <- p_adj_vec
p_adj_mat[lower.tri(p_adj_mat)] <- t(p_adj_mat)[lower.tri(p_adj_mat)]
diag(p_adj_mat) <- 1

# --- Correlation Heatmap (BGI Section 12: Spearman Correlation Heatmap) ---
# Only display |rho| > 0.2 (per BGI report)
cor_display <- cor_mat
cor_display[abs(cor_mat) <= 0.2] <- 0

# Significance star annotation matrix (FDR-corrected)
sig_stars <- matrix("", n_otus, n_otus)
sig_stars[p_adj_mat < 0.05]  <- "*"
sig_stars[p_adj_mat < 0.01]  <- "**"
sig_stars[p_adj_mat < 0.001] <- "***"
sig_stars[abs(cor_mat) <= 0.2] <- ""  # Blank out non-displayed correlations

pheatmap(cor_display,
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = sig_stars,
         number_color = "black",
         fontsize_number = 6,
         main = "Species Spearman Correlation (|rho| > 0.2, FDR-corrected)",
         filename = file.path(output_dir, "Correlation_Heatmap.png"),
         width = 12, height = 10)

# Export correlation results (with FDR-corrected p-values)
cor_results <- data.frame()
for (i in 1:(n_otus - 1)) {
    for (j in (i + 1):n_otus) {
        if (abs(cor_mat[i, j]) > 0.2) {
            cor_results <- rbind(cor_results, data.frame(
                OTU1 = rownames(otu_filt)[i],
                OTU2 = rownames(otu_filt)[j],
                Spearman_rho = round(cor_mat[i, j], 4),
                P_value = p_mat[i, j],
                FDR = p_adj_mat[i, j]
            ))
        }
    }
}
write.table(cor_results, file.path(output_dir, "Correlation_Result.xls"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Network (BGI Section 12: Network Analysis) ---
# Threshold: |rho| > 0.6 and p < 0.05 for network edges
adj_mat <- cor_mat
adj_mat[abs(cor_mat) < 0.6 | p_adj_mat >= 0.05] <- 0
diag(adj_mat) <- 0

g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE,
                                  diag = FALSE)
# Remove isolated nodes
g <- igraph::delete_vertices(g, igraph::degree(g) == 0)

if (vcount(g) > 0) {
    # Edge colors: pink = positive, blue = negative (per BGI report)
    E(g)$color <- ifelse(E(g)$weight > 0, "#FF69B4", "#4169E1")
    E(g)$width <- abs(E(g)$weight) * 3

    # Node size proportional to mean relative abundance
    node_means <- rowMeans(otu_rel[V(g)$name, ])
    V(g)$size <- scales::rescale(node_means, to = c(3, 15))

    png(file.path(output_dir, "CoOccurrence_Network.png"),
        width = 1000, height = 1000, res = 120)
    plot(g, vertex.label.cex = 0.5, vertex.label.color = "black",
         layout = layout_with_fr(g, weights = abs(E(g)$weight)),
         main = "Co-occurrence Network (|rho|>0.6, P<0.05)")
    legend("bottomleft", legend = c("Positive", "Negative"),
           col = c("#FF69B4", "#4169E1"), lwd = 2, bty = "n")
    dev.off()

    # Export edge list
    edge_df <- as_data_frame(g, what = "edges")
    write.table(edge_df, file.path(output_dir, "Network_edges.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
} else {
    print("No significant edges found at |rho|>0.6, P<0.05 threshold.")
}

print("Network and correlation analysis complete.")
