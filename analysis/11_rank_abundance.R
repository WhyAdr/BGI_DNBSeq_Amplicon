# ==============================================================================
# 11_rank_abundance.R
# BGI Amplicon Workflow - Rank-Abundance & Species Accumulation Curves
# ==============================================================================
# Replicates BGI OTU_Rank + Cumulative_Curve directories.
# Software reference: R v3.1.1 per BGI report; R v3.2.1 for accumulation.
# ==============================================================================

library(ggplot2)
library(vegan)

# --- Configuration ---
otu_file <- "../BGI_Result/OTU/OTU_table_for_biom.txt"
meta_file <- "../metadata.tsv"
rank_dir <- "../BGI_Result/OTU_Rank"
cumul_dir <- "../BGI_Result/Cumulative_Curve"
dir.create(rank_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cumul_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "#")
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- Rank-Abundance per sample (BGI Section 6: OTU Rank Curve) ---
# OTU relative abundances ranked in descending order, log-scale Y-axis
rank_data <- do.call(rbind, lapply(common_samples, function(s) {
    counts <- otu[, s]
    counts <- counts[counts > 0]
    counts <- sort(counts, decreasing = TRUE)
    rel_abund <- counts / sum(counts)
    data.frame(Sample = s, Rank = seq_along(rel_abund),
               RelAbundance = rel_abund,
               Group = metadata[s, "Group"], stringsAsFactors = FALSE)
}))

p_rank <- ggplot(rank_data, aes(x = Rank, y = RelAbundance,
                                 color = Group, group = Sample)) +
    geom_line(alpha = 0.6, linewidth = 0.4) +
    scale_y_log10() +
    theme_bw() +
    labs(title = "OTU Rank-Abundance Curves",
         x = "OTU Rank", y = "Relative Abundance (log scale)") +
    theme(legend.position = "right")

ggsave(file.path(rank_dir, "OTU_Rank_Abundance.png"), p_rank, width = 10, height = 6)
ggsave(file.path(rank_dir, "OTU_Rank_Abundance.pdf"), p_rank, width = 10, height = 6)

# Export rank data
rank_export <- do.call(rbind, lapply(common_samples, function(s) {
    counts <- otu[, s]
    counts <- counts[counts > 0]
    counts <- sort(counts, decreasing = TRUE)
    data.frame(Sample = s, Rank = seq_along(counts),
               OTU_count = counts,
               Percent = round(counts / sum(counts) * 100, 4),
               stringsAsFactors = FALSE)
}))
write.table(rank_export, file.path(rank_dir, "OTU_rank_percent.xls"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Species Accumulation Curve (BGI Section 6: Cumulative Curve) ---
set.seed(42)
spec_accum <- specaccum(t(otu), method = "random", permutations = 100)

# Base R plot (matching BGI output style)
png(file.path(cumul_dir, "Cumulative_Curve.png"), width = 800, height = 600, res = 120)
plot(spec_accum, ci.type = "polygon", ci.col = "lightblue", col = "steelblue", lwd = 2,
     main = "Species Accumulation Curve", xlab = "Number of Samples",
     ylab = "Number of OTUs")
dev.off()

pdf(file.path(cumul_dir, "Cumulative_Curve.pdf"), width = 8, height = 6)
plot(spec_accum, ci.type = "polygon", ci.col = "lightblue", col = "steelblue", lwd = 2,
     main = "Species Accumulation Curve", xlab = "Number of Samples",
     ylab = "Number of OTUs")
dev.off()

print("Rank-abundance and species accumulation curve analysis complete.")
