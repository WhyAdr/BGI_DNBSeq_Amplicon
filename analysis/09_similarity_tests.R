# ==============================================================================
# 09_similarity_tests.R
# BGI Amplicon Workflow - ANOSIM + MRPP Hypothesis Tests
# ==============================================================================
# Replicates BGI SimilarityAnalysis directory (MRPP + ANOSIM).
# Software reference: R v3.5.1 (vegan package) per BGI report.
# Uses Bray-Curtis distance on rarefied data.
# ==============================================================================

library(vegan)

# --- Configuration ---
otu_file <- "../BGI_Result/OTU/OTU_table_for_biom.txt"
meta_file <- "../metadata.tsv"
output_dir <- "../BGI_Result/SimilarityAnalysis"
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

# --- Rarefaction (per BGI Section 9 methodology) ---
min_depth <- min(colSums(otu))
set.seed(42)
otu_rare <- as.data.frame(t(rrarefy(t(otu), sample = min_depth)))

# --- Bray-Curtis Distance ---
dist_bray <- vegdist(t(otu_rare), method = "bray")

# --- MRPP (Section 6: MRPP) ---
mrpp_dir <- file.path(output_dir, "1.MRPP")
dir.create(mrpp_dir, showWarnings = FALSE, recursive = TRUE)

mrpp_res <- mrpp(dist_bray, metadata$Group, permutations = 999)

mrpp_out <- data.frame(
    Group = paste(unique(metadata$Group), collapse = "-"),
    Distance = "bray",
    A = round(mrpp_res$A, 4),
    Observe_delta = round(mrpp_res$delta, 4),
    Expect_delta = round(mrpp_res$E.delta, 4),
    P_value = mrpp_res$Pvalue
)
write.table(mrpp_out, file.path(mrpp_dir, "MRPP_result.xls"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- ANOSIM (Section 6: ANOSIM) ---
anosim_dir <- file.path(output_dir, "2.Anosim")
dir.create(anosim_dir, showWarnings = FALSE, recursive = TRUE)

anosim_res <- anosim(dist_bray, metadata$Group, permutations = 999)

anosim_out <- data.frame(
    Group = paste(unique(metadata$Group), collapse = "-"),
    Distance = "bray",
    R_statistic = round(anosim_res$statistic, 4),
    P_value = anosim_res$signif
)
write.table(anosim_out, file.path(anosim_dir, "Anosim_result.xls"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ANOSIM boxplot
png(file.path(anosim_dir, "Anosim_boxplot.png"), width = 800, height = 600, res = 120)
plot(anosim_res, main = paste0("ANOSIM (R=", round(anosim_res$statistic, 3),
                                ", P=", anosim_res$signif, ")"))
dev.off()

print("ANOSIM + MRPP analysis complete.")
