# ==============================================================================
# 00_run_all_groups.R
# BGI Amplicon Workflow - Master Group Iteration Wrapper
# ==============================================================================
# Runs all analysis scripts across each of the 11 BGI group comparisons.
# Each comparison generates a temporary metadata subset, then sources
# individual analysis scripts with the appropriate metadata path.
# ==============================================================================

cat("==================================================================\n")
cat("BGI Amplicon Workflow - Group-wise Analysis Runner\n")
cat("==================================================================\n\n")

# --- Configuration ---
meta_file <- "../metadata.tsv"
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

# --- Define the 11 BGI group comparisons ---
comparisons <- list(
    "A-B" = c("A", "B"),
    "A-C" = c("A", "C"),
    "B-C" = c("B", "C"),
    "A-B-C" = c("A", "B", "C"),
    "A-B-C-D-E-P" = c("A", "B", "C", "D", "E", "P"),
    "F-G-H-I-J-P" = c("F", "G", "H", "I", "J", "P"),
    "K-L-M-N-O-P-Q" = c("K", "L", "M", "N", "O", "P", "Q"),
    "A-B-C-D-E-F-G-H-I-J-P" = c("A","B","C","D","E","F","G","H","I","J","P"),
    "A-B-C-D-E-K-L-M-N-O-P-Q" = c("A","B","C","D","E","K","L","M","N","O","P","Q"),
    "P-Q" = c("P", "Q"),
    "A-B-C-D-E-F-G-H-I-J-K-L-M-N-O-P-Q" = unique(metadata$Group)
)

# --- Scripts to run per comparison ---
# These scripts read meta_file variable from the environment.
# Scripts self-skip gracefully when preconditions (e.g., numeric env vars) aren't met.
analysis_scripts <- c(
    "01_alpha_diversity.R",
    "02_beta_diversity.R",
    "03_taxa_composition.R",
    "04_differential_analysis.R",
    "05_function_prediction.R",
    "06_advanced_analysis.R",
    "07_rarefaction_curves.R",
    "08_pca_analysis.R",
    "09_similarity_tests.R",
    "10_venn_flower.R",
    "11_rank_abundance.R",
    "12_plsda.R",
    "13_network.R",
    "14_enterotypes.R",
    "15_multilevel_taxa.R",
    "16_function_expansion.R",
    "17_unifrac_beta.R",
    "18_nmds.R"
)

# --- Run each comparison ---
for (comp_name in names(comparisons)) {
    cat(sprintf("\n===================================================================\n"))
    cat(sprintf("=== COMPARISON: %s\n", comp_name))
    cat(sprintf("===================================================================\n"))

    groups <- comparisons[[comp_name]]
    meta_sub <- metadata[metadata$Group %in% groups, , drop = FALSE]
    n_samples <- nrow(meta_sub)
    n_groups <- length(unique(meta_sub$Group))

    if (n_samples < 3) {
        cat(sprintf("  SKIPPING: only %d samples.\n", n_samples))
        next
    }

    # Write temporary subset metadata
    tmp_meta <- file.path("..", paste0("metadata_", comp_name, ".tsv"))
    write.table(meta_sub, tmp_meta, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  Subset: %d samples, %d groups -> %s\n", n_samples, n_groups, tmp_meta))

    # Set the output suffix for this comparison
    comp_suffix <- comp_name

    # For each analysis script, modify the meta_file and output_dir in a local environment
    for (script in analysis_scripts) {
        if (!file.exists(script)) {
            cat(sprintf("    [SKIP] %s not found.\n", script))
            next
        }

        cat(sprintf("    [RUN]  %s ...", script))

        tryCatch({
            # Create a clean environment for each script
            env <- new.env(parent = globalenv())

            # Override meta_file to point to the subset
            env$meta_file <- tmp_meta

            # Provide comparison name for scripts that need it (e.g., 17_unifrac_beta.R
            # uses comp_suffix to auto-select comparison-specific phylogenetic trees)
            env$comp_suffix <- comp_suffix

            # --- Output directory interception ---
            # Parse script for all *_dir variable assignments and redirect outputs
            # to BGI_Reproduced/ while keeping input dirs pointing to BGI_Result/
            script_lines <- readLines(script)

            # Input-only dirs — must NOT be redirected
            input_only <- c("otu_dir", "tree_dir_beta", "tree_dir_genus", "base_dir")

            # Flat modules — comp_suffix goes in FILENAME, not directory
            # Match against the last component of the BGI_Result path
            flat_dirs <- c("PLSDA", "NMDS", "Picrust", "Function_Diff")

            # Find all lines like: var_name <- "../BGI_Result/..."
            dir_pattern <- "^if\\s*\\(!exists.*\\)\\s+(\\w+_dir)\\s*<-\\s*\"([^\"]+)\""
            bare_pattern <- "^(\\w+_dir)\\s*<-\\s*\"([^\"]+)\""

            for (line in script_lines) {
                # Try guarded pattern first, then bare assignment
                m <- regmatches(line, regexec(dir_pattern, line))[[1]]
                if (length(m) == 0) m <- regmatches(line, regexec(bare_pattern, line))[[1]]
                if (length(m) < 3) next

                var_name <- m[2]
                var_path <- m[3]

                # Skip input-only dirs
                if (var_name %in% input_only) next

                # Only redirect paths containing BGI_Result
                if (!grepl("BGI_Result", var_path)) next

                new_path <- gsub("BGI_Result", "BGI_Reproduced", var_path)
                # Only nest into comp_suffix subdirectory for non-flat modules
                path_tail <- basename(var_path)
                if (!path_tail %in% flat_dirs) {
                    new_path <- file.path(new_path, comp_suffix)
                }
                env[[var_name]] <- new_path
                dir.create(env[[var_name]], showWarnings = FALSE, recursive = TRUE)
            }

            source(script, local = env)
            cat(" OK\n")
        }, error = function(e) {
            cat(sprintf(" ERROR: %s\n", e$message))
        })
    }
}

cat("\n==================================================================\n")
cat("Group-wise analysis complete.\n")
cat(sprintf("Processed %d comparisons.\n", length(comparisons)))
cat("==================================================================\n")
