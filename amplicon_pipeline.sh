#!/bin/bash

# ==============================================================================
# BGI Amplicon Workflow Pipeline Script (Optimized)
# Based on BGI_Amplicon_Report.md (华大基因)
# ==============================================================================
# This script automates the standard amplicon bioinformatics workflow:
# 1. Data Filtering (Quality trimming, Adapter removal)
# 2. Tag Connection (Paired-end merging)
# 3. OTU Clustering (97% identity, Chimera removal)
# 4. Taxonomy Annotation (RDP Classifier)
# ==============================================================================

# Fail fast on errors, unassigned variables, and pipe failures
set -euo pipefail
# Make globs evaluate to nothing if no files match
shopt -s nullglob

# --- Configuration & Tool Paths ---
THREADS=8
MIN_LEN_PCT=0.75  # 75% of original length
Q_WINDOW=25       # Sliding window size
Q_THRESHOLD=20    # Average Phred score in window
ADAPTER_OVERLAP=15
ADAPTER_ERROR=0.1 # Equivalent to ~3 mismatches in 30bp

# Software Paths
CUTADAPT="cutadapt"
FLASH="flash"
USEARCH="usearch"
RDP_CLASSIFIER="rdp_classifier"

# Reference Databases
GOLD_DB="path/to/gold_v20110519.fasta"
RDP_DB="path/to/rdp_database"

# Input/Output Directories
RAW_DATA_DIR="./raw_data"
CLEAN_DATA_DIR="./BGI_Result/Cleandata"
TAG_DIR="./BGI_Result/Tag"
OTU_DIR="./BGI_Result/OTU"
LOGS_DIR="./logs"

# --- Pre-flight Checks ---
echo "Running pre-flight checks..."
for tool in "$CUTADAPT" "$FLASH" "$USEARCH" "$RDP_CLASSIFIER"; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Error: Required tool '$tool' is not installed or not in PATH."
        exit 1
    fi
done

if [[ ! -d "$RAW_DATA_DIR" ]]; then
    echo "Error: Directory '$RAW_DATA_DIR' does not exist."
    exit 1
fi

mkdir -p "$CLEAN_DATA_DIR" "$TAG_DIR" "$OTU_DIR" "$LOGS_DIR"

# --- Main Pipeline Logic ---

# 1. Data Filtering (Section 4)
echo "--- Step 1: Data Filtering ---"
for f1 in "$RAW_DATA_DIR"/*_1.fq.gz; do
    f2="${f1/_1.fq.gz/_2.fq.gz}"
    sample=$(basename "$f1" _1.fq.gz)
    
    echo "Processing $sample..."
    
    # Cutadapt parameters derived from Figure 2 and Section 4:
    # - Average Q20 within 25bp sliding window (-q)
    # - No N bases allowed (--max-n 0)
    # - Low complexity (ploybase) removal is often handled upstream, but adapter
    #   and general quality are handled here.
    # Optimized cutadapt with multi-threading (-j)
    "$CUTADAPT" \
        -j "$THREADS" \
        -q "$Q_THRESHOLD" \
        -o "$CLEAN_DATA_DIR/${sample}_1.fq.gz" \
        -p "$CLEAN_DATA_DIR/${sample}_2.fq.gz" \
        --max-n 0 \
        --overlap "$ADAPTER_OVERLAP" \
        --error-rate "$ADAPTER_ERROR" \
        "$f1" "$f2" > "$LOGS_DIR/${sample}_cutadapt.log" 2>&1
done

# 2. Tag Connection (Section 5)
echo "--- Step 2: Tag Connection ---"
for f1 in "$CLEAN_DATA_DIR"/*_1.fq.gz; do
    f2="${f1/_1.fq.gz/_2.fq.gz}"
    sample=$(basename "$f1" _1.fq.gz)
    
    # FLASH parameters derived from Table 3:
    # - Min overlapping: 15 bp
    # - Max mismatch ratio: 0.10
    "$FLASH" \
        -m 15 \
        -M 300 \
        -x 0.1 \
        -t "$THREADS" \
        -d "$TAG_DIR" \
        -o "$sample" \
        "$f1" "$f2" > "$LOGS_DIR/${sample}_flash.log" 2>&1

    # Consolidate output as FinalTag.fastq as expected by the next steps
    mv "$TAG_DIR/${sample}.extendedFrags.fastq" "$TAG_DIR/${sample}_FinalTag.fastq"
    gzip -f "$TAG_DIR/${sample}_FinalTag.fastq"
done

# Cleanup unnecessary FLASH output
rm -f "$TAG_DIR"/*.notCombined*.fastq "$TAG_DIR"/*.hist

# 3. OTU Clustering (Section 6)
echo "--- Step 3: OTU Clustering ---"
echo "Merging tags and clustering OTUs..."
# Merge all samples for global clustering
cat "$TAG_DIR"/*_FinalTag.fastq.gz | gunzip > "$OTU_DIR/all_tags.fastq"

# Convert to FASTA and Dereplicate to reduce size
"$USEARCH" -fastq_filter "$OTU_DIR/all_tags.fastq" -fastaout "$OTU_DIR/all_tags.fasta" > "$LOGS_DIR/usearch_filter.log" 2>&1
"$USEARCH" -derep_fulllength "$OTU_DIR/all_tags.fasta" -output "$OTU_DIR/uniques.fasta" -sizeout > "$LOGS_DIR/usearch_derep.log" 2>&1

# Cluster at 97% identity (UPARSE algorithm)
"$USEARCH" -cluster_otus "$OTU_DIR/uniques.fasta" -otus "$OTU_DIR/otus.fasta" -relabel OTU_ > "$LOGS_DIR/usearch_cluster.log" 2>&1

# Chimera removal (UCHIME) against GOLD database
echo "Removing chimeras against GOLD database..."
"$USEARCH" -uchime_ref "$OTU_DIR/otus.fasta" -db "$GOLD_DB" -strand plus -nonchimeras "$OTU_DIR/OTU_final.fasta" > "$LOGS_DIR/usearch_uchime.log" 2>&1

# Create OTU Abundance Table mapping tags back to OTUs
echo "Generating OTU table..."
"$USEARCH" -usearch_global "$OTU_DIR/all_tags.fasta" -db "$OTU_DIR/OTU_final.fasta" -strand plus -id 0.97 -otutabout "$OTU_DIR/OTU_table.txt" > "$LOGS_DIR/usearch_global.log" 2>&1

# Clean intermediate FASTA/FASTQ files to save space
rm -f "$OTU_DIR/all_tags.fastq" "$OTU_DIR/all_tags.fasta" "$OTU_DIR/uniques.fasta" "$OTU_DIR/otus.fasta"

# 4. Taxonomy Annotation (Section 6)
echo "--- Step 4: Taxonomy Annotation ---"
# Using RDP Classifier with confidence threshold 0.6 per BGI report
"$RDP_CLASSIFIER" classify \
    -t "$RDP_DB" \
    -o "$OTU_DIR/OTU_taxonomy.xls" \
    -c 0.6 \
    "$OTU_DIR/OTU_final.fasta" > "$LOGS_DIR/rdp_classifier.log" 2>&1

echo "=============================================================================="
echo "Workflow completed successfully. Results are in BGI_Result/"
echo "=============================================================================="
