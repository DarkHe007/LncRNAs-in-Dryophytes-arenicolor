#!/bin/bash

# Set directories
FINAL_DIR="final_noncoding"
DRYOPHYTES_FASTA="dryophytes_arenicolor.fasta"  # UPDATE THIS PATH
MMSEQS_DIR="mmseqs_dryophytes_results"
FINAL_FILTERED_DIR="final_filtered_noncoding"

# Create output directories
mkdir -p $MMSEQS_DIR $FINAL_FILTERED_DIR

# Check if Dryophytes FASTA exists
if [ ! -f "$DRYOPHYTES_FASTA" ]; then
    echo "Error: Dryophytes arenicolor FASTA file not found at: $DRYOPHYTES_FASTA"
    echo "Please update the DRYOPHYTES_FASTA path in the script"
    exit 1
fi

# Process each stage
for stage in Pre Climax Adult; do
    echo "Processing $stage against Dryophytes arenicolor..."
    
    INPUT_FASTA="$FINAL_DIR/${stage}_high_conf_noncoding.fasta"
    
    if [ ! -f "$INPUT_FASTA" ]; then
        echo "Warning: $INPUT_FASTA not found, skipping..."
        continue
    fi
    
    # Run MMseqs2 easy-search
    mmseqs easy-search \
        $INPUT_FASTA \
        $DRYOPHYTES_FASTA \
        $MMSEQS_DIR/${stage}_hits.tsv \
        $MMSEQS_DIR/tmp_$stage \
        --threads 22 \
        -e 0.01 \
        --format-output "query,target,evalue,pident,alnlen"
    
    # Filter significant hits (E-value < 1e-5, identity > 80%)
    awk '$3 <= 1e-5 && $4 >= 80 {print $1}' $MMSEQS_DIR/${stage}_hits.tsv | \
        sort | uniq > $MMSEQS_DIR/${stage}_significant_hits.txt
    
    # Remove sequences with significant hits to Dryophytes
    seqkit grep -v \
        -f $MMSEQS_DIR/${stage}_significant_hits.txt \
        $INPUT_FASTA \
        -o $FINAL_FILTERED_DIR/${stage}_unique_noncoding.fasta
    
    # Count results
    TOTAL=$(seqkit stats $INPUT_FASTA -T | awk 'NR==2 {print $4}')
    HITS=$(wc -l < $MMSEQS_DIR/${stage}_significant_hits.txt)
    FINAL=$(seqkit stats $FINAL_FILTERED_DIR/${stage}_unique_noncoding.fasta -T | awk 'NR==2 {print $4}')
    
    echo "  $stage: $TOTAL → $FINAL sequences (removed $HITS Dryophytes hits)"
done

echo "Filtering complete!"