#!/bin/bash

# Set directories
CDHIT_DIR="cdhit_results"
MMSEQS_DIR="mmseqs_results"
LNCRNA_DB="lncipedia.fasta"  

# Create output directory
mkdir -p $MMSEQS_DIR

# Process each stage
for stage in Pre Climax Adult; do
    echo "Processing $stage..."
    
    INPUT_FASTA="$CDHIT_DIR/${stage}_nonredundant.fasta"
    OUTPUT_FILE="$MMSEQS_DIR/${stage}_vs_lncipedia.tsv"
    
    # Run easy-search 
    mmseqs easy-search \
        $INPUT_FASTA \
        $LNCRNA_DB \
        $OUTPUT_FILE \
        $MMSEQS_DIR/tmp_$stage \
        --threads 8 \
        -e 0.001 \
        --format-output "query,target,evalue,pident,alnlen,qstart,qend,qlen,tstart,tend,tlen"
    
    echo "Completed $stage → $OUTPUT_FILE"
done

echo "All MMseqs2 easy-search analyses complete!"

# Filter potential lncRNAs

LNCRNA_FASTA="lncipedia.fasta"
mkdir -p mmseqs_results

for stage in Pre Climax Adult; do
    echo "Processing $stage..."
    
    # Run search
    mmseqs easy-search \
        cdhit_results/${stage}_nonredundant.fasta \
        $LNCRNA_FASTA \
        mmseqs_results/${stage}_hits.tsv \
        mmseqs_results/tmp_$stage \
        -e 0.001 \
        --threads 20
    
    # Extract IDs of transcripts with hits
    cut -f1 mmseqs_results/${stage}_hits.tsv | sort | uniq > mmseqs_results/${stage}_lncRNA_ids.txt
    
    # Extract FASTA sequences of potential lncRNAs
    seqtk subseq cdhit_results/${stage}_nonredundant.fasta mmseqs_results/${stage}_lncRNA_ids.txt > mmseqs_results/${stage}_potential_lncRNAs.fasta
    
    # Generate statistics
    TOTAL=$(grep -c ">" cdhit_results/${stage}_nonredundant.fasta)
    HITS=$(wc -l < mmseqs_results/${stage}_lncRNA_ids.txt)
    echo "$stage: $HITS/$TOTAL transcripts hit Lncipedia"
done