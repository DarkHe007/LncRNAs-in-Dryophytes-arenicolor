#!/bin/bash

# Set paths
CDHIT_DIR="cdhit_results"
UNIPROT_FASTA="uniprot_database.fasta"  
MMSEQS_DIR="mmseqs_uniprot_results"

mkdir -p $MMSEQS_DIR

for stage in Pre Climax Adult; do
    echo "Aligning $stage to UniProt..."
    
    mmseqs easy-search \
        $CDHIT_DIR/${stage}_nonredundant.fasta \
        $UNIPROT_FASTA \
        $MMSEQS_DIR/${stage}_uniprot_hits.tsv \
        $MMSEQS_DIR/tmp_$stage \
        --threads 20 \
        -e 0.001
done


FILTERED_DIR="filtered_transcriptomes"

# Create output directory
mkdir -p $FILTERED_DIR

for stage in Pre Climax Adult; do
    echo "Removing UniProt hits from $stage..."
    
    INPUT_FASTA="$CDHIT_DIR/${stage}_nonredundant.fasta"
    HITS_FILE="$MMSEQS_DIR/${stage}_vs_uniprot_significant_ids.txt"
    OUTPUT_FASTA="$FILTERED_DIR/${stage}_noncoding.fasta"
    
    if [ ! -f "$INPUT_FASTA" ]; then
        echo "Warning: $INPUT_FASTA not found, skipping..."
        continue
    fi
    
    if [ ! -f "$HITS_FILE" ]; then
        echo "Warning: $HITS_FILE not found, skipping..."
        continue
    fi
    
    # Use seqkit to remove sequences that have UniProt hits
    seqkit grep -v -f $HITS_FILE $INPUT_FASTA > $OUTPUT_FASTA
    
    # Count sequences before and after
    BEFORE=$(seqkit stats $INPUT_FASTA -T | awk 'NR==2 {print $4}')
    AFTER=$(seqkit stats $OUTPUT_FASTA -T | awk 'NR==2 {print $4}')
    REMOVED=$((BEFORE - AFTER))
    
    echo "  $stage: $BEFORE → $AFTER sequences (removed $REMOVED UniProt hits)"
done

echo "Filtering complete!"