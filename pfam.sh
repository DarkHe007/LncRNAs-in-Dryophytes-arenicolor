#!/bin/bash

# Simplified version - single frame translation
FILTERED_DIR="filtered_transcriptomes_fdr"
TRANSLATION_DIR="translated_sequences"
HMMER_DIR="hmmer_pfam_results"
FINAL_DIR="final_noncoding"

mkdir -p $TRANSLATION_DIR $HMMER_DIR $FINAL_DIR

for stage in Pre Climax Adult; do
    echo "Processing $stage..."
    
    # Translate using longest ORF only (frame 1)
    transeq -sequence $FILTERED_DIR/${stage}_fdr_noncoding.fasta \
            -outseq $TRANSLATION_DIR/${stage}_translated.faa \
            -frame 1 \
            -clean
    
    # Run HMMER
    hmmscan --cpu 20 \
            --tblout $HMMER_DIR/${stage}_pfam.tbl \
            --noali \
            $HMMER_DIR/Pfam-A.hmm \
            $TRANSLATION_DIR/${stage}_translated.faa
    
    # Extract hits and filter
    awk '!/^#/ && $5 < 1e-5 {print $3}' $HMMER_DIR/${stage}_pfam.tbl | sort | uniq > $HMMER_DIR/${stage}_hits.txt
    
    seqkit grep -v -f $HMMER_DIR/${stage}_hits.txt \
        $FILTERED_DIR/${stage}_fdr_noncoding.fasta \
        -o $FINAL_DIR/${stage}_noncoding_no_pfam.fasta
    
    echo "  Completed $stage"
done