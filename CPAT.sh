#!/bin/bash

# Simplified CPAT analysis
CDHIT_DIR="cdhit_results"
CPAT_DIR="cpat_results"
FILTERED_DIR="filtered_transcriptomes"

mkdir -p $CPAT_DIR $FILTERED_DIR

for stage in Pre Climax Adult; do
    echo "Running CPAT for $stage..."
    
    # Run CPAT 
    cpat \
        -g $CDHIT_DIR/${stage}_nonredundant.fasta \
        -x /path/to/xenopus_hexamer.tsv \  
        -m /path/to/xenopus_logit.RData \  
        -o $CPAT_DIR/${stage} \
        --min-orf=50
    
    # Extract non-coding sequences 
    awk '$5 < 0.001 {print $1}' $CPAT_DIR/${stage}.coding_prob > $CPAT_DIR/${stage}_noncoding_ids.txt
    
    seqkit grep -f $CPAT_DIR/${stage}_noncoding_ids.txt \
        $CDHIT_DIR/${stage}_nonredundant.fasta \
        -o $FILTERED_DIR/${stage}_cpat_noncoding.fasta
done

#!/usr/bin/env python3

import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

def main():
    stages = ['Pre', 'Climax', 'Adult']
    fdr_threshold = 0.0001
    method = 'fdr_bh'
    
    print("CPAT FDR Filtering using statsmodels")
    print("====================================\n")
    
    for stage in stages:
        print(f"Processing {stage}...")
        
        # Read CPAT results
        file_path = f"cpat_results/{stage}_cpat.coding_prob"
        cpat_data = pd.read_csv(file_path, sep='\t', header=None)
        cpat_data.columns = ['transcript_id', 'ORF_length', 'ORF_coverage', 
                            'frac_adenine', 'coding_prob', 'coding_label']
        
        # Calculate FDR
        p_values = 1 - cpat_data['coding_prob'].values
        _, fdr_values, _, _ = multipletests(p_values, alpha=fdr_threshold, method=method)
        cpat_data['fdr'] = fdr_values
        
        # Filter significant non-coding transcripts
        significant_nc = cpat_data[
            (cpat_data['fdr'] < fdr_threshold) & 
            (cpat_data['coding_prob'] < 0.5)
        ]['transcript_id'].tolist()
        
        # Filter significant coding transcripts
        significant_coding = cpat_data[
            (cpat_data['fdr'] < fdr_threshold) & 
            (cpat_data['coding_prob'] >= 0.5)
        ]['transcript_id'].tolist()
        
        # Save results
        cpat_data.to_csv(f'cpat_results/{stage}_cpat_with_fdr.csv', index=False)
        
        with open(f'cpat_results/{stage}_significant_noncoding.txt', 'w') as f:
            f.write('\n'.join(significant_nc))
        
        with open(f'cpat_results/{stage}_significant_coding.txt', 'w') as f:
            f.write('\n'.join(significant_coding))
        
        # Print summary
        print(f"  Total transcripts: {len(cpat_data)}")
        print(f"  Significant non-coding: {len(significant_nc)}")
        print(f"  Significant coding: {len(significant_coding)}")
        print()

    print("FDR filtering complete!")

if __name__ == "__main__":
    main()

# Set directories
CDHIT_DIR="cdhit_results"
CPAT_DIR="cpat_results"
FILTERED_DIR="filtered_transcriptomes_fdr"

mkdir -p $FILTERED_DIR

for stage in Pre Climax Adult; do
    echo "Extracting FDR-filtered sequences for $stage..."
    
    # Extract significant non-coding transcripts
    if [ -f "$CPAT_DIR/${stage}_significant_noncoding.txt" ]; then
        seqkit grep -f $CPAT_DIR/${stage}_significant_noncoding.txt \
            $CDHIT_DIR/${stage}_nonredundant.fasta \
            -o $FILTERED_DIR/${stage}_fdr_noncoding.fasta
    fi
    
    # Extract significant coding transcripts
    if [ -f "$CPAT_DIR/${stage}_significant_coding.txt" ]; then
        seqkit grep -f $CPAT_DIR/${stage}_significant_coding.txt \
            $CDHIT_DIR/${stage}_nonredundant.fasta \
            -o $FILTERED_DIR/${stage}_fdr_coding.fasta
    fi
    
    # Count results
    if [ -f "$FILTERED_DIR/${stage}_fdr_noncoding.fasta" ]; then
        NC_COUNT=$(grep -c ">" $FILTERED_DIR/${stage}_fdr_noncoding.fasta)
        CODING_COUNT=$(grep -c ">" $FILTERED_DIR/${stage}_fdr_coding.fasta 2>/dev/null || echo "0")
        echo "  $stage: $NC_COUNT non-coding, $CODING_COUNT coding"
    fi
done

echo "Sequence extraction complete!"