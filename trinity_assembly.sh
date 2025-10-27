#!/bin/bash

# Set directories
TRIMMED_DIR="trimmed_data"
TRINITY_OUTPUT="trinity_assemblies"
ALIGNMENT_DIR="alignment_results"
EXN50_DIR="exn50_stats"

# Create output directories
mkdir -p $TRINITY_OUTPUT $ALIGNMENT_DIR $EXN50_DIR

# Sample configuration
SAMPLES=("Pre" "Climax" "Adult")

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: $SAMPLE"
    
    # Run Trinity 
    mkdir -p $TRINITY_OUTPUT/$SAMPLE
    Trinity \
        --seqType fq \
        --left $TRIMMED_DIR/${SAMPLE}_R1_paired.fastq.gz \
        --right $TRIMMED_DIR/${SAMPLE}_R2_paired.fastq.gz \
        --CPU 8 \
        --max_memory 64G \
        --output $TRINITY_OUTPUT/$SAMPLE \
        --SS_lib_type FR
    
    # Run alignment 
    bowtie2 -p 10 -q --no-unal -k 20 \
        -x $TRINITY_OUTPUT/$SAMPLE/Trinity.fasta \
        -1 $TRIMMED_DIR/${SAMPLE}_R1_paired.fastq.gz \
        -2 $TRIMMED_DIR/${SAMPLE}_R2_paired.fastq.gz \
        2> $ALIGNMENT_DIR/${SAMPLE}_align_stats.txt | \
        samtools view -@10 -Sb -o $ALIGNMENT_DIR/${SAMPLE}_bowtie2.bam
    
    # Generate expression matrix 
    echo "Generating expression matrix for ExN50..."
    
    # Use the alignment-based method for ExN50
    $TRINITY_HOME/util/misc/contig_ExN50_statistic.pl \
        $ALIGNMENT_DIR/${SAMPLE}_bowtie2.bam \
        $TRINITY_OUTPUT/$SAMPLE/Trinity.fasta \
        transcript \
        | tee $EXN50_DIR/${SAMPLE}_ExN50.transcript.stats
    
    echo "Completed ExN50 analysis for $SAMPLE"
done