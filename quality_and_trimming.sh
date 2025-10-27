#!/bin/bash

# Set directories
RAW_DATA="raw_data"
TRIMMED_DATA="trimmed_data"
QC_OUTPUT="qc_results"
ADAPTERS="path/to/adapters/TruSeq3-PE-2.fa"  # Update this path

# Create output directories
mkdir -p $TRIMMED_DATA
mkdir -p $QC_OUTPUT

# Process all samples
for R1_FILE in $RAW_DATA/*_R1*.fastq.gz; do
    # Get base filename
    BASENAME=$(basename $R1_FILE | sed 's/_R1.*//')
    
    # Corresponding R2 file
    R2_FILE=$(echo $R1_FILE | sed 's/_R1/_R2/')
    
    echo "Processing sample: $BASENAME"
    
    # Run Trimmomatic
    java -jar trimmomatic-0.39.jar PE \
        -threads 20 \
        -phred33 \
        $R1_FILE $R2_FILE \
        $TRIMMED_DATA/${BASENAME}_R1_paired.fastq.gz \
        $TRIMMED_DATA/${BASENAME}_R1_unpaired.fastq.gz \
        $TRIMMED_DATA/${BASENAME}_R2_paired.fastq.gz \
        $TRIMMED_DATA/${BASENAME}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:$ADAPTERS:2:30:10:2:keepBothReads \
        SLIDINGWINDOW:4:25 \
        MINLEN:40 \
        

    echo "Completed: $BASENAME"
done

# Generate MultiQC report
multiqc $TRIMMED_DATA -o $QC_OUTPUT