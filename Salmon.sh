#!/bin/bash

# Trinity's full expression analysis pipeline
TRIMMED_DIR="trimmed_data"
TRINITY_DIR="trinity_assemblies"
EXPRESSION_DIR="trinity_expression_matrix"

mkdir -p $EXPRESSION_DIR

# Step 1: Run abundance estimation for each assembly
for stage in Pre Climax Adult; do
    echo "Processing $stage..."
    
    $TRINITY_HOME/util/align_and_estimate_abundance.pl \
        --transcripts $TRINITY_DIR/$stage/Trinity.fasta \
        --seqType fq \
        --left $TRIMMED_DIR/${stage}_R1_paired.fastq.gz \
        --right $TRIMMED_DIR/${stage}_R2_paired.fastq.gz \
        --est_method salmon \        # Trinity recommends salmon
        --trinity_mode \
        --output_dir $EXPRESSION_DIR/${stage}_salmon \
        --thread_count 20
done

# Step 2: Create expression matrix using Trinity's script
$TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
    --est_method salmon \
    --gene_trans_map $TRINITY_DIR/Pre/Trinity.fasta.gene_trans_map \
    --name_sample_by_basedir \
    $EXPRESSION_DIR/Pre_salmon/quant.sf \
    $EXPRESSION_DIR/Climax_salmon/quant.sf \
    $EXPRESSION_DIR/Adult_salmon/quant.sf \
    > $EXPRESSION_DIR/salmon.counts.matrix