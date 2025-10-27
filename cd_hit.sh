#!/bin/bash

# Set directories
TRINITY_OUTPUT="trinity_assemblies"
CONCAT_DIR="concatenated_assemblies"
CDHIT_DIR="cdhit_results"
BUSCO_DIR="busco_results"

# Create output directories
mkdir -p $CONCAT_DIR $CDHIT_DIR

# Simple concatenation by stage
echo "Concatenating transcriptomes..."

# For Pre stage
cat $TRINITY_OUTPUT/Pre/Trinity.fasta > $CONCAT_DIR/Pre_referencia.fasta

# For Climax stage  
cat $TRINITY_OUTPUT/Climax/Trinity.fasta > $CONCAT_DIR/Climax_referencia.fasta

# For Adult stage
cat $TRINITY_OUTPUT/Adult/Trinity.fasta > $CONCAT_DIR/Adult_referencia.fasta

# For Global Reference

cat $CONCAT_DIR/*.fasta > $CONCAT_DIR/Referencia.fasta

echo "Concatenation complete!"

# Run CD-HIT on each concatenated file
echo "Running CD-HIT to remove redundancies..."

cd-hit-est -i $CONCAT_DIR/Pre_referencia.fasta -o $CDHIT_DIR/Pre_nonredundant.fasta -c 0.95 -n 10 -M 32000 -T 20
cd-hit-est -i $CONCAT_DIR/Climax_referencia.fasta -o $CDHIT_DIR/Climax_nonredundant.fasta -c 0.95 -n 10 -M 32000 -T 20 
cd-hit-est -i $CONCAT_DIR/Adult_referencia.fasta -o $CDHIT_DIR/Adult_nonredundant.fasta -c 0.95 -n 10 -M 32000 -T 20
cd-hit-est -i $CONCAT_DIR/Referencia.fasta -o $CDHIT_DIR/Referencia_nonredundant.fasta -c 0.95 -n 10 -M 32000 -T 20

echo "CD-HIT complete!"

# Run BUSCO analysis on non-redundant transcriptomes
echo "Running BUSCO analysis..."
for stage in Pre Climax Adult; do
    echo "Processing $stage..."
    
    INPUT_FASTA="$CDHIT_DIR/${stage}_nonredundant.fasta"
    
    # BUSCO against eukaryota_odb10
    busco -i $INPUT_FASTA \
          -l eukaryota_odb10 \
          -o ${stage}_eukaryota \
          -m transcriptome \
          -c 22 \
          --out_path $BUSCO_DIR
    
    # BUSCO against vertebrata_odb10  
    busco -i $INPUT_FASTA \
          -l vertebrata_odb10 \
          -o ${stage}_vertebrata \
          -m transcriptome \
          -c 22 \
          --out_path $BUSCO_DIR
    
    # BUSCO against tetrapoda_odb10
    busco -i $INPUT_FASTA \
          -l tetrapoda_odb10 \
          -o ${stage}_tetrapoda \
          -m transcriptome \
          -c 22 \
          --out_path $BUSCO_DIR
done

echo "All analyses complete!"