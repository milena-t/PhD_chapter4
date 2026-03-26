#!/bin/bash

# Directories
INPUT_DIR="/proj/naiss2023-6-65/Milena/chapter4/mapping_STAR/star_mapping_res"
OUTPUT_DIR="$INPUT_DIR/picard_marked_indexed"

mkdir -p "$OUTPUT_DIR"


for BAM in "$INPUT_DIR"/*_Aligned.sortedByCoord.out.bam; do
    
    SAMPLE=$(basename "$BAM" "_Aligned.sortedByCoord.out.bam")
    sbatch -o "run_${SAMPLE}" -J "${SAMPLE}" /proj/naiss2023-6-65/Milena/chapter4/PhD_chapter4/RNA_mapping/star_picard_deduplicate.sh $BAM 

done