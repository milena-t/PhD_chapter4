#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -p pelle
#SBATCH -c 8
#SBATCH -t 4:00:00
#SBATCH -J picard_deduplicate
#SBATCH -o picard_deduplicate.out

# Load modules
module load picard/3.4.0-Java-17
module load SAMtools/1.22.1-GCC-13.3.0

# Directories
INPUT_DIR="/proj/naiss2023-6-65/Milena/chapter4/mapping_STAR/star_mapping_res"
OUTPUT_DIR="$INPUT_DIR/picard_marked_indexed"

mkdir -p "$OUTPUT_DIR"

echo "Array task $SLURM_ARRAY_TASK_ID processing $SAMPLE..."

# Loop over samples
for BAM in "$INPUT_DIR"/*_Aligned.sortedByCoord.out.bam; do
    SAMPLE=$(basename "$BAM" "_Aligned.sortedByCoord.out.bam")
    OUTPUT_BAM="$OUTPUT_DIR/${SAMPLE}_marked_duplicates.bam"
    METRICS="$OUTPUT_DIR/${SAMPLE}_markdup_metrics.txt"

    echo " "
    echo " ======================= marking duplicates sample: $SAMPLE ======================= "

    # Mark duplicates
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        -INPUT $BAM \
        -OUTPUT $OUTPUT_BAM \
        -METRICS_FILE $METRICS \
        -VALIDATION_STRINGENCY LENIENT

    # Index the marked BAM
    samtools index "$OUTPUT_BAM"

done

echo " "
echo "All samples mapped successfully!"

