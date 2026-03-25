#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -p pelle
#SBATCH -c 8
#SBATCH -t 4:00:00
#SBATCH -J picard_deduplicate
#SBATCH -o picard_deduplicate.out

# Load modules
module load picard/3.4.0-Java-17 SAMtools/1.22-GCC-13.3.0


# Directories
INPUT_DIR="/proj/naiss2023-6-65/Milena/chapter4/mapping_STAR/star_mapping_res"
OUTPUT_DIR="$INPUT_DIR/picard_marked_indexed"

mkdir -p "$OUTPUT_DIR"


for BAM in "$INPUT_DIR"/*_Aligned.sortedByCoord.out.bam; do

SAMPLE=$(basename "$BAM" "_Aligned.sortedByCoord.out.bam")
OUTPUT_BAM="$OUTPUT_DIR/${SAMPLE}_marked_duplicates.bam"
METRICS="$OUTPUT_DIR/${SAMPLE}_markdup_metrics.txt"

echo "Array task $SLURM_ARRAY_TASK_ID processing $SAMPLE..."

SAMPLE="$(basename "$BAM" "_Aligned.sortedByCoord.out.bam")"

# add readgroups that picard needs
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    INPUT="$BAM" \
    OUTPUT="${BAM%.bam}.rg.bam" \
    RGID=$SAMPLE \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=$SAMPLE \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT


# Mark duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    INPUT="${BAM%.bam}.rg.bam" \
    OUTPUT="$OUTPUT_BAM" \
    METRICS_FILE="$METRICS" \
    VALIDATION_STRINGENCY=LENIENT

# Index the marked BAM
samtools index "$OUTPUT_BAM"

done