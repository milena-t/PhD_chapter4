#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -p pelle
#SBATCH -c 8
#SBATCH -t 4:00:00
#SBATCH -J picard_deduplicate
#SBATCH -o picard_deduplicate.out
#SBATCH --array=1-72

# Load modules
module load picard/3.4.0-Java-17
module load SAMtools/1.22.1-GCC-13.3.0

# Directories
INPUT_DIR="/proj/naiss2023-6-65/Milena/chapter4/mapping_STAR/star_mapping_res"
OUTPUT_DIR="$INPUT_DIR/picard_marked_indexed"

mkdir -p "$OUTPUT_DIR"

#Make an array of all the .bam files 
BAM_FILES=("$INPUT_DIR"/*_Aligned.sortedByCoord.out.bam)

#pick out the specific bam file for each array job 
BAM="${BAM_FILES[$SLURM_ARRAY_TASK_ID - 1]}"  #arrays are 0-based

SAMPLE=$(basename "$BAM" "_Aligned.sortedByCoord.out.bam")
OUTPUT_BAM="$OUTPUT_DIR/${SAMPLE}_marked_duplicates.bam"
METRICS="$OUTPUT_DIR/${SAMPLE}_markdup_metrics.txt"

echo "Array task $SLURM_ARRAY_TASK_ID processing $SAMPLE..."

# Mark duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    -INPUT $BAM \
    -OUTPUT $OUTPUT_BAM \
    -METRICS_FILE $METRICS \
    -VALIDATION_STRINGENCY LENIENT

# Index the marked BAM
samtools index "$OUTPUT_BAM"