#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 16
#SBATCH -t 1:00:00
#SBATCH -J count_reads_Cmac
#SBATCH -o count_reads_Cmac.log
#SBATCH --mail-type=ALL

module load Subread/2.1.1-GCC-13.3.0

#Directories
ANNOTATION="/proj/naiss2023-6-65/Milena/chapter4/annotation/Cmac_Lome_yes_yTor.gtf"
INPUT_DIR="/proj/naiss2023-6-65/Milena/chapter4/mapping_STAR/star_mapping_res/picard_marked_indexed"
OUT_DIR="/proj/naiss2023-6-65/Milena/chapter4/mapping_STAR/gene_counts"

mkdir -p "$OUT_DIR"

# Mode 1: Standard counting (unique mappers only)
# No multimapping, at the gene level, used for differential expression analysis
echo "=== Mode 1: Standard (unique mappers only) ==="
featureCounts -T 16 \
  -a "$ANNOTATION" \
  -o "$OUT_DIR/gene_counts_standard.txt" \
  -p -B -C \
  -g "gene_id" \
  -t "exon" \
  -s 2 \
  "$INPUT_DIR"/*_marked_duplicates.bam
