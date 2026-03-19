#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -p pelle
#SBATCH -c 20
#SBATCH -t 48:00:00
#SBATCH -J salmon_quant
#SBATCH -o salmon_quant.out
#SBATCH --mail-type=ALL
#SBATCH --mem=32G


###############
## Taken directly from sebastian's github!
## https://github.com/sellwe/Master_thesis_sebastian/blob/main/scripts/dataset_1_dominance_kaufmann/salmon_mapping_based/run_salmon_map_consistent.sh
###############

# Load modules
module load Salmon/1.10.3-GCC-13.3.0 

# Paths
WKDIR=/proj/naiss2023-6-65/Milena/chapter4/mapping_salmon

INDEX_DIR="${WKDIR}/salmon_index"
DATA_DIR="${WKDIR}/dedupplicated_reads"
OUT_DIR="${WKDIR}/salmon_quant"

mkdir -p "$OUT_DIR"

# Loop over samples
for R1 in "$DATA_DIR"/*R1_001.fastq.gz; do
    R2="${R1/R1_001.fastq.gz/R2_001.fastq.gz}"
    SAMPLE_NAME=$(basename "$R1" "_R1_001.fastq.gz")
    echo "Quantifying sample: $SAMPLE_NAME"

    salmon quant \
        -i "$INDEX_DIR" \
        -l A \
        -1 "$R1" \
        -2 "$R2" \
        -p 20 \
        --gcBias \
        --seqBias \
        --validateMappings \
        -o "$OUT_DIR/$SAMPLE_NAME"
done

echo "All quantifications completed successfully!"