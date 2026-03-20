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
## https://github.com/sellwe/Master_thesis_sebastian/blob/main/scripts/dataset_1_dominance_kaufmann/star/star_alignment_dominance.sh
###############

# Load modules
module load STAR/2.7.11b-GCC-13.3.0

# Paths
WKDIR=/proj/naiss2023-6-65/Milena/chapter4/mapping_STAR

GENOME_DIR="${WKDIR}/star_index"
RNA_DIR="${WKDIR}/dedupplicated_reads"
OUT_DIR="${WKDIR}/star_mapping_res"

mkdir -p "$OUT_DIR"

# make gtf: 
# gffread Cmac_Lome_yes_yTor.gff -T -o Cmac_Lome_yes_yTor.gtf
GTF="/proj/naiss2023-6-65/Milena/chapter4/annotation/Cmac_Lome_yes_yTor.gtf"
GENOME_FA="/proj/naiss2023-6-65/Milena/chapter4/annotation/Cmac_superscaffolded.fna.masked"

mkdir -p "$GENOME_DIR" "$OUT_DIR"

#Create the STAR genome index
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir "$GENOME_DIR" \
     --genomeFastaFiles "$GENOME_FA" \
     --sjdbGTFfile "$GTF" \
     --sjdbOverhang 149 #max read length 150-1

# Loop over samples
for R1 in "$DATA_DIR"/*R1_001.fastq.gz; do
    R2="${R1/R1_001.fastq.gz/R2_001.fastq.gz}"
    SAMPLE_NAME=$(basename "$R1" "_R1_001.fastq.gz")
    echo " "
    echo " ======================= Quantifying sample: $SAMPLE_NAME ======================= "

    STAR --runThreadN 20 \
       --genomeDir "$GENOME_DIR" \
       --readFilesIn "$R1" "$R2" \
       --readFilesCommand zcat \
       --outSAMtype BAM SortedByCoordinate \
       --quantMode GeneCounts \
       --twopassMode Basic \
       --outFilterMultimapNmax 20 \
       --limitBAMsortRAM 20000000000 \
       --outFileNamePrefix "$OUT_DIR/${SAMPLE}_"
done

echo " "
echo "All samples mapped successfully!"