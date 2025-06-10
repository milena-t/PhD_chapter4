#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 2:00:00
#SBATCH -J sortmerna
#SBATCH -o sortmerna_mRNA2.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools
module load SortMeRNA/4.3.4

sample_path="/proj/naiss2023-6-65/DataDelivery_2023-09-13_15-05-11_snpseq00473/files/WH-3754/230907_M00485_0747_000000000-DM8Y4"
mrna1f="Sample_WH-3754-LS11621-mRNA/WH-3754-LS11621-mRNA_S1_L001_R1_001.fastq.gz"
mrna1r="Sample_WH-3754-LS11621-mRNA/WH-3754-LS11621-mRNA_S1_L001_R2_001.fastq.gz"
mrna2f="Sample_WH-3754-LS31820-mRNA/WH-3754-LS31820-mRNA_S2_L001_R1_001.fastq.gz"
mrna2r="Sample_WH-3754-LS31820-mRNA/WH-3754-LS31820-mRNA_S2_L001_R2_001.fastq.gz"
ribozero1f="Sample_WH-3754-LS11621-RiboZeroPlus/WH-3754-LS11621-RiboZeroPlus_S3_L001_R1_001.fastq.gz"
ribozero1r="Sample_WH-3754-LS11621-RiboZeroPlus/WH-3754-LS11621-RiboZeroPlus_S3_L001_R2_001.fastq.gz"
ribozero2f="Sample_WH-3754-LS31820-RiboZeroPlus/WH-3754-LS31820-RiboZeroPlus_S4_L001_R1_001.fastq.gz"
ribozero2r="Sample_WH-3754-LS31820-RiboZeroPlus/WH-3754-LS31820-RiboZeroPlus_S4_L001_R2_001.fastq.gz"

# echo "$mrna1f $mrna1r  ==================== running  "

# sortmerna -out2 -paired_in -v -threads 20 \
# --workdir ./LS11621-mRNA \
# --ref /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/ncbi_silva_coleoptera_cdhit95.fa \
# -reads "$sample_path/$mrna1f" -reads "$sample_path/$mrna1r" \
# -fastx \
# --no-best \
# --num_alignments 1 \
# -aligned LS11621-mRNA/rRNA_reads \
# -other LS11621-mRNA/non_rRNA_reads

# echo "$mrna1f $mrna1r  ==================== done "
echo "$mrna2f $mrna2r  ==================== running  "

sortmerna -out2 -paired_in -v -threads 20 \
--workdir ./LS31820-mRNA \
--ref /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/ncbi_silva_coleoptera_cdhit95.fa \
-reads "$sample_path/$mrna2f" -reads "$sample_path/$mrna2r" \
-fastx \
--no-best \
--num_alignments 1 \
-aligned LS31820-mRNA/rRNA_reads \
-other LS31820-mRNA/non_rRNA_reads

echo "$mrna2f $mrna2r  ==================== done "
# echo "$ribozero1f $ribozero1r  ==================== running  "

# sortmerna -out2 -paired_in -v -threads 20 \
# --workdir ./LS11621-RiboZeroPlus \
# --ref /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/ncbi_silva_coleoptera_cdhit95.fa \
# -reads "$sample_path/$ribozero1f" -reads "$sample_path/$ribozero1r" \
# -fastx \
# --no-best \
# --num_alignments 1 \
# -aligned LS11621-RiboZeroPlus/rRNA_reads \
# -other LS11621-RiboZeroPlus/non_rRNA_reads

# echo "$ribozero1f $ribozero1r  ==================== done "
# echo "$ribozero2f $ribozero2r  ==================== running  "

# sortmerna -out2 -paired_in -v -threads 20 \
# --workdir ./LS31820-RiboZeroPlus \
# --ref /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/ncbi_silva_coleoptera_cdhit95.fa \
# -reads "$sample_path/$ribozero2f" -reads "$sample_path/$ribozero2r" \
# -fastx \
# --no-best \
# --num_alignments 1 \
# -aligned LS31820-RiboZeroPlus/rRNA_reads \
# -other LS31820-RiboZeroPlus/non_rRNA_reads
