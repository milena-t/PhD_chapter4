#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 2:00:00
#SBATCH -J STAR
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out


module load bioinfo-tools star/2.7.11a

r1="/proj/naiss2024-6-73/Bianca/dedup/paired/reads_for_mapping/R1"
r2="/proj/naiss2024-6-73/Bianca/dedup/paired/reads_for_mapping/R2"

# 1. Generating genome indexes

STAR --runThreadN 12 --runMode genomeGenerate \
--genomeDir /proj/naiss2024-6-73/Bianca/star/third_annot \
--genomeFastaFiles /proj/naiss2023-6-65/douglas/primary1_draft_annotation/primary1.insectalow.masked.fasta \
--sjdbGTFtagExonParentTranscript Parent \   # this flag is because we have gff annotation, not gtf
--sjdbGTFfile /proj/naiss2023-6-65/Milena/coleoptera_sequences/c_maculatus/cmac_Y_s_annotation.gff \
--sjdbOverhang 150 

# Overhang is 150 because the maximum read length is 151. Overhang equation: ovrhng = (max_read_length - 1) 



# 2. Mapping

STAR --runThreadN 12 \
--genomeDir /proj/naiss2024-6-73/Bianca/star/third_annot \
--readFilesIn "$r1/WJ-3841-1-14-16-F_S430_R1_001.fastq.gz","$r1/WJ-3841-1-14-20-M_S426_R1_001.fastq.gz","$r1/WJ-3841-1-16-10-M_S401_R1_001.fastq.gz","$r1/WJ-3841-1-16-17-F_S405_R1_001.fastq.gz","$r1/WJ-3841-1-18-13-F_S420_R1_001.fastq.gz","$r1/WJ-3841-1-18-19-M_S437_R1_001.fastq.gz","$r1/WJ-3841-3-14-23-M_S413_R1_001.fastq.gz","$r1/WJ-3841-3-14-6-F_S363_R1_001.fastq.gz","$r1/WJ-3841-3-16-17-F_S419_R1_001.fastq.gz","$r1/WJ-3841-3-16-26-M_S417_R1_001.fastq.gz","$r1/WJ-3841-3-18-16-M_S397_R1_001.fastq.gz","$r1/WJ-3841-3-18-5-F_S372_R1_001.fastq.gz" "$r2/WJ-3841-1-14-16-F_S430_R2_001.fastq.gz","$r2/WJ-3841-1-14-20-M_S426_R2_001.fastq.gz","$r2/WJ-3841-1-16-10-M_S401_R2_001.fastq.gz","$r2/WJ-3841-1-16-17-F_S405_R2_001.fastq.gz","$r2/WJ-3841-1-18-13-F_S420_R2_001.fastq.gz","$r2/WJ-3841-1-18-19-M_S437_R2_001.fastq.gz","$r2/WJ-3841-3-14-23-M_S413_R2_001.fastq.gz","$r2/WJ-3841-3-14-6-F_S363_R2_001.fastq.gz","$r2/WJ-3841-3-16-17-F_S419_R2_001.fastq.gz","$r2/WJ-3841-3-16-26-M_S417_R2_001.fastq.gz","$r2/WJ-3841-3-18-16-M_S397_R2_001.fastq.gz","$r2/WJ-3841-3-18-5-F_S372_R2_001.fastq.gz" \
--readFilesCommand zcat \                                              # this is because the reads are compressed
--outFileNamePrefix /proj/naiss2024-6-73/Bianca/star/output_third
