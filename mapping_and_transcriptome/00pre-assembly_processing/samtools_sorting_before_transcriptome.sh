#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 2:00:00
#SBATCH -J samtools
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools samtools/1.20

cd /proj/naiss2024-6-73/Bianca/star

# Convert and sort SAM file to BAM
samtools sort -@ 6 -o output_thirdAligned.out.sorted.bam output_thirdAligned.out.sam 
# Index BAM file
samtools index output_thirdAligned.out.sorted.bam
