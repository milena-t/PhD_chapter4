#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -J fastP
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools fastp/0.23.4

# Directory with reads
search_dir="/proj/naiss2024-6-73/Bianca/Cmac_Y_larvae_RNA/"

# Output directory for results
output_dir="/proj/naiss2024-6-73/Bianca/fastp"
output_reports="/proj/naiss2024-6-73/Bianca/fastp/reports"


find -L "$search_dir" -type f -name "*R1_001.fastq.gz" | while read -r r1; do
    # Infer the R2 filename
    r2="${r1/R1_001.fastq.gz/R2_001.fastq.gz}"
    sample=$(basename "$r1" R1_001.fastq.gz)

    # Defining output files
    out_r1="$output_dir/${sample}R1_001.fastq.gz"
    out_r2="$output_dir/${sample}R2_001.fastq.gz"

    echo "Running fastp on $sample ..."
    fastp -w 16 -i "$r1" -I "$r2" -o "$out_r1" -O "$out_r2" -h "$output_reports/${sample}fastp.html"
done
