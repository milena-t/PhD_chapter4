#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 16:00:00
#SBATCH -J fastQC
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools FastQC/0.11.9

# Directory to search
search_dir="/proj/naiss2024-6-73/Bianca/Cmac_Y_larvae_RNA/"

# Output directory for FastQC results
output_dir="/proj/naiss2024-6-73/Bianca/fastqc"


find -L "$search_dir" -type f -name "*L007_R2_001.fastq.gz" | while read -r fastq_file; do
    echo "Running FastQC on $fastq_file ..."
    fastqc "$fastq_file" -t 6 -o "$output_dir"
done

