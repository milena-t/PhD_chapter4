#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 10:00:00
#SBATCH -J second_fastP_for_polyg
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools fastp/0.23.4

# Directory with reads
search_dir="/proj/naiss2024-6-73/Bianca/fastp/polyG_trimming"

# Output directory for results
output_dir="/proj/naiss2024-6-73/Bianca/fastp/second_polyG_trimming"


find -L "$search_dir" -type f -name "*_R1_001.fastq.gz" | while read -r r1; do
    # Infer the R2 filename
    r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    sample=$(basename "$r1" _R1_001.fastq.gz)

    # Defining output files
    out_r1="$output_dir/${sample}_R1_001.fastq.gz"
    out_r2="$output_dir/${sample}_R2_001.fastq.gz"

    echo "Running fastp on $sample ..."
    fastp -w 16 -i "$r1" -I "$r2" -o "$out_r1" -O "$out_r2" -g --poly_g_min_len 7
done 
