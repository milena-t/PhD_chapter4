#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 30:00:00
#SBATCH -J 3_16_dedup
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools 
module load SeqKit/2.4.0

input_dir="/proj/naiss2024-6-73/Bianca/fastp/collapsed_lanes"
output_dir="/proj/naiss2024-6-73/Bianca/dedup"
prefix="WJ-3841-3-16"

for r1 in "$input_dir"/${prefix}*_R1_001.fastq.gz; do
    # Extracting the sample base
    base=$(basename "$r1" _R1_001.fastq.gz)

    # Infering R2
    r2="$input_dir/${base}_R2_001.fastq.gz"

    # Constructing output R1 dedup path
    dedup_r1="$output_dir/${base}_R1_001_dedup.fastq.gz"

    echo " ================ Processing sample: $base ================"
    
    # Deduplicating
    seqkit rmdup -j 12 -s -P "$r1" -o "$dedup_r1"
    
        # -s used because it compares by sequence: length and bases
        # -P is used to only consider the positive strands


    # Pairing R1 with original R2
    seqkit pair -j 12 -1 "$dedup_r1" -2 "$r2" -O "$output_dir" -u

done
