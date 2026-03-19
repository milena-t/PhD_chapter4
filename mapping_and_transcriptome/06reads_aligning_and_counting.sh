# This pipeline uses both Bowtie2 and Salmon. It aligns each RNA sample to the transciptome, then quantifies the reads with salmon.

#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 48:00:00
#SBATCH -J bwt-salmon
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools
module load bowtie2/2.5.2
module load samtools/1.20
module load Salmon/1.10.1

reads="/proj/naiss2024-6-73/Bianca/dedup"
transcriptome="/proj/naiss2024-6-73/Bianca/transcriptome/filtered_transcriptome.fasta"
bams="/proj/naiss2024-6-73/Bianca/transcriptome/alignment_bams"
out="/proj/naiss2024-6-73/Bianca/transcriptome/salmon"


for r1 in "$reads"/*_R1_001.fastq.gz; do
    base=$(basename "$r1" _R1_001.fastq.gz)
    r2="$reads/${base}_R2_001.fastq.gz"

    echo " ================ Processing sample: $base ================"
    
        bowtie2 -p 20 -q --no-unal -k 20 -x "$transcriptome" \
        -1 "$r1" -2 "$r2" | samtools view -@20 -Sb -o "$bams/${base}".bam 
        
        salmon quant -t "$transcriptome" -l A -a "$bams/${base}".bam -o "$out/${base}" -p 20
done
