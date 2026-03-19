#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH -J blast
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools blast/2.15.0+

cd /proj/naiss2024-6-73/Bianca/transcriptome/blast

makeblastdb -in /proj/naiss2023-6-65/Milena/annotation_pipeline/repeatmasking/repeat_libraries/c.mac_HiC_plus_beetle_library.lib.fasta -dbtype nucl -out TE_database

blastn -query /proj/naiss2024-6-73/Bianca/transcriptome/cd-hit/nr95.fasta -db TE_database -outfmt 6 -out TE_results.txt


#The next code selects only the matches with ≥ 90% identity and E-value ≤ 1e-10 and alignment length ≥ 100 bp

awk '$3 >= 90 && $4 >= 100 && $11 <= 1e-10 {print $1}' TE_results.txt | sort | uniq > TE_results_filtered.txt 

