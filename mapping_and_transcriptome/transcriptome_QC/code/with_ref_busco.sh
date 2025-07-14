#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 24:00:00
#SBATCH -J with_ref_busco
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools BUSCO/5.7.1

busco -i /proj/naiss2024-6-73/Bianca/transcriptome/with_ref/trinity_with_ref/Trinity-GG.fasta -m transcriptome -l $BUSCO_LINEAGE_SETS/endopterygota_odb10 -c 12 -o /proj/naiss2024-6-73/Bianca/transcriptome/with_ref/busco
