#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 24:00:00
#SBATCH -J de_novo_busco
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools BUSCO/5.7.1

busco -i /proj/naiss2024-6-73/Bianca/transcriptome/de_novo/trinity_de_novo/trinity.Trinity.fasta -m transcriptome -l $BUSCO_LINEAGE_SETS/endopterygota_odb10 -c 12 -o /proj/naiss2024-6-73/Bianca/transcriptome/de_novo/busco

