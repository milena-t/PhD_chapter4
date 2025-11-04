#!/bin/bash -l

#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J translate_transcriptome
#SBATCH -o translate_transcriptome.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools gffread/0.12.7 samtools/1.20 emboss/6.6.0

CMAC_Bianca_transcriptome=/proj/naiss2024-6-73/Bianca/transcriptome/filtered_transcriptome/filtered_transcriptome.fasta
TRANSLATED=/proj/naiss2024-6-73/Bianca/transcriptome/filtered_transcriptome/filtered_transcriptome_proteinseq.fasta

transeq -sequence $CMAC_Bianca_transcriptome -outseq $TRANSLATED