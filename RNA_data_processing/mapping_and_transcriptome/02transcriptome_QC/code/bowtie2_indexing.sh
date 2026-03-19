#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J de_novo_bwt
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools bowtie2/2.5.2

bowtie2-build trinity.Trinity.fasta trinity.Trinity.fasta
