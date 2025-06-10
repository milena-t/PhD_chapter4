#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00
#SBATCH -J fastqc
#SBATCH -o fastqc.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools 
module load python3
module load FastQC/0.11.9
module load MultiQC/1.12


# takes 1h for all trimmed fasta files of a species
#fastqc -t 10 /proj/naiss2023-6-65/Milena/snpseq00400_WE-3665_reads_trimmed/*paired*

#fastqc -t 10 /proj/naiss2023-6-65/Milena/rnaseq/prototyping/LS*/*fq.gz
#multiqc /proj/naiss2023-6-65/Milena/rnaseq/prototyping/*-mRNA
multiqc /proj/naiss2023-6-65/Milena/rnaseq/prototyping/*-RiboZeroPlus