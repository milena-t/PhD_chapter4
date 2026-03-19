#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 50:00:00
#SBATCH -J salmon_decoy
#SBATCH -o salmon_decoy.log

module load  Salmon/1.10.3-GCC-13.3.0

###############
## Following this tutorial linked in the salmon documentation:
## https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
###############

TRANSCRIPTOME=/proj/naiss2023-6-65/Milena/chapter4/annotation/Cmac_Lome_yes_yTor.faa
ASSEMBLY=/proj/naiss2023-6-65/Milena/chapter4/annotation/Cmac_superscaffolded.fna.masked

## prepare metadata
grep "^>" $ASSEMBLY | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

## make gentrome (first targets then decoys)
cat $TRANSCRIPTOME $ASSEMBLY > gentrome.fa

## run salmon indexing
salmon index -t gentrome.fa -d decoys.txt -p 12 -i salmon_index --gencode