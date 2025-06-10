#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 30:00
#SBATCH -J checksums_RNAseq
#SBATCH -o checksums_RNAseq.sh

md5sum -c /proj/naiss2023-6-65/Cmac_Y_larvae_RNA/checksums.md5
