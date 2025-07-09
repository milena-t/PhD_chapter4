#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH --mem=150G
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 48:00:00
#SBATCH -J trinity
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools trinity/2.14.0

cd /proj/naiss2024-6-73/Bianca/dedup/paired

Trinity --seqType fq --max_memory 148G \
--left WJ-3841-1-14-16-F_S430_R1_001.fastq.gz,WJ-3841-1-14-20-M_S426_R1_001.fastq.gz,WJ-3841-1-16-17-F_S405_R1_001.fastq.gz,WJ-3841-1-16-10-M_S401_R1_001.fastq.gz,WJ-3841-1-18-13-F_S420_R1_001.fastq.gz,WJ-3841-1-18-19-M_S437_R1_001.fastq.gz,WJ-3841-3-14-6-F_S363_R1_001.fastq.gz,WJ-3841-3-14-23-M_S413_R1_001.fastq.gz,WJ-3841-3-16-17-F_S419_R1_001.fastq.gz,WJ-3841-3-16-26-M_S417_R1_001.fastq.gz,WJ-3841-3-18-5-F_S372_R1_001.fastq.gz,WJ-3841-3-18-16-M_S397_R1_001.fastq.gz \
--right WJ-3841-1-14-16-F_S430_R2_001.fastq.gz,WJ-3841-1-14-20-M_S426_R2_001.fastq.gz,WJ-3841-1-16-17-F_S405_R2_001.fastq.gz,WJ-3841-1-16-10-M_S401_R2_001.fastq.gz,WJ-3841-1-18-13-F_S420_R2_001.fastq.gz,WJ-3841-1-18-19-M_S437_R2_001.fastq.gz,WJ-3841-3-14-6-F_S363_R2_001.fastq.gz,WJ-3841-3-14-23-M_S413_R2_001.fastq.gz,WJ-3841-3-16-17-F_S419_R2_001.fastq.gz,WJ-3841-3-16-26-M_S417_R2_001.fastq.gz,WJ-3841-3-18-5-F_S372_R2_001.fastq.gz,WJ-3841-3-18-16-M_S397_R2_001.fastq.gz \
--SS_lib_type FR \
--CPU 20 \
--output /crex/proj/snic2021-6-30/Bianca/trinity
