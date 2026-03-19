#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 72:00:00
#SBATCH -J cd-hit
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools cd-hit/4.8.1

input="/proj/naiss2024-6-73/Bianca/transcriptome/trinity_with_ref/Trinity-GG.fasta"

cd /proj/naiss2024-6-73/Bianca/transcriptome/cd-hit

echo "=========== First tier clurstering ==========="
echo "=========== identity threshold 0.98 ==========="
cd-hit-est -i "$input" -o nr98.fasta -c 0.98 -n 11 -T 16 -M 100000 -l 130 -d 0 -p 1 -g 1

echo "=========== Second tier clurstering ==========="
echo "=========== identity threshold 0.96 ==========="
cd-hit-est -i nr98.fasta -o nr96.fasta -c 0.96 -n 10 -T 16 -M 100000 -l 130 -d 0 -p 1 -g 1

echo "=========== Third tier clurstering ==========="
echo "=========== identity threshold 0.95 ==========="
cd-hit-est -i nr96.fasta -o nr95.fasta -c 0.95 -n 9 -T 16 -M 100000 -l 130 -d 0 -p 1 -g 1


#  --- FLAGS ---

# c = sequence identity threshold, default 0.9. -c 1.0, means 100% identity, is the clustering threshold -c 0.9, means 90% identity, is the clustering threshold.
# n = wordsize
#n 10, 11 for thresholds 0.95 ~ 1.0; -n 8,9 for thresholds 0.90 ~ 0.95; -n 7 for thresholds 0.88 ~ 0.9; -n 6 for thresholds 0.85 ~ 0.88; -n 5 for thresholds 0.80 ~ 0.85; -n 4 for thresholds 0.75 ~ 0.8
# T = threads
# M = memory limit
# l = length of throwaway seq. If a transcript is less than 120 it gets thrown away
# d = only takes the clean IDs from each transcript
# p = print alignment overlap in .clstr file
# g = 0 (default) means a sequence is clustered to the first cluster that meet the threshold (fast cluster). If set to 1, the program will cluster it into the most similar cluster that meet the threshold (accurate but slow mode). We choose accurate and slow

