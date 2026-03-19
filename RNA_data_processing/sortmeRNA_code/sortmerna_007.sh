#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 56:00:00
#SBATCH -J 007_3_18
#SBATCH --mail-type=ALL
#SBATCH --mail-user biancasarcani@gmail.com
#SBATCH --output=%x.%j.out

module load bioinfo-tools SortMeRNA/4.3.4

search_dir="/proj/naiss2024-6-73/Bianca/fastp/3/18"
working_dir="/proj/naiss2024-6-73/Bianca/sortmerna"
reference="/proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/ncbi_silva_coleoptera_cdhit95.fa"


find -L "$search_dir" -type f -name "*L007_R1_001.fastq.gz" | while read -r r1; do
    # Infer the R2 filename
    r2="${r1/L007_R1_001.fastq.gz/L007_R2_001.fastq.gz}"
    sample=$(basename "$r1" _L007_R1_001.fastq.gz)

    echo "===================== Running sortmeRNA on $sample ... ====================="

	sortmerna --threads 4 \
        --ref "$reference" \
	--workdir "$working_dir/007/3/18" \
	--reads "$r1" \
	--reads "$r2" \
	--fastx \
	--paired_in \
        --out2 --v \
	--no-best --num_alignments 1 \
	--aligned "$working_dir/aligned/${sample}_L007_rRNA" \
	--other "$working_dir/other/${sample}_L007_nonrRNA"

    # Emptying the intermediary files used for the analysis (we only need the output)
	rm -r "$working_dir/007/3/18/kvdb"
        rm -r "$working_dir/007/3/18/readb" 
      
    echo "===================== sortmeRNA on $sample is done ====================="

done

