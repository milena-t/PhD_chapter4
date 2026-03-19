for sample in $(ls *_L006_R1_001.fastq.gz | sed 's/_L006_R1_001.fastq.gz//' | sort -u); do
    echo "Processing $sample..."

    # Merging R1 files
    cat ${sample}_L00*_R1_001.fastq.gz > ${sample}_R1_001.fastq.gz

    # Merging R2 files
    cat ${sample}_L00*_R2_001.fastq.gz > ${sample}_R2_001.fastq.gz
done
