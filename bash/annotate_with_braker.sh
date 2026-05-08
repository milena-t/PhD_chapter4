#!/bin/sh

ORTHODB_ARTHROPODA=/proj/naiss2023-6-65/Milena/annotation_pipeline/annotation_protein_data/OrthoDB_Arthropoda_v11.fa

LOME_DIR=/proj/naiss2023-6-65/Milena/chapter4/annotation/Lome_annot_new_RNAseq
cd $LOME_DIR

sbatch --job-name="Cmac_LOME_DE" --output="Cmac_LOME_diverse.out" -t 5-00:00:00 \
/proj/naiss2023-6-65/Milena/chapter4/PhD_chapter4/bash/braker3_run.sh Cmac_Lome_diverse ${LOME_DIR}/assembly_genomic.fna.masked \
$ORTHODB_ARTHROPODA \
ERR12383247_trimmed,ERR12383253_trimmed,ERR12383273_trimmed,ERR12383281_trimmed,ERR12383299_trimmed,ERR12383251_trimmed,ERR12383316_trimmed,ERR12383270_trimmed,ERR12383276_trimmed,ERR12383296_trimmed,ERR12383303_trimmed,Sample_WJ-3841-1-14-9-M_S376,Sample_WJ-3841-1-16-20-M_S431,Sample_WJ-3841-3-14-17-M_S411,Sample_WJ-3841-3-16-11-M_S365 \
/proj/naiss2023-6-65/Milena/chapter4/annotation/Lome_annot_new_RNAseq/Cmac_Lome_RNA/

# for SAMPLE in "Sample_WJ-3841-1-14-9-M" "Sample_WJ-3841-1-16-20-M" "Sample_WJ-3841-3-14-17-M" "Sample_WJ-3841-3-16-11-M"; do READS="${SAMPLE#Sample_}"; ll /proj/naiss2023-6-65/Milena/chapter4/Bianca/dedupplicated_reads/${READS}*  ;  done
# for SAMPLE in $(ls) ; do echo $SAMPLE ; S1="${SAMPLE/_001/}" ; S2="${S1/R/}"; mv $SAMPLE $S2 ;done

# for SAMPLE in "ERR12383247" "ERR12383253" "ERR12383273" "ERR12383281" "ERR12383299" "ERR12383251" "ERR12383316" "ERR12383270" "ERR12383276" "ERR12383296" "ERR12383303"; do ln -s /proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/fastp_multiqc/${SAMPLE}_trimmed* . ;  done

bash /proj/naiss2023-6-65/Milena/chapter4/PhD_chapter4/bash/braker3_run.sh Cmac_Lome_diverse /proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lome_superscaffolded_comparison/Cmac_Lome_diverse/assembly_genomic.fna.masked $ORTHODB_ARTHROPODA ERR12383247_trimmed,ERR12383253_trimmed,ERR12383273_trimmed,ERR12383281_trimmed,ERR12383299_trimmed,ERR12383251_trimmed,ERR12383316_trimmed,ERR12383270_trimmed,ERR12383276_trimmed,ERR12383296_trimmed,ERR12383303_trimmed,Sample_WJ-3841-1-14-9-M_S376,Sample_WJ-3841-1-16-20-M_S431,Sample_WJ-3841-3-14-17-M_S411,Sample_WJ-3841-3-16-11-M_S365 /proj/naiss2023-6-65/Milena/chapter4/annotation/Lome_annot_new_RNAseq/Cmac_Lome_RNA/