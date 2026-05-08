#!/bin/sh

ORTHODB_ARTHROPODA=/proj/naiss2023-6-65/Milena/annotation_pipeline/annotation_protein_data/OrthoDB_Arthropoda_v11.fa

# T. freemani annotation
# TFRE_DIR=/proj/naiss2023-6-65/Milena/chapter3/species_assemblies
# cd $TFRE_DIR
# sbatch --job-name="T_freemani_annotation" --output="T_freemani_annotation.out" -t 5-00:00:00 \
# /proj/naiss2023-6-65/Milena/annotation_pipeline/braker3_singularity_with_RNAseq_in_SNIC_TMP.sh \
# T_freemani ${TFRE_DIR}/Tfre_GCA_022388455.1.fasta.masked \
# $ORTHODB_ARTHROPODA SRR15965976 SRR15965980 SRR15965983 SRR15965985 SRR14070854 SRR14070855 SRR14070871 SRR14070870 

# SRR15965976 SRR15965980 SRR15965983 SRR15965985 \ # larvae RNA different stages and population densities
# SRR14070854 SRR14070855 \ # male adult
# SRR14070871 SRR14070870 \ # female adult

# sbatch --job-name="C_magnifica_annotation" --output="C_magnifica_annotation_no_RNAseq.out" -t 5-00:00:00 \
# /proj/naiss2023-6-65/Milena/annotation_pipeline/braker3_singularity_with_RNAseq_in_SNIC_TMP.sh \
# C_magnifica ${TFRE_DIR}/Cmag_GCA_965644565.1.fasta.masked \
# $ORTHODB_ARTHROPODA

LOME_DIR=/proj/naiss2023-6-65/Milena/chapter4/annotation/Lome_annot_new_RNAseq
cd $LOME_DIR
sbatch --job-name="Cmac_LOME_DE" --output="Cmac_LOME_diverse.out" -t 5-00:00:00 \
/proj/naiss2023-6-65/Milena/chapter4/PhD_chapter4/bash/braker3_run.sh Cmac_Lome_diverse ${LOME_DIR}/assembly_genomic.fna.masked \
$ORTHODB_ARTHROPODA ERR12383247,ERR12383253,ERR12383273,ERR12383281,ERR12383299,ERR12383251,ERR12383316,ERR12383270,ERR12383276,ERR12383296,ERR12383303,/proj/naiss2023-6-65/Milena/chapter4/annotation/Lome_annot_new_RNAseq/larval_RNA/Sample_WJ-3841-1-14-9-M,/proj/naiss2023-6-65/Milena/chapter4/annotation/Lome_annot_new_RNAseq/larval_RNA/Sample_WJ-3841-1-16-20-M,/proj/naiss2023-6-65/Milena/chapter4/annotation/Lome_annot_new_RNAseq/larval_RNA/Sample_WJ-3841-3-14-17-M,/proj/naiss2023-6-65/Milena/chapter4/annotation/Lome_annot_new_RNAseq/larval_RNA/Sample_WJ-3841-3-16-11-M,
# for SAMPLE in $(ls) ;do cd  $SAMPLE ; READS="${SAMPLE#Sample_}";  pwd ; ln -s /proj/naiss2023-6-65/Milena/chapter4/Bianca/dedupplicated_reads/${READS}* . ;  cd .. ; done