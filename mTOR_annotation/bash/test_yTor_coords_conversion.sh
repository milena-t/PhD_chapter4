#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -n 5
#SBATCH -t 1:00:00
#SBATCH -J read_transcripts_from_gff
#SBATCH -o read_transcripts_from_gff.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# use gffread to extract the protein coding sequences
# -M :  cluster the input transcripts into loci, discarding "duplicated" transcripts (those with the same exact introns and fully contained or equal boundaries)
# -x :  write a FASTA file with spliced CDS for each GFF transcript

module load gffread/0.12.7-GCCcore-13.3.0 SAMtools/1.22-GCC-13.3.0 AGAT/1.6.1-GCCcore-13.3.0

ANNOT_GFF_RAW=$1
ASSEMBLY=/proj/naiss2023-6-65/Milena/chapter4/annotation/Cmac_superscaffolded.fna.masked
ASSEMBLY_K=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_maculatus/assembly_genomic.fna.masked

ANNOT_TRANSCRIPTS=${ANNOT_GFF_RAW}_isoform_filtered_transcripts.fna
ANNOT_PROTEINS=${ANNOT_GFF_RAW}_isoform_filtered_proteins.faa
TRANSEQ_PATH=/proj/naiss2023-6-65/Milena/software_install/emboss/EMBOSS-6.6.0/EMBOSS-6.6.0/bin/transeq

echo $(pwd)
echo $(ls -lh $ASSEMBLY)

## fix IDs and cds reading frame with AGAT
ANNOT_GFF_ID=${ANNOT_GFF_RAW}_AGAT_ID.gff
agat_sp_manage_IDs.pl --gff $ANNOT_GFF_RAW -o $ANNOT_GFF_ID
ANNOT_GFF=${ANNOT_GFF_ID}_CDS.gff
agat_sp_fix_cds_phases.pl --gff $ANNOT_GFF_ID --fasta $ASSEMBLY -o $ANNOT_GFF
rm $ANNOT_GFF_ID

# index assemblies (greatly decreases computing time, and won't work for the more fragmented callosobruchus assemblies otherwise)
# samtools faidx $ASSEMBLY
# extract transcript sequences

echo "gffread $ANNOT_GFF -M -x $ANNOT_TRANSCRIPTS -g $ASSEMBLY"
gffread $ANNOT_GFF -M -x $ANNOT_TRANSCRIPTS -g $ASSEMBLY

# gffread $ANNOT_GFF -g $ASSEMBLY -y $ANNOT_PROTEINS

# change fasta headers to include species names
# sed -i "s/>/>${SPECIES_NAME}_/g" $ANNOT_TRANSCRIPTS
# translate transcript sequences
$TRANSEQ_PATH -sequence $ANNOT_TRANSCRIPTS -outseq $ANNOT_PROTEINS
ls -lh $ANNOT_TRANSCRIPTS
echo "###########################################"
