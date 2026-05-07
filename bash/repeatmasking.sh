#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH -t 5-10:00:00
#SBATCH --mail-type=ALL
#SBATCH -J repeatmasking
#SBATCH -o repeatmasking.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# Philipp thesis pdf page 154 repeatmasker settings: -xsmall -s -u -engine ncbi -gff
# /proj/naiss2023-6-65/Milena/annotation_pipeline

module load RepeatModeler/2.0.7-foss-2024a
module load RepeatMasker/4.2.1-foss-2024a
module load SAMtools/1.22.1-GCC-13.3.0


"""
NOTES ON USEAGE:

Usage: $0 path_to_assembly repeat_library_dir species_identifier_string

if the repeat_library_dir exists, then the script will assume that there is a working repeat library in it and skip the library creation step and go straight to masking.
repeat_library_dir should contain identifying information, like a species name.

Good blog post with more detailed explanations: https://darencard.net/blog/2022-07-09-genome-repeat-annotation/

====================================
"""


if [ $# -lt 2 ]; then
    echo "Usage: $0 path_to_assembly repeat_library_dir species_identifier_string"
    echo "you have $#"
    exit 1
fi



ASSEMBLY=$1
LIBRARIES_DIR=$2
SPECIES_IDENT=$3

# sbatch /proj/naiss2023-6-65/Milena/chapter4/PhD_chapter4/bash/repeatmasking.sh /proj/naiss2023-6-65/Milena/Y_L_assembly_annotation/Cmac_y_L_assembly.fna /proj/naiss2023-6-65/Milena/Y_L_assembly_annotation/Cmac_y_L_repeat_families Cmac_y_L 
# sbatch /proj/naiss2023-6-65/Milena/chapter4/PhD_chapter4/bash/repeatmasking.sh Cmac_y_L_assembly.fna Cmac_y_L_repeat_families Cmac_y_L 

## make custom repeat library based on the species assembly
# RepeatModeler uses a NCBI BLASTDB as input to the repeat modeling pipeline, BuildDatabase is a wrapper to make this database for all future steps

if [ -d $LIBRARIES_DIR ]; then
  echo "Directory '$LIBRARIES_DIR' already exists, assume it has a repeat library in it: ${SPECIES_IDENT}_repeats-families.fa"
else
  mkdir -p "$LIBRARIES_DIR"
  BuildDatabase -name "${LIBRARIES_DIR}/${SPECIES_IDENT}_repeats" $ASSEMBLY  # this takes like 15 mins for Cmac
  echo "=====================> build database done"
  RepeatModeler -database "${LIBRARIES_DIR}/${SPECIES_IDENT}_repeats" -threads 20 -LTRStruct  # this takes over a day for Cmac
  echo "=====================> repeatmodeller done"
fi


### combine repeatmodeler library and custom curated library from the Zooeco people
# cat $CMAC_CURATED_REPEATS "$LIBRARIES_DIR/A_obt_repeats-families.fa" > $LIBRARIES_COMBINED

echo "REPEAT LIBRARY: ${LIBRARIES_DIR}/${SPECIES_IDENT}_repeats-families.fa"
echo "ASSEMBLY: $ASSEMBLY"

echo "  --> truncate assembly fasta headers to max 50 characters because repeatmasker wants that"
head -n 1 $ASSEMBLY

# cut to 49 characters
# echo "sed -i 's/^(>.{48}).*/\1/' $ASSEMBLY"
# sed -E -i 's/^(>.{48}).*/\1/' $ASSEMBLY

# cut after first underscore
echo "'s/^>([^_]*)_.*/>\1/' $ASSEMBLY"
sed -r -i 's/^>([^_]*)_.*/>\1/' $ASSEMBLY
head -n 1 $ASSEMBLY

echo "  --> index the assembly to match the new fasta headers"
samtools faidx $ASSEMBLY # indexed file will be in the same directory as $ASSEMBLY, not working directory
head -n 1 "${ASSEMBLY}.fai"


echo " "
echo "  --> start repeatmasking"
RepeatMasker -lib "${LIBRARIES_DIR}/${SPECIES_IDENT}_repeats-families.fa" -xsmall -s -u -engine ncbi -gff -pa 20 $ASSEMBLY



# CMAC_ENA="/proj/naiss2023-6-65/Milena/coleoptera_sequences/c_maculatus/Cmac_from_ENA_GCA_951848785.1.fasta" # 938 contigs
# RepeatMasker -lib $CMAC_CURATED_REPEATS -xsmall -s -u -engine ncbi -gff -pa 20 $CMAC_ENA

# CMAC_SUPERSCAFFOLDED="/proj/naiss2023-6-65/Milena/coleoptera_sequences/c_maculatus/Cmac_from_ENA_GCA_951848785.1.fasta" # 753 contigs
# RepeatMasker -lib $CMAC_CURATED_REPEATS -xsmall -s -u -engine ncbi -gff -pa 20 $CMAC_SUPERSCAFFOLDED


