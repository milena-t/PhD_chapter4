#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J gff_formatting
#SBATCH -o gff_formatting.log
#SBATCH --mail-type=ALL

# interactive -A uppmax2026-1-8 -t 5:00:00

# module load gffread/0.12.7-GCCcore-13.3.0 AGAT/1.6.1-GCCcore-13.3.0
mamba activate mamba_agat
cd /proj/naiss2023-6-65/Milena/chapter4/annotation

GTF_ANNOT=/proj/naiss2023-6-65/Milena/chapter4/annotation/Cmac_Lome_no_yTor.gtf
GFF_ANNOT=/proj/naiss2023-6-65/Milena/chapter4/annotation/Cmac_Lome_no_yTor.gff

agat_convert_sp_gxf2gxf.pl -gtf $GTF_ANNOT -o $GFF_ANNOT


## -T when converting gff to gtf for salmon/star mapping
# gffread $GFF_ANNOT -o $GTF_ANNOT -T