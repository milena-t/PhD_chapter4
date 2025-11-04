#!/bin/bash -l


YTOR_QUERY=/Users/miltr339/work/c_maculatus/annotation_comparison/superscaffolded_annotation/yTOR_annotation/yTor_Cmac_kaufmann2023.faa

CMAC_Lome_PROTEINS=/Users/miltr339/work/c_maculatus/annotation_comparison/superscaffolded_annotation/Cmac_Lome_diverse_filtered.faa
CMAC_Lu_PROTEINS=/Users/miltr339/work/c_maculatus/annotation_comparison/superscaffolded_annotation/Cmac_Lu_simple_filtered.faa
CMAC_SI_PROTEINS=/Users/miltr339/work/c_maculatus/annotation_comparison/superscaffolded_annotation/Cmac_SI_diverse_filtered.faa
CMAC_noRNA_PROTEINS=/Users/miltr339/work/c_maculatus/annotation_comparison/superscaffolded_annotation/yTOR_annotation/C_maculatus_filtered_proteinfasta_TE_filtered.fa
CMAC_nonSCAFFOLDED_noRNA_PROTEINS=/Users/miltr339/work/c_maculatus/annotation_comparison/non_superscaffolded/Cmac_isoform_filtered_no_bad_proteins.faa
CMAC_nonSCAFFOLDED_yesRNA_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/Kaufmann2023_updated_RNAseq_annotation/Cmac_not_superscaffolded/Cmac_isoform_filtered_no_bad_proteins.faa

CMAC_Bianca_transcriptome=/proj/naiss2024-6-73/Bianca/transcriptome/filtered_transcriptome/filtered_transcriptome.fasta

# makeblastdb -in $CMAC_Lome_PROTEINS -dbtype prot
# makeblastdb -in $CMAC_Lu_PROTEINS -dbtype prot
# makeblastdb -in $CMAC_SI_PROTEINS -dbtype prot

# blastp -query $YTOR_QUERY -db $CMAC_Lome_PROTEINS -out CMAC_Lome_ytor_blast.out -outfmt 6 -num_threads 4 -evalue 1e-10
# blastp -query $YTOR_QUERY -db $CMAC_Lu_PROTEINS -out CMAC_Lu_ytor_blast.out -outfmt 6 -num_threads 4 -evalue 1e-10
# blastp -query $YTOR_QUERY -db $CMAC_SI_PROTEINS -out CMAC_SI_ytor_blast.out -outfmt 6 -num_threads 4 -evalue 1e-10

# makeblastdb -in $CMAC_noRNA_PROTEINS -dbtype prot
# blastp -query $YTOR_QUERY -db $CMAC_noRNA_PROTEINS -out CMAC_noRNA_ytor_blast.out -outfmt 6 -num_threads 4 -evalue 1e-10

# makeblastdb -in $CMAC_nonSCAFFOLDED_noRNA_PROTEINS -dbtype prot
# blastp -query $YTOR_QUERY -db $CMAC_nonSCAFFOLDED_noRNA_PROTEINS -out CMAC_nonSCAFFOLDED_noRNA_ytor_blast.out -outfmt 6 -num_threads 4 -evalue 1e-10

# makeblastdb -in $CMAC_nonSCAFFOLDED_yesRNA_PROTEINS -dbtype prot
# blastp -query $YTOR_QUERY -db $CMAC_nonSCAFFOLDED_yesRNA_PROTEINS -out CMAC_nonSCAFFOLDED_noRNA_ytor_blast.out -outfmt 6 -num_threads 4 -evalue 1e-10

makeblastdb -in $CMAC_Bianca_transcriptome -dbtype prot
blastp -query $YTOR_QUERY -db $CMAC_Bianca_transcriptome -out yTOR_against_Bianca_transcriptome_blast.out -outfmt 6 -num_threads 4 -evalue 1e-10