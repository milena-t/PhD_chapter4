#!/bin/bash -l


YTOR_QUERY=/Users/miltr339/work/c_maculatus/annotation_comparison/superscaffolded_annotation/yTOR_annotation/yTor_Cmac_kaufmann2023.faa

CMAC_Lome_PROTEINS=/Users/miltr339/work/c_maculatus/annotation_comparison/superscaffolded_annotation/Cmac_Lome_diverse_filtered.faa
CMAC_Lu_PROTEINS=/Users/miltr339/work/c_maculatus/annotation_comparison/superscaffolded_annotation/Cmac_Lu_simple_filtered.faa
CMAC_SI_PROTEINS=/Users/miltr339/work/c_maculatus/annotation_comparison/superscaffolded_annotation/Cmac_SI_diverse_filtered.faa
CMAC_noRNA_PROTEINS=/Users/miltr339/work/c_maculatus/annotation_comparison/superscaffolded_annotation/yTOR_annotation/C_maculatus_filtered_proteinfasta_TE_filtered.fa
CMAC_nonSCAFFOLDED_noRNA_PROTEINS=/Users/miltr339/work/c_maculatus/annotation_comparison/non_superscaffolded/Cmac_isoform_filtered_no_bad_proteins.faa

# makeblastdb -in $CMAC_Lome_PROTEINS -dbtype prot
# makeblastdb -in $CMAC_Lu_PROTEINS -dbtype prot
# makeblastdb -in $CMAC_SI_PROTEINS -dbtype prot

# blastp -query $YTOR_QUERY -db $CMAC_Lome_PROTEINS -out CMAC_Lome_ytor_blast.out -outfmt 6 -num_threads 4 -evalue 1e-10
# blastp -query $YTOR_QUERY -db $CMAC_Lu_PROTEINS -out CMAC_Lu_ytor_blast.out -outfmt 6 -num_threads 4 -evalue 1e-10
# blastp -query $YTOR_QUERY -db $CMAC_SI_PROTEINS -out CMAC_SI_ytor_blast.out -outfmt 6 -num_threads 4 -evalue 1e-10

# makeblastdb -in $CMAC_noRNA_PROTEINS -dbtype prot
# blastp -query $YTOR_QUERY -db $CMAC_noRNA_PROTEINS -out CMAC_noRNA_ytor_blast.out -outfmt 6 -num_threads 4 -evalue 1e-10

makeblastdb -in $CMAC_nonSCAFFOLDED_noRNA_PROTEINS -dbtype prot
blastp -query $YTOR_QUERY -db $CMAC_nonSCAFFOLDED_noRNA_PROTEINS -out CMAC_nonSCAFFOLDED_noRNA_ytor_blast.out -outfmt 6 -num_threads 4 -evalue 1e-10