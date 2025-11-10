# Differential expression analysis on Larval RNAseq data

## literature notes

* [Shaw 2022](https://www.cell.com/trends/genetics/fulltext/S0168-9525(22)00084-1) Evolution of gene regulation on sex chromosomes. 
  * Mostly for young sex chromosomes, so potentially not super relevant for us, but check the methods of how they identified regulatory mutations.


## Notes from NBIS bioinformatics adisory program meeting Nov 2024

### relevant and actionable

* [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard) marks reads that are e.g. PCR duplicates in a SAM/BAM file. They are *not* removed
* manual curation of the TOR in the BRAKER3 annotation: re-run annotation
  * check the optional flag to run braker3 with exonerate
  * add larval RNA data
  * Manually curate TOR gene models in IGV based on RNA read alignment (bam files)
  * [Jbrowse](https://jbrowse.org/jb2/) is a GUI that keeps track of edits and exports them as gff

### other 

* **RNA from whole bodies may include contamination from the gut microbiome.** This can be filtered pretty easily with a bacterial reference library, but since our 80% mapping rate is already acceptable and these reads shouldn't map to the genome anyways I am not implementing it
* **Transcriptome.** If we really wanted to make the transcriptome work, this is what they recommend:
  * does translate (emboss) automate the selection of the translation table? It's probably fine because I also use it for the annotation stuff where it works but might be worth it to check.
  * The internal stop codons may be some weird artefact of the scaffolding or how CDhit picks its representative sequence.
  * "Manual" clustering that involves looking at the genomic coordinates and not the nucleotide sequences may be better.
  * AGAT has a tool that can make gene models based on genomic position of reads