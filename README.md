# Differential expression analysis on Larval RNAseq data

## 1. Preprocessing and testing

This analysis is based on a lot of work to properly filter the data and compare mapping pipelines to handle the high rRNA contamination rate and multimapping issues.

### 1.1 rRNA filtering and other reads preprocessing

This was done by Bianca Sarcani, all data here: `PhD_chapter4/RNA_data_processing`. Basically:

* filter rRNA reads with sortmeRNA and custom rRNA library (see back to Axel's work about generating the reference library)
* collapse lanes
* deduplicate reads (removePCR duplicates and poly-G)

Bianca continued here and assembled a transcriptome but I will not use it here

### 1.2 Evaluation and comparison of mapping pipelines

This was done by Sebastian Ellwe, see repository here: https://github.com/sellwe/Master_thesis_sebastian. He compared:

* STAR
* Salmon
  * salmon-mapping (map RNA to proteinfasta with decoy dataset)
  * salmon-align (refine existing mapped bam file)

From his results I made the choice to use Salmon-mapping, see `PhD_chapter4/mapping_and_transcriptome`.

### 1.3 yTor annotation manual curation

This work is based on the selection lines evaluated by Kaufmann *et al.* ([link](https://doi.org/10.1093/molbev/msad167)), but I am using a different version of the *C. maculatus* annotation that does not have the gene structures of the Y TOR copy number variation. Therefore I lift these gene structures from the original annotation and transform the coordinates to match my superscaffolded version of the assembly/annotation, see `PhD_chapter4/mTOR_annotation`.

## 2. STAR mapping

See `PhD_chapter4/RNA_mapping`.

## 3. DE analysis

### PCA plots

<details>
<summary>Toggle down for plots</summary>

<p float="left">
  <img src="data/DE_figures/PCA_sex_line_all_counts.png" width="49%" />
  <img src="data/DE_figures/PCA_sex_day_all_counts.png" width="49%" />
</p>

When plotting all samples at once, the line is a clear separator, but the day not as clearly. Day 14 seems to be mostly to the left, but 16 and 18 are across the entire range.

</details>

