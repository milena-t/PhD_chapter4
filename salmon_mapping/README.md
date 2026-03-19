# Mapping the RNAseq data with salmon

I will use Salmon mapping, which maps the raw reads against a reference transcriptome (proteinfasta from genome annotation) with a decoy dataset. This is not the same as salmon alignment mode.

## Make decoy

There are two options according to the [Salmon wiki](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode), and I will aim for the more comprehensive method that uses the whole genome, and not the other one even though it uses less memory.


