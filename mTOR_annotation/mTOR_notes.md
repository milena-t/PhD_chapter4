# mTOR annotation

## Old annotation: why am I doing this

There is an older version of the annotation done with BRAKER2 and TSEBRA separately, published in Kaufmann 2023 ([here](https://academic.oup.com/mbe/article/40/8/msad167/7227908)) which has like 35k genes, a lot for beetles. When I re-annotate a superscaffolded version of this assembly with BRAKER3 and the same RNAseq data, I find a much more reasonable 15k genes. I think this is because of the way that TSEBRA sets it's default values when combining evidence from different sources. the TSEBRA internal default weight value for RNAseq is only half of the value that BRAKER3 sets for TSEBRA. I assume this is why the old BRAKER2 annotation (which used TSEBRA-specific defaults) finds about the same number of genes as my annotation with no RNAseq, because the RNA evidence is valued down. When it is weighted higher like in BRAKER3, then more false positives are filtered out, reducing the number of genes overall, especially for shorter transcripts. 

However, this old annotation contains some manual curation for the copy number variation of TOR on the Y contigs. To be able to use the updated annotation for this project, I need to include the manual curation and also do a functional annotation.

## Functional annotation



## Manual curation 

This is the workflow that Doug used and thankfully documented really well in `/proj/naiss2023-6-65/douglas/nobackup/Callosobruchus_maculatus/mTor`

### mTOR consensus sequence

They already have candidate proteins from Cmac that they think are the duplicates on the Y, so Doug starts with making an alignment between them, the autosomal one, and several other species (`sequences/mTor_sequences.faa`):
<details>
  <summary>list</summary>
* [Homo sapiens]:\t NP_004949.1 serine/threonine kinas mTOR isoform 1 
* [Callosobruchus maculatus]:\t VEN43112.1 unnamed protein product 
* [Callosobruchus maculatus]:\t VEN51984.1 unnamed protein product 
* [Anoplophora glabripennis]:\t XP_018572076.1 serine/threonine-protein kinase Tor 
*  virgifera virgifera]:\t XP_028145210.1 serine/threonine-protein kinase Tor [Diabrotica
* [Leptinotarsa decemlineata]:\t XP_023015005.1 target of rapamycin 
* [Leptinotarsa decemlineata]:\t ALE20544.1 mTOR 
* [Brassicogethes aeneus]:\t CAH0563318.1 unnamed protein product 
* [Brassicogethes aeneus]:\t CAH0563403.1 unnamed protein product 
* [Rhynchophorus ferrugineus]:\t KAF7282696.1 hypothetical protein GWI33_002162 
* [Aethina tumida]:\t XP_019880827.1 PREDICTED: LOW QUALITY PROTEIN: target of rapamycin-like 
* [Tribolium castaneum]:\t XP_971819.1 PREDICTED: serine/threonine-protein kinase mTOR 
* [Sitophilus oryzae]:\t XP_030750054.1 serine/threonine-protein kinase Tor 
* [Asbolus verrucosus]:\t RZC37432.1 serine/threonine-protein kinase mTOR 
* [Tenebrio molitor]:\t CAH1377105.1 unnamed protein product 
* [Tenebrio molitor]:\t AKB11618.1 target of rapamycin 
* [Ignelater luminosus]:\t KAF2880605.1 hypothetical protein ILUMI_25569 
* [Lamprigera yunnana]:\t KAF5285925.1 hypothetical protein FQA39_LY04386 
* [Nicrophorus vespilloides]:\t XP_017768823.1 PREDICTED: target of rapamycin 
* [Photinus pyralis]:\t XP_031352545.1 serine/threonine-protein kinase Tor 
* [Onthophagus taurus]:\t XP_022907797.1 target of rapamycin 
* [Agrilus planipennis]:\t XP_025831250.1 serine/threonine-protein kinase Tor 
* [Coccinella septempunctata]:\t XP_044759281.1 serine/threonine-protein kinase Tor 
* [Harmonia axyridis]:\t XP_045478375.1 serine/threonine-protein kinase Tor 
* [Abscondita terminalis]:\t KAF5280820.1 hypothetical protein FQR65_LT14927 
* [Propylea japonica]:\t UIB01653.1 serine/threonine-protein kinase mTOR 

</details>

Since it has been identified in the other Cmac annotation already, I will not do it from scratch again in this one, I will just blast the Cmac yTor transcripts against the proteinfasta files of all the RNA annotations I did for the comparison in chapter 1.