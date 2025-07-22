# mTOR annotation

## Old annotation: why am I doing this

There is an older version of the annotation done with BRAKER2 and TSEBRA separately, published in Kaufmann 2023 ([here](https://academic.oup.com/mbe/article/40/8/msad167/7227908)) which has like 35k genes, a lot for beetles. When I re-annotate a superscaffolded version of this assembly with BRAKER3 and the same RNAseq data, I find a much more reasonable 15k genes. I think this is because of the way that TSEBRA sets it's default values when combining evidence from different sources. the TSEBRA internal default weight value for RNAseq is only half of the value that BRAKER3 sets for TSEBRA. I assume this is why the old BRAKER2 annotation (which used TSEBRA-specific defaults) finds about the same number of genes as my annotation with no RNAseq, because the RNA evidence is valued down. When it is weighted higher like in BRAKER3, then more false positives are filtered out, reducing the number of genes overall, especially for shorter transcripts. 

However, this old annotation contains some manual curation for the copy number variation of TOR on the Y contigs. To be able to use the updated annotation for this project, I need to include the manual curation and also do a functional annotation.

## Functional annotation

If I want Bianca to use this for the differential expression analysis, I should do a functional annotation, Doug used eggnogg, I can probably figure it out.

## Manual curation of Y-Tor

This is the workflow that Doug used and thankfully documented really well in `/proj/naiss2023-6-65/douglas/nobackup/Callosobruchus_maculatus/mTor`

### mTor consensus sequence 

They already have candidate proteins from Cmac that they think are the duplicates on the Y, so Doug starts with making an alignment between them, the autosomal one, and several other species (`sequences/mTor_sequences.faa`):
<details>
  <summary>list</summary>
<li>[Homo sapiens]: NP_004949.1 serine/threonine kinas mTOR isoform 1 </li>
<li>[Callosobruchus maculatus]: VEN43112.1 unnamed protein product </li>
<li>[Callosobruchus maculatus]: VEN51984.1 unnamed protein product </li>
<li>[Anoplophora glabripennis]: XP_018572076.1 serine/threonine-protein kinase Tor </li>
<li>[Diabrotica virgifera virgifera]: XP_028145210.1 serine/threonine-protein kinase Tor </li>
<li>[Leptinotarsa decemlineata]: XP_023015005.1 target of rapamycin </li>
<li>[Leptinotarsa decemlineata]: ALE20544.1 mTOR </li>
<li>[Brassicogethes aeneus]: CAH0563318.1 unnamed protein product </li>
<li>[Brassicogethes aeneus]: CAH0563403.1 unnamed protein product </li>
<li>[Rhynchophorus ferrugineus]: KAF7282696.1 hypothetical protein GWI33_002162 </li>
<li>[Aethina tumida]: XP_019880827.1 PREDICTED: LOW QUALITY PROTEIN: target of rapamycin-like </li>
<li>[Tribolium castaneum]: XP_971819.1 PREDICTED: serine/threonine-protein kinase mTOR </li>
<li>[Sitophilus oryzae]: XP_030750054.1 serine/threonine-protein kinase Tor </li>
<li>[Asbolus verrucosus]: RZC37432.1 serine/threonine-protein kinase mTOR </li>
<li>[Tenebrio molitor]: CAH1377105.1 unnamed protein product </li>
<li>[Tenebrio molitor]: AKB11618.1 target of rapamycin </li>
<li>[Ignelater luminosus]: KAF2880605.1 hypothetical protein ILUMI_25569 </li>
<li>[Lamprigera yunnana]: KAF5285925.1 hypothetical protein FQA39_LY04386 </li>
<li>[Nicrophorus vespilloides]: XP_017768823.1 PREDICTED: target of rapamycin </li>
<li>[Photinus pyralis]: XP_031352545.1 serine/threonine-protein kinase Tor </li>
<li>[Onthophagus taurus]: XP_022907797.1 target of rapamycin </li>
<li>[Agrilus planipennis]: XP_025831250.1 serine/threonine-protein kinase Tor </li>
<li>[Coccinella septempunctata]: XP_044759281.1 serine/threonine-protein kinase Tor </li>
<li>[Harmonia axyridis]: XP_045478375.1 serine/threonine-protein kinase Tor </li>
<li>[Abscondita terminalis]: KAF5280820.1 hypothetical protein FQR65_LT14927 </li>
<li>[Propylea japonica]: UIB01653.1 serine/threonine-protein kinase mTOR </li>

</details>

### blastp hits of the yTor in the other annotations

Since it has been identified in the other Cmac annotation already, I will not do it from scratch again in this one, I will just blast the Cmac yTor transcripts against the proteinfasta files of all the RNA annotations I did for the comparison in chapter 1. I expect them to be on contig `utg000322l_1`in the old assembly, and therefore `scaffold_26` or `scaffold_48` in the superscaffolded assembly (that the annotations are based on.)

* **Lome** results: There is only one hit and it is `Cmac_Lome_diverse_g1010.t1_1` for both the queries. the other ones have high e-values but the sequence identity is only 35% or lower. It is on `scafold_1`
*  **Lu** results: Same as above, only `Cmac_Lu2024_simple_g1006.t1_1` has a sequence identity above 35%. It is on `scafold_1`
*  **SI** results: Same as the other two, only `Cmac_SI_diverse_g963.t1_1`. It is on `scafold_1`

None of the RNA-based annotations detect the y-TOR, so I'm going to try with the uniform annotation that does not use RNAseq. it has much more hits with above99% sequence identity:
* VEN43112.1 (longer query)
    * C_maculatus_g11558.t1_1 : `scaffold_271`
    * C_maculatus_g23887.t1_1 : `scaffold_6`
    * C_maculatus_g23872.t1_1 : `scaffold_6`
    * C_maculatus_g23878.t1_1 : `scaffold_6`
    * C_maculatus_g23876.t1_1 : `scaffold_6`
    * C_maculatus_g23870.t1_1 : `scaffold_6`
    * C_maculatus_g23885.t1_1 : `scaffold_6`
    * C_maculatus_g11556.t1_1 : `scaffold_271`
  
* VEN51984.1 (shorter query)
    * C_maculatus_g23876.t1_1 : `scaffold_6`
    * C_maculatus_g23870.t1_1 : `scaffold_6`
    * C_maculatus_g23885.t1_1 : `scaffold_6`
    * C_maculatus_g11556.t1_1 : `scaffold_271`
    * other hits with low sequence identity
  
This does also not find the Y-Tor, but a bunch of stuff on Scaffold 6? None of the Y contigs in the old assembly are placed on scaffold 6. 

### blastp hits in non-superscaffolded annotation

Since no TOR copy is found on any Y contig in the superscaffolded annotation, I will check the uniform annotation i have for the non-superscaffolded one. These are the blast results for the same two query proteins as above. Two results with 100% sequence identity are highlighted, the rest are above 99%.
* VEN43112.1 (longer query)
  * C_maculatus_g11558.t1_1 (100% seq ident) : `utg000092` (`scaffold_1`)
  * C_maculatus_g23887.t1_1 : **`utg000322` (Y)**
  * C_maculatus_g23872.t1_1 : **`utg000322` (Y)**
  * C_maculatus_g23878.t1_1 : **`utg000322` (Y)**
  * C_maculatus_g23876.t1_1 : **`utg000322` (Y)** (100% seq ident for shorter query)
  * C_maculatus_g23870.t1_1 : **`utg000322` (Y)**
  * C_maculatus_g23885.t1_1 : **`utg000322` (Y)**
  * C_maculatus_g11556.t1_1 : `utg000092` (`scaffold_1`)

* VEN51984.1 (shorter query), all hits are also hits with longer query
  * C_maculatus_g23876.t1_1 (100% seq ident) : **`utg000322` (Y)**
  * C_maculatus_g23870.t1_1 : **`utg000322` (Y)**
  * C_maculatus_g23885.t1_1 : **`utg000322` (Y)**
  * C_maculatus_g11556.t1_1 : `utg000092` (`scaffold_1`)
