

<!-- Tab links -->
<script>
function openTab(evt, cityName) {
  // Declare all variables
  var i, tabcontent, tablinks;

  // Get all elements with class="tabcontent" and hide them
  tabcontent = document.getElementsByClassName("tabcontent");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  tablinks = document.getElementsByClassName("tablinks");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(cityName).style.display = "block";
  evt.currentTarget.className += " active";
}

</script>

<style>

/* Style the tab */
.tab {
  overflow: hidden;
  border: 1px solid #ccc;
  background-color: #f1f1f1;
}

/* Style the buttons that are used to open the tab content */
.tab button {
  background-color: inherit;
  float: left;
  border: none;
  outline: none;
  cursor: pointer;
  padding: 14px 16px;
  transition: 0.3s;
}

/* Change background color of buttons on hover */
.tab button:hover {
  background-color: #ddd;
}

/* Create an active/current tablink class */
.tab button.active {
  background-color: #ccc;
}

/* Style the tab content */
.tabcontent {
  display: none;
  padding: 6px 12px;
  border: 1px solid #ccc;
  border-top: none;
}

</style>


# Differential expression analysis on Larval RNAseq data

## 1. Preprocessing and testing

This analysis is based on a lot of work to properly filter the data and compare mapping pipelines to handle the high rRNA contamination rate and multimapping issues.

### 1.1 rRNA filtering and other reads preprocessing

This was done by Bianca Sarcani, all data here: `RNA_data_processing`. Basically:

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

From his results I made the choice to use STAR, see `mapping_and_transcriptome`. In short, with the new annotation the mapping rate is good enough with just STAR and not the complicated optimization algorithm that salmon uses.

### 1.3 yTor annotation manual curation

This work is based on the selection lines evaluated by Kaufmann *et al.* ([link](https://doi.org/10.1093/molbev/msad167)), but I am using a different version of the *C. maculatus* annotation that does not have the gene structures of the Y TOR copy number variation. Therefore I lift these gene structures from the original annotation and transform the coordinates to match my superscaffolded version of the assembly/annotation, see `mTOR_annotation`.

The genes are named `yTor-A`, `yTor-B`, and `yTor-C`, and the autosomal Tor is `transcript-30111`, or `gene-30110`.

## 2. STAR mapping

See the `RNA_mapping` directory in this repository for by-sample mapping rate information. Mean uniquely mapped reads: 83.76% and mean multimapped reads: 9.74%. 

## 3. DE analysis

I will use edgeR to get the logFC of all comparisons of interest below, but the plotting is done in python. 

### 3.1 yTor expression

After reading the data and filtering for minimum expression thresholds, yTor-A and yTor-C are expressed, but not yTor-B. Since we know that yTor-C (orange) in SL1 is the closest to the SL3 yTor, it makes sense that it shows expression in SL1 and SL3, while yTor-A, which is more diverged from the SL3 yTor, does not map reads from these samples. Furthermore, yTor-C is also the closes to aTor, which likely explains the female samples, which are mismapped reads that only map to yTor-C and not yTor-A.

The autosomal Tor `gene-30110` is expressed much higher than any y-linked copy, which makes sense since it is a highly conserved gene in the insulin signalling pathway. There seems to be no difference between SL1 and SL3.

<details>
<summary>normalized counts plots</summary> 

`(aTor,(yTor-SL3,(yTor-C,(yTor-A,yTor-B))))`

<p float="left">
  <img src="data/yTor_analysis/yTor_counts.png" width="49%" />
  <img src="data/yTor_analysis/all_Tor_counts.png" width="49%" />
</p>

The tick labels are according to this naming scheme: `SL`-`day`-`sample_ID`-`sex`.

</details>

### 3.2 PCA plots

PCAs are based on log-transformed normalized counts. Lines are SL1 and SL3 which are the small (1) and large (3) males respectively. Intuitively, the tor copy number is reversed, with the small males in SL1 having three copies, and the large males in SL3 having only one. The days are day 14, 16, or 18 of larval development. 


<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'All samples')">All samples</button>
  <button class="tablinks" onclick="openTab(event, 'Only one sex')">Only one sex</button>
</div>

<div id="All samples" class="tabcontent">


<p float="left">
  <img src="data/DE_figures/PCA_sex_line_all_counts.png" width="32%" />
  <img src="data/DE_figures/PCA_sex_day_all_counts.png" width="32%" />
</p>

Line is a clear separator (left), but not day (right). Day 14 seems to be mostly to the left, but 16 and 18 are across the entire range. There is a trend of females being more to the left and males more to the right but it is not as nice of a separation as line.

</div>

<div id="Only one sex" class="tabcontent">

<p float="left">
  <img src="data/DE_figures/PCA_M_day_line.png" width="32%" />
  <img src="data/DE_figures/PCA_F_day_line.png" width="32%" />
</p>

Since we are interested in the male variation and the females are mostly control, we have much fewer female than male samples.

line is the lagest difference, which is the same as when all samples are plotted. The days show a trend where day 18 is more to the left, 16 is intermediate, and 14 is to the right, but it is more of a gradual transition and not a clear separation. This makes sense since age is continuous, and we can only sample with limited precision, so it is likely that there is some age variation present within all the day categories that is represented here. (Except that one outlier in females (right) for SL3 day18 all the way on the right.)

</div>

<details>
<summary>I have also generated a MDS plot based on the edgeR data structure using plotMDS(), but they mostly show the same results as the PCA plots. Toggle down for MDS plots</summary>

<p float="left">
  <img src="data/DE_figures/MDS_males_only.png" width="32%" />
  <img src="data/DE_figures/MDS_males_and_females.png" width="32%" />
</p>

Males and females are kind of but not super clearly separated, but for only male samples, the SL1 and SL3 border is relatively clear.

</details>

### 3.3 DE analysis of sex-separated samples 

I started the differential expression analysis with only samples from one sex at a time. The contrasts are within each day (14, 16, 18), and always `SL1 - SL3`. `SL1` are the small males (three Tor copies), and all genes identified as "upregulated" are higher expressed in `SL1`. I am trying both `glmLRT` and `glmQLFTest` to test for differential expression, both fit negative binomial GLMs with the first one being more simple but having a higher false-positive error, while the second takes more variation in dispersion into account. I show both here and for the lines comparison, but I will only plot the results from `glmQLFTest`.

In all contrasts, genes that are differentially expressed in both females and males between the lines are excluded from the smear plots (see venn-diagram overlaps). 

#### Differential expression between SL1 and SL3

Additionally, in day 14 and day 16, there is a larger number of upregulated (higher in `SL1`) genes, while day 18 about the same number as up- and downregulated genes. This looks like there is a stronger line-difference in day 14 and 16, which becomes reduced in day 18.

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'males_lines_DE')">males (Smear plots)</button>
  <button class="tablinks" onclick="openTab(event, 'females_lines_DE')">females (Smear plots)</button>
  <button class="tablinks" onclick="openTab(event, 'males_lines_DE_volcano')">males (Volcano plots)</button>
  <button class="tablinks" onclick="openTab(event, 'females_lines_DE_volcano')">females (Volcano plots)</button>
</div>

<div id="males_lines_DE" class="tabcontent">

<table>

| `glmQLFTest`  | Day 14        | Day 16        | Day 18        |
| ------------- | ------------- | ------------- | ------------- |
| Downregulated | 93            | 199           | 16            |
| no difference | 10326         | 10164         | 10563         |
| Upregulated   | 148           | 186           | 36            |

</table>

<p float="left">
  <img src="data/DE_figures_python/smear_M_1-3_day14.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_1-3_day16.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_1-3_day18.png" width="32%" />
</p>

Upregulation here means that genes are expressed higher in SL1 than SL3. Day 14 and 16 have more significantly differentially expressed genes than day 18. In day 18, the larvae are close to pupation, which could mean that they are switching from gene expression related to grwoth and digestion to what they need for pupation instead, which is not related to the Y-haplotype difference any more, resulting in less DE between the lines on day 18.

<p float="left">
  <img src="data/DE_figures_python/Venn_males_age_by_line_bias.png" width="25%" />
</p>

</div>

<div id="females_lines_DE" class="tabcontent">

<table>

| `glmQLFTest`  | Day 14        | Day 16        | Day 18        |
| ------------- | ------------- | ------------- | ------------- |
| Downregulated | 5             | 8             | 29            |
| no difference | 9581          | 9553          | 9559          |
| Upregulated   | 1             | 8             | 47            |

</table>


<p float="left">
  <img src="data/DE_figures_python/smear_F_1-3_day14.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_1-3_day16.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_1-3_day18.png" width="32%" />
</p>

Fewer DE genes than males in day 14 and 16, and similar number in day 18. Within females, similar amounts of DE genes on day 14, 16, and 18, which is different from the male samples where day 18 is a clear outlier. This supports the hypothesis that this is caused by the line difference in growth which impacts day 14 and 16 more than 18, and therefore the DE genes here are not related to the growth differences between the lines that impact the males.

<p float="left">
  <img src="data/DE_figures_python/Venn_females_age_by_line_bias.png" width="25%" />
  <img src="data/DE_figures_python/Venn_day14_f_vs_m.png" width="24%" />
  <img src="data/DE_figures_python/Venn_day16_f_vs_m.png" width="24%" />
  <img src="data/DE_figures_python/Venn_day18_f_vs_m.png" width="24%" />
</p>

When comparing the genes that are DE between lines in the females to the ones in the males, we find that they are mostly a subset of the male ones in day 14 and day 16, and only start to diverge in day 18. These are probably not differences mediated by the Y-chromosome haplotype, and therefore should be excluded when comparing lines. Lists are below (toggle):

<details>
<summary>Day 14 (69 genes)</summary>

```
'gene-428738', 'gene-224697', 'gene-222350', 'gene-428765', 'gene-222600', 'gene-224875', 'gene-241001', 'gene-430032', 'gene-220028', 'gene-222486', 'gene-241055', 'gene-224357', 'gene-226245', 'gene-225738', 'gene-224968', 'gene-222531', 'gene-430263', 'gene-224201', 'gene-225107', 'gene-225236', 'gene-225140', 'gene-224227', 'gene-390616', 'gene-225709', 'gene-225325', 'gene-222332', 'gene-222519', 'gene-430314', 'gene-120952', 'gene-240871', 'gene-224860', 'gene-326873', 'gene-240929', 'gene-80359', 'gene-84970', 'gene-322912', 'gene-326849', 'gene-81427', 'gene-323148', 'gene-322927', 'gene-224782', 'gene-218529', 'gene-224743', 'gene-240623', 'gene-222383', 'gene-225173', 'gene-222365', 'gene-222344', 'gene-237881', 'gene-430068', 'gene-224956', 'gene-225720', 'gene-224682', 'gene-431701', 'gene-222555', 'gene-224896', 'gene-403809', 'gene-240910', 'gene-323803', 'gene-390956', 'gene-430080', 'gene-225635', 'gene-240833', 'gene-224593', 'gene-241126', 'gene-225030', 'gene-240691', 'gene-391222', 'gene-90157'
```

</details>

<details>
<summary>Day 16 (87 genes)</summary>

```
'gene-224697', 'gene-222600', 'gene-224875', 'gene-241001', 'gene-220028', 'gene-222486', 'gene-224357', 'gene-224968', 'gene-222159', 'gene-323148', 'gene-223773', 'gene-224782', 'gene-240623', 'gene-225173', 'gene-222344', 'gene-225720', 'gene-431701', 'gene-222555', 'gene-323803', 'gene-225635', 'gene-430080', 'gene-87700', 'gene-330102', 'gene-225030', 'gene-223419', 'gene-90157', 'gene-241262', 'gene-428738', 'gene-222350', 'gene-428765', 'gene-224079', 'gene-225325', 'gene-222332', 'gene-430314', 'gene-120952', 'gene-223491', 'gene-84970', 'gene-322927', 'gene-237881', 'gene-430068', 'gene-224956', 'gene-224682', 'gene-224896', 'g14784', 'gene-240833', 'gene-240691', 'gene-286545', 'gene-223318', 'gene-124877', 'gene-225738', 'gene-222531', 'gene-430263', 'gene-407280', 'gene-225140', 'gene-224227', 'gene-225709', 'gene-224890', 'gene-80359', 'gene-322912', 'gene-227370', 'gene-224743', 'gene-406796', 'gene-240910', 'gene-390956', 'gene-391222', 'gene-430032', 'gene-229506', 'gene-241055', 'gene-226245', 'gene-225107', 'gene-224201', 'gene-225236', 'gene-390616', 'gene-282853', 'gene-222519', 'gene-240871', 'gene-224860', 'gene-326873', 'gene-240929', 'gene-326849', 'gene-222383', 'gene-222365', 'gene-403809', 'gene-224593', 'gene-241126', 'gene-222746', 'gene-238407'
```

</details>

<details>
<summary>Day 18 (21 genes)</summary>

```
'gene-428738', 'gene-224697', 'gene-224875', 'gene-241055', 'gene-223758', 'gene-430263', 'gene-225236', 'gene-225325', 'gene-301479', 'gene-120952', 'gene-223491', 'gene-224860', 'gene-223773', 'gene-240623', 'gene-406796', 'gene-225720', 'gene-238849', 'gene-227308', 'gene-240833', 'gene-224593', 'gene-240691'
```

</details>

</div>




<div id="males_lines_DE_volcano" class="tabcontent">



<p float="left">
  <img src="data/DE_figures_python/smear_M_1-3_day14_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_1-3_day16_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_1-3_day18_volcano.png" width="32%" />
</p>
For venn-diagrams see smear plots

</div>

<div id="females_lines_DE_volcano" class="tabcontent">



<p float="left">
  <img src="data/DE_figures_python/smear_F_1-3_day14_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_1-3_day16_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_1-3_day18_volcano.png" width="32%" />
</p>

For venn-diagrams see smear plots

</div>




#### Scatterplots of pairwise comparisons

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'males_day_DE')">males</button>
  <button class="tablinks" onclick="openTab(event, 'females_day_DE')">females</button>
</div>

<!-- Tab content -->
<div id="males_day_DE" class="tabcontent">


<p float="left">
  <img src="data/DE_figures_python/LFC_scatter_males_day_14_vs_day_16_line_bias.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_males_day_14_vs_day_18_line_bias.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_males_day_16_vs_day_18_line_bias.png" width="32%" />
</p>

<table>
<tr><th>Day 14 and day 16 </th><th>Day 14 and day 18</th><th>Day 16 and day 18</th></tr>
<tr><td>

| category      | num DE genes  |
| ------------- | ------------- |
| both          | 176           |
| d14 exclusive | 134           |
| d16 exclusive | 296           |
| shared with f | 89            |

</td><td>

| category      | num DE genes  |
| ------------- | ------------- |
| both          | 47            |
| d14 exclusive | 263           |
| d18 exclusive | 26            |
| shared with f | 76            |

</td><td>

| category      | num DE genes  |
| ------------- | ------------- |
| both          | 48            |
| d16 exclusive | 424           |
| d18 exclusive | 25            |
| shared with f | 91            |

</td></tr> </table>

('shared with f' genes are genes that are shared with f in one or both time points in the comparison, so the numbers are a little higher than the lists above)

</div>

<div id="females_day_DE" class="tabcontent">


<p float="left">
  <img src="data/DE_figures_python/LFC_scatter_females_day_14_vs_day_16_line_bias.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_females_day_14_vs_day_18_line_bias.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_females_day_16_vs_day_18_line_bias.png" width="32%" />
</p>

<table>
<tr><th>Day 14 and day 16 </th><th>Day 14 and day 18</th><th>Day 16 and day 18</th></tr>
<tr><td>

| category      | num DE genes  |
| ------------- | ------------- |
| both          | 67            |
| d14 exclusive | 36            |
| d16 exclusive | 8             |

</td><td>

| category      | num DE genes  |
| ------------- | ------------- |
| both          | 60            |
| d14 exclusive | 37            |
| d18 exclusive | 15            |

</td><td>

| category      | num DE genes  |
| ------------- | ------------- |
| both          | 75            |
| d16 exclusive | 22            |
| d18 exclusive | 28            |

</td></tr> </table>

</div>


#### Interaction of line by day effects

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'males_line_day_DE')">males (Smear plots)</button>
  <button class="tablinks" onclick="openTab(event, 'females_line_day_DE')">females (Smear plots)</button>
  <button class="tablinks" onclick="openTab(event, 'males_line_day_DE_volcano')">males (Volcano plots)</button>
<button class="tablinks" onclick="openTab(event, 'females_line_day_DE_volcano')">females (Volcano plots)</button>
</div>

<div id="males_line_day_DE" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_M_SL_1_3_14_16.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL_1_3_18_14.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL_1_3_18_16.png" width="32%" />
</p>

Contrasts specified like this: (either option is fine, gives the exact same results)
```
(SL1_18-SL1_14)-(SL3_18-SL3_14), (SL1_14-SL1_16)-(SL3_14-SL3_16), (SL1_18-SL1_16)-(SL3_18-SL3_16)
(SL1_14-SL3_14)-(SL1_16-SL3_16), (SL1_18-SL3_18)-(SL1_16-SL3_16), (SL1_18-SL3_18)-(SL1_14-SL3_14)
```

</div>

<div id="females_line_day_DE" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_F_SL_1_3_14_16.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL_1_3_18_14.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL_1_3_18_16.png" width="32%" />
</p>
Contrasts specified like this: (either option is fine, gives the exact same results)
```
(SL1_18-SL1_14)-(SL3_18-SL3_14), (SL1_14-SL1_16)-(SL3_14-SL3_16), (SL1_18-SL1_16)-(SL3_18-SL3_16)
(SL1_14-SL3_14)-(SL1_16-SL3_16), (SL1_18-SL3_18)-(SL1_16-SL3_16), (SL1_18-SL3_18)-(SL1_14-SL3_14)
```

</div>


<div id="males_line_day_DE_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_M_SL_1_3_14_16_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL_1_3_18_14_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL_1_3_18_16_volcano.png" width="32%" />
</p>

Contrasts specified like this: (either option is fine, gives the exact same results)
```
(SL1_18-SL1_14)-(SL3_18-SL3_14), (SL1_14-SL1_16)-(SL3_14-SL3_16), (SL1_18-SL1_16)-(SL3_18-SL3_16)
(SL1_14-SL3_14)-(SL1_16-SL3_16), (SL1_18-SL3_18)-(SL1_16-SL3_16), (SL1_18-SL3_18)-(SL1_14-SL3_14)
```

</div>

<div id="females_line_day_DE_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_F_SL_1_3_14_16_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL_1_3_18_14_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL_1_3_18_16_volcano.png" width="32%" />
</p>

Contrasts specified like this: (either option is fine, gives the exact same results)
```
(SL1_18-SL1_14)-(SL3_18-SL3_14), (SL1_14-SL1_16)-(SL3_14-SL3_16), (SL1_18-SL1_16)-(SL3_18-SL3_16)
(SL1_14-SL3_14)-(SL1_16-SL3_16), (SL1_18-SL3_18)-(SL1_16-SL3_16), (SL1_18-SL3_18)-(SL1_14-SL3_14)
```

</div>

<details>
<summary>Interpretation of line:day interaction</summary>

Interaction line by day: when genes are upregulated (in the line contrast) in one day and then downregulated in another day. These cases exist, see scatterplots above, but there is apparently none that overcome the significance threshold, only one gene is significantly DE when the day contrast is 18-14. 

The original plan was to then plot the LFC of the day contrast for both lines for the genes that are significant in the interaction, so `LFC(SL1d18-SL1d14)` by `LFC(SL3d18-SL3d14)`.


</details>



#### Differential expression within lines between day18 and mean(day14,day16)

I hypothesize that day 14 and 16 are where a lot of growth happens and SL1 and SL3 differ, while day 18 is the transition to pupation where the line differences become less substantial. I will therefore see what genes are involved in growth specifically by looking at the contrast between day 18 and the mean of day 14 and day 16. 

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'males_lines_DE_tables')">males (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'females_lines_DE_tables')">females (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'males_lines_DE_tables_volcano')">males (Volcano plot)</button>
  <button class="tablinks" onclick="openTab(event, 'females_lines_DE_tables_volcano')">females (Volcano plot)</button>
</div>

<div id="males_lines_DE_tables" class="tabcontent">

<table>

| `glmQLFTest`  | Line 1        | Line 3        |
| ------------- | ------------- | ------------- |
| Downregulated | 639           | 1271          |
| no difference | 8150          | 6653          |
| Upregulated   | 1850          | 2766          |

</table>

<p float="left">
  <img src="data/DE_figures_python/smear_M_SL1.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL3.png" width="32%" />
</p>

<p float="left">
  <img src="data/DE_figures_python/Venn_males_line_by_age_bias.png" width="25%" />
  <img src="data/DE_figures_python/LFC_scatter_males_line_by_age_bias.png" width="25%" />
</p>

Upregulation means that genes get transcribed more for day 18. In general, SL3 has more change than SL1

| category      | num DE genes  |
| ------------- | ------------- |
| both          | 1984          |
| SL1 exclusive | 502           |
| SL3 exclusive | 1999          |
| shared with f | 93            |

</div>

<div id="females_lines_DE_tables" class="tabcontent">

<table>

| `glmQLFTest`  | Line 1        | Line 3        |
| ------------- | ------------- | ------------- |
| Downregulated | 0             | 11            |
| no difference | 9656          | 9569          |
| Upregulated   | 0             | 76            |

</table>

<p float="left">
  <img src="data/DE_figures_python/smear_F_SL1.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL3.png" width="32%" />
</p>

<p float="left">
  <img src="data/DE_figures_python/Venn_females_line_by_age_bias.png" width="25%" />
  <img src="data/DE_figures_python/LFC_scatter_females_line_by_age_bias.png" width="25%" />
  <img src="data/DE_figures_python/Venn_SL3_f_vs_m.png" width="25%" />
</p>

In SL3, 74 of 87 DE genes in females are also DE in males, 13 are unique to females. 

| category      | num DE genes  |
| ------------- | ------------- |
| both          | 0             |
| SL1 exclusive | 0             |
| SL3 exclusive | 87            |

<details>
<summary>74 genes DE in SL3 in m and f</summary>

```
'gene-99775', 'gene-40274', 'gene-304827', 'gene-92346', 'gene-306335', 'gene-285669', 'gene-120763', 'gene-2286', 'gene-97407', 'gene-232392', 'gene-328941', 'gene-166511', 'gene-39692', 'gene-384091', 'gene-74686', 'gene-122220', 'gene-218723', 'gene-414353', 'gene-312890', 'gene-153482', 'gene-39770', 'gene-132340', 'gene-253632', 'gene-378608', 'gene-206556', 'gene-336703', 'gene-21229', 'gene-166391', 'gene-120784', 'gene-87502', 'gene-317372', 'gene-73253', 'gene-211196', 'gene-9548', 'gene-60190', 'gene-234650', 'gene-410057', 'gene-121262', 'gene-100036', 'gene-227137', 'gene-75744', 'gene-279912', 'gene-343203', 'gene-233901', 'gene-163028', 'gene-238407', 'gene-39680', 'gene-198700', 'gene-231228', 'gene-410209', 'gene-143368', 'gene-388261', 'gene-30328', 'gene-226944', 'gene-334263', 'gene-130081', 'gene-47823', 'gene-228519', 'gene-350792', 'gene-277340', 'gene-182683', 'gene-60151', 'gene-206576', 'gene-202718', 'gene-130096', 'gene-377275', 'gene-244780', 'gene-288834', 'gene-189246', 'gene-48535', 'gene-205011', 'gene-253157', 'gene-62891', 'gene-69698'
```

</details>

</div>


<div id="males_lines_DE_tables_volcano" class="tabcontent">


<p float="left">
  <img src="data/DE_figures_python/smear_M_SL1_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL3_volcano.png" width="32%" />
</p>


</div>

<div id="females_lines_DE_tables_volcano" class="tabcontent">


<p float="left">
  <img src="data/DE_figures_python/smear_F_SL1_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL3_volcano.png" width="32%" />
</p>


</div>





#### Differential expression within lines with all pairwise day comparisons

Males have a massive change in gene regulation, both up- and down between day 18 and the earlier stages (larger difference in SL3) but females change not at all in SL1 and only a little bit in SL3. To be sure that this is real I also do individual comparisons between all the days even though that is more difficult to interpret.

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'males_lines_DE_tables_single_days')">males (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'females_lines_DE_tables_single_days')">females (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'males_lines_DE_tables_single_days_volcano')">males (Volcano plot)</button>
  <button class="tablinks" onclick="openTab(event, 'females_lines_DE_tables_single_days_volcano')">females (Volcano plot)</button>
</div>

<div id="males_lines_DE_tables_single_days" class="tabcontent">

<table>

<tr><th>SL1</th><th>SL3</th></tr>
<tr><td>

| `glmQLFTest`  | d18-d14       | d18-d16       | d14-d16       |
| ------------- | ------------- | ------------- | ------------- |
| Downregulated | 1213          | 15            | 1487          |
| no difference | 7086          | 10559         | 8257          |
| Upregulated   | 2337          | 62            | 892           |

</td><td>

| `glmQLFTest`  | d18-d14       | d18-d16       | d14-d16       |
| ------------- | ------------- | ------------- | ------------- |
| Downregulated | 1159          | 606           | 607           |
| no difference | 7125          | 8956          | 9972          |
| Upregulated   | 2352          | 1074          | 57            |

</td><td> </table>

<p float="left">
  <img src="data/DE_figures_python/smear_M_SL1_18_14.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL1_18_16.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL1_14_16.png" width="32%" />
</p>

<p float="left">
  <img src="data/DE_figures_python/smear_M_SL3_18_14.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL3_18_16.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL3_14_16.png" width="32%" />
</p>

in SL1, the difference is more between day 14 and day 16/18, since the day 18 - day 16 plot has only very few DE genes. In SL3, the difference is more between day 18 and day14/16, since the day14/16 plot has the least DE genes. Maybe this means that SL1 starts preparation for pupation between day 14 and 16? There is no clear indication for SL3, which has larger individuals, so maybe it is even earlier? Note that this does not show a morphological change yet, the earliest individuals show a start of pupation only at day 18 (did we sequence those then?).

</div>

<div id="females_lines_DE_tables_single_days" class="tabcontent">

<table>

<tr><th>SL1</th><th>SL3</th></tr>
<tr><td>

| `glmQLFTest`  | d18-d14       | d18-d16       | d14-d16       |
| ------------- | ------------- | ------------- | ------------- |
| Downregulated | 55            | 0             | 0             |
| no difference | 9595          | 9656          | 9656          |
| Upregulated   | 6             | 0             | 0             |

</td><td>

| `glmQLFTest`  | d18-d14       | d18-d16       | d14-d16       |
| ------------- | ------------- | ------------- | ------------- |
| Downregulated | 30            | 1             | 2             |
| no difference | 9558          | 9655          | 9653          |
| Upregulated   | 68            | 0             | 1             |

</td><td> </table>

<p float="left">
  <img src="data/DE_figures_python/smear_F_SL1_18_14.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL1_14_16.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL1_18_16.png" width="32%" />
</p>

<p float="left">
  <img src="data/DE_figures_python/smear_F_SL3_18_14.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL3_14_16.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL3_18_16.png" width="32%" />
</p>

Only difference between day 14 and 18, no other developmental cues? Do only the males do special stuff before pupation?

</div>


<div id="males_lines_DE_tables_single_days_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_M_SL1_18_14_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL1_18_16_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL1_14_16_volcano.png" width="32%" />
</p>

<p float="left">
  <img src="data/DE_figures_python/smear_M_SL3_18_14_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL3_18_16_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_M_SL3_14_16_volcano.png" width="32%" />
</p>

</div>

<div id="females_lines_DE_tables_single_days_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_F_SL1_18_14_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL1_14_16_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL1_18_16_volcano.png" width="32%" />
</p>

<p float="left">
  <img src="data/DE_figures_python/smear_F_SL3_18_14_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL3_14_16_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_F_SL3_18_16_volcano.png" width="32%" />
</p>


</div>



### 3.4 DE analysis of line-separated samples

I will now split the data by line to see sex differences in expression during the development stages. 

<details>
<summary>MDS plots for SL1 and SL3</summary>

<p float="left">
  <img src="data/DE_figures/MDS_SL1_only.png" width="32%" />
  <img src="data/DE_figures/MDS_SL3_only.png" width="32%" />
</p>

</details>


#### Differential expression between males and females

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'SL1_smear')">SL1 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'SL3_smear')">SL3 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'SL1_volcano')">SL1 (Volcano plot)</button>
  <button class="tablinks" onclick="openTab(event, 'SL3_volcano')">SL3 (Volcano plot)</button>
</div>

<div id="SL1_smear" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_SL1_day14_F-M.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL1_day16_F-M.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL1_day18_F-M.png" width="32%" />
</p>

| `glmQLFTest`  | Day 14        | Day 16        | Day 18        |
| ------------- | ------------- | ------------- | ------------- |
| Downregulated | 265           | 767           | 1039          |
| no difference | 10724         | 10155         | 9889          |
| Upregulated   | 4             | 71            | 65            |

<p float="left">
  <img src="data/DE_figures_python/Venn_SL1_age_by_sex_bias.png" width="25%" />
</p>

Day 14 has the least sex differences, which is when they are for sure still in larval development. Above, I speculate that SL1 starts pupation (preparation) at day 16, which is when sex differences can first be introduced, so it makes sense that the sex-bias (mostly only male-bias) increases during those days. The venn-diagramm shows that day 16 and 18 mostly have the same DE genes.

</div>

<div id="SL3_smear" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_SL3_day14_F-M.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL3_day16_F-M.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL3_day18_F-M.png" width="32%" />
</p>

| `glmQLFTest`  | Day 14        | Day 16        | Day 18        |
| ------------- | ------------- | ------------- | ------------- |
| Downregulated | 569           | 875           | 861           |
| no difference | 10015         | 9676          | 9721          |
| Upregulated   | 3             | 36            | 5             |

<p float="left">
  <img src="data/DE_figures_python/Venn_SL3_age_by_sex_bias.png" width="25%" />
</p>

A similar trend as SL1 but not as strong, more genes are already sex-biased on day 14. The male-bias is more severe. This supports the speculation that the larger SL3 individuals start preparation for pupation even earlier (above shows that in SL1, day 16/18 is clearly very similar, while 14/16 and 16/18 differ a lot, maybe indicating a major developmental milestone before day 16. In SL3 all pairwise day comparisons are very different, and we do not see the same milestone. So this means it either happens earlier or later, and these results indicate earlier). 

</div>



<div id="SL1_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_SL1_day14_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL1_day16_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL1_day18_F-M_volcano.png" width="32%" />
</p>


</div>

<div id="SL3_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_SL3_day14_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL3_day16_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL3_day18_F-M_volcano.png" width="32%" />
</p>

</div>


#### Scatterplots of pairwise comparisons

Since we are explicitly checking for sex-bias, it makes no sense to highlight genes with shared sex-bias in SL1 and SL3 for example, so I exclude the highlight here.

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'SL1_overlap')">SL1</button>
  <button class="tablinks" onclick="openTab(event, 'SL3_overlap')">SL3</button>
</div>

<div id="SL1_overlap" class="tabcontent">


<p float="left">
  <img src="data/DE_figures_python/LFC_scatter_SL1_day_14_vs_day_18_sex_bias.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_SL1_day_14_vs_day_16_sex_bias.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_SL1_day_16_vs_day_18_sex_bias.png" width="32%" />
</p>


|               | Day 14 - Day 16 | Day 14 - Day 18 | Day 16 - Day 18 |
| ------------- | --------------- | --------------- | --------------- |
| shared        | 257             | 259             | 746             |
| younger only  | 12              | 10              | 92              |
| older only    | 581             | 845             | 358             |


</div>

<div id="SL3_overlap" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/LFC_scatter_SL3_day_14_vs_day_18_sex_bias.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_SL3_day_14_vs_day_16_sex_bias.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_SL3_day_16_vs_day_18_sex_bias.png" width="32%" />
</p>


|               | Day 14 - Day 16 | Day 14 - Day 18 | Day 16 - Day 18 |
| ------------- | --------------- | --------------- | --------------- |
| shared        | 535             | 516             | 769             |
| younger only  | 37              | 56              | 142             |
| older only    | 376             | 350             | 97              |

</div>

#### interaction of sex by day effects

<span style="color:#BD351E;"><b>QUESTION: </b></span>What kind of interaction is sensible here? This is the only one I could come up with: `(F_18-M_18)-(F_14-M_14)` so the interaction highlights genes that are female biased in one time point and then male biased in the other. 

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'SL1_interaction_smear')">SL1 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'SL3_interaction_smear')">SL3 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'SL1_interaction_volcano')">SL1 (Volcano plot)</button>
  <button class="tablinks" onclick="openTab(event, 'SL3_interaction_volcano')">SL3 (Volcano plot)</button>
</div>

<div id="SL1_interaction_smear" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_SL1_day14_16_F-M.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL1_day18_14_F-M.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL1_day18_16_F-M.png" width="32%" />
</p>

`(F_14-M_14)-(F_16-M_16), (F_18-M_18)-(F_14-M_14), (F_18-M_18)-(F_16-M_16)`

| `glmQLFTest`  | Day 14-16     | Day 18-14     | Day 18-16     |
| ------------- | ------------- | ------------- | ------------- |
| Downregulated | 0             | 66            | 0             |
| no difference | 10993         | 10926         | 10993         |
| Upregulated   | 0             | 1             | 0             |

</div>


<div id="SL3_interaction_smear" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_SL3_day14_16_F-M.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL3_day18_14_F-M.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL3_day18_16_F-M.png" width="32%" />
</p>

`(F_14-M_14)-(F_16-M_16), (F_18-M_18)-(F_14-M_14), (F_18-M_18)-(F_16-M_16)`

| `glmQLFTest`  | Day 14-16     | Day 18-14     | Day 18-16     |
| ------------- | ------------- | ------------- | ------------- |
| Downregulated | 0             | 28            | 0             |
| no difference | 10587         | 10553         | 10587         |
| Upregulated   | 0             | 6             | 0             |

</div>


<div id="SL1_interaction_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_SL1_day14_16_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL1_day18_14_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL1_day18_16_F-M_volcano.png" width="32%" />
</p>

`(F_14-M_14)-(F_16-M_16), (F_18-M_18)-(F_14-M_14), (F_18-M_18)-(F_16-M_16)`

| `glmQLFTest`  | Day 14-16     | Day 18-14     | Day 18-16     |
| ------------- | ------------- | ------------- | ------------- |
| Downregulated | 0             | 66            | 0             |
| no difference | 10993         | 10926         | 10993         |
| Upregulated   | 0             | 1             | 0             |

</div>


<div id="SL3_interaction_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_SL3_day14_16_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL3_day18_14_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_SL3_day18_16_F-M_volcano.png" width="32%" />
</p>

`(F_14-M_14)-(F_16-M_16), (F_18-M_18)-(F_14-M_14), (F_18-M_18)-(F_16-M_16)`

| `glmQLFTest`  | Day 14-16     | Day 18-14     | Day 18-16     |
| ------------- | ------------- | ------------- | ------------- |
| Downregulated | 0             | 28            | 0             |
| no difference | 10587         | 10553         | 10587         |
| Upregulated   | 0             | 6             | 0             |

</div>

### 3.5 DE analysis of day-separated samples

Analyze line and sex differences.



#### Differential expression between SL1 and SL3

In Females, there should be no line bias, so all genes that show line-bias in females (see lists in tabs) are excluded for the day contrasts.

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'day14_SL1-3_smear')">Day 14 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day14_SL1-3_volcano')">Day 14 (Volcano plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day16_SL1-3_smear')">Day 16 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day16_SL1-3_volcano')">Day 16 (Volcano plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day18_SL1-3_smear')">Day 18 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day18_SL1-3_volcano')">Day 18 (Volcano plot)</button>
</div>

<div id="day14_SL1-3_smear" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_day14_F_1-3.png" width="32%" />
  <img src="data/DE_figures_python/smear_day14_M_1-3.png" width="32%" />
  <img src="data/DE_figures_python/Venn_day14_line_bias_by_sex.png" width="25%" />
</p>

<details>
<summary>130 DE genes shared between males and females (excluded from plots, toggle down to see list)</summary>

```
['gene-327441', 'gene-241001', 'gene-403652', 'gene-403851', 'gene-224860', 'gene-218086', 'gene-390956', 'gene-224682', 'gene-240910', 'gene-81551', 'gene-392224', 'gene-240983', 'gene-326825', 'gene-222519', 'gene-237881', 'gene-224968', 'gene-220544', 'gene-226245', 'gene-391198', 'gene-403818', 'gene-239553', 'gene-90157', 'gene-84970', 'gene-219019', 'gene-407253', 'gene-224845', 'gene-81640', 'gene-222383', 'gene-220249', 'gene-430080', 'gene-430044', 'gene-391222', 'gene-222365', 'gene-240833', 'gene-428765', 'gene-240691', 'gene-80359', 'gene-392248', 'gene-390678', 'gene-225629', 'gene-120952', 'gene-428756', 'gene-224956', 'gene-403809', 'gene-224227', 'gene-240935', 'gene-225720', 'gene-225140', 'gene-224614', 'gene-224357', 'gene-392159', 'gene-222531', 'gene-225325', 'gene-117712', 'gene-84949', 'gene-222344', 'gene-406796', 'gene-240929', 'gene-392290', 'gene-81599', 'gene-225107', 'gene-231925', 'gene-214979', 'gene-220028', 'gene-260693', 'gene-224277', 'gene-224250', 'gene-224782', 'gene-390616', 'gene-224593', 'gene-322912', 'gene-431701', 'gene-80466', 'gene-403706', 'gene-323148', 'gene-323803', 'gene-283443', 'gene-326909', 'gene-225173', 'gene-243308', 'gene-430068', 'gene-222430', 'gene-240871', 'gene-237318', 'gene-224697', 'gene-222332', 'gene-241126', 'gene-89234', 'gene-406468', 'gene-328764', 'gene-326810', 'gene-88715', 'gene-222486', 'gene-224743', 'gene-221953', 'gene-225635', 'gene-222501', 'gene-428738', 'gene-430032', 'gene-241108', 'gene-81572', 'gene-392186', 'gene-222555', 'gene-81427', 'gene-225738', 'gene-403700', 'gene-234575', 'gene-240623', 'gene-222350', 'gene-395080', 'gene-225236', 'gene-224201', 'gene-241055', 'gene-225709', 'gene-326873', 'gene-326849', 'gene-322927', 'gene-80484', 'gene-395143', 'gene-225030', 'gene-222600', 'gene-224875', 'gene-390637', 'gene-224307', 'gene-403583', 'gene-224896', 'gene-430263', 'gene-428747', 'gene-430314', 'gene-403902']
```

</details>

</div>

<div id="day14_SL1-3_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_day14_F_1-3_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_day14_M_1-3_volcano.png" width="32%" />
  <img src="data/DE_figures_python/Venn_day14_line_bias_by_sex.png" width="25%" />
</p>

<details>
<summary>130 DE genes shared between males and females (excluded from plots, toggle down to see list)</summary>

```
['gene-327441', 'gene-241001', 'gene-403652', 'gene-403851', 'gene-224860', 'gene-218086', 'gene-390956', 'gene-224682', 'gene-240910', 'gene-81551', 'gene-392224', 'gene-240983', 'gene-326825', 'gene-222519', 'gene-237881', 'gene-224968', 'gene-220544', 'gene-226245', 'gene-391198', 'gene-403818', 'gene-239553', 'gene-90157', 'gene-84970', 'gene-219019', 'gene-407253', 'gene-224845', 'gene-81640', 'gene-222383', 'gene-220249', 'gene-430080', 'gene-430044', 'gene-391222', 'gene-222365', 'gene-240833', 'gene-428765', 'gene-240691', 'gene-80359', 'gene-392248', 'gene-390678', 'gene-225629', 'gene-120952', 'gene-428756', 'gene-224956', 'gene-403809', 'gene-224227', 'gene-240935', 'gene-225720', 'gene-225140', 'gene-224614', 'gene-224357', 'gene-392159', 'gene-222531', 'gene-225325', 'gene-117712', 'gene-84949', 'gene-222344', 'gene-406796', 'gene-240929', 'gene-392290', 'gene-81599', 'gene-225107', 'gene-231925', 'gene-214979', 'gene-220028', 'gene-260693', 'gene-224277', 'gene-224250', 'gene-224782', 'gene-390616', 'gene-224593', 'gene-322912', 'gene-431701', 'gene-80466', 'gene-403706', 'gene-323148', 'gene-323803', 'gene-283443', 'gene-326909', 'gene-225173', 'gene-243308', 'gene-430068', 'gene-222430', 'gene-240871', 'gene-237318', 'gene-224697', 'gene-222332', 'gene-241126', 'gene-89234', 'gene-406468', 'gene-328764', 'gene-326810', 'gene-88715', 'gene-222486', 'gene-224743', 'gene-221953', 'gene-225635', 'gene-222501', 'gene-428738', 'gene-430032', 'gene-241108', 'gene-81572', 'gene-392186', 'gene-222555', 'gene-81427', 'gene-225738', 'gene-403700', 'gene-234575', 'gene-240623', 'gene-222350', 'gene-395080', 'gene-225236', 'gene-224201', 'gene-241055', 'gene-225709', 'gene-326873', 'gene-326849', 'gene-322927', 'gene-80484', 'gene-395143', 'gene-225030', 'gene-222600', 'gene-224875', 'gene-390637', 'gene-224307', 'gene-403583', 'gene-224896', 'gene-430263', 'gene-428747', 'gene-430314', 'gene-403902']
```

</details>

</div>

<div id="day16_SL1-3_smear" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_day16_F_1-3.png" width="32%" />
  <img src="data/DE_figures_python/smear_day16_M_1-3.png" width="32%" />
  <img src="data/DE_figures_python/Venn_day16_line_bias_by_sex.png" width="25%" />
  
</p>

<details>
<summary>142 DE genes shared between males and females (excluded from plots, toggle down to see list)</summary>

```
['gene-241001', 'gene-282746', 'gene-224860', 'gene-227308', 'gene-390956', 'gene-224682', 'gene-240910', 'gene-399475', 'gene-240983', 'gene-222519', 'gene-282458', 'gene-237881', 'gene-224968', 'gene-220544', 'gene-226245', 'gene-240860', 'gene-391198', 'gene-282524', 'gene-403818', 'gene-90157', 'gene-84970', 'gene-224845', 'gene-81640', 'gene-222383', 'gene-220249', 'gene-430080', 'gene-430044', 'gene-391222', 'gene-222365', 'gene-240833', 'gene-428765', 'gene-240691', 'gene-80359', 'gene-390678', 'gene-225629', 'gene-238849', 'gene-282620', 'gene-120952', 'gene-428756', 'gene-224956', 'gene-403809', 'gene-282853', 'gene-224227', 'gene-225720', 'gene-225140', 'gene-238407', 'gene-239506', 'gene-224357', 'gene-222531', 'gene-241682', 'gene-225325', 'gene-84949', 'gene-222344', 'gene-406796', 'gene-240929', 'gene-282551', 'gene-240638', 'gene-392290', 'gene-225158', 'gene-225107', 'gene-220028', 'gene-260693', 'gene-224277', 'gene-224782', 'gene-390616', 'gene-224593', 'gene-322912', 'gene-431701', 'gene-428113', 'gene-380466', 'gene-323148', 'gene-323803', 'gene-282665', 'gene-326909', 'gene-428104', 'gene-225173', 'gene-227370', 'gene-399317', 'gene-224890', 'gene-282347', 'gene-282784', 'gene-430068', 'gene-282491', 'gene-400393', 'gene-222430', 'gene-399484', 'gene-240871', 'gene-224697', 'gene-222332', 'gene-241126', 'gene-89234', 'gene-399424', 'gene-282362', 'g14784', 'gene-223491', 'gene-400384', 'gene-222486', 'gene-224743', 'gene-243299', 'gene-221953', 'gene-225635', 'gene-222501', 'gene-428738', 'gene-400402', 'gene-430032', 'gene-241108', 'gene-282398', 'gene-81572', 'gene-392186', 'gene-222555', 'gene-81427', 'gene-225738', 'gene-282701', 'gene-282886', 'gene-240623', 'gene-222159', 'gene-222350', 'gene-395080', 'gene-225236', 'gene-224201', 'gene-241055', 'gene-286545', 'gene-225709', 'gene-431030', 'gene-326849', 'gene-322927', 'gene-395143', 'gene-225030', 'gene-222600', 'gene-224875', 'gene-223773', 'gene-390637', 'gene-224307', 'gene-231854', 'gene-224896', 'gene-430263', 'gene-428747', 'gene-399270', 'gene-241262', 'gene-430314', 'gene-242512', 'gene-282590']
```

</details>

</div>

<div id="day16_SL1-3_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_day16_F_1-3_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_day16_M_1-3_volcano.png" width="32%" />
  <img src="data/DE_figures_python/Venn_day16_line_bias_by_sex.png" width="25%" />
  
</p>

<details>
<summary>142 DE genes shared between males and females (excluded from plots, toggle down to see list)</summary>

```
['gene-241001', 'gene-282746', 'gene-224860', 'gene-227308', 'gene-390956', 'gene-224682', 'gene-240910', 'gene-399475', 'gene-240983', 'gene-222519', 'gene-282458', 'gene-237881', 'gene-224968', 'gene-220544', 'gene-226245', 'gene-240860', 'gene-391198', 'gene-282524', 'gene-403818', 'gene-90157', 'gene-84970', 'gene-224845', 'gene-81640', 'gene-222383', 'gene-220249', 'gene-430080', 'gene-430044', 'gene-391222', 'gene-222365', 'gene-240833', 'gene-428765', 'gene-240691', 'gene-80359', 'gene-390678', 'gene-225629', 'gene-238849', 'gene-282620', 'gene-120952', 'gene-428756', 'gene-224956', 'gene-403809', 'gene-282853', 'gene-224227', 'gene-225720', 'gene-225140', 'gene-238407', 'gene-239506', 'gene-224357', 'gene-222531', 'gene-241682', 'gene-225325', 'gene-84949', 'gene-222344', 'gene-406796', 'gene-240929', 'gene-282551', 'gene-240638', 'gene-392290', 'gene-225158', 'gene-225107', 'gene-220028', 'gene-260693', 'gene-224277', 'gene-224782', 'gene-390616', 'gene-224593', 'gene-322912', 'gene-431701', 'gene-428113', 'gene-380466', 'gene-323148', 'gene-323803', 'gene-282665', 'gene-326909', 'gene-428104', 'gene-225173', 'gene-227370', 'gene-399317', 'gene-224890', 'gene-282347', 'gene-282784', 'gene-430068', 'gene-282491', 'gene-400393', 'gene-222430', 'gene-399484', 'gene-240871', 'gene-224697', 'gene-222332', 'gene-241126', 'gene-89234', 'gene-399424', 'gene-282362', 'g14784', 'gene-223491', 'gene-400384', 'gene-222486', 'gene-224743', 'gene-243299', 'gene-221953', 'gene-225635', 'gene-222501', 'gene-428738', 'gene-400402', 'gene-430032', 'gene-241108', 'gene-282398', 'gene-81572', 'gene-392186', 'gene-222555', 'gene-81427', 'gene-225738', 'gene-282701', 'gene-282886', 'gene-240623', 'gene-222159', 'gene-222350', 'gene-395080', 'gene-225236', 'gene-224201', 'gene-241055', 'gene-286545', 'gene-225709', 'gene-431030', 'gene-326849', 'gene-322927', 'gene-395143', 'gene-225030', 'gene-222600', 'gene-224875', 'gene-223773', 'gene-390637', 'gene-224307', 'gene-231854', 'gene-224896', 'gene-430263', 'gene-428747', 'gene-399270', 'gene-241262', 'gene-430314', 'gene-242512', 'gene-282590']
```

</details>

</div>

<div id="day18_SL1-3_smear" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_day18_F_1-3.png" width="32%" />
  <img src="data/DE_figures_python/smear_day18_M_1-3.png" width="32%" />
  <img src="data/DE_figures_python/Venn_day18_line_bias_by_sex.png" width="25%" />
  
</p>

</div>

<div id="day18_SL1-3_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_day18_F_1-3_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_day18_M_1-3_volcano.png" width="32%" />
  <img src="data/DE_figures_python/Venn_day18_line_bias_by_sex.png" width="25%" />
  
</p>

</div>






#### Differential expression between males and females

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'day14_FM_smear')">Day 14 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day14_FM_volcano')">Day 14 (Volcano plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day16_FM_smear')">Day 16 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day16_FM_volcano')">Day 16 (Volcano plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day18_FM_smear')">Day 18 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day18_FM_volcano')">Day 18 (Volcano plot)</button>
</div>

<div id="day14_FM_smear" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_day14_SL1_F-M.png" width="32%" />
  <img src="data/DE_figures_python/smear_day14_SL3_F-M.png" width="32%" />
  <img src="data/DE_figures_python/Venn_day14_sex_bias_by_line.png" width="25%" />
</p>

</div>

<div id="day14_FM_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_day14_SL1_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_day14_SL3_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/Venn_day14_sex_bias_by_line.png" width="25%" />
</p>

</div>

<div id="day16_FM_smear" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_day16_SL1_F-M.png" width="32%" />
  <img src="data/DE_figures_python/smear_day16_SL3_F-M.png" width="32%" />
  <img src="data/DE_figures_python/Venn_day16_sex_bias_by_line.png" width="25%" />
</p>

</div>

<div id="day16_FM_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_day16_SL1_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_day16_SL3_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/Venn_day16_sex_bias_by_line.png" width="25%" />
</p>

</div>

<div id="day18_FM_smear" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_day18_SL1_F-M.png" width="32%" />
  <img src="data/DE_figures_python/smear_day18_SL3_F-M.png" width="32%" />
  <img src="data/DE_figures_python/Venn_day16_sex_bias_by_line.png" width="25%" />
</p>

</div>

<div id="day18_FM_volcano" class="tabcontent">

<p float="left">
  <img src="data/DE_figures_python/smear_day18_SL1_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/smear_day18_SL3_F-M_volcano.png" width="32%" />
  <img src="data/DE_figures_python/Venn_day16_sex_bias_by_line.png" width="25%" />
</p>

</div>



#### Interaction sex-bias by line-bias

`(F_1 - M_1) - (F_3 - M_3)`: Female-biased in one line and male-biased in the other. Also plot a scatter of the genes significant in the interaction in the same scatter plot as above if they are significant.

<span style="color:#BD351E;"><b>QUESTION: </b></span>What woud the biological interpretation of this be? Does it still make sense even though `F1-F3` should be 0 in all cases?

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'day14_interaction_smear')">Day 14 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day14_interaction_volcano')">Day 14 (Volcano plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day16_interaction_smear')">Day 16 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day16_interaction_volcano')">Day 16 (Volcano plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day18_interaction_smear')">Day 18 (Smear plot)</button>
  <button class="tablinks" onclick="openTab(event, 'day18_interaction_volcano')">Day 18 (Volcano plot)</button>
</div>

<div id="day14_interaction_smear" class="tabcontent">
<p float="left">
  <img src="data/DE_figures_python/smear_day14_F-M_by_1-3.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_interaction_day14_line_bias.png" width="25%" />
  <img src="data/DE_figures_python/LFC_scatter_interaction_day14_sex_bias.png" width="25%" />  
</p>

<details>
<summary>gene IDs list</summary>

* Downregulated (160)

```
['gene-124877', 'gene-124865', 'gene-239506', 'gene-223773', 'gene-224079', 'gene-223758', 'gene-424750', 'gene-223491', 'gene-371957', 'gene-127740', 'gene-429522', 'gene-90190', 'gene-87554', 'gene-426056', 'gene-231854', 'gene-230270', 'gene-288660', 'gene-232048', 'gene-220055', 'gene-77070', 'gene-290138', 'gene-120832', 'gene-39545', 'gene-286732', 'gene-130096', 'gene-241328', 'gene-218301', 'gene-76019', 'gene-166391', 'gene-130081', 'gene-38977', 'gene-231604', 'gene-223318', 'gene-39692', 'gene-158197', 'gene-414353', 'gene-241268', 'gene-266875', 'gene-222737', 'gene-88032', 'gene-115312', 'gene-372264', 'gene-122220', 'gene-288356', 'gene-70385', 'gene-204669', 'gene-268146', 'gene-374892', 'gene-375015', 'gene-406510', 'gene-326329', 'gene-407355', 'gene-327787', 'gene-210277', 'gene-125867', 'gene-223419', 'gene-346228', 'gene-143368', 'gene-39770', 'gene-414838', 'gene-161821', 'gene-125638', 'gene-39680', 'gene-75424', 'gene-78574', 'gene-149137', 'gene-217258', 'gene-210286', 'gene-53501', 'gene-240397', 'gene-336688', 'gene-95461', 'gene-100036', 'gene-71572', 'gene-227137', 'gene-236155', 'gene-241262', 'gene-336088', 'gene-351784', 'gene-263126', 'gene-188501', 'gene-127607', 'gene-6017', 'gene-223599', 'gene-24347', 'gene-73885', 'gene-302994', 'gene-215608', 'gene-62772', 'gene-266887', 'gene-391847', 'gene-223554', 'gene-246615', 'gene-289849', 'gene-60429', 'gene-229506', 'gene-233901', 'gene-124680', 'gene-244661', 'gene-69775', 'gene-272401', 'gene-89057', 'gene-286545', 'gene-124766', 'gene-346500', 'gene-2467', 'gene-83083', 'gene-89798', 'gene-2286', 'gene-333750', 'gene-311581', 'gene-406519', 'gene-67440', 'gene-241193', 'gene-161848', 'gene-283197', 'gene-223791', 'gene-72046', 'gene-85294', 'g11517', 'gene-241506', 'gene-224028', 'gene-129852', 'gene-335551', 'gene-223088', 'gene-73754', 'gene-39533', 'gene-266297', 'gene-218529', 'gene-323550', 'gene-250391', 'gene-377275', 'gene-83830', 'gene-425122', 'gene-68813', 'gene-57335', 'gene-269058', 'gene-238320', 'gene-347681', 'gene-215430', 'gene-53269', 'gene-31110', 'gene-68612', 'gene-212020', 'gene-55557', 'gene-223285', 'gene-280630', 'gene-211847', 'gene-271655', 'gene-406603', 'gene-277078', 'gene-203445', 'gene-77109', 'gene-120964', 'gene-7268', 'gene-179860', 'gene-62891', 'gene-77835', 'gene-30328', 'gene-127773']
```

* Upregulated (140)

```
['gene-227308', 'gene-224369', 'gene-225158', 'gene-221012', 'g14784', 'gene-125210', 'gene-224890', 'gene-238849', 'gene-403878', 'gene-120660', 'gene-237378', 'gene-119161', 'gene-241682', 'gene-243630', 'gene-218813', 'gene-421265', 'gene-84577', 'gene-370643', 'gene-240602', 'gene-227370', 'gene-240860', 'gene-5731', 'gene-330102', 'gene-217099', 'gene-48714', 'gene-81750', 'gene-370487', 'gene-90428', 'gene-395158', 'gene-217443', 'gene-371230', 'gene-303243', 'gene-124377', 'gene-68558', 'gene-370842', 'gene-215103', 'gene-62927', 'gene-301479', 'gene-370323', 'gene-64491', 'gene-227284', 'gene-118620', 'gene-221469', 'gene-417051', 'gene-250269', 'gene-370562', 'gene-333008', 'gene-392350', 'gene-82722', 'gene-264249', 'gene-81581', 'gene-407735', 'gene-73288', 'gene-272192', 'gene-262796', 'gene-81518', 'gene-226254', 'gene-57617', 'gene-408949', 'gene-120070', 'gene-81458', 'gene-398993', 'gene-81476', 'gene-334677', 'gene-421549', 'gene-370595', 'gene-425532', 'gene-241841', 'gene-81675', 'gene-118985', 'gene-81418', 'gene-65889', 'gene-88458', 'gene-240638', 'g11957', 'gene-80415', 'gene-87700', 'gene-225355', 'gene-76284', 'gene-401486', 'gene-227071', 'gene-399223', 'gene-263313', 'gene-313589', 'gene-288413', 'gene-166814', 'gene-237857', 'gene-65163', 'gene-87487', 'gene-229515', 'gene-219157', 'gene-350813', 'gene-234256', 'gene-64470', 'gene-57689', 'gene-263588', 'gene-242691', 'gene-294364', 'gene-400426', 'gene-222159', 'gene-268137', 'gene-88092', 'gene-324340', 'gene-262046', 'gene-90524', 'g5814', 'gene-410993', 'gene-399203', 'gene-215313', 'gene-358397', 'gene-361221', 'gene-417093', 'gene-261982', 'gene-378262', 'gene-131338', 'gene-77588', 'gene-152204', 'gene-26634', 'gene-58495', 'gene-221504', 'g1591', 'gene-392536', 'gene-407714', 'gene-370967', 'gene-266345', 'gene-153689', 'gene-125849', 'gene-5187', 'gene-331395', 'gene-55887', 'gene-130632', 'gene-73160', 'gene-263597', 'gene-331896', 'gene-271117', 'gene-59093', 'gene-431030', 'gene-397251', 'gene-283696', 'gene-119684']
```

</details>

</div>

<div id="day14_interaction_volcano" class="tabcontent">
<p float="left">
  <img src="data/DE_figures_python/smear_day14_F-M_by_1-3_volcano.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_interaction_day14_line_bias.png" width="25%" />
  <img src="data/DE_figures_python/LFC_scatter_interaction_day14_sex_bias.png" width="25%" />
</p>
</div>

<div id="day16_interaction_smear" class="tabcontent">
<p float="left">
  <img src="data/DE_figures_python/smear_day16_F-M_by_1-3.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_interaction_day16_line_bias.png" width="25%" />
  <img src="data/DE_figures_python/LFC_scatter_interaction_day16_sex_bias.png" width="25%" />  
</p>

<details>
<summary>gene IDs list</summary>

 * Downregulated (173)

```
['gene-424768', 'gene-263126', 'gene-240397', 'gene-424750', 'gene-211137', 'gene-243753', 'gene-215596', 'gene-63245', 'gene-223722', 'gene-428089', 'gene-54363', 'gene-424726', 'gene-87554', 'gene-223758', 'gene-55869', 'gene-90190', 'gene-223419', 'gene-223734', 'gene-367523', 'gene-42302', 'gene-113411', 'gene-14276', 'gene-217473', 'gene-7220', 'gene-63617', 'gene-424863', 'gene-229506', 'gene-6223', 'gene-34632', 'gene-119402', 'gene-312890', 'gene-69401', 'gene-241506', 'gene-23181', 'gene-127607', 'gene-336178', 'gene-8365', 'gene-38885', 'gene-328746', 'gene-388769', 'gene-86738', 'gene-14252', 'gene-231493', 'gene-125638', 'gene-393138', 'gene-279676', 'gene-221980', 'gene-407280', 'gene-68612', 'gene-152989', 'gene-104371', 'gene-234686', 'gene-271655', 'gene-336688', 'gene-241856', 'gene-306977', 'gene-127740', 'gene-233883', 'gene-57335', 'gene-38977', 'gene-85401', 'gene-353013', 'gene-269365', 'gene-346228', 'gene-201605', 'gene-127849', 'gene-205849', 'gene-152165', 'gene-431362', 'gene-313040', 'gene-174594', 'gene-294618', 'gene-128145', 'gene-228118', 'gene-88032', 'gene-32436', 'gene-424896', 'gene-124877', 'gene-231604', 'gene-420308', 'gene-149557', 'gene-327074', 'gene-426056', 'gene-223318', 'gene-283260', 'gene-223512', 'gene-230270', 'gene-80062', 'gene-279912', 'gene-231540', 'gene-333603', 'gene-116902', 'gene-410561', 'gene-259467', 'gene-183665', 'gene-240184', 'gene-367095', 'gene-424914', 'gene-195335', 'gene-345135', 'gene-351334', 'gene-254128', 'gene-222737', 'gene-277218', 'gene-23163', 'gene-23042', 'gene-48598', 'gene-68000', 'gene-185170', 'gene-402875', 'gene-197114', 'gene-421566', 'gene-233485', 'gene-268996', 'gene-272072', 'gene-30322', 'gene-55240', 'gene-39933', 'gene-405355', 'gene-402536', 'gene-241238', 'gene-310012', 'gene-333433', 'gene-282008', 'gene-255088', 'gene-324223', 'gene-17262', 'gene-30595', 'gene-90918', 'gene-279975', 'gene-321078', 'gene-327616', 'gene-254639', 'gene-372264', 'gene-101822', 'gene-266887', 'gene-174585', 'gene-246732', 'gene-286289', 'gene-21635', 'gene-287310', 'gene-275750', 'gene-362311', 'gene-64360', 'gene-425173', 'gene-232048', 'gene-196966', 'gene-428194', 'gene-89219', 'gene-268975', 'gene-55252', 'gene-183393', 'gene-182960', 'gene-329366', 'gene-6199', 'gene-269228', 'gene-222746', 'gene-336348', 'gene-289528', 'gene-14156', 'gene-369289', 'gene-228465', 'gene-117360', 'gene-16790', 'gene-218421', 'gene-188158', 'gene-286744', 'gene-270896', 'gene-428071', 'gene-360503', 'gene-330422', 'gene-58752', 'gene-220788']
```
	 
  * Upregulated (108)

```
['gene-240602', 'gene-24185', 'gene-24290', 'gene-23840', 'gene-23597', 'gene-24203', 'gene-24120', 'gene-24132', 'gene-23538', 'gene-423321', 'gene-87700', 'gene-24088', 'gene-24221', 'gene-390687', 'gene-24278', 'gene-24167', 'gene-13404', 'gene-23514', 'gene-407253', 'gene-326882', 'gene-23365', 'gene-15763', 'gene-23689', 'gene-428729', 'gene-23413', 'gene-23834', 'gene-327441', 'gene-23893', 'gene-392224', 'gene-23884', 'gene-240935', 'gene-80466', 'gene-421265', 'gene-120660', 'gene-428774', 'gene-24052', 'gene-90307', 'gene-24079', 'gene-122692', 'gene-417051', 'gene-282611', 'gene-403700', 'gene-224250', 'gene-282641', 'gene-328764', 'gene-27466', 'gene-370842', 'gene-84577', 'gene-219019', 'gene-410366', 'gene-400426', 'gene-214979', 'gene-403652', 'gene-227164', 'gene-326810', 'gene-130665', 'gene-9713', 'gene-17601', 'gene-12075', 'gene-218086', 'gene-84224', 'gene-237372', 'gene-219157', 'gene-326873', 'gene-329410', 'gene-403706', 'gene-81551', 'gene-245743', 'gene-81599', 'gene-411056', 'gene-127707', 'gene-395158', 'gene-431788', 'gene-217099', 'gene-80484', 'gene-403851', 'gene-23911', 'gene-370550', 'gene-81675', 'g2779', 'gene-215563', 'gene-392159', 'gene-224369', 'gene-90503', 'gene-24914', 'gene-216914', 'gene-81418', 'gene-370643', 'gene-15918', 'gene-237351', 'gene-81476', 'gene-81581', 'gene-285127', 'gene-222242', 'gene-125849', 'gene-349163', 'gene-227137', 'gene-282437', 'gene-13124', 'gene-17046', 'gene-10767', 'gene-225355', 'gene-237857', 'gene-119684', 'gene-16667', 'gene-81458', 'gene-80415', 'gene-326825']
```

</details>

</div>

<div id="day16_interaction_volcano" class="tabcontent">
<p float="left">
  <img src="data/DE_figures_python/smear_day16_F-M_by_1-3_volcano.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_interaction_day16_line_bias.png" width="25%" />
  <img src="data/DE_figures_python/LFC_scatter_interaction_day16_sex_bias.png" width="25%" />
</p>
</div>

<div id="day18_interaction_smear" class="tabcontent">
<p float="left">
  <img src="data/DE_figures_python/smear_day18_F-M_by_1-3.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_interaction_day18_line_bias.png" width="25%" />
  <img src="data/DE_figures_python/LFC_scatter_interaction_day18_sex_bias.png" width="25%" />  
</p>

<details>
<summary>gene IDs list</summary>

* Downregulated (3)

```
['gene-428071', 'gene-233410', 'gene-428104']
```

* Upregulated (11)

```
['gene-396259', 'gene-395158', 'gene-224875', 'gene-395143', 'gene-90157', 'gene-224860', 'gene-428738', 'gene-395080', 'gene-224697', 'gene-301479', 'gene-242595']
```

</details>

</div>

<div id="day18_interaction_volcano" class="tabcontent">
<p float="left">
  <img src="data/DE_figures_python/smear_day18_F-M_by_1-3_volcano.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_interaction_day18_line_bias.png" width="25%" />
  <img src="data/DE_figures_python/LFC_scatter_interaction_day18_sex_bias.png" width="25%" />
</p>
</div>




#### Overlap of DE scatterplots


<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'day14_scatter')">Day 14</button>
  <button class="tablinks" onclick="openTab(event, 'day16_scatter')">Day 16</button>
  <button class="tablinks" onclick="openTab(event, 'day18_scatter')">Day 18</button>
</div>

<div id="day14_scatter" class="tabcontent">
<p float="left">
  <img src="data/DE_figures_python/LFC_scatter_day14_line_bias_by_sex.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_day14_sex_bias_by_line.png" width="32%" />
</p>
</div>

<div id="day16_scatter" class="tabcontent">
<p float="left">
  <img src="data/DE_figures_python/LFC_scatter_day16_line_bias_by_sex.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_day16_sex_bias_by_line.png" width="32%" />
</p>
</div>

<div id="day18_scatter" class="tabcontent">
<p float="left">
  <img src="data/DE_figures_python/LFC_scatter_day18_line_bias_by_sex.png" width="32%" />
  <img src="data/DE_figures_python/LFC_scatter_day18_sex_bias_by_line.png" width="32%" />
</p>
</div>