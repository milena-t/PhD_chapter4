"""
basic plotting functions for the RNAseq analysis
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats


def plot_counts(counts_table:str, geneIDs_list, outfile_name:str, y_label="normalized counts", mean_per_sample=False, remove_females = False):
    """
    plot the counts of the geneIDs in all samples
    except samples where all samples have zero counts
    if mean_per_sample, plot the mean counts for every sample with standard error
    if geneIDs_list is a dict with {geneID : legend_name} then the legend will be printed with the proper legend name and not the geneID 
    """
    assert len(geneIDs_list)>0
    counts_df = pd.read_csv(counts_table, sep="\t", index_col=0)
    # counts_present = counts_df.filter(items = geneIDs_list, axis = 1)
    headers = counts_df.columns.tolist()
    gene_counts = {}
    nonexpressed = []
    if type(geneIDs_list) == dict:
        usedict_labels = True
        geneIDs_dict = geneIDs_list
        geneIDs_list = list(geneIDs_list.keys())
    else:
        usedict_labels = False

    for geneID in geneIDs_list:
        try:
            counts_dict = counts_df.loc[geneID].to_dict()
        except:
            nonexpressed.append(geneID)
            continue
        gene_counts[geneID] = counts_dict
    print(f"{len(nonexpressed)} out of {len(geneIDs_list)} genes not expressed")

    nonzero_samples = []
    for sample in headers:
        try:
            sample_counts = [gene_counts[geneID][sample] for geneID in gene_counts.keys()]
        except: 
            raise RuntimeError(f"{sample} not found in {gene_counts[geneIDs_list[0]].keys()}")
        if sum(sample_counts) > 0:
            nonzero_samples.append(sample)
    print(f"out of {len(headers)} there are {len(nonzero_samples)} samples that have at least one count in one gene!")# \n{nonzero_samples}")
    if remove_females:
        nonzero_samples = [sample for sample in nonzero_samples if "-F" not in sample]
        print(f"remove female samples! {len(nonzero_samples)} samples left.")

    ### plotting
    fig, ax = plt.subplots(1,1, figsize=(20, 10)) # for more than three rows

    fs = 25
    ps = fs*15 # point size
    lw=2
    linest = ":"

    # sort nonzero samples
    females = [sample for sample in nonzero_samples if "-F_" in sample]
    f1 = [sample for sample in females if "-1-" in sample]
    f3 = [sample for sample in females if "-3-" in sample]
    males = [sample for sample in nonzero_samples if "-M_" in sample]
    m1 = [sample for sample in males if "-1-" in sample]
    m3 = [sample for sample in males if "-3-" in sample]
    nonzero_samples_sorted = f1+f3+m1+m3
    assert len(nonzero_samples) == len(nonzero_samples_sorted)
    
    if mean_per_sample:

        ## make medians and standard errors
        medians_dict = [0.0 for sample in nonzero_samples_sorted]
        errors_dict = [0.0 for sample in nonzero_samples_sorted]
        tick_labels = ["" for sample in nonzero_samples_sorted]
        tick_pos = [i for i, sample in enumerate(nonzero_samples_sorted)]
        
        for i, sample in enumerate(nonzero_samples_sorted):
            curr_counts = []
            count_expressed = 0
            for geneID,sample_counts in gene_counts.items():
                if sample_counts[sample] >0:
                    curr_counts.append(sample_counts[sample])
                    count_expressed+=1
            
            if len(curr_counts)>0:
                medians_dict[i] = np.median(curr_counts)
                errors_dict[i] = stats.sem(curr_counts)
            else:
                medians_dict[i] = np.nan
                errors_dict[i] = np.nan
            sample_ = sample.replace("WJ-3841-","").split("_")[0]
            tick_labels[i] = f"SL{sample_} ({count_expressed})"


        ax.errorbar(tick_pos, medians_dict, xerr = 0, yerr = errors_dict, color="#683257", linewidth =lw,
                    marker = ".", markersize=20, linestyle = linest)
        if "log" not in y_label:
            ymin, ymax = ax.get_ylim()
            ax.set_ylim([-0.1,ymax])
        tick_cols = ["#000000" if "M" in sample else "#8A8A8A" for sample in nonzero_samples_sorted ]
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_labels)

    
    else:
        
        if usedict_labels:
            for geneID,sample_counts in gene_counts.items():
                y_vec = [sample_counts[sample] for sample in nonzero_samples_sorted]
                ax.plot(nonzero_samples_sorted, y_vec,"-o", label = geneIDs_dict[geneID])
        else:
            for geneID,sample_counts in gene_counts.items():
                y_vec = [sample_counts[sample] for sample in nonzero_samples_sorted]
                ax.plot(nonzero_samples_sorted, y_vec,"-o",label=geneID)
    
        ax.set_xticklabels([sample.replace("WJ-3841-","").split("_")[0] for sample in nonzero_samples_sorted])
        plt.legend(fontsize=fs)

        tick_cols = ["#000000" if "M" in sample else "#8A8A8A" for sample in nonzero_samples_sorted ]
        for tick_label, color in zip(ax.get_xticklabels(), tick_cols):
            tick_label.set_color(color)
    
    ax.tick_params(axis='x', labelsize=fs*0.75,labelrotation=90)#, colors)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_ylabel(f"{y_label}", fontsize = fs)
    plt.tight_layout()
    plt.savefig(outfile_name, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {outfile_name}")


if __name__ == "__main__":
    
    username = "miltr339"
    counts_file = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/gene_counts_normalized_nolog.tsv"

    if True:
        yTor_IDs = {"yTor-all":"yTor","gene-30110":"aTor"}
        yTor_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/merged_yTor_aTor_counts.png"
        plot_counts(counts_table=counts_file, geneIDs_list=yTor_IDs, outfile_name=yTor_plot, remove_females=True)

        yTor_IDs = ["yTor-all"]
        yTor_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/merged_yTor_counts.png"
        plot_counts(counts_table=counts_file, geneIDs_list=yTor_IDs, outfile_name=yTor_plot, remove_females=True)
        if False:
            ## split yTor lists
            # all TOR
            yTor_IDs = ["yTor-A", "yTor-B", "yTor-C","gene-30110"]
            yTor_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/all_Tor_counts.png"
            plot_counts(counts_table=counts_file, geneIDs_list=yTor_IDs, outfile_name=yTor_plot)
            # only yTOR
            yTor_IDs = ["yTor-A", "yTor-B", "yTor-C"]
            yTor_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/yTor_counts.png"
            plot_counts(counts_table=counts_file, geneIDs_list=yTor_IDs, outfile_name=yTor_plot)
    if False:
        MSL2_IDs = {"gene-371922" : "Y-copy","gene-343165": "A-copy"} # Y-copy,A-copy
        MSL2_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/MSL2_counts.png"
        plot_counts(counts_table=counts_file, geneIDs_list=MSL2_IDs, outfile_name=MSL2_plot)

    if False:
        # all Y expressed
        y_expr = ["gene-371805","yTor-A","gene-371844","yTor-C","gene-371913","gene-371922","gene-371957","gene-372053","gene-372068","gene-372216","gene-372264"]
        y_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/y_expr_counts.png"
        plot_counts(counts_table=counts_file, geneIDs_list=y_expr, outfile_name=y_plot)
        y_expr = ["gene-371805","yTor-A","gene-371844","yTor-C","gene-371913","gene-371922","gene-371957","gene-372053","gene-372216"]
        y_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/y_low_expr_counts.png"
        plot_counts(counts_table=counts_file, geneIDs_list=y_expr, outfile_name=y_plot)