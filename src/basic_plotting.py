"""
basic plotting functions for the RNAseq analysis
"""

import matplotlib.pyplot as plt
import pandas as pd


def plot_counts(counts_table:str, geneIDs_list:list, outfile_name:str):
    """
    plot the counts of the geneIDs in all samples
    except samples where all samples have zero counts
    """
    assert len(geneIDs_list)>0
    counts_df = pd.read_csv(counts_table, sep="\t", index_col=0)
    # counts_present = counts_df.filter(items = geneIDs_list, axis = 1)
    headers = counts_df.columns.tolist()
    gene_counts = {}
    for geneID in geneIDs_list:
        try:
            counts_dict = counts_df.loc[geneID].to_dict()
        except:
            print(f"{geneID} not expressed!")
            continue
        gene_counts[geneID] = counts_dict

    nonzero_samples = []

    for sample in headers:
        try:
            sample_counts = [gene_counts[geneID][sample] for geneID in gene_counts.keys()]
        except: 
            raise RuntimeError(f"{sample} not found in {gene_counts[geneIDs_list[0]].keys()}")
        if sum(sample_counts) > 0:
            nonzero_samples.append(sample)
    print(f"out of {len(headers)} there are {len(nonzero_samples)} samples that have at least one count in one gene! \n{nonzero_samples}")

    ### plotting
    fig, ax = plt.subplots(1,1, figsize=(15, 10)) # for more than three rows

    fs = 25
    ps = fs*15 # point size
    for geneID,sample_counts in gene_counts.items():
        y_vec = [sample_counts[sample] for sample in nonzero_samples]
        ax.plot(nonzero_samples, y_vec,"-o",label=geneID)
    
    ax.set_xticklabels([sample.replace("WJ-3841-","").split("_")[0] for sample in nonzero_samples])
    tick_cols = ["#000000" if "M" in sample else "#8A8A8A" for sample in nonzero_samples ]
    ax.tick_params(axis='x', labelsize=fs*0.75,labelrotation=90)#, colors)
    for tick_label, color in zip(ax.get_xticklabels(), tick_cols):
        tick_label.set_color(color)
    ax.tick_params(axis='y', labelsize=fs)
    
    ax.set_ylabel(f"normalized counts", fontsize = fs)
    plt.legend(fontsize=fs)
    plt.tight_layout()
    plt.savefig(outfile_name, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as: {outfile_name}")


if __name__ == "__main__":
    
    username = "miltr339"
    counts_file = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/gene_counts_normalized_nolog.tsv"

    yTor_IDs = ["yTor-A", "yTor-B", "yTor-C"]
    yTor_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_counts.png"
    plot_counts(counts_table=counts_file, geneIDs_list=yTor_IDs, outfile_name=yTor_plot)