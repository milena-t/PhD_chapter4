"""
basic plotting functions for the RNAseq analysis
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from matplotlib.ticker import FuncFormatter


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
    plt.savefig(outfile_name, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {outfile_name}")


def plot_exon_intron_structure(plot_filename, mTor_gff, yL_gff, yS_gff):
    TORs = ["yTor-A", "yTor-B", "yTor-C"]
    yTor_A_start = []
    yTor_A_end = []
    yTor_A_categories = []
    yTor_B_start = []
    yTor_B_end = []
    yTor_B_categories = []
    yTor_C_start = []
    yTor_C_end = []
    yTor_C_categories = []
    categories = ["exon", "cds", "intron"]
    category_colors = {categories[0]: "mediumseagreen",
                    categories[2]: "indianred",
                    categories[1]: "cornflowerblue"}

    genes_range = [] # 0:1 = mTor; 2:3 = A ; 4:5 = B ; 6:7 = C

    mTor_start = []
    mTor_end = []
    mTor_categories = []

    ## read intron/exon/cds locations from the gff files for all TOR copies

    with open(mTor_gff, "r") as annotation_file:
        lines = annotation_file.readlines()
        for line in lines:
            line = line.split("\t")
            if line[2] in categories:
                if "mTor" in line[-1]:
                    mTor_start.append(int(line[3]))
                    mTor_end.append(int(line[4]))
                    mTor_categories.append(line[2])     
            if line[2] == "gene":
                genes_range.append(int(line[3]))
                genes_range.append(int(line[4]))
                pass
            else:
                pass

    with open(yS_gff, "r") as annotation_file, open(yL_gff, "w") as torc_outfile:
        torc_outfile.write("type,beginning,end\n")
        lines = annotation_file.readlines()
        for line in lines:
            line = line.split("\t")
            if line[2] in categories:
                if "yTor-A" in line[-1]:
                    yTor_A_start.append(int(line[3]))
                    yTor_A_end.append(int(line[4]))
                    yTor_A_categories.append(line[2])     
                elif "yTor-B" in line[-1]:
                    yTor_B_start.append(int(line[3]))
                    yTor_B_end.append(int(line[4]))
                    yTor_B_categories.append(line[2])
                elif "yTor-C" in line[-1]:
                    yTor_C_start.append(int(line[3]))
                    yTor_C_end.append(int(line[4]))
                    yTor_C_categories.append(line[2])
                    if line[2] not in " cds ":
                        torc_outfile.write(",".join(line[2:5])+"\n")
            if line[2] == "gene":
                genes_range.append(int(line[3]))
                genes_range.append(int(line[4]))
            else:
                pass

    ## PLOTTING ##
    if True:

        plt.rcParams['text.usetex'] = True

        lw = 1
        offset_factor = 0.1 # moving the exon/intron/coding sections on top of each other for the TOR plots
        y_factor_A = genes_range[2] - genes_range[0] # transform all TOR copies onto the same bp numbers as mTor
        y_factor_B = genes_range[4] - genes_range[0]
        y_factor_C = genes_range[6] - genes_range[0]

        fig, ax = plt.subplots(1,1, figsize=(15, 5))
        # plt.figure(figsize=(15, 5))
        fs=20

        print("\naTor")
        print(len(mTor_start), len(mTor_end))
        for i in range(len(mTor_start)):
            if mTor_categories[i] == "exon":
                pass
            else:
                y_ax_offset = (categories.index(mTor_categories[i])+15)*offset_factor # +15 is the top gene in the plot, the other ones are at 10, 5 and 0
                plt.fill([mTor_start[i], mTor_end[i], mTor_end[i],  mTor_start[i]], [y_ax_offset+offset_factor/2, y_ax_offset+offset_factor/2, y_ax_offset-offset_factor/2, y_ax_offset-offset_factor/2], color = category_colors[mTor_categories[i]], linewidth = lw, alpha = 0.5)


        print("\nyTor-A")
        print(len(yTor_A_start), len(yTor_A_end))
        for i in range(len(yTor_A_start)):
            if mTor_categories[i] == "exon":
                pass
            else:
                y_ax_offset = (categories.index(yTor_A_categories[i])+5)*offset_factor
                plt.fill([yTor_A_start[i]-y_factor_A, yTor_A_end[i]-y_factor_A, yTor_A_end[i]-y_factor_A,  yTor_A_start[i]-y_factor_A], [y_ax_offset+offset_factor/2, y_ax_offset+offset_factor/2, y_ax_offset-offset_factor/2, y_ax_offset-offset_factor/2], color = category_colors[yTor_A_categories[i]], linewidth = lw, alpha = 0.5)


        print("yTor-B")
        print(len(yTor_B_start), len(yTor_B_end))
        for i in range(len(yTor_B_start)):
            if mTor_categories[i] == "exon":
                pass
            else:
                y_ax_offset = (categories.index(yTor_B_categories[i])+0)*offset_factor
                plt.fill([yTor_B_start[i]-y_factor_B, yTor_B_end[i]-y_factor_B, yTor_B_end[i]-y_factor_B,  yTor_B_start[i]-y_factor_B], [y_ax_offset+offset_factor/2, y_ax_offset+offset_factor/2, y_ax_offset-offset_factor/2, y_ax_offset-offset_factor/2], color = category_colors[yTor_A_categories[i]], linewidth = lw, alpha = 0.5)

        print("yTor-C")
        print(len(yTor_C_start), len(yTor_C_end))
        for i in range(len(yTor_C_start)):
            if mTor_categories[i] == "exon":
                pass
            else:
                y_ax_offset = (categories.index(yTor_C_categories[i])+10)*offset_factor
                plt.fill([yTor_C_start[i]-y_factor_C, yTor_C_end[i]-y_factor_C, yTor_C_end[i]-y_factor_C,  yTor_C_start[i]-y_factor_C], [y_ax_offset+offset_factor/2, y_ax_offset+offset_factor/2, y_ax_offset-offset_factor/2, y_ax_offset-offset_factor/2], color = category_colors[yTor_A_categories[i]], linewidth = lw, alpha = 0.5)



        for i in range(len(genes_range)):
            plt.plot([genes_range[i], genes_range[i]],[-0.5, 2.5], color = "silver", linestyle = "--")


        all_coordinates = yTor_A_start+yTor_A_end+yTor_B_start+yTor_B_end+yTor_C_start+yTor_C_end
        #print(len(all_coordinates), str(59*6))

        # adjust plot labes horizontal or veritcal (0 = horizontal, 90 = vertical)
        rotate_factor = 90

        plt.axis([min(genes_range)-min(genes_range)*0.0005, genes_range[1]+genes_range[1]*0.001, -0.2, 2])
        # plt.xlabel("length in $10^3$ bp", rotation = 0, fontsize = fs)

        # custom y ticks to label the tor copies
        ticks = [0.1, 0.6, 1.1, 1.6]
        tick_names = [ "y-TOR B", "y-TOR A","y-TOR C", "a-Tor"]
        plt.yticks(ticks, tick_names)#, rotation = rotate_factor)
        ax.tick_params(axis='y', labelsize=fs)

        # custom x ticks 
        tick_intervals = int(((genes_range[1] - genes_range[0])/6)/1000)*1000
        print("tick intervals: ", tick_intervals)
        ticks = [genes_range[0], genes_range[0]+tick_intervals, genes_range[0]+2*tick_intervals, genes_range[0]+3*tick_intervals, genes_range[0]+4*tick_intervals, genes_range[0]+5*tick_intervals, genes_range[1]]
        tick_names = [0, tick_intervals, 2*tick_intervals, 3*tick_intervals, 4*tick_intervals, 5*tick_intervals, genes_range[1]-genes_range[0]]
        tick_names = [f"{i/1000:.0f}kb" for i in tick_names]
        plt.xticks(ticks, tick_names, rotation = 0, fontsize = fs)
        # ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x < 1 else f'{x / 1e3:.0f} kb'))
        ax.tick_params(axis='x', labelsize=fs)#, rotation=rotate_factor)#, colors)
        
        plt.title("Intron/exon structure of Y-Tor copies", fontsize = fs*1.25)

        plt.legend(labelcolor = [ "cornflowerblue", "indianred"], loc = "lower right",  labels = [ "exon", "intron"], fontsize = fs)
        #plt.legend(labelcolor = [ "cornflowerblue", "mediumseagreen",  "indianred"], loc = "lower right",  labels = [ "cds", "exon", "intron"], fontsize = "large")

        plt.tight_layout()
    
        plt.savefig(plot_filename, format='png', dpi=300, transparent = True) 


if __name__ == "__main__":
    
    username = "milena"
    counts_file = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/gene_counts_normalized_nolog.tsv"

    if False:
        yTor_IDs = {"yTor-all":"yTor","gene-30110":"aTor"}
        yTor_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/merged_yTor_aTor_counts.png"
        plot_counts(counts_table=counts_file, geneIDs_list=yTor_IDs, outfile_name=yTor_plot, remove_females=True)

        yTor_IDs = ["yTor-all"]
        yTor_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/merged_yTor_counts.png"
        plot_counts(counts_table=counts_file, geneIDs_list=yTor_IDs, outfile_name=yTor_plot, remove_females=True)
        if False:
            ## split yTor lists -> HAS TO USE OLD FILE WITH ACTUAL SPLIT GENEIDS!
            ## since I generate this with edgeR from the raw counts, the current file is the merged one used for all other analysis.
            ## to re-run these plots I need to re-create the old file by running the split yTor raw counts one through edgeR again

            # all TOR
            yTor_IDs = ["yTor-A", "yTor-B", "yTor-C","gene-30110"]
            yTor_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/all_Tor_counts.png"
            plot_counts(counts_table=counts_file, geneIDs_list=yTor_IDs, outfile_name=yTor_plot)
            # only yTOR
            yTor_IDs = ["yTor-A", "yTor-B", "yTor-C"]
            yTor_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/yTor_counts.png"
            plot_counts(counts_table=counts_file, geneIDs_list=yTor_IDs, outfile_name=yTor_plot)
    if False:
        MSL2_IDs = {"gene-371922" : "Y-MSL2","gene-343165": "A-MSL2"} # Y-copy,A-copy
        MSL2_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/MSL2_counts.png"
        plot_counts(counts_table=counts_file, geneIDs_list=MSL2_IDs, outfile_name=MSL2_plot)

    if False:
        # all Y expressed
        y_expr = ["gene-371805","gene-371844","gene-371889","gene-371913","gene-371922","gene-371957","gene-372053","gene-372068","gene-372216","gene-372264","yTor-all"]
        y_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/y_expr_counts.png"
        plot_counts(counts_table=counts_file, geneIDs_list=y_expr, outfile_name=y_plot)
        # remove highly expressed genes
        y_expr = ["gene-371805","gene-371844","gene-371889","gene-371913","gene-371922","gene-371957","gene-372053","gene-372216","yTor-all"]
        y_plot = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/y_low_expr_counts.png"
        plot_counts(counts_table=counts_file, geneIDs_list=y_expr, outfile_name=y_plot)

    if True:
        ## plot exon intron structure
        plot_exon_intron_structure(plot_filename = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/y_s_tor_copies.png",
                                   mTor_gff=f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/mTOR_annotation.gff",
                                   yL_gff=f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/y_l_TorC_exon_coordinates.csv",
                                   yS_gff=f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis/Y_L_TORs_annotation_with_introns.gff")