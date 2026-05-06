"""
reproduce plots from R since there are contradictions. 
I will use the tables in data made with the topTags() function with no filtering so it just shows all expressed genes in the comparison
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3
import math

def get_tables(username = "miltr339"):
    """
    tables are either split by sex, so that the line and day contrasts are made on a subset that is only males or only females,
    or split by line so that the day and sex contrasts are made on a subset of only one line at a time.
    """

    tables_dir = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/"

    out_dict = {
        
        "sex_separated" : {
            "females" :{
                "SL1_14 - SL3_14" : f"{tables_dir}DE_genes_F_1-3_day14.txt",
                "SL1_16 - SL3_16" : f"{tables_dir}DE_genes_F_1-3_day16.txt",
                "SL1_18 - SL3_18" : f"{tables_dir}DE_genes_F_1-3_day18.txt",
                "SL1_18 - (SL1_14+SL1_16)/2" : f"{tables_dir}DE_genes_F_SL1.txt", 
                "SL3_18 - (SL3_14+SL3_16)/2" : f"{tables_dir}DE_genes_F_SL3.txt", 
                "SL1_18 - SL1_14" : f"{tables_dir}DE_genes_F_SL1_18_14.txt",
                "SL1_18 - SL1_16" : f"{tables_dir}DE_genes_F_SL1_18_16.txt",
                "SL1_14 - SL1_16" : f"{tables_dir}DE_genes_F_SL1_14_16.txt",
                "SL3_18 - SL3_14" : f"{tables_dir}DE_genes_F_SL3_18_14.txt",
                "SL3_18 - SL3_16" : f"{tables_dir}DE_genes_F_SL3_18_16.txt",
                "SL3_14 - SL3_16" : f"{tables_dir}DE_genes_F_SL3_14_16.txt",
                "(SL1_18 - SL1_14) - (SL3_18 - SL3_14)" : f"{tables_dir}DE_genes_F_SL_1_3_18_14.txt",
                "(SL1_14 - SL1_16) - (SL3_14 - SL3_16)" : f"{tables_dir}DE_genes_F_SL_1_3_14_16.txt",
                "(SL1_18 - SL1_16) - (SL3_18 - SL3_16)" : f"{tables_dir}DE_genes_F_SL_1_3_18_16.txt",
            },
            "males" : {
                "SL1_14 - SL3_14" : f"{tables_dir}DE_genes_M_1-3_day14.txt",
                "SL1_16 - SL3_16" : f"{tables_dir}DE_genes_M_1-3_day16.txt",
                "SL1_18 - SL3_18" : f"{tables_dir}DE_genes_M_1-3_day18.txt",
                "SL1_18 - (SL1_14+SL1_16)/2" : f"{tables_dir}DE_genes_M_SL1.txt",
                "SL3_18 - (SL3_14+SL3_16)/2" : f"{tables_dir}DE_genes_M_SL3.txt",
                "SL1_18 - SL1_14" : f"{tables_dir}DE_genes_M_SL1_18_14.txt",
                "SL1_18 - SL1_16" : f"{tables_dir}DE_genes_M_SL1_18_16.txt",
                "SL1_14 - SL1_16" : f"{tables_dir}DE_genes_M_SL1_14_16.txt",
                "SL3_18 - SL3_14" : f"{tables_dir}DE_genes_M_SL3_18_14.txt",
                "SL3_18 - SL3_16" : f"{tables_dir}DE_genes_M_SL3_18_16.txt",
                "SL3_14 - SL3_16" : f"{tables_dir}DE_genes_M_SL3_14_16.txt",
                "(SL1_18 - SL1_14) - (SL3_18 - SL3_14)" : f"{tables_dir}DE_genes_M_SL_1_3_18_14.txt",
                "(SL1_14 - SL1_16) - (SL3_14 - SL3_16)" : f"{tables_dir}DE_genes_M_SL_1_3_14_16.txt",
                "(SL1_18 - SL1_16) - (SL3_18 - SL3_16)" : f"{tables_dir}DE_genes_M_SL_1_3_18_16.txt",
            }
        },
        "line_separated" : {
            "SL1" : {
                "F_14 - M_14" : f"{tables_dir}DE_genes_SL1_day14_F-M.txt",
                "F_16 - M_16" : f"{tables_dir}DE_genes_SL1_day16_F-M.txt",
                "F_18 - M_18" : f"{tables_dir}DE_genes_SL1_day18_F-M.txt",
                "(F_14 - M_14) - (F_16 - M_16)" : f"{tables_dir}DE_genes_SL1_day14_16_F-M.txt",
                "(F_18 - M_18) - (F_14 - M_14)" : f"{tables_dir}DE_genes_SL1_day18_14_F-M.txt",
                "(F_18 - M_18) - (F_16 - M_16)" : f"{tables_dir}DE_genes_SL1_day18_16_F-M.txt",
            },
            "SL3" : {
                "F_14 - M_14" : f"{tables_dir}DE_genes_SL3_day14_F-M.txt",
                "F_16 - M_16" : f"{tables_dir}DE_genes_SL3_day16_F-M.txt",
                "F_18 - M_18" : f"{tables_dir}DE_genes_SL3_day18_F-M.txt",
                "(F_14 - M_14) - (F_16 - M_16)" : f"{tables_dir}DE_genes_SL3_day14_16_F-M.txt",
                "(F_18 - M_18) - (F_14 - M_14)" : f"{tables_dir}DE_genes_SL3_day18_14_F-M.txt",
                "(F_18 - M_18) - (F_16 - M_16)" : f"{tables_dir}DE_genes_SL3_day18_16_F-M.txt",
            }
        },
        "day_separated" : {
            "day14" : {
                "F_1 - M_1" : f"{tables_dir}DE_genes_day14_SL1_F-M.txt",
                "F_3 - M_3" : f"{tables_dir}DE_genes_day14_SL3_F-M.txt",
                "F_1 - F_3" : f"{tables_dir}DE_genes_day14_F_1-3.txt",
                "M_1 - M_3" : f"{tables_dir}DE_genes_day14_M_1-3.txt",
                "(F_1 - M_1) - (F_3 - M_3)" : f"{tables_dir}DE_genes_day14_F-M_by_1-3.txt",
            },
            "day16" : {
                "F_1 - M_1" : f"{tables_dir}DE_genes_day16_SL1_F-M.txt",
                "F_3 - M_3" : f"{tables_dir}DE_genes_day16_SL3_F-M.txt",
                "F_1 - F_3" : f"{tables_dir}DE_genes_day16_F_1-3.txt",
                "M_1 - M_3" : f"{tables_dir}DE_genes_day16_M_1-3.txt",
                "(F_1 - M_1) - (F_3 - M_3)" : f"{tables_dir}DE_genes_day16_F-M_by_1-3.txt",
            },
            "day18" : {
                "F_1 - M_1" : f"{tables_dir}DE_genes_day18_SL1_F-M.txt",
                "F_3 - M_3" : f"{tables_dir}DE_genes_day18_SL3_F-M.txt",
                "F_1 - F_3" : f"{tables_dir}DE_genes_day18_F_1-3.txt",
                "M_1 - M_3" : f"{tables_dir}DE_genes_day18_M_1-3.txt",
                "(F_1 - M_1) - (F_3 - M_3)" : f"{tables_dir}DE_genes_day18_F-M_by_1-3.txt",
            },
        }
    }

    contrast_plot_titles = {
        "SL1_14 - SL3_14" : f"day 14",
        "SL1_16 - SL3_16" : f"day 16",
        "SL1_18 - SL3_18" : f"day 18",
        "SL1_18 - (SL1_14+SL1_16)/2" : f"SL1",
        "SL3_18 - (SL3_14+SL3_16)/2" : f"SL3",
        "F_18 - M_18" : f"day 18",
        "F_14 - M_14" : f"day 14",
        "F_16 - M_16" : f"day 16",
        "SL1_18 - SL1_14" : f"SL1",
        "SL1_18 - SL1_16" : f"SL1",
        "SL1_14 - SL1_16" : f"SL1",
        "SL3_18 - SL3_14" : f"SL3",
        "SL3_18 - SL3_16" : f"SL3",
        "SL3_14 - SL3_16" : f"SL3",
        "(SL1_18 - SL1_14) - (SL3_18 - SL3_14)" : f"SL1-SL3 by day 18-14",
        "(SL1_14 - SL1_16) - (SL3_14 - SL3_16)" : f"SL1-SL3 by day 14-16",
        "(SL1_18 - SL1_16) - (SL3_18 - SL3_16)" : f"SL1-SL3 by day 18-16",
        "(F_14 - M_14) - (F_16 - M_16)" : f"F-M by day 14-16",
        "(F_18 - M_18) - (F_14 - M_14)" : f"F-M by day 18-14",
        "(F_18 - M_18) - (F_16 - M_16)" : f"F-M by day 18-16",
        "F_1 - M_1" : f"SL1 F-M",
        "F_3 - M_3" : f"SL3 F-M",
        "F_1 - F_3" : f"females SL1-3",
        "M_1 - M_3" : f"males SL1-3",
        "(F_1 - M_1) - (F_3 - M_3)" : f"F-M by SL1-SL3",
    }
    return out_dict,contrast_plot_titles


def plot_smear(table_path, contrast, smear_plot_name = "smear_plot.png", p_sig = 0.05, min_LFC = 0, title = "significant logFC", excl_genes_list = [], x_axis = "logcpm"):
    """
    replicate smear-plot from R, logCPM by logFC, significance highlighted in red
    return the number of up/downregulated and no difference genes
    """
    df = pd.read_csv(table_path, sep="\t", skiprows=0)
    if len(excl_genes_list) >0:
        old_len = df.shape[0]
        df = df.drop(index=excl_genes_list, errors='ignore') # only drop existing labels, ignore the rest
        new_len = df.shape[0]
        print(f"\t{len(excl_genes_list)} genes dropped from excl_genes_list (gene number {old_len} -> {new_len})")

    cols = {"nonsig" : "#243742", "sig" : "#BD351E"}

    if  x_axis == "logcpm":
        fig, ax = plt.subplots(1,1, figsize=(18, 10)) 
    elif x_axis == "fdr_p":
        fig, ax = plt.subplots(1,1, figsize=(18, 12)) 
    fs = 35
    point_size_factor = 5
    ps = fs*point_size_factor # point size
    
    ### plot the unsignificant first so they are below and the significant ones are above
    df_nonsig = df.loc[df['FDR'] >= p_sig]
    df_sig = df.loc[df['FDR'] < p_sig]
    print(f"\tnum sig genes: {df_sig.shape[0]}")
    if min_LFC>0:
        df_sig = df_sig.loc[abs(df_sig['logFC']) >= min_LFC]
        print(f"\tnum sig genes with LFC > {min_LFC}: {df_sig.shape[0]}")

    if  x_axis == "logcpm":
        ax.scatter(df_nonsig["logCPM"], df_nonsig["logFC"], color = cols["nonsig"], alpha=0.5, s=ps)
        ax.scatter(df_sig["logCPM"], df_sig["logFC"], color = cols["sig"], alpha=1, s=ps)
    elif x_axis == "fdr_p":
        ax.scatter(df_nonsig["logFC"], df_nonsig["FDR"], color = cols["nonsig"], alpha=0.5, s=ps)
        ax.scatter(df_sig["logFC"], df_sig["FDR"], color = cols["sig"], alpha=1, s=ps)
    
    if "SL" in contrast and "/2" not in contrast:
        if "SL1" in contrast and "SL3" in contrast:
            label_contrast = contrast.replace("_14", "").replace("_16", "").replace("_18", "")
        else:
            label_contrast = contrast.replace("SL1_", "day ").replace("SL3_", "day ")    
    elif "F" in contrast and "(" not in contrast:
        label_contrast = contrast.replace("_14", "").replace("_16", "").replace("_18", "")
    else:
        label_contrast = contrast.replace("SL1", "day").replace("SL3", "day")
    if  x_axis == "logcpm":
        ax.set_ylabel(f"logFC ({label_contrast})", fontsize = fs)
        ax.set_xlabel(f"log CPM", fontsize = fs)
    elif x_axis == "fdr_p":
        ax.set_xlabel(f"logFC ({label_contrast})", fontsize = fs)
        ax.set_ylabel(f"FDR p-value", fontsize = fs)
        ax.yaxis.set_inverted(True) 
    ax.tick_params(axis='x', labelsize=fs*0.9)
    ax.tick_params(axis='y', labelsize=fs*0.9)
    ax.set_title(title, fontsize = fs*1.25)

    if  x_axis == "logcpm":
        min_line = min(df["logCPM"])-0.25
        max_line = max(df["logCPM"])+0.25
        if min_LFC > 0:
            ax.hlines(y=min_LFC, xmin=min_line, xmax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
            ax.hlines(y=-1*min_LFC, xmin=min_line, xmax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
        else:
            ax.hlines(y=1, xmin=min_line, xmax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
            ax.hlines(y=-1, xmin=min_line, xmax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
        ax.set_xlim([min_line,max_line])
    elif x_axis == "fdr_p":
        pass
        # min_line = min(df["logFC"])
        # max_line = max(df["logFC"])
        # ax.vlines(x=0, ymin=min_line, ymax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
        # ax.set_xlim([min_line,max_line])
    

    if x_axis == "fdr_p":
        smear_plot_name = smear_plot_name.replace(".png", "_volcano.png")
    plt.tight_layout()
    fig.subplots_adjust(left=0.125) # or whatever
    plt.savefig(smear_plot_name, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {smear_plot_name}")
    plt.clf()
    plt.cla()
    plt.close()

    # numbers for table
    upreg = df_sig.loc[df_sig['logFC'] > 0]
    # upreg = upreg.shape[0]
    upreg = upreg["Gene"].tolist()
    downreg = df_sig.loc[df_sig['logFC'] < 0]
    # downreg = downreg.shape[0]
    downreg = downreg["Gene"].tolist()
    # all_rows = df.shape[0]
    # nodiff = all_rows - upreg - downreg
    df_nonsig = df.drop(index=downreg, errors='ignore')
    df_nonsig = df_nonsig.drop(index=upreg, errors='ignore')
    nodiff = df_nonsig["Gene"].tolist()

    out_dict = {"Downregulated" : downreg, "no difference" : nodiff,  "Upregulated" : upreg}

    return out_dict


def plot_venn_DE_genes(tables_dict:dict, p_sig = 0.05, min_LFC = 0, venn_filename = "venn_diagram.png", venn_title = "", get_shared_list = False, plot=True):
    """
    plot a venn diagram of two or three lists of DE genes from R 
    showing the number of genes that are shared in all cases
    """
    if len(tables_dict) !=2 and len(tables_dict) != 3:
        raise RuntimeError(f"list should have only 2 or 3 elements but it has has length {len(tables_dict)}! \n{tables_dict}")
    
    sig_geneIDs_lists = {table_title : [] for table_title in tables_dict.keys()}
    for table_title,table_path in tables_dict.items():
        df = pd.read_csv(table_path, sep="\t", skiprows=0)
        
        df_sig = df.loc[df['FDR'] < p_sig]
        if min_LFC>0:
            df_sig = df_sig.loc[abs(df_sig['logFC']) >= min_LFC]
        sig_geneIDs_lists[table_title] = set(df_sig["Gene"].tolist())
        
    subsets = [id_list for id_list in sig_geneIDs_lists.values()]
    labels = [id_lab for id_lab in sig_geneIDs_lists.keys()]

    if plot: 
        fs = 15 # fontsize
        if len(tables_dict)==2:
            v = venn2(subsets=subsets, set_labels=labels)
        elif len(tables_dict)==3:
            v = venn3(subsets=subsets, set_labels=labels)
        for text in v.set_labels:
            try:
                text.set_fontsize(fs)
            except:
                pass
        for text in v.subset_labels:
            try:
                text.set_fontsize(fs)
            except:
                pass

        plt.title(venn_title, fontsize = fs)
        plt.tight_layout()
        plt.savefig(venn_filename, dpi = 300, transparent = True)
        print(f"plot saved in current working directory as: {venn_filename}")
        
        plt.clf()
        plt.cla()
        plt.close()
    
    if get_shared_list:
        if len(tables_dict)==2:
            intersection = set(subsets[0]) & set(subsets[1])
            return intersection
        elif len(tables_dict)==3:
            print(f"only return overlap list for two-way comparison! TODO implement if you want three")
            return None


def plot_sig_LFC_overlap(tables_dict:dict, p_sig = 0.05, min_LFC = 0, LFC_filename = "sig_LFC_scatter.png", LFC_title = "", excl_geneIDs =[], incl_geneIDs = []):
    """ 
    plot a scatterplot of sig. DE genes with LFC values
    """
    if len(tables_dict) !=2 :
        raise RuntimeError(f"list should have only 2 elements but it has has length {len(tables_dict)}! \n{tables_dict}")
    
    sig_geneIDs_lists = {table_title : [] for table_title in tables_dict.keys()}
    tables_df = {table_title : [] for table_title in tables_dict.keys()}
    tables_df_all = {table_title : [] for table_title in tables_dict.keys()}

    for table_title,table_path in tables_dict.items():

        df = pd.read_csv(table_path, sep="\t", skiprows=0)

        if len(incl_geneIDs) >0:
            full_size = df.shape[0]
            print(f"\t{table_title} : only including {len(incl_geneIDs)} from {full_size}")
            df = df[df["Gene"].isin(incl_geneIDs)]
            # if including only a subset of genes don't filter for only significant ones, include everything
        
        df_sig = df.loc[df['FDR'] < p_sig]
        tables_df_all[table_title] = df
        
        if min_LFC>0:
            df_sig = df_sig.loc[abs(df_sig['logFC']) >= min_LFC]
        
        tables_df[table_title] = df_sig
        sig_geneIDs_lists[table_title] = df_sig["Gene"].tolist()

    ## make the overlap of significanlty DE genes
    table_a,table_b = tables_dict.keys()
    set_a = set(sig_geneIDs_lists[table_a])
    set_b = set(sig_geneIDs_lists[table_b])

    lists = {
        "shared" : list(set_a & set_b),
        table_a : list(set_a - set_b),
        table_b : list(set_b - set_a)
        }

    fig, ax = plt.subplots(1,1, figsize=(13, 13)) 
    fs = 35
    point_size_factor = 8
    ps = fs*point_size_factor # point size

    sig_colors = ["#BD351E","#EA882C"] # red , orange
    colors_dict = {table_title : sig_colors[i] for i,table_title in enumerate(tables_dict.keys())}
    colors_dict["shared"] = "#3C7FA7" # blue
    colors_dict["nonsig"] = "#4B3B47" # mauve shadow

    if len(incl_geneIDs) == 0:
        excl_counter = { cat : 0 for cat in [table_a,table_b,"shared"]}
    else:
        excl_counter = { cat : 0 for cat in ["nonsig",table_a,table_b,"shared"]}
        all_sig = lists['shared']+lists[table_a]+lists[table_b]
        lists["nonsig"] = [geneID for geneID in incl_geneIDs if geneID not in all_sig]
    for cat in excl_counter.keys():
        for geneID in lists[cat]:
            try:
                x = tables_df_all[table_a].loc[geneID,"logFC"]
            except:
                x = 0
                print(geneID)
            try:
                y = tables_df_all[table_b].loc[geneID,"logFC"]
            except:
                y = 0
                print(geneID)
            if geneID in excl_geneIDs:
                excl_counter[cat]+=1
                ax.scatter(x,y,color = colors_dict[cat], s=ps*1.5, alpha = 1, marker="1")
            else:
                ax.scatter(x,y,color = colors_dict[cat], s=ps, alpha = 0.75)

    ax.set_xlabel(f"logFC {table_a}", fontsize = fs)
    ax.set_ylabel(f"logFC {table_b}", fontsize = fs)
    ax.tick_params(axis='x', labelsize=fs*0.9)
    ax.tick_params(axis='y', labelsize=fs*0.9)
    ax.set_title(LFC_title, fontsize = fs)

    min_yline,max_yline = ax.get_ylim()
    min_xline,max_xline = ax.get_xlim()
    min_yline = min_yline-0.5
    max_yline = max_yline+0.5
    min_xline = min_xline-0.5
    max_xline = max_xline+0.5
    if min_LFC == 0:
        ax.hlines(y=1, xmin=min_xline, xmax=max_xline, linewidth=point_size_factor*0.5, linestyle = ":", color="#818181",zorder=0)
        ax.hlines(y=-1, xmin=min_xline, xmax=max_xline, linewidth=point_size_factor*0.5, linestyle = ":", color="#818181",zorder=0)
        ax.vlines(x=1, ymin=min_yline, ymax=max_yline, linewidth=point_size_factor*0.5, linestyle = ":", color="#818181",zorder=0)
        ax.vlines(x=-1, ymin=min_yline, ymax=max_yline, linewidth=point_size_factor*0.5, linestyle = ":", color="#818181",zorder=0)
    else:
        ax.hlines(y=min_LFC, xmin=min_xline, xmax=max_xline, linewidth=point_size_factor*0.5, linestyle = ":", color="#818181",zorder=0)
        ax.hlines(y=-1*min_LFC, xmin=min_xline, xmax=max_xline, linewidth=point_size_factor*0.5, linestyle = ":", color="#818181",zorder=0)
        ax.vlines(x=min_LFC, ymin=min_yline, ymax=max_yline, linewidth=point_size_factor*0.5, linestyle = ":", color="#818181",zorder=0)
        ax.vlines(x=-1*min_LFC, ymin=min_yline, ymax=max_yline, linewidth=point_size_factor*0.5, linestyle = ":", color="#818181",zorder=0)
    
    ax.hlines(y=0, xmin=min_xline, xmax=max_xline, linewidth=point_size_factor, linestyle = ":", color="#2F3E1D",zorder=0)
    ax.vlines(x=0, ymin=min_yline, ymax=max_yline, linewidth=point_size_factor, linestyle = ":", color="#2F3E1D",zorder=0)
    
    ## make legend points
    for cat in reversed(list(excl_counter.keys())):
        ax.scatter(1000,1000,color = colors_dict[cat], s=ps, alpha = 0.75, label = cat)
    
    if excl_geneIDs != []:
        if "females" in LFC_filename:
            shared_lab = "also DE\nin males"
        else:
            shared_lab = "also DE\nin females"
        ax.scatter(1000,1000,color = "#666666", s=ps*1.5, alpha = 1, label = shared_lab, marker="1")
    
    ax.set_ylim([min_yline,max_yline])
    ax.set_xlim([min_xline,max_xline])

    plt.legend(fontsize = fs*0.8, title ="gene sig. in", title_fontsize = fs*0.8)
    plt.tight_layout()
    plt.savefig(LFC_filename, dpi = 300, transparent = True)
    if excl_geneIDs != []:
        print(f"genes excluded due to shared sex-bias: {excl_counter}")
    print(f"plot saved in current working directory as: {LFC_filename}")
    plt.clf()
    plt.cla()
    plt.close()

    lengths = { key : len(val) for key,val in lists.items()}
    return(lengths)


if __name__ == "__main__":

    username = "milena"
    table_paths,contrast_plot_titles = get_tables(username=username)
    out_path_figs = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures_python"
    
    ### genes that are line-biased in both males and females (should be excluded from male analysis since they can't be related to the Y-haplotype)
    excl_line_bias_lists = {
        "sex_separated" : {
            "day14" : ['gene-428738', 'gene-224697', 'gene-222350', 'gene-428765', 'gene-222600', 'gene-224875', 'gene-241001', 'gene-430032', 'gene-220028', 'gene-222486', 'gene-241055', 'gene-224357', 'gene-226245', 'gene-225738', 'gene-224968', 'gene-222531', 'gene-430263', 'gene-224201', 'gene-225107', 'gene-225236', 'gene-225140', 'gene-224227', 'gene-390616', 'gene-225709', 'gene-225325', 'gene-222332', 'gene-222519', 'gene-430314', 'gene-120952', 'gene-240871', 'gene-224860', 'gene-326873', 'gene-240929', 'gene-80359', 'gene-84970', 'gene-322912', 'gene-326849', 'gene-81427', 'gene-323148', 'gene-322927', 'gene-224782', 'gene-218529', 'gene-224743', 'gene-240623', 'gene-222383', 'gene-225173', 'gene-222365', 'gene-222344', 'gene-237881', 'gene-430068', 'gene-224956', 'gene-225720', 'gene-224682', 'gene-431701', 'gene-222555', 'gene-224896', 'gene-403809', 'gene-240910', 'gene-323803', 'gene-390956', 'gene-430080', 'gene-225635', 'gene-240833', 'gene-224593', 'gene-241126', 'gene-225030', 'gene-240691', 'gene-391222', 'gene-90157'],
            "day16" : ['gene-224697', 'gene-222600', 'gene-224875', 'gene-241001', 'gene-220028', 'gene-222486', 'gene-224357', 'gene-224968', 'gene-222159', 'gene-323148', 'gene-223773', 'gene-224782', 'gene-240623', 'gene-225173', 'gene-222344', 'gene-225720', 'gene-431701', 'gene-222555', 'gene-323803', 'gene-225635', 'gene-430080', 'gene-87700', 'gene-330102', 'gene-225030', 'gene-223419', 'gene-90157', 'gene-241262', 'gene-428738', 'gene-222350', 'gene-428765', 'gene-224079', 'gene-225325', 'gene-222332', 'gene-430314', 'gene-120952', 'gene-223491', 'gene-84970', 'gene-322927', 'gene-237881', 'gene-430068', 'gene-224956', 'gene-224682', 'gene-224896', 'g14784', 'gene-240833', 'gene-240691', 'gene-286545', 'gene-223318', 'gene-124877', 'gene-225738', 'gene-222531', 'gene-430263', 'gene-407280', 'gene-225140', 'gene-224227', 'gene-225709', 'gene-224890', 'gene-80359', 'gene-322912', 'gene-227370', 'gene-224743', 'gene-406796', 'gene-240910', 'gene-390956', 'gene-391222', 'gene-430032', 'gene-229506', 'gene-241055', 'gene-226245', 'gene-225107', 'gene-224201', 'gene-225236', 'gene-390616', 'gene-282853', 'gene-222519', 'gene-240871', 'gene-224860', 'gene-326873', 'gene-240929', 'gene-326849', 'gene-222383', 'gene-222365', 'gene-403809', 'gene-224593', 'gene-241126', 'gene-222746', 'gene-238407'],
            "day18" : ['gene-428738', 'gene-224697', 'gene-224875', 'gene-241055', 'gene-223758', 'gene-430263', 'gene-225236', 'gene-225325', 'gene-301479', 'gene-120952', 'gene-223491', 'gene-224860', 'gene-223773', 'gene-240623', 'gene-406796', 'gene-225720', 'gene-238849', 'gene-227308', 'gene-240833', 'gene-224593', 'gene-240691'],
            },
        "line_separated" : {
            "SL1" : [],
            "SL3" : ['gene-99775', 'gene-40274', 'gene-304827', 'gene-92346', 'gene-306335', 'gene-285669', 'gene-120763', 'gene-2286', 'gene-97407', 'gene-232392', 'gene-328941', 'gene-166511', 'gene-39692', 'gene-384091', 'gene-74686', 'gene-122220', 'gene-218723', 'gene-414353', 'gene-312890', 'gene-153482', 'gene-39770', 'gene-132340', 'gene-253632', 'gene-378608', 'gene-206556', 'gene-336703', 'gene-21229', 'gene-166391', 'gene-120784', 'gene-87502', 'gene-317372', 'gene-73253', 'gene-211196', 'gene-9548', 'gene-60190', 'gene-234650', 'gene-410057', 'gene-121262', 'gene-100036', 'gene-227137', 'gene-75744', 'gene-279912', 'gene-343203', 'gene-233901', 'gene-163028', 'gene-238407', 'gene-39680', 'gene-198700', 'gene-231228', 'gene-410209', 'gene-143368', 'gene-388261', 'gene-30328', 'gene-226944', 'gene-334263', 'gene-130081', 'gene-47823', 'gene-228519', 'gene-350792', 'gene-277340', 'gene-182683', 'gene-60151', 'gene-206576', 'gene-202718', 'gene-130096', 'gene-377275', 'gene-244780', 'gene-288834', 'gene-189246', 'gene-48535', 'gene-205011', 'gene-253157', 'gene-62891', 'gene-69698']
            },
        "day_separated" : {
            "day14" : ['gene-327441', 'gene-241001', 'gene-403652', 'gene-403851', 'gene-224860', 'gene-218086', 'gene-390956', 'gene-224682', 'gene-240910', 'gene-81551', 'gene-392224', 'gene-240983', 'gene-326825', 'gene-222519', 'gene-237881', 'gene-224968', 'gene-220544', 'gene-226245', 'gene-391198', 'gene-403818', 'gene-239553', 'gene-90157', 'gene-84970', 'gene-219019', 'gene-407253', 'gene-224845', 'gene-81640', 'gene-222383', 'gene-220249', 'gene-430080', 'gene-430044', 'gene-391222', 'gene-222365', 'gene-240833', 'gene-428765', 'gene-240691', 'gene-80359', 'gene-392248', 'gene-390678', 'gene-225629', 'gene-120952', 'gene-428756', 'gene-224956', 'gene-403809', 'gene-224227', 'gene-240935', 'gene-225720', 'gene-225140', 'gene-224614', 'gene-224357', 'gene-392159', 'gene-222531', 'gene-225325', 'gene-117712', 'gene-84949', 'gene-222344', 'gene-406796', 'gene-240929', 'gene-392290', 'gene-81599', 'gene-225107', 'gene-231925', 'gene-214979', 'gene-220028', 'gene-260693', 'gene-224277', 'gene-224250', 'gene-224782', 'gene-390616', 'gene-224593', 'gene-322912', 'gene-431701', 'gene-80466', 'gene-403706', 'gene-323148', 'gene-323803', 'gene-283443', 'gene-326909', 'gene-225173', 'gene-243308', 'gene-430068', 'gene-222430', 'gene-240871', 'gene-237318', 'gene-224697', 'gene-222332', 'gene-241126', 'gene-89234', 'gene-406468', 'gene-328764', 'gene-326810', 'gene-88715', 'gene-222486', 'gene-224743', 'gene-221953', 'gene-225635', 'gene-222501', 'gene-428738', 'gene-430032', 'gene-241108', 'gene-81572', 'gene-392186', 'gene-222555', 'gene-81427', 'gene-225738', 'gene-403700', 'gene-234575', 'gene-240623', 'gene-222350', 'gene-395080', 'gene-225236', 'gene-224201', 'gene-241055', 'gene-225709', 'gene-326873', 'gene-326849', 'gene-322927', 'gene-80484', 'gene-395143', 'gene-225030', 'gene-222600', 'gene-224875', 'gene-390637', 'gene-224307', 'gene-403583', 'gene-224896', 'gene-430263', 'gene-428747', 'gene-430314', 'gene-403902'],
            "day16" : ['gene-241001', 'gene-282746', 'gene-224860', 'gene-227308', 'gene-390956', 'gene-224682', 'gene-240910', 'gene-399475', 'gene-240983', 'gene-222519', 'gene-282458', 'gene-237881', 'gene-224968', 'gene-220544', 'gene-226245', 'gene-240860', 'gene-391198', 'gene-282524', 'gene-403818', 'gene-90157', 'gene-84970', 'gene-224845', 'gene-81640', 'gene-222383', 'gene-220249', 'gene-430080', 'gene-430044', 'gene-391222', 'gene-222365', 'gene-240833', 'gene-428765', 'gene-240691', 'gene-80359', 'gene-390678', 'gene-225629', 'gene-238849', 'gene-282620', 'gene-120952', 'gene-428756', 'gene-224956', 'gene-403809', 'gene-282853', 'gene-224227', 'gene-225720', 'gene-225140', 'gene-238407', 'gene-239506', 'gene-224357', 'gene-222531', 'gene-241682', 'gene-225325', 'gene-84949', 'gene-222344', 'gene-406796', 'gene-240929', 'gene-282551', 'gene-240638', 'gene-392290', 'gene-225158', 'gene-225107', 'gene-220028', 'gene-260693', 'gene-224277', 'gene-224782', 'gene-390616', 'gene-224593', 'gene-322912', 'gene-431701', 'gene-428113', 'gene-380466', 'gene-323148', 'gene-323803', 'gene-282665', 'gene-326909', 'gene-428104', 'gene-225173', 'gene-227370', 'gene-399317', 'gene-224890', 'gene-282347', 'gene-282784', 'gene-430068', 'gene-282491', 'gene-400393', 'gene-222430', 'gene-399484', 'gene-240871', 'gene-224697', 'gene-222332', 'gene-241126', 'gene-89234', 'gene-399424', 'gene-282362', 'g14784', 'gene-223491', 'gene-400384', 'gene-222486', 'gene-224743', 'gene-243299', 'gene-221953', 'gene-225635', 'gene-222501', 'gene-428738', 'gene-400402', 'gene-430032', 'gene-241108', 'gene-282398', 'gene-81572', 'gene-392186', 'gene-222555', 'gene-81427', 'gene-225738', 'gene-282701', 'gene-282886', 'gene-240623', 'gene-222159', 'gene-222350', 'gene-395080', 'gene-225236', 'gene-224201', 'gene-241055', 'gene-286545', 'gene-225709', 'gene-431030', 'gene-326849', 'gene-322927', 'gene-395143', 'gene-225030', 'gene-222600', 'gene-224875', 'gene-223773', 'gene-390637', 'gene-224307', 'gene-231854', 'gene-224896', 'gene-430263', 'gene-428747', 'gene-399270', 'gene-241262', 'gene-430314', 'gene-242512', 'gene-282590'],
            "day18" : [],
            },
    }

            
    ############################################
    ######### MAKE ALL THE SMEAR PLOTS #########
    ############################################
    if False:
        for separation, seps_dict in table_paths.items():
            print(f"\n=========================== {separation} ===========================")
            lists = {}
            
            if "day" not in separation:
                continue

            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")

                numbers = {}
                lists[category] = {}
                for contrast, table_path in paths_dict.items():
                    
                    if "(" not in contrast:
                        continue
                    
                    print(f"{separation}:{category} --> contrast: {contrast}")
                    table_name = table_path.split("/")[-1].replace(".txt", "").replace("DE_genes_", "")
                    smear_name = f"{out_path_figs}/smear_{table_name}.png"
                    smear_title = contrast_plot_titles[contrast]

                    excl_line_bias = excl_line_bias_lists[separation]
                    excl_list = []
                    excl_list_name = smear_title.replace(" ","")
                    if excl_list_name in excl_line_bias:
                        # only do the exclusion list when it's actually in a relevant contrast
                        excl_list = excl_line_bias[excl_list_name]
                        print(f"\texcluding genes from list '{excl_list_name}'")
                    elif "day" in separation and category in excl_line_bias:
                        excl_list = excl_line_bias[category]
                        print(f"\texcluding genes from list '{excl_list_name}'")
        
                    smear_lists = plot_smear(table_path=table_path, contrast=contrast, smear_plot_name=smear_name, title = smear_title, excl_genes_list=excl_list, x_axis="fdr_p")
                    smear_lists = plot_smear(table_path=table_path, contrast=contrast, smear_plot_name=smear_name, title = smear_title, excl_genes_list=excl_list, x_axis="logcpm")
                    if True:
                        Downlist = smear_lists["Downregulated"]
                        Uplist = smear_lists["Upregulated"]
                        print(f"\t * Downregulated : {Downlist}\n\t * Upregulated : {Uplist}")

                    smear_nums = {gene_set : len(gene_list) for gene_set,gene_list in smear_lists.items()}
                    numbers[smear_title] = smear_nums
                    lists[category][smear_title] = smear_lists

                print(numbers)

    ############################################
    ######## MAKE ALL THE VENN DIAGRAMS ########
    ############################################
    
    ## standard sets matching the tabs in the html
    if False:
        venn_sets = {
            "sex_separated" : {
                "females" : {
                    "age by line bias" :["SL1_14 - SL3_14","SL1_16 - SL3_16","SL1_18 - SL3_18"],
                    "line by age bias" : ["SL1_18 - (SL1_14+SL1_16)/2","SL3_18 - (SL3_14+SL3_16)/2"]
                },
                "males" : {
                    "age by line bias" : ["SL1_14 - SL3_14","SL1_16 - SL3_16","SL1_18 - SL3_18"],
                    "line by age bias" : ["SL1_18 - (SL1_14+SL1_16)/2","SL3_18 - (SL3_14+SL3_16)/2"]
                }
            },
            "line_separated" : {
                "SL1" : {"age by sex bias" : ["F_14 - M_14","F_16 - M_16","F_18 - M_18"]},
                "SL3" : {"age by sex bias" : ["F_14 - M_14","F_16 - M_16","F_18 - M_18"]}
            },
            "day_separated" : {
                "day14" : {
                    "sex bias by line" : ["F_1 - M_1","F_3 - M_3"],
                    "line bias by sex" : ["F_1 - F_3","M_1 - M_3"],
                },
                "day16" : {
                    "sex bias by line" : ["F_1 - M_1","F_3 - M_3"],
                    "line bias by sex" : ["F_1 - F_3","M_1 - M_3"],
                    },
                "day18" : {
                    "sex bias by line" : ["F_1 - M_1","F_3 - M_3"],
                    "line bias by sex" : ["F_1 - F_3","M_1 - M_3"],
                    },
            }
        }
        for separation, seps_dict in table_paths.items():
            print(f"\n=========================== {separation} ===========================")

            if "day" not in separation:
                continue

            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")

                for venn_cat, venn_contrasts_list in venn_sets[separation][category].items():
                    print(f"{venn_cat} : {venn_contrasts_list}")
                    venn_paths_dict = {contrast_plot_titles[contrast] : paths_dict[contrast] for contrast in venn_contrasts_list}
                    venn_filename_ = venn_cat.replace(" ", "_")
                    venn_filename = f"{out_path_figs}/Venn_{category}_{venn_filename_}.png"
                    venn_title = f"sig. DE genes overlap ({category})\n{venn_cat}"
                    plot_venn_DE_genes(venn_paths_dict, venn_filename=venn_filename, venn_title=venn_title)

    ##### make all the lists of genes to exclude because they are line-biased in females

    ## compare if the same genes are DE between lines within days in males as in females
    if False:
        venn_sets_day = {
            "sex_separated" : {
                "day14" : {
                    "females" : ["SL1_14 - SL3_14"],
                    "males" : ["SL1_14 - SL3_14"]
                },
                "day16" : {
                    "females" : ["SL1_16 - SL3_16"],
                    "males" : ["SL1_16 - SL3_16"]
                },
                "day18" : {
                    "females" : ["SL1_18 - SL3_18"],
                    "males" : ["SL1_18 - SL3_18"]
                }
            }
        }
        
        for separation, days_dict in venn_sets_day.items():
            print(f"\n=========================== {separation} ===========================")
            for day, sexes_contrasts_dict in days_dict.items():
                print(f"\n ------------------- {day} -------------------")

                venn_paths_dict = {}
                for sex, venn_contrasts_list in sexes_contrasts_dict.items():
                    print(f"{sex} : {venn_contrasts_list}")

                    venn_paths_dict[sex] = table_paths[separation][sex][venn_contrasts_list[0]]
                
                venn_filename = f"{out_path_figs}/Venn_{day}_f_vs_m.png"
                day_ = day.replace("day", "day ")
                venn_title = f"sig. DE genes overlap ({day_})\nfemales and males"
                shared_list = plot_venn_DE_genes(venn_paths_dict, venn_filename=venn_filename, venn_title=venn_title, get_shared_list=True)
                print(shared_list)
    
    ## compare if the same genes are DE between days within lines in males as in females
    if False:
        venn_sets_line = {
            "sex_separated" : {
                "SL1" : {
                    "females" : ["SL1_18 - (SL1_14+SL1_16)/2"],
                    "males" : ["SL1_18 - (SL1_14+SL1_16)/2"]
                },
                "SL3" : {
                    "females" : ["SL3_18 - (SL3_14+SL3_16)/2"],
                    "males" : ["SL3_18 - (SL3_14+SL3_16)/2"]
                },
            }
        }
        for separation, days_dict in venn_sets_line.items():
            print(f"\n=========================== {separation} ===========================")
            for line, sexes_contrasts_dict in days_dict.items():
                print(f"\n ------------------- {line} -------------------")

                venn_paths_dict = {}
                for sex, venn_contrasts_list in sexes_contrasts_dict.items():
                    print(f"{sex} : {venn_contrasts_list}")

                    venn_paths_dict[sex] = table_paths[separation][sex][venn_contrasts_list[0]]
                
                venn_filename = f"{out_path_figs}/Venn_{line}_f_vs_m.png"
                venn_title = f"sig. DE genes overlap ({line})\nfemales and males"
                shared_list = plot_venn_DE_genes(venn_paths_dict, venn_filename=venn_filename, venn_title=venn_title, get_shared_list=True)
                print(f"{len(shared_list)} genes : \n{shared_list}")

    ## compare if the same genes are DE between lines within sexes in days
    if False:
        venn_sets_day = {
            "day_separated" : {
                "day14" : {
                    "males" : ["M_1 - M_3"],
                    "females" : ["F_1 - F_3"],
                },
                "day16" : {
                    "males" : ["M_1 - M_3"],
                    "females" : ["F_1 - F_3"],
                    },
                "day18" : {
                    "males" : ["M_1 - M_3"],
                    "females" : ["F_1 - F_3"],
                    },
            }
        }
        for separation, days_dict in venn_sets_day.items():
            print(f"\n=========================== {separation} ===========================")
            for day, sexes_contrasts_dict in days_dict.items():
                print(f"\n ------------------- {day} -------------------")

                venn_paths_dict = {}
                for sex, venn_contrasts_list in sexes_contrasts_dict.items():
                    print(f"{sex} : {venn_contrasts_list}")

                    venn_paths_dict[sex] = table_paths[separation][day][venn_contrasts_list[0]]
                
                venn_filename = f"{out_path_figs}/Venn_{day}_f_vs_m.png"
                venn_title = f"sig. DE genes overlap ({day})\nfemales and males"
                shared_list = plot_venn_DE_genes(venn_paths_dict, venn_filename=venn_filename, venn_title=venn_title, get_shared_list=True, plot=False)
                print(f"{len(shared_list)} genes : \n{shared_list}")

    ############################################
    ######## MAKE ALL THE SB COMPARISONS #######
    ############################################

    if False:
        LFC_comp_sets = {
            "sex_separated" : {
                "females" : {
                    "day 14 vs. day 16 line bias" : ["SL1_14 - SL3_14","SL1_16 - SL3_16"],
                    "day 14 vs. day 18 line bias" : ["SL1_14 - SL3_14","SL1_18 - SL3_18"],
                    "day 16 vs. day 18 line bias" : ["SL1_16 - SL3_16","SL1_18 - SL3_18"],
                    "line by age bias" : ["SL1_18 - (SL1_14+SL1_16)/2","SL3_18 - (SL3_14+SL3_16)/2"]
                },
                "males" : {
                    "day 14 vs. day 16 line bias" : ["SL1_14 - SL3_14","SL1_16 - SL3_16"],
                    "day 14 vs. day 18 line bias" : ["SL1_14 - SL3_14","SL1_18 - SL3_18"],
                    "day 16 vs. day 18 line bias" : ["SL1_16 - SL3_16","SL1_18 - SL3_18"],
                    "line by age bias" : ["SL1_18 - (SL1_14+SL1_16)/2","SL3_18 - (SL3_14+SL3_16)/2"]
                }
            },
            "line_separated" : {
                "SL1" : {
                    "day 14 vs. day 16 sex bias" : ["F_14 - M_14","F_16 - M_16"],
                    "day 14 vs. day 18 sex bias" : ["F_14 - M_14","F_18 - M_18"],
                    "day 16 vs. day 18 sex bias" : ["F_16 - M_16","F_18 - M_18"],
                },
                "SL3" : {
                    "day 14 vs. day 16 sex bias" : ["F_14 - M_14","F_16 - M_16"],
                    "day 14 vs. day 18 sex bias" : ["F_14 - M_14","F_18 - M_18"],
                    "day 16 vs. day 18 sex bias" : ["F_16 - M_16","F_18 - M_18"],
                }
            },
            "day_separated" : {
                "day14" : {
                    "sex bias by line" : ["F_1 - M_1","F_3 - M_3"],
                    "line bias by sex" : ["F_1 - F_3","M_1 - M_3"],
                },
                "day16" : {
                    "sex bias by line" : ["F_1 - M_1","F_3 - M_3"],
                    "line bias by sex" : ["F_1 - F_3","M_1 - M_3"],
                    },
                "day18" : {
                    "sex bias by line" : ["F_1 - M_1","F_3 - M_3"],
                    "line bias by sex" : ["F_1 - F_3","M_1 - M_3"],
                    },
            }
        }

        def get_excl_genes_list(contrasts:list, contrast_plot_titles:dict, excl_line_bias_lists:dict):
            excl_list = []
            for contrast in contrasts:
                day_ = contrast_plot_titles[contrast]
                day = day_.replace(" ", "")
                if day in  excl_line_bias_lists:
                    excl_list.extend(excl_line_bias_lists[day])
            # make unique
            return list(set(excl_list))

            

        for separation, seps_dict in table_paths.items():
            print(f"\n=========================== {separation} ===========================")

            if "day" not in separation:
                continue

            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")

                for LFC_cat, LFC_contrasts_list in LFC_comp_sets[separation][category].items():
                    print(f"{LFC_cat} : {LFC_contrasts_list}")
                    LFC_filename_ = LFC_cat.replace(" ", "_").replace(".", "")
                    LFC_filename = f"{out_path_figs}/LFC_scatter_{category}_{LFC_filename_}.png"

                    LFC_paths_dict = {contrast_plot_titles[contrast] : paths_dict[contrast] for contrast in LFC_contrasts_list}
                    excl_line_bias = excl_line_bias_lists[separation]
                    if "day" not in separation:
                        excl_geneIDs = get_excl_genes_list(contrasts = LFC_contrasts_list, contrast_plot_titles=contrast_plot_titles, excl_line_bias_lists=excl_line_bias)
                    else:
                        excl_geneIDs = excl_line_bias[category]
                    
                    if separation == "line_separated":
                        # when looking at sex bias explicitly (within each line) it makes no sense to exclude genes based on shared sex-bias
                        excl_geneIDs = []

                    if "line bias" in LFC_cat:
                        title = f"{category}: contrast SL1 - SL3"
                    elif "sex bias" in LFC_cat:
                        title = f"{category}: female - male"
                    elif "age bias" in LFC_cat:
                        title = f"{category}: day 18 - mean(day 14, day 16)"
                    else:
                        title = LFC_cat

                    if excl_geneIDs == []:
                        print(f"no excluded genes")
                    else:
                        print(f"{len(excl_geneIDs)} excluded geneIDs")

                    numbers = plot_sig_LFC_overlap(LFC_paths_dict, LFC_filename = LFC_filename, LFC_title = title, excl_geneIDs = excl_geneIDs)
                    print(f"\t{numbers}")
    
    ### plot only the sig DE genes in the interactions
    if True:
        sig_IDs_list = {
            "day_separated" : {
                "day14" : ['gene-124877', 'gene-124865', 'gene-239506', 'gene-223773', 'gene-224079', 'gene-223758', 'gene-424750', 'gene-223491', 'gene-371957', 'gene-127740', 'gene-429522', 'gene-90190', 'gene-87554', 'gene-426056', 'gene-231854', 'gene-230270', 'gene-288660', 'gene-232048', 'gene-220055', 'gene-77070', 'gene-290138', 'gene-120832', 'gene-39545', 'gene-286732', 'gene-130096', 'gene-241328', 'gene-218301', 'gene-76019', 'gene-166391', 'gene-130081', 'gene-38977', 'gene-231604', 'gene-223318', 'gene-39692', 'gene-158197', 'gene-414353', 'gene-241268', 'gene-266875', 'gene-222737', 'gene-88032', 'gene-115312', 'gene-372264', 'gene-122220', 'gene-288356', 'gene-70385', 'gene-204669', 'gene-268146', 'gene-374892', 'gene-375015', 'gene-406510', 'gene-326329', 'gene-407355', 'gene-327787', 'gene-210277', 'gene-125867', 'gene-223419', 'gene-346228', 'gene-143368', 'gene-39770', 'gene-414838', 'gene-161821', 'gene-125638', 'gene-39680', 'gene-75424', 'gene-78574', 'gene-149137', 'gene-217258', 'gene-210286', 'gene-53501', 'gene-240397', 'gene-336688', 'gene-95461', 'gene-100036', 'gene-71572', 'gene-227137', 'gene-236155', 'gene-241262', 'gene-336088', 'gene-351784', 'gene-263126', 'gene-188501', 'gene-127607', 'gene-6017', 'gene-223599', 'gene-24347', 'gene-73885', 'gene-302994', 'gene-215608', 'gene-62772', 'gene-266887', 'gene-391847', 'gene-223554', 'gene-246615', 'gene-289849', 'gene-60429', 'gene-229506', 'gene-233901', 'gene-124680', 'gene-244661', 'gene-69775', 'gene-272401', 'gene-89057', 'gene-286545', 'gene-124766', 'gene-346500', 'gene-2467', 'gene-83083', 'gene-89798', 'gene-2286', 'gene-333750', 'gene-311581', 'gene-406519', 'gene-67440', 'gene-241193', 'gene-161848', 'gene-283197', 'gene-223791', 'gene-72046', 'gene-85294', 'g11517', 'gene-241506', 'gene-224028', 'gene-129852', 'gene-335551', 'gene-223088', 'gene-73754', 'gene-39533', 'gene-266297', 'gene-218529', 'gene-323550', 'gene-250391', 'gene-377275', 'gene-83830', 'gene-425122', 'gene-68813', 'gene-57335', 'gene-269058', 'gene-238320', 'gene-347681', 'gene-215430', 'gene-53269', 'gene-31110', 'gene-68612', 'gene-212020', 'gene-55557', 'gene-223285', 'gene-280630', 'gene-211847', 'gene-271655', 'gene-406603', 'gene-277078', 'gene-203445', 'gene-77109', 'gene-120964', 'gene-7268', 'gene-179860', 'gene-62891', 'gene-77835', 'gene-30328', 'gene-127773','gene-227308', 'gene-224369', 'gene-225158', 'gene-221012', 'g14784', 'gene-125210', 'gene-224890', 'gene-238849', 'gene-403878', 'gene-120660', 'gene-237378', 'gene-119161', 'gene-241682', 'gene-243630', 'gene-218813', 'gene-421265', 'gene-84577', 'gene-370643', 'gene-240602', 'gene-227370', 'gene-240860', 'gene-5731', 'gene-330102', 'gene-217099', 'gene-48714', 'gene-81750', 'gene-370487', 'gene-90428', 'gene-395158', 'gene-217443', 'gene-371230', 'gene-303243', 'gene-124377', 'gene-68558', 'gene-370842', 'gene-215103', 'gene-62927', 'gene-301479', 'gene-370323', 'gene-64491', 'gene-227284', 'gene-118620', 'gene-221469', 'gene-417051', 'gene-250269', 'gene-370562', 'gene-333008', 'gene-392350', 'gene-82722', 'gene-264249', 'gene-81581', 'gene-407735', 'gene-73288', 'gene-272192', 'gene-262796', 'gene-81518', 'gene-226254', 'gene-57617', 'gene-408949', 'gene-120070', 'gene-81458', 'gene-398993', 'gene-81476', 'gene-334677', 'gene-421549', 'gene-370595', 'gene-425532', 'gene-241841', 'gene-81675', 'gene-118985', 'gene-81418', 'gene-65889', 'gene-88458', 'gene-240638', 'g11957', 'gene-80415', 'gene-87700', 'gene-225355', 'gene-76284', 'gene-401486', 'gene-227071', 'gene-399223', 'gene-263313', 'gene-313589', 'gene-288413', 'gene-166814', 'gene-237857', 'gene-65163', 'gene-87487', 'gene-229515', 'gene-219157', 'gene-350813', 'gene-234256', 'gene-64470', 'gene-57689', 'gene-263588', 'gene-242691', 'gene-294364', 'gene-400426', 'gene-222159', 'gene-268137', 'gene-88092', 'gene-324340', 'gene-262046', 'gene-90524', 'g5814', 'gene-410993', 'gene-399203', 'gene-215313', 'gene-358397', 'gene-361221', 'gene-417093', 'gene-261982', 'gene-378262', 'gene-131338', 'gene-77588', 'gene-152204', 'gene-26634', 'gene-58495', 'gene-221504', 'g1591', 'gene-392536', 'gene-407714', 'gene-370967', 'gene-266345', 'gene-153689', 'gene-125849', 'gene-5187', 'gene-331395', 'gene-55887', 'gene-130632', 'gene-73160', 'gene-263597', 'gene-331896', 'gene-271117', 'gene-59093', 'gene-431030', 'gene-397251', 'gene-283696', 'gene-119684'],
                "day16" : ['gene-424768', 'gene-263126', 'gene-240397', 'gene-424750', 'gene-211137', 'gene-243753', 'gene-215596', 'gene-63245', 'gene-223722', 'gene-428089', 'gene-54363', 'gene-424726', 'gene-87554', 'gene-223758', 'gene-55869', 'gene-90190', 'gene-223419', 'gene-223734', 'gene-367523', 'gene-42302', 'gene-113411', 'gene-14276', 'gene-217473', 'gene-7220', 'gene-63617', 'gene-424863', 'gene-229506', 'gene-6223', 'gene-34632', 'gene-119402', 'gene-312890', 'gene-69401', 'gene-241506', 'gene-23181', 'gene-127607', 'gene-336178', 'gene-8365', 'gene-38885', 'gene-328746', 'gene-388769', 'gene-86738', 'gene-14252', 'gene-231493', 'gene-125638', 'gene-393138', 'gene-279676', 'gene-221980', 'gene-407280', 'gene-68612', 'gene-152989', 'gene-104371', 'gene-234686', 'gene-271655', 'gene-336688', 'gene-241856', 'gene-306977', 'gene-127740', 'gene-233883', 'gene-57335', 'gene-38977', 'gene-85401', 'gene-353013', 'gene-269365', 'gene-346228', 'gene-201605', 'gene-127849', 'gene-205849', 'gene-152165', 'gene-431362', 'gene-313040', 'gene-174594', 'gene-294618', 'gene-128145', 'gene-228118', 'gene-88032', 'gene-32436', 'gene-424896', 'gene-124877', 'gene-231604', 'gene-420308', 'gene-149557', 'gene-327074', 'gene-426056', 'gene-223318', 'gene-283260', 'gene-223512', 'gene-230270', 'gene-80062', 'gene-279912', 'gene-231540', 'gene-333603', 'gene-116902', 'gene-410561', 'gene-259467', 'gene-183665', 'gene-240184', 'gene-367095', 'gene-424914', 'gene-195335', 'gene-345135', 'gene-351334', 'gene-254128', 'gene-222737', 'gene-277218', 'gene-23163', 'gene-23042', 'gene-48598', 'gene-68000', 'gene-185170', 'gene-402875', 'gene-197114', 'gene-421566', 'gene-233485', 'gene-268996', 'gene-272072', 'gene-30322', 'gene-55240', 'gene-39933', 'gene-405355', 'gene-402536', 'gene-241238', 'gene-310012', 'gene-333433', 'gene-282008', 'gene-255088', 'gene-324223', 'gene-17262', 'gene-30595', 'gene-90918', 'gene-279975', 'gene-321078', 'gene-327616', 'gene-254639', 'gene-372264', 'gene-101822', 'gene-266887', 'gene-174585', 'gene-246732', 'gene-286289', 'gene-21635', 'gene-287310', 'gene-275750', 'gene-362311', 'gene-64360', 'gene-425173', 'gene-232048', 'gene-196966', 'gene-428194', 'gene-89219', 'gene-268975', 'gene-55252', 'gene-183393', 'gene-182960', 'gene-329366', 'gene-6199', 'gene-269228', 'gene-222746', 'gene-336348', 'gene-289528', 'gene-14156', 'gene-369289', 'gene-228465', 'gene-117360', 'gene-16790', 'gene-218421', 'gene-188158', 'gene-286744', 'gene-270896', 'gene-428071', 'gene-360503', 'gene-330422', 'gene-58752', 'gene-220788','gene-240602', 'gene-24185', 'gene-24290', 'gene-23840', 'gene-23597', 'gene-24203', 'gene-24120', 'gene-24132', 'gene-23538', 'gene-423321', 'gene-87700', 'gene-24088', 'gene-24221', 'gene-390687', 'gene-24278', 'gene-24167', 'gene-13404', 'gene-23514', 'gene-407253', 'gene-326882', 'gene-23365', 'gene-15763', 'gene-23689', 'gene-428729', 'gene-23413', 'gene-23834', 'gene-327441', 'gene-23893', 'gene-392224', 'gene-23884', 'gene-240935', 'gene-80466', 'gene-421265', 'gene-120660', 'gene-428774', 'gene-24052', 'gene-90307', 'gene-24079', 'gene-122692', 'gene-417051', 'gene-282611', 'gene-403700', 'gene-224250', 'gene-282641', 'gene-328764', 'gene-27466', 'gene-370842', 'gene-84577', 'gene-219019', 'gene-410366', 'gene-400426', 'gene-214979', 'gene-403652', 'gene-227164', 'gene-326810', 'gene-130665', 'gene-9713', 'gene-17601', 'gene-12075', 'gene-218086', 'gene-84224', 'gene-237372', 'gene-219157', 'gene-326873', 'gene-329410', 'gene-403706', 'gene-81551', 'gene-245743', 'gene-81599', 'gene-411056', 'gene-127707', 'gene-395158', 'gene-431788', 'gene-217099', 'gene-80484', 'gene-403851', 'gene-23911', 'gene-370550', 'gene-81675', 'g2779', 'gene-215563', 'gene-392159', 'gene-224369', 'gene-90503', 'gene-24914', 'gene-216914', 'gene-81418', 'gene-370643', 'gene-15918', 'gene-237351', 'gene-81476', 'gene-81581', 'gene-285127', 'gene-222242', 'gene-125849', 'gene-349163', 'gene-227137', 'gene-282437', 'gene-13124', 'gene-17046', 'gene-10767', 'gene-225355', 'gene-237857', 'gene-119684', 'gene-16667', 'gene-81458', 'gene-80415', 'gene-326825'],
                "day18" : ['gene-428071', 'gene-233410', 'gene-428104','gene-396259', 'gene-395158', 'gene-224875', 'gene-395143', 'gene-90157', 'gene-224860', 'gene-428738', 'gene-395080', 'gene-224697', 'gene-301479', 'gene-242595'],
            }
        }
        contrasts_interaction_list = {
            "day_separated" : {
                "sex bias" : ["F_1 - M_1","F_3 - M_3"],
                "line bias" : ["F_1 - F_3","M_1 - M_3"],
            }
        }

        for separation, seps_dict in table_paths.items():
            print(f"\n=========================== {separation} ===========================")

            if "day" not in separation:
                continue

            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")

                for bias_cat, contrasts_list in contrasts_interaction_list[separation].items():

                    incl_geneIDs = sig_IDs_list[separation][category]
                    print(f"{bias_cat} : {len(incl_geneIDs)} genes")

                    LFC_paths_dict = {contrast_plot_titles[contrast] : paths_dict[contrast] for contrast in contrasts_list}
                    LFC_filename_ = bias_cat.replace(" ", "_")
                    LFC_filename = f"{out_path_figs}/LFC_scatter_interaction_{category}_{LFC_filename_}.png"

                    numbers = plot_sig_LFC_overlap(LFC_paths_dict, LFC_filename = LFC_filename, LFC_title = bias_cat, incl_geneIDs=incl_geneIDs)
                    print(f"\t{numbers}")

