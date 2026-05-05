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
            },
            "SL3" : {
                "F_14 - M_14" : f"{tables_dir}DE_genes_SL3_day14_F-M.txt",
                "F_16 - M_16" : f"{tables_dir}DE_genes_SL3_day16_F-M.txt",
                "F_18 - M_18" : f"{tables_dir}DE_genes_SL3_day18_F-M.txt"
            }
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
    elif "F" in contrast:
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



def plot_sig_LFC_overlap(tables_dict:dict, p_sig = 0.05, min_LFC = 0, LFC_filename = "sig_LFC_scatter.png", LFC_title = "", excl_geneIDs =[]):
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
        tables_df_all[table_title] = df

        df_sig = df.loc[df['FDR'] < p_sig]
        if min_LFC>0:
            df_sig = df_sig.loc[abs(df_sig['logFC']) >= min_LFC]
        
        tables_df[table_title] = df_sig
        sig_geneIDs_lists[table_title] = df_sig["Gene"].tolist()

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

    excl_counter = { cat : 0 for cat in [table_a,table_b,"shared"]}
    for cat in [table_a,table_b,"shared"]:
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
    for cat in [table_a,table_b,"shared"]:
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
    print(f"genes excluded due to shared sex-bias: {excl_counter}")
    print(f"plot saved in current working directory as: {LFC_filename}")
    plt.clf()
    plt.cla()
    plt.close()

    lengths = { key : len(val) for key,val in lists.items()}
    return(lengths)


if __name__ == "__main__":

    username = "miltr339"
    table_paths,contrast_plot_titles = get_tables(username=username)
    out_path_figs = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures_python"
    
    ### genes that are line-biased in both males and females (should be excluded from male analysis since they can't be related to the Y-haplotype)
    excl_line_bias_lists = {
        "day14" : ['gene-428738', 'gene-224697', 'gene-222350', 'gene-428765', 'gene-222600', 'gene-224875', 'gene-241001', 'gene-430032', 'gene-220028', 'gene-222486', 'gene-241055', 'gene-224357', 'gene-226245', 'gene-225738', 'gene-224968', 'gene-222531', 'gene-430263', 'gene-224201', 'gene-225107', 'gene-225236', 'gene-225140', 'gene-224227', 'gene-390616', 'gene-225709', 'gene-225325', 'gene-222332', 'gene-222519', 'gene-430314', 'gene-120952', 'gene-240871', 'gene-224860', 'gene-326873', 'gene-240929', 'gene-80359', 'gene-84970', 'gene-322912', 'gene-326849', 'gene-81427', 'gene-323148', 'gene-322927', 'gene-224782', 'gene-218529', 'gene-224743', 'gene-240623', 'gene-222383', 'gene-225173', 'gene-222365', 'gene-222344', 'gene-237881', 'gene-430068', 'gene-224956', 'gene-225720', 'gene-224682', 'gene-431701', 'gene-222555', 'gene-224896', 'gene-403809', 'gene-240910', 'gene-323803', 'gene-390956', 'gene-430080', 'gene-225635', 'gene-240833', 'gene-224593', 'gene-241126', 'gene-225030', 'gene-240691', 'gene-391222', 'gene-90157'],
        "day16" : ['gene-224697', 'gene-222600', 'gene-224875', 'gene-241001', 'gene-220028', 'gene-222486', 'gene-224357', 'gene-224968', 'gene-222159', 'gene-323148', 'gene-223773', 'gene-224782', 'gene-240623', 'gene-225173', 'gene-222344', 'gene-225720', 'gene-431701', 'gene-222555', 'gene-323803', 'gene-225635', 'gene-430080', 'gene-87700', 'gene-330102', 'gene-225030', 'gene-223419', 'gene-90157', 'gene-241262', 'gene-428738', 'gene-222350', 'gene-428765', 'gene-224079', 'gene-225325', 'gene-222332', 'gene-430314', 'gene-120952', 'gene-223491', 'gene-84970', 'gene-322927', 'gene-237881', 'gene-430068', 'gene-224956', 'gene-224682', 'gene-224896', 'g14784', 'gene-240833', 'gene-240691', 'gene-286545', 'gene-223318', 'gene-124877', 'gene-225738', 'gene-222531', 'gene-430263', 'gene-407280', 'gene-225140', 'gene-224227', 'gene-225709', 'gene-224890', 'gene-80359', 'gene-322912', 'gene-227370', 'gene-224743', 'gene-406796', 'gene-240910', 'gene-390956', 'gene-391222', 'gene-430032', 'gene-229506', 'gene-241055', 'gene-226245', 'gene-225107', 'gene-224201', 'gene-225236', 'gene-390616', 'gene-282853', 'gene-222519', 'gene-240871', 'gene-224860', 'gene-326873', 'gene-240929', 'gene-326849', 'gene-222383', 'gene-222365', 'gene-403809', 'gene-224593', 'gene-241126', 'gene-222746', 'gene-238407'],
        "day18" : ['gene-428738', 'gene-224697', 'gene-224875', 'gene-241055', 'gene-223758', 'gene-430263', 'gene-225236', 'gene-225325', 'gene-301479', 'gene-120952', 'gene-223491', 'gene-224860', 'gene-223773', 'gene-240623', 'gene-406796', 'gene-225720', 'gene-238849', 'gene-227308', 'gene-240833', 'gene-224593', 'gene-240691'],
        "SL1" : [],
        "SL3" : ['gene-99775', 'gene-40274', 'gene-304827', 'gene-92346', 'gene-306335', 'gene-285669', 'gene-120763', 'gene-2286', 'gene-97407', 'gene-232392', 'gene-328941', 'gene-166511', 'gene-39692', 'gene-384091', 'gene-74686', 'gene-122220', 'gene-218723', 'gene-414353', 'gene-312890', 'gene-153482', 'gene-39770', 'gene-132340', 'gene-253632', 'gene-378608', 'gene-206556', 'gene-336703', 'gene-21229', 'gene-166391', 'gene-120784', 'gene-87502', 'gene-317372', 'gene-73253', 'gene-211196', 'gene-9548', 'gene-60190', 'gene-234650', 'gene-410057', 'gene-121262', 'gene-100036', 'gene-227137', 'gene-75744', 'gene-279912', 'gene-343203', 'gene-233901', 'gene-163028', 'gene-238407', 'gene-39680', 'gene-198700', 'gene-231228', 'gene-410209', 'gene-143368', 'gene-388261', 'gene-30328', 'gene-226944', 'gene-334263', 'gene-130081', 'gene-47823', 'gene-228519', 'gene-350792', 'gene-277340', 'gene-182683', 'gene-60151', 'gene-206576', 'gene-202718', 'gene-130096', 'gene-377275', 'gene-244780', 'gene-288834', 'gene-189246', 'gene-48535', 'gene-205011', 'gene-253157', 'gene-62891', 'gene-69698']
    }

    ############################################
    ######### MAKE ALL THE SMEAR PLOTS #########
    ############################################
    if True:
        for separation, seps_dict in table_paths.items():
            print(f"\n=========================== {separation} ===========================")
            lists = {}
            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")
                numbers = {}
                lists[category] = {}
                for contrast, table_path in paths_dict.items():
                    print(f"{separation}:{category} --> contrast: {contrast}")
                    table_name = table_path.split("/")[-1].replace(".txt", "").replace("DE_genes_", "")
                    smear_name = f"{out_path_figs}/smear_{table_name}.png"
                    smear_title = contrast_plot_titles[contrast]

                    excl_list = []
                    excl_list_name = smear_title.replace(" ","")
                    if excl_list_name in excl_line_bias_lists:
                        excl_list = excl_line_bias_lists[excl_list_name]
                        print(f"\texcluding genes from list '{excl_list_name}'")
        
                    smear_lists = plot_smear(table_path=table_path, contrast=contrast, smear_plot_name=smear_name, title = smear_title, excl_genes_list=excl_list, x_axis="fdr_p")
                    # Downlist = smear_lists["Downregulated"]
                    # Uplist = smear_lists["Upregulated"]
                    # print(f"\t * Downregulated : {Downlist}\n\t * Upregulated : {Uplist}")
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
            }
        }
        for separation, seps_dict in table_paths.items():
            print(f"\n=========================== {separation} ===========================")
            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")

                for venn_cat, venn_contrasts_list in venn_sets[separation][category].items():
                    print(f"{venn_cat} : {venn_contrasts_list}")
                    venn_paths_dict = {contrast_plot_titles[contrast] : paths_dict[contrast] for contrast in venn_contrasts_list}
                    venn_filename_ = venn_cat.replace(" ", "_")
                    venn_filename = f"{out_path_figs}/Venn_{category}_{venn_filename_}.png"
                    venn_title = f"sig. DE genes overlap ({category})\n{venn_cat}"
                    plot_venn_DE_genes(venn_paths_dict, venn_filename=venn_filename, venn_title=venn_title)

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

            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")

                for LFC_cat, LFC_contrasts_list in LFC_comp_sets[separation][category].items():
                    print(f"{LFC_cat} : {LFC_contrasts_list}")
                    LFC_filename_ = LFC_cat.replace(" ", "_").replace(".", "")
                    LFC_filename = f"{out_path_figs}/LFC_scatter_{category}_{LFC_filename_}.png"

                    LFC_paths_dict = {contrast_plot_titles[contrast] : paths_dict[contrast] for contrast in LFC_contrasts_list}
                    excl_geneIDs = get_excl_genes_list(contrasts = LFC_contrasts_list, contrast_plot_titles=contrast_plot_titles, excl_line_bias_lists=excl_line_bias_lists)
                    
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