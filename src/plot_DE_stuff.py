"""
reproduce plots from R since there are contradictions. 
I will use the tables in data made with the topTags() function with no filtering so it just shows all expressed genes in the comparison
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3


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
            },
            "males" : {
                "SL1_14 - SL3_14" : f"{tables_dir}DE_genes_M_1-3_day14.txt",
                "SL1_16 - SL3_16" : f"{tables_dir}DE_genes_M_1-3_day16.txt",
                "SL1_18 - SL3_18" : f"{tables_dir}DE_genes_M_1-3_day18.txt",
                "SL1_18 - (SL1_14+SL1_16)/2" : f"{tables_dir}DE_genes_M_SL1.txt",
                "SL3_18 - (SL3_14+SL3_16)/2" : f"{tables_dir}DE_genes_M_SL3.txt",
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
        "SL1_14 - SL3_14" : f"day 14",
        "SL1_16 - SL3_16" : f"day 16",
        "SL1_18 - SL3_18" : f"day 18",
        "SL1_18 - (SL1_14+SL1_16)/2" : f"SL1",
        "SL3_18 - (SL3_14+SL3_16)/2" : f"SL3",
        "F_14 - M_14" : f"day 14",
        "F_16 - M_16" : f"day 16",
        "F_18 - M_18" : f"day 18",
        "F_14 - M_14" : f"day 14",
        "F_16 - M_16" : f"day 16",
    }
    return out_dict,contrast_plot_titles


def plot_smear(table_path, contrast, smear_plot_name = "smear_plot.png", p_sig = 0.05, min_LFC = 0, title = "significant logFC"):
    """
    replicate smear-plot from R, logCPM by logFC, significance highlighted in red
    return the number of up/downregulated and no difference genes
    """
    df = pd.read_csv(table_path, sep="\t", skiprows=0)
    
    cols = {"nonsig" : "#243742", "sig" : "#BD351E"}

    fig, ax = plt.subplots(1,1, figsize=(18, 10)) 
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

    ax.scatter(df_nonsig["logCPM"], df_nonsig["logFC"], color = cols["nonsig"], alpha=0.5, s=ps)
    ax.scatter(df_sig["logCPM"], df_sig["logFC"], color = cols["sig"], alpha=1, s=ps)
    
    if "SL" in contrast and "/2" not in contrast:
        label_contrast = contrast.replace("_14", "").replace("_16", "").replace("_18", "")
    elif "F" in contrast:
        label_contrast = contrast.replace("_14", "").replace("_16", "").replace("_18", "")
    else:
        label_contrast = contrast.replace("SL1", "day").replace("SL3", "day")
    ax.set_ylabel(f"logFC ({label_contrast})", fontsize = fs)
    ax.set_xlabel(f"log CPM", fontsize = fs)
    ax.tick_params(axis='x', labelsize=fs*0.9)
    ax.tick_params(axis='y', labelsize=fs*0.9)
    ax.set_title(title, fontsize = fs*1.25)

    min_line = min(df["logCPM"])
    max_line = max(df["logCPM"])
    if min_LFC > 0:
        ax.hlines(y=min_LFC, xmin=min_line, xmax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
        ax.hlines(y=-1*min_LFC, xmin=min_line, xmax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
    else:
        ax.hlines(y=1, xmin=min_line, xmax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
        ax.hlines(y=-1, xmin=min_line, xmax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
    ax.set_xlim([min_line-0.25,max_line+0.25])
    
    plt.tight_layout()
    plt.savefig(smear_plot_name, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {smear_plot_name}")

    # numbers for table
    upreg = df_sig.loc[df_sig['logFC'] > 0]
    upreg = upreg.shape[0]
    downreg = df_sig.loc[df_sig['logFC'] < 0]
    downreg = downreg.shape[0]
    all_rows = df.shape[0]
    nodiff = all_rows - upreg - downreg

    out_dict = {"Downregulated" : downreg, "no difference" : nodiff,  "Upregulated" : upreg}

    return out_dict


def plot_venn_DE_genes(tables_dict:dict, p_sig = 0.05, min_LFC = 0, venn_filename = "venn_diagram.png", venn_title = ""):
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


def plot_sig_LFC_overlap(tables_dict:dict, p_sig = 0.05, min_LFC = 0, LFC_filename = "sig_LFC_scatter.png", LFC_title = ""):
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
            ax.scatter(x,y,color = colors_dict[cat], s=ps, alpha = 0.75)

    ax.set_xlabel(f"logFC {table_a}", fontsize = fs)
    ax.set_ylabel(f"logFC {table_b}", fontsize = fs)
    ax.tick_params(axis='x', labelsize=fs*0.9)
    ax.tick_params(axis='y', labelsize=fs*0.9)
    ax.set_title(LFC_title, fontsize = fs*1.25)

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
    
    ax.set_ylim([min_yline,max_yline])
    ax.set_xlim([min_xline,max_xline])

    plt.legend(fontsize = fs*0.8, title ="gene sig. in", title_fontsize = fs*0.8)
    plt.tight_layout()
    plt.savefig(LFC_filename, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {LFC_filename}")

    lengths = { key : len(val) for key,val in lists.items()}
    return(lengths)


if __name__ == "__main__":

    username = "miltr339"
    table_paths,contrast_plot_titles = get_tables(username=username)
    out_path_figs = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures_python"
    

    ############################################
    ######### MAKE ALL THE SMEAR PLOTS #########
    ############################################
    if False:
        for separation, seps_dict in table_paths.items():
            print(f"\n=========================== {separation} ===========================")
            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")

                numbers = {}
                for contrast, table_path in paths_dict.items():
                    print(f"{separation}:{category} --> contrast: {contrast}")
                    table_name = table_path.split("/")[-1].replace(".txt", "").replace("DE_genes_", "")
                    smear_name = f"{out_path_figs}/smear_{table_name}.png"
                    smear_title = contrast_plot_titles[contrast]
        
                    smear_nums = plot_smear(table_path=table_path, contrast=contrast, smear_plot_name=smear_name, title = smear_title)
                    numbers[smear_title] = smear_nums
        
                print(numbers)

    ############################################
    ######## MAKE ALL THE VENN DIAGRAMS ########
    ############################################
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

    ############################################
    ######## MAKE ALL THE SB COMPARISONS #######
    ############################################

    if True:
        LFC_comp_sets = {
            "sex_separated" : {
                "females" : {
                    "day 14 vs. day 16 line bias" : ["SL1_14 - SL3_14","SL1_16 - SL3_16"],
                    "day 14 vs. day 18 line bias" : ["SL1_14 - SL3_14","SL1_18 - SL3_18"],
                    "day 16 vs. day 16 line bias" : ["SL1_16 - SL3_16","SL1_18 - SL3_18"],
                    "line by age bias" : ["SL1_18 - (SL1_14+SL1_16)/2","SL3_18 - (SL3_14+SL3_16)/2"]
                },
                "males" : {
                    "day 14 vs. day 16 line bias" : ["SL1_14 - SL3_14","SL1_16 - SL3_16"],
                    "day 14 vs. day 18 line bias" : ["SL1_14 - SL3_14","SL1_18 - SL3_18"],
                    "day 16 vs. day 16 line bias" : ["SL1_16 - SL3_16","SL1_18 - SL3_18"],
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

        for separation, seps_dict in table_paths.items():
            print(f"\n=========================== {separation} ===========================")
            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")

                for LFC_cat, LFC_contrasts_list in LFC_comp_sets[separation][category].items():
                    print(f"{LFC_cat} : {LFC_contrasts_list}")
                    LFC_filename_ = LFC_cat.replace(" ", "_").replace(".", "")
                    LFC_filename = f"{out_path_figs}/LFC_scatter_{category}_{LFC_filename_}.png"
                    LFC_paths_dict = {contrast_plot_titles[contrast] : paths_dict[contrast] for contrast in LFC_contrasts_list}

                    if "line bias" in LFC_cat:
                        title = f"{category}: contrast SL1 - SL3"
                    elif "sex bias" in LFC_cat:
                        title = f"{category}: female - male"
                    elif "age bias" in LFC_cat:
                        title = f"{category}: day 18 - mean(day 14, day 16)"
                    else:
                        title = LFC_cat

                    numbers = plot_sig_LFC_overlap(LFC_paths_dict, LFC_filename = LFC_filename, LFC_title = title)
                    print(f"\t{numbers}")