"""
reproduce plots from R since there are contradictions. 
I will use the tables in data made with the topTags() function with no filtering so it just shows all expressed genes in the comparison
"""

import pandas as pd
import matplotlib.pyplot as plt


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




if __name__ == "__main__":

    username = "miltr339"
    table_paths,contrast_plot_titles = get_tables(username=username)
    out_path_figs = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures_python"
    
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