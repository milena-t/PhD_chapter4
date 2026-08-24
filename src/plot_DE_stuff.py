"""
reproduce plots from R since there are contradictions. 
I will use the tables in data made with the topTags() function with no filtering so it just shows all expressed genes in the comparison
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib_venn import venn2,venn3
import math
import scipy.stats as sts
import warnings
import upsetplot
import Y_expression_quantification as time_series_plots

def get_tables(username = "miltr339"):
    """
    tables are either split by sex, so that the line and day contrasts are made on a subset that is only males or only females,
    or split by line so that the day and sex contrasts are made on a subset of only one line at a time.
    """

    tables_dir = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/"

    out_dict = {

        "no_separation" : {
            "day_random" : {
                "SL1_F - SL1_M" : f"{tables_dir}full_dataset_day_ignored_SL1_F-M.txt",
                "SL3_F - SL3_M" : f"{tables_dir}full_dataset_day_ignored_SL3_F-M.txt",
                "SL1_F - SL3_F" : f"{tables_dir}full_dataset_day_ignored_F_SL1-SL3.txt",
                "SL1_M - SL3_M" : f"{tables_dir}full_dataset_day_ignored_M_SL1-SL3.txt",
                "(SL1_F - SL3_F) - (SL1_M - SL3_M)" : f"{tables_dir}full_dataset_day_ignored_line_sex_interaction.txt",
                "day16,day18" : f"{tables_dir}full_dataset_day_random_factor_days_diff.txt",
            },
            "line_random" : {
                "day14_F - day14_M" : f"{tables_dir}full_dataset_line_ignored_day14_F-M.txt",
                "day16_F - day16_M" : f"{tables_dir}full_dataset_line_ignored_day16_F-M.txt",
                "day18_F - day18_M" : f"{tables_dir}full_dataset_line_ignored_day18_F-M.txt",
                "day14_M - day16_M" : f"{tables_dir}full_dataset_line_ignored_M_14_16.txt",
                "day16_M - day18_M" : f"{tables_dir}full_dataset_line_ignored_M_16_18.txt",
                "day14_F - day16_F" : f"{tables_dir}full_dataset_line_ignored_F_14_16.txt",
                "day16_F - day18_F" : f"{tables_dir}full_dataset_line_ignored_F_16_18.txt",
                "(day14_F - day14_M) - (day16_F - day16_M) - (day18_F - day18_M)" : f"{tables_dir}full_dataset_line_ignored_day_sex_interaction.txt",
                "line3" : f"{tables_dir}full_dataset_line_random_factor_line_diff.txt",
            },
            "sex_random" : {
                "(SL1_d14 - SL3_d14) - (SL1_d16 - SL3_d16) - (SL1_d18 - SL3_d18)" : f"{tables_dir}full_dataset_line_ignored_day_line_interaction.txt",
            }
        },
        
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
                "SL1 - SL3" : f"{tables_dir}only_F_no_time_DE_genes_1-3.txt",
                "day14 - day16" : f"{tables_dir}only_F_no_line_DE_genes_day_14-16.txt",
                "day16 - day18" : f"{tables_dir}only_F_no_line_DE_genes_day_16-18.txt",
                "day14 - day18" : f"{tables_dir}only_F_no_line_DE_genes_day_14-18.txt",
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
                "SL1 - SL3" : f"{tables_dir}only_M_no_time_DE_genes_1-3.txt",
                "day14 - day16" : f"{tables_dir}only_M_no_line_DE_genes_day_14-16.txt",
                "day16 - day18" : f"{tables_dir}only_M_no_line_DE_genes_day_16-18.txt",
                "day14 - day18" : f"{tables_dir}only_M_no_line_DE_genes_day_14-18.txt",
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
        "SL1 - SL3" : "SL1-SL3 days merged",
        "SL1_F - SL1_M" : "SL1 F-M",
        "SL3_F - SL3_M" : "SL3 F-M",
        "SL1_F - SL3_F" : "SL1-SL3 F",
        "SL1_M - SL3_M" : "SL1-SL3 M",
        "(SL1_F - SL3_F) - (SL1_M - SL3_M)" : " line by sex interaction",
        "day16,day18" : "day effect",
        "day14_F - day14_M" : "day 14 sex bias",
        "day16_F - day16_M" : "day 16 sex bias",
        "day18_F - day18_M" : "day 18 sex bias",
        "(day14_F - day14_M) - (day16_F - day16_M) - (day18_F - day18_M)" : "day by sex interaction",
        "(SL1_d14 - SL3_d14) - (SL1_d16 - SL3_d16) - (SL1_d18 - SL3_d18)" : "day by line interaction",
        "line3" : "line effect",
        "day14 - day16" : "day 14 - 16",
        "day14 - day18" : "day 14 - 18",
        "day16 - day18" : "day 16 - 18",
        "day14_M - day16_M" : "males day 14 - 16",
        "day16_M - day18_M" : "males day 16 - 18",
        "day14_F - day16_F" : "females day 14 - 16",
        "day16_F - day18_F" : "females day 16 - 18",
    }
    return out_dict,contrast_plot_titles


def plot_smear(table_path, contrast, smear_plot_name = "smear_plot.png", p_sig = 0.05, min_LFC = 0, title = "significant logFC", excl_genes_list = [], x_axis = "logcpm", plot = True, coefficients=False, highlight_yTOR=False):
    """
    replicate smear-plot from R, logCPM by logFC, significance highlighted in red
    return the number of up/downregulated and no difference genes
    the x_axis can be set to logCPM which is a smear plot or to fdr_p which is a volcano plot
    if no plot is generated (plot=False) only the list of sig genes is returned
    """
    df = pd.read_csv(table_path, sep="\t", skiprows=0)
    if len(excl_genes_list) >0:
        old_len = df.shape[0]
        df = df.drop(index=excl_genes_list, errors='ignore') # only drop existing labels, ignore the rest
        new_len = df.shape[0]
        print(f"\t{len(excl_genes_list)} genes dropped from excl_genes_list (gene number {old_len} -> {new_len})")

    cols = {"nonsig" : "#243742", "sig" : "#BD351E"}

    fs = 55
    fig_width=18
    fig_height=15
    point_size_factor = 5
    ps = fs*point_size_factor # point size
    
    if x_axis == "fdr_p":
        df['FDR_volc'] = df['FDR'].transform(lambda x: math.log10(x))
    ### plot the unsignificant first so they are below and the significant ones are above
    # if no sex bias so no minLFC threshold
    df_nonsig = df.loc[df['FDR'] >= p_sig] 
    df_sig = df.loc[df['FDR'] < p_sig]
    print(f"\tnum sig genes: {df_sig.shape[0]}")
    
    # if sex bias and minLFC threshold
    if min_LFC>0 and coefficients==False:
        df_nonsig = df.loc[(df['FDR'] >= p_sig) | (abs(df['logFC']) < min_LFC)]
        df_sig = df_sig.loc[abs(df_sig['logFC']) >= min_LFC]
        print(f"\tnum sig genes with LFC > {min_LFC}: {df_sig.shape[0]}")
    
    if coefficients == False:
        num_up = df_sig.loc[df_sig['logFC'] > min_LFC].shape[0]
        num_down = df_sig.loc[df_sig['logFC'] < min_LFC].shape[0]
        label_sig = f"up:         {num_up}\ndown:    {num_down}"

    num_nonsig = df_nonsig.shape[0]
    label_nonsig = f"nonsig.: {num_nonsig}"
    
    if coefficients==True:
        # make individual plots for each coefficient
        logfc_colnames = [ coln for coln in df.columns if "logFC" in coln]
        sig_dfs_coeff = {logfc_colname : df_sig.loc[abs(df_sig[logfc_colname]) >= min_LFC] for logfc_colname in logfc_colnames}
        out_dict = {logfc_colname.split(".")[-1] : {} for logfc_colname in logfc_colnames}

        for logfc_colname, sig_df_coeff in sig_dfs_coeff.items():
            num_up = sig_df_coeff.loc[sig_df_coeff[logfc_colname] > min_LFC].shape[0]
            num_down = sig_df_coeff.loc[sig_df_coeff[logfc_colname] < min_LFC].shape[0]
            label_sig = f"up:         {num_up}\ndown:    {num_down}"

            coeff_name = logfc_colname.split(".")[-1]
            print(f"\tcoef:{coeff_name} num sig genes with LFC > {min_LFC}: {sig_df_coeff.shape[0]}")
            if  x_axis == "logcpm":
                fig, ax = plt.subplots(1,1, figsize=(fig_width, fig_height-5)) 
            elif x_axis == "fdr_p":
                fig, ax = plt.subplots(1,1, figsize=(fig_width, fig_height)) 

            if  x_axis == "logcpm":
                ax.scatter(df_nonsig["logCPM"], df_nonsig[logfc_colname], color = cols["nonsig"], alpha=0.4, s=ps)
                ax.scatter(df_sig["logCPM"], df_sig[logfc_colname], color = cols["sig"], alpha=0.4, s=ps)
            elif x_axis == "fdr_p":
                ax.scatter(df_nonsig[logfc_colname], df_nonsig["FDR_volc"], color = cols["nonsig"], alpha=0.4, s=ps, label=label_nonsig)
                ax.scatter(df_sig[logfc_colname], df_sig["FDR_volc"], color = cols["sig"], alpha=0.4, s=ps, label=label_sig)
                # if Y-TOR is significant highlight it
                if highlight_yTOR:
                    y_tor_sig_test = df_sig[df_sig["Gene"] == "yTor-all"]
                    if not y_tor_sig_test.empty:
                        y_tor_sig = y_tor_sig_test.iloc[[0]]   # double brackets keep it as a DataFrame, not a Series
                        ax.scatter(y_tor_sig["logFC"], y_tor_sig["FDR_volc"], color = "#599BB7", alpha=1, s=ps, label="y-TOR")  
        
            if  x_axis == "logcpm":
                ax.set_ylabel(f"{logfc_colname}", fontsize = fs)
                ax.set_xlabel(f"log CPM", fontsize = fs)
            elif x_axis == "fdr_p":
                ax.set_xlabel(f"{logfc_colname}", fontsize = fs)
                ax.set_ylabel(f"log10  FDR p-value", fontsize = fs)
                ax.yaxis.set_inverted(True) 
            ax.tick_params(axis='x', labelsize=fs*0.9)
            ax.tick_params(axis='y', labelsize=fs*0.9)
            coeff_title = f"{title}:{coeff_name}"
            ax.set_title(coeff_title, fontsize = fs*1.25)

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
                min_line = min(df["FDR_volc"])-0.25
                max_line = max(df["FDR_volc"])+0.25
                if min_LFC > 0:
                    ax.vlines(x=min_LFC, ymin=min_line, ymax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
                    ax.vlines(x=-1*min_LFC, ymin=min_line, ymax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
                else:
                    ax.vlines(x=1, ymin=min_line, ymax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
                    ax.vlines(x=-1, ymin=min_line, ymax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")

            smear_plot_name_coeff = smear_plot_name.replace(".png", f"_{coeff_name}.png")
            if x_axis == "fdr_p":
                smear_plot_name_coeff = smear_plot_name.replace(".png", f"_{coeff_name}_volcano.png")
            
            if x_axis == "fdr_p":
                plt.legend(fontsize = fs, loc='upper center') #, title ="gene counts", title_fontsize = fs*0.7
            # layout rect=(left, bottom, right, top)
            plt.tight_layout()# rect=[0.0, 0.05, 1, 1])
            fig.subplots_adjust(left=0.18) # or whatever
            plt.savefig(smear_plot_name_coeff, dpi = 300, transparent = True)
            print(f"plot saved in current working directory as: {smear_plot_name_coeff}")
            plt.clf()
            plt.cla()
            plt.close()

            # numbers for table
            upreg = df_sig.loc[df_sig[logfc_colname] > 0]
            # upreg = upreg.shape[0]
            upreg = upreg["Gene"].tolist()
            downreg = df_sig.loc[df_sig[logfc_colname] < 0]
            # downreg = downreg.shape[0]
            downreg = downreg["Gene"].tolist()
            # all_rows = df.shape[0]
            # nodiff = all_rows - upreg - downreg
            df_nonsig = df.drop(index=downreg, errors='ignore')
            df_nonsig = df_nonsig.drop(index=upreg, errors='ignore')
            nodiff = df_nonsig["Gene"].tolist()

            out_dict[coeff_name] = {"Downregulated" : downreg, "no difference" : nodiff,  "Upregulated" : upreg}

    else:    
        if plot:
            if  x_axis == "logcpm":
                fig, ax = plt.subplots(1,1, figsize=(fig_width, fig_height-5))
            elif x_axis == "fdr_p":
                fig, ax = plt.subplots(1,1, figsize=(fig_width, fig_height))

            if  x_axis == "logcpm":
                ax.scatter(df_nonsig["logCPM"], df_nonsig["logFC"], color = cols["nonsig"], alpha=0.5, s=ps)
                ax.scatter(df_sig["logCPM"], df_sig["logFC"], color = cols["sig"], alpha=1, s=ps)
            elif x_axis == "fdr_p":
                ax.scatter(df_nonsig["logFC"], df_nonsig["FDR_volc"], color = cols["nonsig"], alpha=0.5, s=ps, label=label_nonsig)
                ax.scatter(df_sig["logFC"], df_sig["FDR_volc"], color = cols["sig"], alpha=1, s=ps, label=label_sig)
                # if Y-TOR is significant highlight it
                if highlight_yTOR:
                    y_tor_sig_test = df_sig[df_sig["Gene"] == "yTor-all"]
                    if not y_tor_sig_test.empty:
                        y_tor_sig = y_tor_sig_test.iloc[[0]]   # double brackets keep it as a DataFrame, not a Series
                        ax.scatter(y_tor_sig["logFC"], y_tor_sig["FDR_volc"], color = "#599BB7", alpha=1, s=ps, label="y-TOR")  


            
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
                ax.set_ylabel(f"log2FC ({label_contrast})", fontsize = fs)
                ax.set_xlabel(f"log CPM", fontsize = fs)
            elif x_axis == "fdr_p":
                ax.set_xlabel(f"logFC", fontsize = fs)
                ax.set_ylabel(f"log10  FDR p-value", fontsize = fs)
                ax.yaxis.set_inverted(True) 
            ax.tick_params(axis='x', labelsize=fs*0.9)
            ax.tick_params(axis='y', labelsize=fs*0.9)

            title_ = title.replace("SL1-3", "line-bias")
            title_ = title_.replace("SL1-SL3", "line-bias")
            title_ = title_.replace("SL1", "small-Y").replace("SL3", "large-Y")
            print(f"\t plot title: {title_}")
            ax.set_title(title_, fontsize = fs*1.25)

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
                min_line = min(df["FDR_volc"])-0.25
                max_line = max(df["FDR_volc"])+0.1
                if min_LFC > 0:
                    ax.vlines(x=min_LFC, ymin=min_line, ymax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
                    ax.vlines(x=-1*min_LFC, ymin=min_line, ymax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
                else:
                    ax.vlines(x=1, ymin=min_line, ymax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")
                    ax.vlines(x=-1, ymin=min_line, ymax=max_line, linewidth=point_size_factor, linestyle = ":", color="#818181")            
                ax.set_ylim([max_line, min_line])

            if x_axis == "fdr_p":
                smear_plot_name = smear_plot_name.replace(".png", "_volcano.png")
                ax.yaxis.set_major_locator(MaxNLocator(integer=True)) # force y axis as integers to make the y axis label visible and not outside of bounds
            if highlight_yTOR:
                smear_plot_name = smear_plot_name.replace(".png", "_y_TOR.png")

            if x_axis == "fdr_p":
                plt.legend(fontsize = fs, loc='upper center') #, title ="gene counts", title_fontsize = fs*0.7
            plt.tight_layout()
            fig.subplots_adjust(left=0.14) # or whatever
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


def plot_sig_LFC_overlap(tables_dict:dict, p_sig = 0.05, min_LFC = 0, LFC_filename = "sig_LFC_scatter.png", LFC_title = "", excl_geneIDs =[], incl_geneIDs = [], intersection_nums = False):
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
            excl_geneIDs = []
            # full_size = df.shape[0]
            # print(f"\t{table_title} : only including {len(incl_geneIDs)} from {full_size}")
            # df = df[df["Gene"].isin(incl_geneIDs)]
            # # if including only a subset of genes don't filter for only significant ones, include everything
        
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
        table_a : list(set_a - set_b),
        table_b : list(set_b - set_a),
        "shared" : list(set_a & set_b),
        }

    fig, ax = plt.subplots(1,1, figsize=(13, 13)) 
    fs = 45
    point_size_factor = 8
    ps = fs*point_size_factor # point size

    sig_colors = ["#BD351E","#EA882C"] # red , orange
    colors_dict = {table_title : sig_colors[i] for i,table_title in enumerate(tables_dict.keys())}
    colors_dict["shared"] = "#3C7FA7" # blue
    colors_dict["neither"] = "#4B3B47" # mauve shadow

    if len(incl_geneIDs) == 0:
        excl_counter = { cat : 0 for cat in [table_a,table_b,"shared"]}
        all_sig = lists['shared']+lists[table_a]+lists[table_b]
    else:
        excl_counter = { cat : 0 for cat in [table_a,table_b,"shared","neither"]}
        all_sig = lists['shared']+lists[table_a]+lists[table_b]
        nonsig_cat_list = [geneID for geneID in incl_geneIDs if geneID not in all_sig]
        if len(nonsig_cat_list) > 0 and incl_geneIDs != [None]:
            lists["neither"] = nonsig_cat_list
        elif incl_geneIDs == [None]:
            lists["neither"] = []

        if intersection_nums:
            # make a list of the genes that are sig. in both in the format required for the GO enrichment
            sig_list_outfile = LFC_filename.split("/")[-1].replace(".png", "_shared_sig_DE_list.txt")
            sig_list_outfile = f"{lists_outdir}/{sig_list_outfile}"
            with open(sig_list_outfile, "w") as sig_list_file:
                DE_list_outfile = [f"{geneID},1" for geneID in lists['shared']]
                DE_string = "\n".join(DE_list_outfile)
                sig_list_file.write(f"geneID,sig_DE\n{DE_string}\n") 
                singleDE_list_outfile = [f"{geneID},0" for geneID in lists[table_a]+lists[table_b]+nonsig_cat_list]
                singleDE_string = "\n".join(singleDE_list_outfile)
                sig_list_file.write(f"{singleDE_string}\n") # needs the newline character so that R can read the list right

            print(f" * list of sig DE genes written to: {sig_list_outfile}")

    nonsig_incl = 0
    sig_incl = 0
    count_points = 0
    for cat in lists.keys():
        if len(lists[cat])>50 or cat!="shared":
            print(f"\t* {cat} ({len(lists[cat])})")
        else:
            print(f"\t* {cat} ({len(lists[cat])})\n\t{lists[cat]}")

        for geneID in lists[cat]:
            try:
                y = tables_df_all[table_a].loc[geneID,"logFC"]
            except:
                y = 0
                print(geneID)
            try:
                x = tables_df_all[table_b].loc[geneID,"logFC"]
            except:
                x = 0
                print(geneID)
            if len(excl_geneIDs)>0:
                if geneID in excl_geneIDs:
                    excl_counter[cat]+=1
                    ax.scatter(x,y,color = colors_dict[cat], s=ps*1.5, alpha = 1, marker="1")
                else:
                    ax.scatter(x,y,color = colors_dict[cat], s=ps, alpha = 0.75)
                count_points += 1
            elif len(incl_geneIDs)>0:
                if geneID not in incl_geneIDs:
                    excl_counter[cat]+=1
                    ax.scatter(x,y,color = colors_dict[cat], s=ps*1.5, alpha = 1, marker="1")
                    nonsig_incl+=1
                else:
                    ax.scatter(x,y,color = colors_dict[cat], s=ps, alpha = 0.75)
                    sig_incl+=1
                count_points += 1
            else:
                ax.scatter(x,y,color = colors_dict[cat], s=ps, alpha = 0.75)
                count_points += 1
            

    print(f" points with coords: {count_points}")
    label_a = table_a.replace("SL1-3", "line-bias")
    label_a = label_a.replace("SL1-SL3", "line-bias")
    label_a = label_a.replace("F-M", "sex-bias")
    label_a = label_a.replace("SL1", "small-Y").replace("SL3", "large-Y")
    ax.set_ylabel(f"logFC {label_a}", fontsize = fs)
    label_b = table_b.replace("SL1-3", "line-bias")
    label_b = label_b.replace("SL1-SL3", "line-bias")
    label_b = label_b.replace("F-M", "sex-bias")
    label_b = label_b.replace("SL1", "small-Y").replace("SL3", "large-Y")
    ax.set_xlabel(f"logFC {label_b}", fontsize = fs)
    ax.tick_params(axis='x', labelsize=fs*0.9)
    ax.tick_params(axis='y', labelsize=fs*0.9)
    ax.set_title(LFC_title, fontsize = fs*0.95)

    min_yline,max_yline = ax.get_ylim()
    min_xline,max_xline = ax.get_xlim()
    min_yline = min_yline-0.5
    max_yline = max_yline+0.5
    min_xline = min_xline-0.5
    max_xline = max_xline+0.5
    if incl_geneIDs != []:
        # min_yline = min_yline*1.2
        max_xline = max_xline*1.2
        min_xline = min_xline*1.4

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

    if excl_geneIDs != []:
        if "females" in LFC_filename:
            shared_lab = "also DE\nin males"
        else:
            shared_lab = "also DE\nin females"
        int_sig = ax.scatter(1000,1000,color = "#666666", s=ps*1.5, alpha = 1, label = shared_lab, marker="1")
        int_legend="No"

    elif incl_geneIDs != []:
        int_nonsig = ax.scatter(1000,1000,color = "#666666", s=ps*1.5, alpha = 1, label = f"nonsig. ({nonsig_incl})", marker="1")
        int_sig = ax.scatter(1000,1000,color = "#666666", s=ps*1.5, alpha = 1, label = f"significant ({sig_incl})", marker="o")
        int_legend = plt.legend(handles = [int_sig, int_nonsig], fontsize = fs*0.75, title ="line-by-sex\ninteraction", title_fontsize = fs*0.7, loc='upper left')
    else:
        int_legend="No"

    ## make legend points
    main_sig = []
    if int_legend!="No":
        markertype="s"
    if int_legend=="No":
        markertype="o"
    # for cat in reversed(list(lists.keys())):
    for cat in list(lists.keys()):
        legend_label = cat.replace("SL1-3", "line-bias")
        legend_label = legend_label.replace("SL1-SL3", "line-bias")
        legend_label = legend_label.replace("females", "F")
        legend_label = legend_label.replace("males", "M")
        legend_label = legend_label.replace("shared", "both")
        legend_label = legend_label.replace("SL1", "small-Y").replace("SL3", "large-Y")
        if intersection_nums and legend_label=="both":
            legend_label = f"both ({len(lists['shared'])})"
        main_sig_ind = ax.scatter(1000,1000,color = colors_dict[cat], s=ps, alpha = 0.75, label = legend_label, marker=markertype)
        main_sig.append(main_sig_ind)

    ax.set_ylim([min_yline,max_yline])
    ax.set_xlim([min_xline,max_xline])

    if int_legend!="No":
        main_sig_legend = plt.legend(handles = main_sig, fontsize = fs*0.75, title ="gene sig. in\nmain effect", title_fontsize = fs*0.7, loc='lower right')
        ax.add_artist(int_legend)
    if int_legend=="No":
        main_sig_legend = plt.legend(handles = main_sig, fontsize = fs*0.75, title ="gene sig. in\nmain effect", title_fontsize = fs*0.7, loc='lower right')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True)) # force y axis as integers to make the y axis label visible and not outside of bounds

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


def plot_sig_gene_overlap(geneIDs_dict:dict, plot_filename:str, plot_title  = "", min_overlap = 0):
    """
    make an upseetplot of sig DE genes within each data separation
    """
    data = upsetplot.from_contents(geneIDs_dict)
    if min_overlap>0:
        upsetplot.UpSet(data, subset_size="count", show_counts=True, min_subset_size=min_overlap, sort_by="cardinality", sort_categories_by="input").plot()
    else:
        upsetplot.UpSet(data, subset_size="count", sort_by="cardinality", show_counts=True, min_subset_size=min_overlap).plot()

    if len(plot_title)>0:
        plt.title(plot_title)
    plt.tight_layout()
    plt.savefig(plot_filename, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {plot_filename}")
    plt.clf()
    plt.cla()
    plt.close()
    return data


def plot_logFC_boxplots(infiles_dict, p_sig = 0.05, min_LFC = -1, only_all_intersection = False, plot_filename  = "LogFC_boxplot.png", plot_title = ""):
    """
    Make boxplot of abs LogFC of sig DE genes for all files specified in infiles_dict with the dict keys as axis labels
    """
    tables_list = []
    xtick_labels = []
    # only use genes present in all tables
    all_genes = []
    if only_all_intersection:
        for name,table_path in infiles_dict.items():
            df = pd.read_csv(table_path, sep="\t", skiprows=0)
            df_sig = df.loc[df['FDR'] < p_sig]
            if min_LFC!=0:
                # df_sig = df_sig.loc[abs(df_sig['logFC']) >= min_LFC]
                df_sig = df_sig.loc[df_sig['logFC'] < min_LFC] ### only male-bias, min is -1 and no abs()
            genes = df_sig["Gene"].tolist()
            print(f"\t{name} --> {len(genes)} sig. male-biased genes")
            all_genes.append(set(genes))
    
    intersection_genes = all_genes[0] & all_genes[1] & all_genes[2] & all_genes[3] & all_genes[4] & all_genes[5]
    print(f"{len(intersection_genes)} shared sig genes")

    for name,table_path in infiles_dict.items():
        df = pd.read_csv(table_path, sep="\t", skiprows=0)
        df_sig = df.loc[df['FDR'] < p_sig]
        if min_LFC!=0:
            # df_sig = df_sig.loc[abs(df_sig['logFC']) >= min_LFC]
            df_sig = df_sig.loc[df_sig['logFC'] < min_LFC]
        if len(intersection_genes)>0:
            before_downsample = df_sig.shape[0]
            df_sig = df_sig[df_sig['Gene'].isin(intersection_genes)]
            print(f"filter {before_downsample} -> {df_sig.shape[0]}")
        data = df_sig["logFC"].tolist()
        tables_list.append(data)
        name_ = name.replace("SL1", "small-Y").replace("SL3", "large-Y").replace(":", f"\n")
        if only_all_intersection:
            xtick_labels.append(f"{name_}")
        else:
            xtick_labels.append(f"{name_}\n({len(data)})")
        # sig_geneIDs_lists[name] = df_sig["Gene"].tolist()
    tick_pos = range(len(xtick_labels))

    colors_dict = {
        "fill" : "#F2933A", # uniform_filtered orange
        "edge" : "#C36711", # darker orange
        "medians" : "#FFBB7C", # lighter orange
        "X_fill" : "#b82946", # native red
        "X_edge" : "#861D32", #dark red
        "X_medians" : "#D86A80" # light red
    }
    
    fs = 50 # font size
    lw = 6
    width = 0.7
    # set figure aspect ratio
    aspect_ratio = 14 / 8
    height_pixels = 1500  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels
    fig, ax = plt.subplots(1,1,figsize=(width_pixels / 100, height_pixels / 100), dpi=100)

    bp = ax.boxplot(tables_list, positions=tick_pos, widths=width, patch_artist=True)
    if min_LFC>0:
        ax.axhline(y=0, color='#B78F85', linestyle='--', linewidth=lw)

    ax.tick_params(axis='y', labelsize=fs*0.9)
    ax.tick_params(axis='x', labelsize=fs*0.9)
    ax.set_xticklabels(xtick_labels, fontsize=fs)
    ax.set_ylabel("log2FC (F - M)", fontsize = fs)#, x=0.0, y=0.625)
    if len(plot_title)>0:
        ax.set_title(plot_title, fontsize=fs*1.2)#, fontstyle='italic')

    ## modify boxplot colors
    if True:
        for i, box in enumerate(bp['boxes']):
            if i%2==0:
                box.set(facecolor=colors_dict["fill"], edgecolor=colors_dict["edge"], linewidth=lw)
            else:
                box.set(facecolor=colors_dict["X_fill"], edgecolor=colors_dict["X_edge"], linewidth=lw)
        for i, median in enumerate(bp['medians']):
            if i%2==0:
                median.set(color=colors_dict['medians'], linewidth=lw)
            else:
                median.set(color=colors_dict['X_medians'], linewidth=lw)
        for i, whisker in enumerate(bp['whiskers']):
            # print(f"whisker: {i}")
            if i//2 % 2==0:
                whisker.set(color=colors_dict['edge'], linestyle='-',linewidth=lw)
            else:
                whisker.set(color=colors_dict['X_edge'], linestyle='-',linewidth=lw)
        for i, cap in enumerate(bp['caps']):
            if i//2 % 2==0:
                cap.set(color=colors_dict['edge'],linewidth=lw)
            else:
                cap.set(color=colors_dict['X_edge'],linewidth=lw)
        for i, flier in enumerate(bp['fliers']):
            if i%2==0:
                flier.set(marker='.', markersize = lw*6, markerfacecolor=colors_dict['edge'], markeredgecolor=colors_dict['edge'])
            else:
                flier.set(marker='.', markersize = lw*6, markerfacecolor=colors_dict['X_edge'], markeredgecolor=colors_dict['X_edge'])

    ## make statistical annotation
    if True:
        def add_significance_bar_log(ax, x1, x2, data, y, factor=1.1, color='black', lw=lw, fs=fs):
            """
            This function was modified from one created by claude code
            """
            # run the test
            data1 = data[x1]
            data2 = data[x2]
            stat, p = sts.mannwhitneyu(data1, data2, alternative='two-sided')

            # convert p-value to stars
            text = 'ns'
            color="#8e8e8e" #light grey
            tick_top = y + factor*0.7
            if p < 0.05:
                text = '*'
                color="#343434"# darker grey
                ax.text((x1+x2)/2, tick_top*0.975, text, ha='center', va='bottom', color=color, fontsize=fs)
                print(f"  * Mann-Whitney {xtick_labels[x1]}-{xtick_labels[x2]} : \t U statistic: {stat}, p-value: {p}")
            else:
                print(f"    Mann-Whitney {xtick_labels[x1]}-{xtick_labels[x2]} : \t U statistic: {stat}, p-value: {p}")
                ax.text((x1+x2)/2, tick_top*0.975, text, ha='center', va='bottom', color=color, fontsize=fs)
            if y>0:
                bar_raise = 0
            else:
                bar_raise = 0
            ax.plot([x1, x1, x2, x2], [y+bar_raise, tick_top, tick_top, y+bar_raise], lw=lw, color=color)

        ymax = max([max(box) for box in tables_list])

        add_significance_bar_log(ax=ax, x1=tick_pos[0], x2=tick_pos[1], data=tables_list, y=ymax+2, lw=lw, fs=fs)
        add_significance_bar_log(ax=ax, x1=tick_pos[2], x2=tick_pos[3], data=tables_list, y=ymax+2, lw=lw, fs=fs)
        add_significance_bar_log(ax=ax, x1=tick_pos[4], x2=tick_pos[5], data=tables_list, y=ymax+2, lw=lw, fs=fs)
        add_significance_bar_log(ax=ax, x1=tick_pos[0], x2=tick_pos[4], data=tables_list, y=ymax+4
        , lw=lw, fs=fs)
        add_significance_bar_log(ax=ax, x1=tick_pos[1], x2=tick_pos[5], data=tables_list, y=ymax+6, lw=lw, fs=fs)

    xmin_,xmax_ = ax.get_xlim()
    print(f"--------> {xmin_} to {xmax_}")
    ax.plot([xmin_, xmax_], [0, 0], color = "#7A6266", linestyle="dashed", linewidth=lw)

    plt.tight_layout()
    plt.savefig(plot_filename, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {plot_filename}")


def plot_sig_LFC_diff(tables_diff:list, table_LB:str, p_sig = 0.05, min_LFC = 0, LFC_filename = "sig_LFC_scatter.png", LFC_title = "", excl_geneIDs = [], intersection_nums = False):
    """ 
    plot a scatterplot of sig. DE genes with LFC values
    """
    if len(tables_diff) !=2 :
        raise RuntimeError(f"list should have only 2 elements but it has has length {len(tables_diff)}! \n{tables_diff}")
    
    df_SB1 = pd.read_csv(tables_diff[0], sep="\t", skiprows=0)
    df_SB1_sig = df_SB1.loc[df_SB1['FDR'] < p_sig]
    SB1_sig = df_SB1_sig["Gene"].tolist()
    df_SB1 = df_SB1.drop(columns=["logCPM","F","PValue","FDR"])

    df_SB3 = pd.read_csv(tables_diff[1], sep="\t", skiprows=0)
    df_SB3_sig = df_SB3.loc[df_SB3['FDR'] < p_sig]
    SB3_sig = df_SB3_sig["Gene"].tolist()
    df_SB3 = df_SB3.drop(columns=["logCPM","F","PValue","FDR"])
    
    df_sb = pd.merge(df_SB1,df_SB3, on = "Gene") # keeps only genes present in both, default inner join
    # df_sb["logFC_diff"] = abs(df_sb["logFC_x"]-df_sb["logFC_y"])
    df_sb["logFC_diff"] = df_sb["logFC_x"] - df_sb["logFC_y"]
    df_sb = df_sb.drop(columns=["logFC_x","logFC_y"])

    df_lb = pd.read_csv(table_LB, sep="\t", skiprows=0)
    df_lb = df_lb.loc[df_lb['FDR'] < p_sig] # only genes that are sig. line biased in males
    df_lb = df_lb.drop(columns=["logCPM","F","PValue","FDR"])
    
    df = pd.merge(df_sb,df_lb, on = "Gene")
    print(f"{df.shape[0]} geneIDs included (significant line bias in males)")
    print(f"{len(excl_geneIDs)} geneIDs excluded (significant line bias also in females)")
    
    fig, ax = plt.subplots(1,1, figsize=(13, 13)) 
    fs = 45
    point_size_factor = 8
    ps = fs*point_size_factor # point size

    colors_dict = {"small-Y SB" : "#BD351E" , "large-Y SB" : "#EA882C"}
    colors_dict["both SB"] = "#3C7FA7" # blue
    shared_IDs = {label : [] for label in colors_dict.keys()}
    colors_dict["neither"] = "#4B3B47" # mauve shadow
    colors_count = {label : 0 for label in colors_dict.keys()}

    excl_count = 0
    for geneID in df_lb["Gene"].tolist():
        if geneID in excl_geneIDs:
            excl_count+=1
        else:
            try:
                y = df.loc[df["Gene"]==geneID, "logFC_diff"]
            except:
                y = 0
                print(geneID)
            try:
                x = df.loc[df["Gene"]==geneID, "logFC"]
            except:
                x = 0
                print(geneID)

            if geneID in SB1_sig and geneID not in SB3_sig:
                c = colors_dict["small-Y SB"]
                colors_count["small-Y SB"]+=1
                shared_IDs["small-Y SB"].append(geneID)
            elif geneID not in SB1_sig and geneID in SB3_sig:
                c = colors_dict["large-Y SB"]
                colors_count["large-Y SB"]+=1
                shared_IDs["large-Y SB"].append(geneID)
            elif geneID in SB1_sig and geneID in SB3_sig:
                c = colors_dict["both SB"]
                colors_count["both SB"]+=1
                shared_IDs["both SB"].append(geneID)
            elif geneID not in SB1_sig and geneID not in SB3_sig:
                continue
                c = colors_dict["neither"]
                colors_count["neither"]+=1
            
            ax.scatter(x,y,color = c, s=ps, alpha = 0.75)

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
    main_sig = []
    markertype="o"
    # for cat in reversed(list(lists.keys())):
    if intersection_nums:
        for legend_label, count in colors_count.items():
            main_sig_ind = ax.scatter(1000,1000,color = colors_dict[legend_label], s=ps, alpha = 0.75, label = f"{legend_label} ({count})", marker=markertype)
            main_sig.append(main_sig_ind)
    else:
        for legend_label, count in colors_count.items():
            main_sig_ind = ax.scatter(1000,1000,color = colors_dict[legend_label], s=ps, alpha = 0.75, label = f"{legend_label}", marker=markertype)
            main_sig.append(main_sig_ind)

    ax.set_ylim([min_yline,max_yline])
    ax.set_xlim([min_xline,max_xline])

    main_sig_legend = plt.legend(handles = main_sig, fontsize = fs*0.75, title ="gene sig. in\nmain effect", title_fontsize = fs*0.7)#, loc='lower right')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True)) # force y axis as integers to make the y axis label visible and not outside of bounds

    ax.set_ylabel(f"small-Y SB - large-Y SB", fontsize = fs)
    ax.set_xlabel(f"logFC male line-bias", fontsize = fs)
    ax.tick_params(axis='x', labelsize=fs*0.9)
    ax.tick_params(axis='y', labelsize=fs*0.9)
    ax.set_title(LFC_title, fontsize = fs*0.95)

    plt.tight_layout()
    plt.savefig(LFC_filename, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {LFC_filename}")
    plt.clf()
    plt.cla()
    plt.close()

    return(shared_IDs)


if __name__ == "__main__":

    warnings.filterwarnings("ignore")
    username = "miltr339"
    table_paths,contrast_plot_titles = get_tables(username=username)
    out_path_figs = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures_python"
    
    ### genes that are line-biased in both males and females (should be excluded from male analysis in the sex-separated data since they can't be related to the Y-haplotype)
    excl_line_bias_lists = {
        "no_separation" : {
            "all_samples" : [],
        },
        "sex_separated" : {
            "day14" : ['gene-225709', 'gene-240623', 'gene-431701', 'gene-224956', 'gene-84970', 'gene-240871', 'gene-222486', 'gene-224697', 'gene-224743', 'gene-241001', 'gene-225236', 'gene-430032', 'gene-428738', 'gene-323803', 'gene-224682', 'gene-390616', 'gene-241126', 'gene-222600', 'gene-240910', 'gene-391222', 'gene-237881', 'gene-222519', 'gene-224227', 'gene-430314', 'gene-326873', 'gene-222365', 'gene-323148', 'gene-225738', 'gene-224896', 'gene-225173', 'gene-226245', 'gene-224357', 'gene-222383', 'gene-222555', 'gene-225030', 'gene-220028', 'gene-241055', 'gene-430080', 'gene-225720', 'gene-322927', 'gene-240929', 'gene-240833', 'gene-225140', 'gene-224968', 'gene-225635', 'gene-225107', 'gene-428765', 'gene-430068', 'gene-224593', 'gene-218529', 'gene-224782', 'gene-81427', 'gene-222332', 'gene-120952', 'gene-224860', 'gene-390956', 'gene-224201', 'gene-222350', 'gene-326849', 'gene-240691', 'gene-322912', 'gene-225325', 'gene-222344', 'gene-224875', 'gene-80359', 'gene-403809', 'gene-430263', 'gene-222531', 'gene-90157'],
            "day16" : ['gene-431701', 'gene-407280', 'gene-323803', 'gene-223491', 'gene-227370', 'gene-237881', 'gene-430314', 'gene-224896', 'gene-222555', 'gene-225030', 'gene-430080', 'gene-225140', 'gene-224593', 'gene-120952', 'gene-222531', 'gene-224201', 'gene-406796', 'gene-240691', 'gene-403809', 'gene-222746', 'gene-225709', 'gene-240623', 'gene-240871', 'gene-241001', 'gene-225236', 'gene-224079', 'gene-428738', 'gene-224682', 'gene-87700', 'gene-222519', 'gene-224227', 'gene-224890', 'gene-222365', 'gene-225738', 'gene-241262', 'gene-223773', 'gene-225720', 'gene-322927', 'gene-430068', 'gene-238407', 'gene-229506', 'gene-223318', 'gene-80359', 'gene-222486', 'gene-224697', 'gene-391222', 'gene-222159', 'gene-323148', 'gene-282853', 'gene-222383', 'gene-286545', 'gene-241055', 'g14784', 'gene-240929', 'gene-240833', 'gene-224968', 'gene-225107', 'gene-224860', 'gene-390956', 'gene-326849', 'gene-322912', 'gene-222344', 'gene-225325', 'gene-224956', 'gene-84970', 'gene-224743', 'gene-430032', 'gene-241126', 'gene-222600', 'gene-240910', 'gene-326873', 'gene-225173', 'gene-226245', 'gene-224357', 'gene-330102', 'gene-220028', 'gene-124877', 'gene-225635', 'gene-428765', 'gene-224782', 'gene-222332', 'gene-222350', 'gene-224875', 'gene-223419', 'gene-430263', 'gene-390616', 'gene-90157'],
            "day18" : ['gene-240623', 'gene-224697', 'gene-225236', 'gene-428738', 'gene-223491', 'gene-223758', 'gene-223773', 'gene-241055', 'gene-225720', 'gene-240833', 'gene-224593', 'gene-120952', 'gene-224860', 'gene-406796', 'gene-225325', 'gene-240691', 'gene-224875', 'gene-227308', 'gene-301479', 'gene-430263', 'gene-238849'],
            "SL1-SL3 days merged" : ['gene-5731', 'gene-391222', 'gene-124877', 'gene-370643', 'gene-222159', 'g5034', 'gene-328746', 'gene-215103', 'gene-224896', 'gene-80359', 'gene-225140', 'gene-241055', 'gene-222350', 'gene-224875', 'gene-88032', 'gene-327441', 'gene-241682', 'gene-222365', 'gene-370323', 'gene-210286', 'gene-214979', 'gene-221101', 'gene-225355', 'gene-241886', 'gene-84949', 'gene-237857', 'gene-428738', 'gene-224860', 'gene-238320', 'gene-240833', 'gene-84224', 'gene-424914', 'gene-224890', 'gene-224227', 'gene-370595', 'gene-326849', 'gene-90157', 'gene-240910', 'gene-269365', 'gene-38977', 'gene-430080', 'gene-120964', 'gene-243299', 'gene-282008', 'gene-131074', 'gene-120660', 'gene-84970', 'gene-241506', 'gene-428765', 'gene-224682', 'gene-225738', 'gene-223773', 'gene-302994', 'gene-240397', 'gene-225236', 'gene-87700', 'gene-224956', 'gene-218529', 'gene-224201', 'gene-89234', 'gene-223491', 'gene-118887', 'gene-224593', 'gene-323803', 'gene-81427', 'gene-231854', 'gene-227308', 'gene-430032', 'gene-222519', 'gene-395158', 'gene-231493', 'g14784', 'gene-225325', 'gene-55252', 'gene-238849', 'gene-278356', 'gene-406510', 'gene-407280', 'gene-215357', 'g422', 'gene-216986', 'gene-240056', 'gene-424896', 'gene-120832', 'gene-237378', 'gene-327074', 'gene-431030', 'gene-224782', 'gene-223419', 'gene-322912', 'gene-225030', 'gene-431701', 'gene-237881', 'gene-60429', 'gene-407355', 'gene-84577', 'gene-231604', 'gene-227164', 'gene-326873', 'gene-75424', 'gene-240691', 'gene-270971', 'gene-227370', 'gene-241841', 'gene-224968', 'gene-225107', 'gene-224079', 'gene-424863', 'gene-425532', 'gene-219157', 'gene-222531', 'gene-407253', 'gene-64665', 'gene-426056', 'gene-222344', 'gene-421265', 'gene-223758', 'gene-241268', 'gene-430263', 'gene-240623', 'gene-406796', 'gene-222600', 'gene-324176', 'gene-230270', 'gene-222486', 'gene-368309', 'gene-220028', 'gene-241328', 'gene-225635', 'gene-241238', 'gene-73288', 'gene-216914', 'gene-330102', 'gene-268146', 'gene-224357', 'gene-219019', 'gene-303243', 'gene-301479', 'gene-390616', 'gene-222555', 'gene-241262', 'gene-263597', 'gene-55240', 'gene-322927', 'gene-403809', 'gene-223318', 'gene-89057', 'gene-224697', 'gene-224118', 'gene-210277', 'gene-71572', 'gene-241856', 'gene-410330', 'gene-223213', 'gene-225921', 'gene-240929', 'gene-403902', 'gene-323991', 'gene-223599', 'gene-392186', 'gene-120952', 'gene-225709', 'gene-60610', 'gene-5642', 'gene-68612', 'gene-218086', 'gene-221980', 'gene-239553', 'gene-130722', 'gene-224743', 'gene-222332', 'gene-299793', 'gene-222383', 'gene-70017', 'gene-225173', 'gene-430068', 'gene-225720', 'gene-241126', 'gene-328764', 'gene-222746', 'gene-421566', 'gene-423321', 'gene-226245', 'gene-421549', 'gene-327787', 'gene-428089', 'gene-429522', 'gene-241001']
            },
        "line_separated" : {
            "SL1" : [],
            "SL3" : ['gene-122220', 'gene-312890', 'gene-47823', 'gene-92346', 'gene-206576', 'gene-74686', 'gene-2286', 'gene-334263', 'gene-100036', 'gene-163028', 'gene-226944', 'gene-189246', 'gene-384091', 'gene-202718', 'gene-238407', 'gene-414353', 'gene-205011', 'gene-75744', 'gene-343203', 'gene-121262', 'gene-60151', 'gene-232392', 'gene-182683', 'gene-206556', 'gene-166511', 'gene-73253', 'gene-48535', 'gene-388261', 'gene-228519', 'gene-253157', 'gene-211196', 'gene-227137', 'gene-317372', 'gene-132340', 'gene-60190', 'gene-306335', 'gene-21229', 'gene-234650', 'gene-39692', 'gene-288834', 'gene-62891', 'gene-410057', 'gene-253632', 'gene-198700', 'gene-9548', 'gene-233901', 'gene-69698', 'gene-378608', 'gene-218723', 'gene-285669', 'gene-277340', 'gene-336703', 'gene-30328', 'gene-40274', 'gene-97407', 'gene-39680', 'gene-39770', 'gene-244780', 'gene-328941', 'gene-130081', 'gene-410209', 'gene-130096', 'gene-120784', 'gene-120763', 'gene-143368', 'gene-87502', 'gene-166391', 'gene-99775', 'gene-304827', 'gene-350792', 'gene-231228', 'gene-153482', 'gene-279912', 'gene-377275']
            },
        "day_separated" : {
            "day14" : ['gene-390678', 'gene-326825', 'gene-428756', 'gene-392159', 'gene-237318', 'gene-80484', 'gene-222486', 'gene-326909', 'gene-391198', 'gene-220544', 'gene-225236', 'gene-323148', 'gene-90157', 'gene-403706', 'gene-428738', 'gene-328764', 'gene-241055', 'gene-390956', 'gene-241108', 'gene-403583', 'gene-323803', 'gene-80466', 'gene-89234', 'gene-237881', 'gene-224357', 'gene-226245', 'gene-395080', 'gene-220249', 'gene-430044', 'gene-222365', 'gene-240929', 'gene-392186', 'gene-430032', 'gene-224782', 'gene-239553', 'gene-240623', 'gene-224743', 'gene-403700', 'gene-117712', 'gene-224860', 'gene-407253', 'gene-225635', 'gene-403902', 'gene-225709', 'gene-225030', 'gene-406468', 'gene-395143', 'gene-322927', 'gene-222600', 'gene-240871', 'gene-84949', 'gene-81551', 'gene-224845', 'gene-431701', 'gene-240833', 'gene-222383', 'gene-225738', 'gene-224227', 'gene-243308', 'gene-224277', 'gene-241126', 'gene-240910', 'gene-240935', 'gene-392248', 'gene-283443', 'gene-403809', 'gene-222344', 'gene-218086', 'gene-390637', 'gene-392224', 'gene-225629', 'gene-224697', 'gene-224682', 'gene-406796', 'gene-214979', 'gene-222531', 'gene-231925', 'gene-430068', 'gene-225173', 'gene-80359', 'gene-241001', 'gene-81640', 'gene-81599', 'gene-225140', 'gene-391222', 'gene-225720', 'gene-219019', 'gene-222332', 'gene-430263', 'gene-390616', 'gene-120952', 'gene-222555', 'gene-224968', 'gene-224614', 'gene-224250', 'gene-428747', 'gene-224307', 'gene-222501', 'gene-81427', 'gene-222430', 'gene-221953', 'gene-81572', 'gene-225107', 'gene-430314', 'gene-327441', 'gene-222350', 'gene-240983', 'gene-240691', 'gene-84970', 'gene-322912', 'gene-224201', 'gene-224956', 'gene-428765', 'gene-224593', 'gene-392290', 'gene-326873', 'gene-224896', 'gene-326849', 'gene-403818', 'gene-326810', 'gene-403851', 'gene-88715', 'gene-234575', 'gene-220028', 'gene-260693', 'gene-225325', 'gene-224875', 'gene-222519', 'gene-430080', 'gene-403652'],
            "day16" : ['gene-390678', 'gene-431030', 'gene-231854', 'gene-428756', 'gene-282886', 'gene-282746', 'gene-222486', 'gene-222159', 'gene-220544', 'gene-391198', 'gene-326909', 'gene-225236', 'gene-223491', 'gene-323148', 'gene-90157', 'gene-400402', 'gene-428738', 'gene-224890', 'gene-241055', 'gene-390956', 'gene-282347', 'gene-241108', 'gene-238407', 'gene-323803', 'gene-89234', 'gene-282491', 'gene-237881', 'gene-224357', 'gene-226245', 'gene-282853', 'gene-395080', 'gene-220249', 'gene-380466', 'gene-282398', 'gene-399475', 'gene-430044', 'gene-286545', 'gene-222365', 'gene-240929', 'gene-392186', 'gene-430032', 'gene-224782', 'gene-428113', 'gene-241682', 'gene-240623', 'gene-224743', 'gene-222519', 'gene-224860', 'gene-242512', 'gene-225635', 'gene-225709', 'gene-225030', 'g14784', 'gene-223773', 'gene-395143', 'gene-322927', 'gene-222600', 'gene-227308', 'gene-240871', 'gene-84949', 'gene-224845', 'gene-282590', 'gene-431701', 'gene-240833', 'gene-222383', 'gene-225738', 'gene-224227', 'gene-282551', 'gene-399317', 'gene-224277', 'gene-241126', 'gene-240638', 'gene-240910', 'gene-225158', 'gene-403809', 'gene-222344', 'gene-390637', 'gene-225629', 'gene-399424', 'gene-224697', 'gene-400393', 'gene-224682', 'gene-243299', 'gene-406796', 'gene-241262', 'gene-282620', 'gene-222531', 'gene-430068', 'gene-225173', 'gene-400384', 'gene-80359', 'gene-227370', 'gene-241001', 'gene-81640', 'gene-225140', 'gene-391222', 'gene-428104', 'gene-225720', 'gene-222332', 'gene-430263', 'gene-282458', 'gene-390616', 'gene-120952', 'gene-399484', 'gene-222555', 'gene-224968', 'gene-282362', 'gene-428747', 'gene-224307', 'gene-222501', 'gene-81427', 'gene-239506', 'gene-222430', 'gene-221953', 'gene-282524', 'gene-81572', 'gene-225107', 'gene-399270', 'gene-430314', 'gene-222350', 'gene-240983', 'gene-240691', 'gene-84970', 'gene-238849', 'gene-322912', 'gene-224201', 'gene-224956', 'gene-428765', 'gene-224593', 'gene-224896', 'gene-240860', 'gene-326849', 'gene-282784', 'gene-403818', 'gene-282665', 'gene-220028', 'gene-260693', 'gene-225325', 'gene-224875', 'gene-392290', 'gene-430080', 'gene-282701'],
            "day18" : [],
            },
    }

            
    ############################################
    ######### MAKE ALL THE SMEAR PLOTS ######### # or upset plots 
    ############################################
    
    # only one of the below ones can be true at the same time! if both are false, smear/volcano plots are created by default
    ############
    make_upset = False # don't plot the smear/volcano plots but insetad make category-wise upset plots of DE genes
    ############
    make_list_outfiles = False # don't plot anything, instead make output files with lists of significant geneIDs for each contrast
    lists_outdir = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists"
    ############
    highlight_yTOR = False

    if False:
        
        if make_upset or make_list_outfiles:
            plot=False
            ## when making the sex-separated upset plot for the days, also include the no-separation day-sex interaction
            no_sep_day_by_sex = []
        else:
            plot=True

        for separation, seps_dict in table_paths.items():
            if separation != "day_separated":
            # if (separation == "sex_separated" or separation == "no_separation") == False:
                print(f"ignore {separation}")
                continue

            print(f"\n=========================== {separation} ===========================")
            
            # for upsetplot
            sep_DE_genes = {}
            # for intersetction geneIDs to use in GO enrichment
            gneIDs_all = []

            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")

                upset_all_cat=True
                # if "line_random" not in category:
                #     print(f"ignore {category}")
                #     upset_all_cat=False
                #     continue
    
                numbers = {}
                sig_DE_genes = {}
                
                for contrast, table_path in paths_dict.items():
                    
                    # if "(" not in contrast:
                    #     print(f"ignore {category}:{contrast}")
                    #     continue


                    # if "day" not in contrast:
                    #     print(f"ignore {contrast}")
                    #     continue

                    if "F" in contrast and "M" in contrast and "(" not in contrast:
                        # sex-biased contrast -> use |LFC| = 1 as min
                        # don't do it for the interaction
                        min_LFC = 1
                    else:
                        min_LFC = 0

                    if "line_random" in category and "(" in contrast:
                        min_LFC = 1

                    if contrast == "day16,day18" or contrast == "line3":
                        # this is not a contrasts but a test for the effect of coefficients, therefore the output table is structured a little differently and the plotting function needs to know
                        coefficients = True
                        print(f" ---<>--> {separation}:{category} --> coefficient(s): {contrast}")
                    else:
                        coefficients = False
                        print(f" ---<>--> {separation}:{category} --> contrast: {contrast}")

                    table_name = table_path.split("/")[-1].replace(".txt", "").replace("DE_genes_", "")
                    smear_name = f"{out_path_figs}/smear_{table_name}.png"

                    if "random" in category:
                        smear_title = f"{contrast_plot_titles[contrast]}"
                    elif "day" in category or "day" in contrast:
                        if coefficients==False:
                            smear_title = f"{category} {contrast_plot_titles[contrast]}"
                        elif coefficients == True:
                            smear_title = f"{contrast_plot_titles[contrast]}"
                    else:
                        smear_title = contrast_plot_titles[contrast]
                        smear_title = smear_title.replace(" days merged", "")

                    excl_line_bias = excl_line_bias_lists[separation]
                    excl_list = []
                    excl_list_name = smear_title.replace(" ","")
                    if excl_list_name in excl_line_bias:
                        # only do the exclusion list when it's actually in a relevant contrast
                        excl_list = excl_line_bias[excl_list_name]
                        print(f"\texcluding {len(excl_list)} genes from list '{excl_list_name}'")
                    elif "day" in separation and category in excl_line_bias:
                        excl_list = excl_line_bias[category]
                        print(f"\texcluding {len(excl_list)} genes from list '{excl_list_name}'")
        
                    # smear plot
                    smear_lists = plot_smear(table_path=table_path, contrast=contrast, smear_plot_name=smear_name, min_LFC=min_LFC, title = smear_title, excl_genes_list=excl_list, x_axis="logcpm", plot=plot, coefficients=coefficients)
                    # volcano plot
                    smear_lists = plot_smear(table_path=table_path, contrast=contrast, smear_plot_name=smear_name, min_LFC=min_LFC, title = smear_title, excl_genes_list=excl_list, x_axis="fdr_p", plot=plot, coefficients=coefficients, highlight_yTOR=highlight_yTOR)
                    
                    if False:
                        Downlist = smear_lists["Downregulated"]
                        Uplist = smear_lists["Upregulated"]
                        print(f"\t * Downregulated : {Downlist}\n\t * Upregulated : {Uplist}")
                    elif coefficients==False:
                        DE_list = smear_lists["Downregulated"]+smear_lists["Upregulated"]
                        nonDE_list = smear_lists["no difference"]
                        all_geneIDs = gneIDs_all + smear_lists["Downregulated"]+smear_lists["Upregulated"]+smear_lists["no difference"]

                        if make_list_outfiles:
                            sig_list_outfile = table_path.split("/")[-1].replace("DE_genes_", "sig_DE_list_")
                            sig_list_outfile = f"{lists_outdir}/{sig_list_outfile}"
                            with open(sig_list_outfile, "w") as sig_list_file:
                                DE_list_outfile = [f"{geneID},1" for geneID in DE_list]
                                DE_string = "\n".join(DE_list_outfile)
                                sig_list_file.write(f"geneID,sig_DE\n{DE_string}\n") 
                                nonDE_list_outfile = [f"{geneID},0" for geneID in nonDE_list]
                                nonDE_string = "\n".join(nonDE_list_outfile)
                                sig_list_file.write(f"{nonDE_string}\n") # needs the newline character so that R can read the list right

                            print(f" * list of sig DE genes written to: {sig_list_outfile}")

                        smear_nums = {gene_set : len(gene_list) for gene_set,gene_list in smear_lists.items()}
                        numbers[contrast] = smear_nums

                        if category =="line_random":
                            # for full dataset line random effect DO include interactions
                            sig_DE_genes[contrast] = DE_list
                            sep_DE_genes[f"{category}: {contrast}"] = DE_list
                        elif "(" not in contrast:
                            # if not full dataset line random effect DON'T include interactions in upset plot
                            sig_DE_genes[contrast] = DE_list
                            sep_DE_genes[f"{category}: {contrast}"] = DE_list

                        if "(" in contrast and category =="line_random":
                            ## save the day-sex interaction for the full dataset separately to use in an upsetplot for the sex separated samples
                            no_sep_day_by_sex = DE_list
                            print(f"!!!!!!!!!! interaction day-by-sex full data: {len(no_sep_day_by_sex)} genes")
                    
                    elif coefficients:
                        for coeff, coeff_list in smear_lists.items():
                            numbers[f"coef:{coeff}"] = {gene_set : len(gene_list) for gene_set,gene_list in coeff_list.items()}

                if plot:
                    for contrast, number in numbers.items():
                        print(f"\tNUMBERS SUMMARY:")
                        print(f"\t{contrast} : {number}")
                if make_upset:
                    ## upsetplot within each separation
                    min_overlap = 30
                    plot_title = f"{separation}: all categories\nsig. DE genes overlap\n(min. overlap size: {min_overlap})"
                    if separation == "no_separation":
                        # fix the contrast names when the interaction is included in the no_separation data
                        sig_DE_genes = { contrast_plot_titles[key] : val for key,val in sig_DE_genes.items()}
                        plot_title = f"min. overlap size: {min_overlap}"
                    plot_sig_gene_overlap(sig_DE_genes, plot_filename= f"{out_path_figs}/upsetplot_{separation}_{category}.png", plot_title = plot_title, min_overlap = min_overlap)

            if make_upset and upset_all_cat:
                ## upsetplot for each contrast between categories
                min_overlap = 15
                plot_title = f"{separation}: all categories\nsig. DE genes overlap (min. overlap size: {min_overlap})"
                plot_sig_gene_overlap(sep_DE_genes, plot_filename= f"{out_path_figs}/upsetplot_{separation}_all_categories.png", plot_title = plot_title, min_overlap = min_overlap)

                if separation == "day_separated":
                    # plot only male line bias in the day separated data
                    sep_DE_genes_line = {sep.replace(": M_1 - M_3", ": males") : l_ for sep,l_ in sep_DE_genes.items() if "M_1 - M_3" in sep}
                    min_overlap_ = 0
                    plot_title_ = f"male line-biased genes\n(small-large)"
                    upset_data_line = plot_sig_gene_overlap(sep_DE_genes_line, plot_filename= f"{out_path_figs}/upsetplot_{separation}_male_line_bias.png", plot_title = plot_title_, min_overlap = min_overlap_)
                    # plot only sex bias for each line in the day separated data
                    sep_DE_genes_1 = {sep.replace(": F_1 - M_1", ": small-Y") : l_ for sep,l_ in sep_DE_genes.items() if "F_1 - M_1" in sep}
                    sep_DE_genes_3 = {sep.replace(": F_3 - M_3", ": large-Y") : l_ for sep,l_ in sep_DE_genes.items() if "F_3 - M_3" in sep}
                    sep_DE_genes__ = sep_DE_genes_1|sep_DE_genes_3
                    sep_DE_genes_sex = dict(reversed(sep_DE_genes__.items()))
                    min_overlap_ = 15
                    plot_title_ = f"Sex biased genes\n(min. overlap size: {min_overlap})"
                    upset_data_sex = plot_sig_gene_overlap(sep_DE_genes_sex, plot_filename= f"{out_path_figs}/upsetplot_{separation}_all_lines_sex_bias.png", plot_title = plot_title_, min_overlap = min_overlap_)

                if separation == "sex_separated":
                    # plot only male line bias in the day separated data
                    sep_DE_genes_ = {sep : l_ for sep,l_ in sep_DE_genes.items() if "day" in sep and "females" not in sep}
                    sep_DE_genes_["full dataset: day-by-sex"] = no_sep_day_by_sex
                    min_overlap_ = 20
                    plot_title_ = f"min. overlap size: {min_overlap_}"
                    upset_data = plot_sig_gene_overlap(sep_DE_genes_, plot_filename= f"{out_path_figs}/upsetplot_{separation}_day_bias.png", plot_title = plot_title_, min_overlap = min_overlap_)

                #######################################
                #### make list of all the significantly line-biased genes from Fig 1 to do the GO enrichment                
                #######################################
                if True and separation == "day_separated":

                    print(f"\n\n\n\n ---<>--> male line bias upset data")

                    ######### male samples line bias
                    if True:
                        # filter to only include genes that are sig. in at least two days
                        mask = (
                            upset_data_line.index.get_level_values('day14: males').astype(int) +
                            upset_data_line.index.get_level_values('day16: males').astype(int) +
                            upset_data_line.index.get_level_values('day18: males').astype(int)
                        ) >= 2
                        filt = upset_data_line[mask]
                        print(f"filtered for all geneIDs that are sig in at least two:")
                        print(f"{len(filt)}")

                        filt.to_csv(f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/day_separated_male_line_bias_overlap_sigIDs.txt", sep="\t")

                        # filter to only include genes that are sig. only day 14 and 16
                        filt = upset_data_line[upset_data_line.index.get_level_values('day14: males') & upset_data_line.index.get_level_values('day16: males')]
                        filt = filt["id"].tolist()
                        nonsig_geneIDs = [id for id in list(set(all_geneIDs)) if id not in filt]
                        nonsig = ",0\n".join(nonsig_geneIDs)+",0\n"

                        filt = ",1\n".join(filt)+",1\n"
                        outfile_path = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/sig_DE_list_day_sep_M1-M3_day14and16.txt"
                        with open(outfile_path, "w") as outfile:
                            outfile.write("geneID,sig_DE\n")
                            outfile.write(filt)
                            outfile.write(nonsig)

                    ######### all lines sex bias
                    if True:
                        # filter to only include genes that are sig. in at least two days
                        mask_all = (
                            upset_data_sex.index.get_level_values('day14: small-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day16: small-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day18: small-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day14: large-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day16: large-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day18: large-Y').astype(int)
                        ) ==6
                        filt_all = upset_data_sex[mask_all]
                        print(f"filtered for all geneIDs that are sig in all:")
                        print(f"{len(filt_all)}")

                        # filt_all.to_csv(f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/day_separated_all_lines_sex_bias_overlap_all_days_sigIDs.txt", sep="\t")
                        filt_all = filt_all["id"].tolist()
                        nonsig_geneIDs = [id for id in list(set(all_geneIDs)) if id not in filt_all]
                        nonsig = ",0\n".join(nonsig_geneIDs)+",0\n"

                        filt_all = ",1\n".join(filt_all)+",1\n"
                        outfile_path = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/day_separated_all_lines_sex_bias_overlap_all_days_sigIDs.txt"
                        with open(outfile_path, "w") as outfile:
                            outfile.write("geneID,sig_DE\n")
                            outfile.write(filt_all)
                            outfile.write(nonsig)


                        # filter to only include genes that are sig. only day 16 and 18 in both lines
                        mask_late = (
                            upset_data_sex.index.get_level_values('day14: small-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day16: small-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day18: small-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day14: large-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day16: large-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day18: large-Y').astype(int)
                        ) ==4
                        filt_late = upset_data_sex[mask_late]
                        filt_late = filt_late[filt_late.index.get_level_values('day16: small-Y') & filt_late.index.get_level_values('day18: small-Y') & filt_late.index.get_level_values('day16: large-Y') & filt_late.index.get_level_values('day18: large-Y')]
                        print(f"filtered for all geneIDs that are sig day 16 and 18 both lines:")
                        print(f"{len(filt_late)}")
                        filt_late = filt_late["id"].tolist()
                        nonsig_geneIDs = [id for id in list(set(all_geneIDs)) if id not in filt_late]
                        nonsig = ",0\n".join(nonsig_geneIDs)+",0\n"

                        filt_late = ",1\n".join(filt_late)+",1\n"
                        outfile_path = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/day_separated_all_lines_sex_bias_overlap_day1618_sigIDs.txt"
                        with open(outfile_path, "w") as outfile:
                            outfile.write("geneID,sig_DE\n")
                            outfile.write(filt_late)
                            outfile.write(nonsig)

                    ######### sex bias in all exc. day14:smallY
                    if True:    
                        # filter only genes that are sig. not in day14 small-Y, just everything else
                        mask_early = (
                            upset_data_sex.index.get_level_values('day14: small-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day16: small-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day18: small-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day14: large-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day16: large-Y').astype(int) +
                            upset_data_sex.index.get_level_values('day18: large-Y').astype(int)
                        ) ==5

                        filt_early = upset_data_sex[mask_early]
                        filt_early = filt_early[~filt_early.index.get_level_values('day14: small-Y')]
                        print(f"filtered for all geneIDs that are sig. in everything except day14:smallY")
                        print(f"{len(filt_early)}")
                        filt_early = filt_early["id"].tolist()
                        nonsig_geneIDs = [id for id in list(set(all_geneIDs)) if id not in filt_early]
                        nonsig = ",0\n".join(nonsig_geneIDs)+",0\n"

                        filt_early = ",1\n".join(filt_early)+",1\n"
                        outfile_path = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/day_separated_all_lines_sex_bias_overlap_all_excl_day14smallY_sigIDs.txt"
                        with open(outfile_path, "w") as outfile:
                            outfile.write("geneID,sig_DE\n")
                            outfile.write(filt_early)
                            outfile.write(nonsig)


    ############################################
    ######## MAKE ALL THE VENN DIAGRAMS ########
    ############################################
    
    ## standard sets matching the tabs in the html
    if False:
        venn_sets = {
            "no_separation" : { "all_samples" : {"" : []}},
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
            if separation != "day_separated":
                continue

            print(f"\n=========================== {separation} ===========================")

            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")

                for venn_cat, venn_contrasts_list in venn_sets[separation][category].items():
                    if venn_contrasts_list == []:
                        continue
                    print(f"{venn_cat} : {venn_contrasts_list}")
                    venn_paths_dict = {contrast_plot_titles[contrast] : paths_dict[contrast] for contrast in venn_contrasts_list}
                    venn_filename_ = venn_cat.replace(" ", "_")
                    venn_filename = f"{out_path_figs}/Venn_{category}_{venn_filename_}.png"
                    venn_title = f"sig. DE genes overlap ({category})\n{venn_cat}"
                    plot_venn_DE_genes(venn_paths_dict, venn_filename=venn_filename, venn_title=venn_title)

    ###########################################
    ##### make all the lists of genes to exclude because they are line-biased in females
    ## THESE ARE HARDCODED ABOVE! COPY-PASTE FROM THE TERMINAL, NO ACTUAL VARIABLES ARE SAVED, EVERYTHING IS IN 'excl_line_bias_lists'
    ###########################################

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
                },
                "days_merged" : {
                    "females" : ["SL1 - SL3"],
                    "males" : ["SL1 - SL3"]
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
                if day == "days_merged":
                    day_ = day.replace("_", " ")
                    venn_filename = venn_filename.replace("f_vs_m", "SL1-SL3")
                    venn_title = f"sig. line biased genes overlap\nfemales and males"
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
                    "interaction" : ["(F_1 - M_1) - (F_3 - M_3)"],
                },
                "day16" : {
                    "males" : ["M_1 - M_3"],
                    "females" : ["F_1 - F_3"],
                    "interaction" : ["(F_1 - M_1) - (F_3 - M_3)"],
                    },
                "day18" : {
                    "males" : ["M_1 - M_3"],
                    "females" : ["F_1 - F_3"],
                    "interaction" : ["(F_1 - M_1) - (F_3 - M_3)"],
                    },
            }
        }
        for separation, days_dict in venn_sets_day.items():
            print(f"\n=========================== {separation} ===========================")
            interaction_paths = {}
            for day, sexes_contrasts_dict in days_dict.items():
                print(f"\n ------------------- {day} -------------------")

                venn_paths_dict = {}
                for sex, venn_contrasts_list in sexes_contrasts_dict.items():

                    if sex=="interaction":
                        interaction_paths[day] = table_paths[separation][day][venn_contrasts_list[0]]
                    else:
                        print(f"{sex} : {venn_contrasts_list}")
                        venn_paths_dict[sex] = table_paths[separation][day][venn_contrasts_list[0]]

                venn_filename = f"{out_path_figs}/Venn_{day}_f_vs_m.png"
                venn_title = f"sig. DE genes overlap ({day})\nfemales and males"
                shared_list = plot_venn_DE_genes(venn_paths_dict, venn_filename=venn_filename, venn_title=venn_title, get_shared_list=True, plot=False)
                print(f"{len(shared_list)} genes : \n{shared_list}")

            print(f"\n ------------------- all day interactions -------------------")
            for day,path in interaction_paths.items():
                print(f"{day} : {path}")
            venn_filename = f"{out_path_figs}/Venn_day_sep_all_interactions.png"
            venn_title = f"sig. DE genes overlap of all days\nline-by-sex interaction"
            shared_list = plot_venn_DE_genes(interaction_paths, venn_filename=venn_filename, venn_title=venn_title, get_shared_list=True, plot=True)
            # print(f"{len(shared_list)} genes : \n{shared_list}")

    ## compare sex bias and interaction on full dataset with days merged/ignored
    if False:
        venn_sets_all = {
            "no_separation" : {
                "all_samples:sex_bias" : {
                    "SL1" : ["SL1_F - SL1_M"],
                    "SL3" : ["SL3_F - SL3_M"],
                    "line by sex" : ["(SL1_F - SL3_F) - (SL1_M - SL3_M)"]
                },
                "all_samples:line_bias" : {
                    "females" : ["SL1_F - SL3_F"],
                    "males" : ["SL1_M - SL3_M"],
                    "line by sex" : ["(SL1_F - SL3_F) - (SL1_M - SL3_M)"]
                }
            }
        }
        
        for separation, samples_dict in venn_sets_all.items():
            print(f"\n=========================== {separation} ===========================")
            for sample, sexes_contrasts_dict in samples_dict.items():
                print(f"\n ------------------- {sample} -------------------")

                venn_paths_dict = {}
                sample_ds = sample.split(":")[0]
                for line, venn_contrasts_list in sexes_contrasts_dict.items():
                    print(f"{line} : {venn_contrasts_list}")
                    venn_paths_dict[line] = table_paths[separation][sample_ds][venn_contrasts_list[0]]
                
                sample_=sample.replace(":", "_")
                venn_filename = f"{out_path_figs}/Venn_{sample_}_no_age.png"
                day_ = sample.replace("day", "day ")
                if sample == "all_samples":
                    day_ = sample.replace("_", " ")
                venn_title = f"sig. DE genes overlap ({day_})"
                shared_list = plot_venn_DE_genes(venn_paths_dict, venn_filename=venn_filename, venn_title=venn_title, get_shared_list=True)
                print(shared_list)

    ############################################
    ######## MAKE ALL THE SB COMPARISONS #######
    ############################################
    ##
    ## plot the scatterplots to compare LogFC of two different contrasts
    if False:
        LFC_comp_sets = {
            "no_separation" : {
                "all_samples" : {
                    "sex bias" : ["SL1_F - SL1_M", "SL3_F - SL3_M"],
                    "line bias" : ["SL1_F - SL3_F", "SL1_M - SL3_M"],
                },
            },
        
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

            if separation != "day_separated":
                print(f"ignore {separation}")
                continue
            
            print(f"\n=========================== {separation} ===========================")

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
    if False:
        sig_IDs_list = {
            ## geneIDs that are significant in the day separated line-by-sex interaction
            ## from PhD_chapter4/data/sig_DE_genes_lists/sig_DE_list_day14_F-M_by_1-3.txt and other days
            "day_separated" : {
                "day14" : ["gene-237342","gene-181562","gene-58400"],
                "day16" : [None], # no sig. interaction
                "day18" : ["gene-426041","gene-48523","gene-73742","gene-2355","gene-237494","gene-279406","gene-237342","gene-276797","gene-367071"],
            }
        }
        contrasts_interaction_list = {
            "day_separated" : {
                "sex bias" : ["F_1 - M_1","F_3 - M_3"],
                "line bias" : ["F_1 - F_3","M_1 - M_3"],
                "sb_SL1 lb_F" : ["F_1 - M_1","F_1 - F_3"], # sex bias SL1, line bias Females
                "sb_SL3 lb_F" : ["F_3 - M_3","F_1 - F_3"],
                "sb_SL1 lb_M" : ["F_1 - M_1","M_1 - M_3"],
                "sb_SL3 lb_M" : ["F_3 - M_3","M_1 - M_3"],
                "sb_diff lb_M" : {"sb_diff" : ["F_1 - M_1","F_3 - M_3"], "lb_M" : ["M_1 - M_3"]} # do line biased genes have a stronger or weaker sex bias
            }
        }

        for separation, seps_dict in table_paths.items():
            print(f"\n=========================== {separation} ===========================")

            # only relevant interactions are in the day-separated data
            if "day" not in separation:
                continue
            
            shared_SB_LB_for_time_series = {"small-Y SB" : [] , "large-Y SB" : [], "both SB" : []}

            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")

                for bias_cat, contrasts_list in contrasts_interaction_list[separation].items():
                    
                    if "bias" in bias_cat or "lb_F" in bias_cat:
                        print(f"skip {bias_cat}:({contrasts_list})")
                        continue

                    LFC_filename_ = bias_cat.replace(" ", "_")
                    LFC_filename = f"{out_path_figs}/LFC_scatter_{category}_{LFC_filename_}.png"
                    category_ = category.replace("day", "day ")
                    bias_cat_ = bias_cat.replace("lb_", ", line-bias ")

                    if "diff" in bias_cat:
                        bias_cat_ = bias_cat_.replace("sb_diff ", "sex-bias")
                        plot_title = f"{category_}: {bias_cat_}"
                        
                        tables_diff = [paths_dict[contrast] for contrast in contrasts_list["sb_diff"]]
                        table_SB = paths_dict[contrasts_list["lb_M"][0]]
                        excl_geneIDs = excl_line_bias_lists[separation][category]

                        shared_IDs = plot_sig_LFC_diff(tables_diff=tables_diff, table_LB=table_SB, LFC_filename = LFC_filename, excl_geneIDs=excl_geneIDs, LFC_title = plot_title, intersection_nums = True )

                        if False:
                            shared_SB_LB_for_time_series["small-Y SB"].extend(shared_IDs["small-Y SB"])
                            shared_SB_LB_for_time_series["large-Y SB"].extend(shared_IDs["large-Y SB"])
                            shared_SB_LB_for_time_series["both SB"].extend(shared_IDs["both SB"])
                            ## plot time series of expression of shared line- and sex-biased genes
                            count_files = time_series_plots.get_counts_paths(username=username)
                            samples_group_dict = time_series_plots.samples_group()
                            
                            plot_dir = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures_python/counts_time_series"
                            plot_file = f"{plot_dir}/time_series_day_and_sex_bias_{category}.png"

                            plot_title = f"{category_}: genes that are line-biased in males and sex-biased"
                            
                            time_series_plots.plot_counts_sum_sets(counts_table=count_files["no_log"], geneIDs_lists_dict = shared_IDs, 
                                                outfile_name = plot_file, y_label= "normalized counts", errorbars=True, samples_group_dict = samples_group_dict, plot_title=plot_title)

                    else:
                        continue

                        try:
                            incl_geneIDs = sig_IDs_list[separation][category]
                        except:
                            print(f"no interaction genes included for category '{category}'") 
                            incl_geneIDs = []

                        if bias_cat != "line bias":
                            incl_geneIDs = []
                            print(f"{bias_cat} : no interaction included")
                        else:
                            print(f"{bias_cat} : interaction has {len(incl_geneIDs)} genes")

                        LFC_paths_dict = {contrast_plot_titles[contrast] : paths_dict[contrast] for contrast in contrasts_list}
                        # plot_title = f"{bias_cat} in males and females"
                        bias_cat_ = bias_cat.replace("sb_", "sex-bias ").replace("lb_", ", line-bias ")
                        bias_cat_ = bias_cat_.replace("SL1", "small-Y").replace("SL3", "large-Y")

                        if len(incl_geneIDs)>0:
                            plot_title = f"{category_}: {bias_cat_} and interaction"
                        else:
                            plot_title = f"{category_}: {bias_cat_}"
                        if len(plot_title)>33:
                            plot_title = plot_title.replace(" and", "\nand")
                        
                        numbers = plot_sig_LFC_overlap(LFC_paths_dict, LFC_filename = LFC_filename, LFC_title = plot_title , incl_geneIDs=incl_geneIDs, intersection_nums = True )
                        print(f"\t{numbers}")

            if False:
                ## plot time series of expression of shared line- and sex-biased genes
                count_files = time_series_plots.get_counts_paths(username=username)
                samples_group_dict = time_series_plots.samples_group()
                sex_chromosomes_superscaffolded,y_contigs,x_contigs,tor_related = time_series_plots.get_y_information()

                plot_dir = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures_python/counts_time_series"
                plot_file = f"{plot_dir}/time_series_day_and_sex_bias_all_days.png"

                plot_title = f"genes that are line-biased in males and sex-biased on at least one day"
                unique_shared_SB_LB_for_time_series = {cat : list(set(geneids)) for cat,geneids in shared_SB_LB_for_time_series.items()}

                print(f"\n\tY-PROPORTION:")
                ## check how many of the sex- and line-biased genes are on the Y
                for cat,geneIDs_cat in unique_shared_SB_LB_for_time_series.items():
                    num_genes_cat = len(geneIDs_cat)
                    intersection = list(set(y_contigs["expressed"]) & set(geneIDs_cat))
                    Y_count = len(intersection)
                    frac = 100*Y_count/num_genes_cat
                    print(f"\t - {cat}: {frac:.2f}% ({Y_count}/{num_genes_cat}) : {intersection}")
                
                print(f"\n")
                time_series_plots.plot_counts_sum_sets(counts_table=count_files["no_log"], geneIDs_lists_dict = unique_shared_SB_LB_for_time_series, 
                                    outfile_name = plot_file, y_label= "normalized counts", errorbars=True, samples_group_dict = samples_group_dict, plot_title=plot_title)
                print(unique_shared_SB_LB_for_time_series)

    #############################################
    ####### PLOT LOGFC MAGNITUDE BOXPLOTS #######
    #############################################

    ## plot LogFC boxplots of male-biased genes of several contrasts within each separation
    if True:
        LFC_comp_sets = {
            "day_separated" : {
                "day14" : {
                    "SL1" : "F_1 - M_1",
                    "SL3" : "F_3 - M_3",
                },
                "day16" : {
                    "SL1" : "F_1 - M_1",
                    "SL3" : "F_3 - M_3",
                    },
                "day18" : {
                    "SL1" : "F_1 - M_1",
                    "SL3" : "F_3 - M_3",
                    },
            }
        }
        for separation, seps_dict in table_paths.items():
            if separation not in LFC_comp_sets.keys():
                continue
            print(f"\n=========================== {separation} ===========================")

            boxplot_contrasts = {}
            for category, paths_dict in seps_dict.items():
                # print(f"\n ------------------- {category} -------------------")
                for name, contrast in LFC_comp_sets[separation][category].items():
                    boxplot_contrasts[f"{category}:{name}"] = paths_dict[LFC_comp_sets[separation][category][name]] 

            # set minLFC to only include male-biased genes
            # plot_logFC_boxplots(infiles_dict= boxplot_contrasts, min_LFC=-1, plot_filename=f"{out_path_figs}/{separation}_sex_bias_LFC_boxplot.png")
            plot_logFC_boxplots(infiles_dict= boxplot_contrasts, min_LFC=-1, only_all_intersection=True, plot_filename=f"{out_path_figs}/{separation}_sex_bias_LFC_boxplot_shared_sig_only.png")
