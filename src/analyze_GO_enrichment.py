"""
Analyze the significantly enriched GO terms in all the contrasts
"""

import pandas as pd
import upsetplot
import matplotlib.pyplot as plt


def get_tables(username="miltr339"):
    """
    Tables from GO enrichment analysis in PhD_chapter4/src/GO_enrichment.Rmd
    """
    tables_dir = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/"

    out_dict = {
        
        "sex_separated" : {
            "females" :{
                "SL1_14 - SL3_14" : f"{tables_dir}GO_enrichment_F_1-3_day14.txt",
                "SL1_16 - SL3_16" : f"{tables_dir}GO_enrichment_F_1-3_day16.txt",
                "SL1_18 - SL3_18" : f"{tables_dir}GO_enrichment_F_1-3_day18.txt",
                "SL3_18 - SL3_14" : f"{tables_dir}GO_enrichment_F_SL3_18_14.txt",
            },
            "males" : {
                "SL1_14 - SL3_14" : f"{tables_dir}GO_enrichment_M_1-3_day14.txt",
                "SL1_16 - SL3_16" : f"{tables_dir}GO_enrichment_M_1-3_day16.txt",
                "SL1_18 - SL3_18" : f"{tables_dir}GO_enrichment_M_1-3_day18.txt",
                "SL1_18 - SL1_14" : f"{tables_dir}GO_enrichment_M_SL1_18_14.txt",
                "SL1_18 - SL1_16" : f"{tables_dir}GO_enrichment_M_SL1_18_16.txt",
                "SL1_14 - SL1_16" : f"{tables_dir}GO_enrichment_M_SL1_14_16.txt",
                "SL3_18 - SL3_14" : f"{tables_dir}GO_enrichment_M_SL3_18_14.txt",
                "SL3_18 - SL3_16" : f"{tables_dir}GO_enrichment_M_SL3_18_16.txt",
                "SL3_14 - SL3_16" : f"{tables_dir}GO_enrichment_M_SL3_14_16.txt",
            }
        },
        "line_separated" : {
            "SL1" : {
                "F_14 - M_14" : f"{tables_dir}GO_enrichment_SL1_day14_F-M.txt",
                "F_16 - M_16" : f"{tables_dir}GO_enrichment_SL1_day16_F-M.txt",
                "F_18 - M_18" : f"{tables_dir}GO_enrichment_SL1_day18_F-M.txt",
                "(F_18 - M_18) - (F_14 - M_14)" : f"{tables_dir}GO_enrichment_SL1_day18_14_F-M.txt",
            },
            "SL3" : {
                "F_14 - M_14" : f"{tables_dir}GO_enrichment_SL3_day14_F-M.txt",
                "F_16 - M_16" : f"{tables_dir}GO_enrichment_SL3_day16_F-M.txt",
                "F_18 - M_18" : f"{tables_dir}GO_enrichment_SL3_day18_F-M.txt",
                "(F_18 - M_18) - (F_14 - M_14)" : f"{tables_dir}GO_enrichment_SL3_day18_14_F-M.txt",
            }
        },
        "day_separated" : {
            "day14" : {
                "F_1 - M_1" : f"{tables_dir}GO_enrichment_day14_SL1_F-M.txt",
                "F_3 - M_3" : f"{tables_dir}GO_enrichment_day14_SL3_F-M.txt",
                "F_1 - F_3" : f"{tables_dir}GO_enrichment_day14_F_1-3.txt",
                "M_1 - M_3" : f"{tables_dir}GO_enrichment_day14_M_1-3.txt",
                "(F_1 - M_1) - (F_3 - M_3)" : f"{tables_dir}GO_enrichment_day14_F-M_by_1-3.txt",
            },
            "day16" : {
                "F_1 - M_1" : f"{tables_dir}GO_enrichment_day16_SL1_F-M.txt",
                "F_3 - M_3" : f"{tables_dir}GO_enrichment_day16_SL3_F-M.txt",
                "F_1 - F_3" : f"{tables_dir}GO_enrichment_day16_F_1-3.txt",
                "M_1 - M_3" : f"{tables_dir}GO_enrichment_day16_M_1-3.txt",
                "(F_1 - M_1) - (F_3 - M_3)" : f"{tables_dir}GO_enrichment_day16_F-M_by_1-3.txt",
            },
            "day18" : {
                "F_1 - M_1" : f"{tables_dir}GO_enrichment_day18_SL1_F-M.txt",
                "F_3 - M_3" : f"{tables_dir}GO_enrichment_day18_SL3_F-M.txt",
                "F_1 - F_3" : f"{tables_dir}GO_enrichment_day18_F_1-3.txt",
                "M_1 - M_3" : f"{tables_dir}GO_enrichment_day18_M_1-3.txt",
                "(F_1 - M_1) - (F_3 - M_3)" : f"{tables_dir}GO_enrichment_day18_F-M_by_1-3.txt",
            },
        }
    }

    return out_dict


def read_R_GO_table(table_filepath:str, p_val = 0.05):
    """
    reads the output from the enrichment analysis with TopGO and returns a list of GO terms that are enriched with a p-value below p_val
    """
    df = pd.read_csv(table_filepath)
    df = df.loc[df['classicFisher'] <= p_val]
    sig_terms = df.shape[0]
    GO_list = df["GO.ID"].tolist()
    print(f" * {sig_terms} sig. enriched GO terms")
    return GO_list


def plot_enriched_GO_overlap(GO_Terms_dict:dict, plot_filename:str, plot_title  = "", min_overlap = 0):
    """
    make an upseetplot of sig DE genes within each data separation
    """
    data = upsetplot.from_contents(GO_Terms_dict)
    if min_overlap>0:
        upsetplot.UpSet(data, subset_size="count", show_counts=True, min_subset_size=min_overlap, sort_by="cardinality", sort_categories_by="input").plot()
    else:
        upsetplot.UpSet(data, subset_size="count", sort_by="cardinality", show_counts=True, min_subset_size=min_overlap).plot()

    if min_overlap >0:
        plt.subplots_adjust(top=0.85)    
    plt.title(plot_title)
    plt.tight_layout()
    plt.savefig(plot_filename, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {plot_filename}")
    plt.clf()
    plt.cla()
    plt.close()


if __name__ == "__main__":
    username = "miltr339"
    go_tables_paths = get_tables(username=username)
    out_path_figs = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures_python/GO_enrichment"

    if True:
        for separation, seps_dict in go_tables_paths.items():
            print(f"\n=========================== {separation} ===========================")

            sep_GO_terms = {}
            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")
                
                sig_GO_terms = {}
                for contrast, table_path in paths_dict.items():
                    
                    print(f"{table_path}")
                    if "(" not in contrast:
                        GO_list = read_R_GO_table(table_filepath=table_path)
                        sig_GO_terms[contrast] = GO_list
                        sep_GO_terms[f"{category}: {contrast}"] = GO_list

                min_overlap = 5
                separation_  = separation.replace("_", "-")
                if min_overlap>0:
                    plot_title = f"{separation_}:{category}\nsig. GO:terms overlap\n(min. overlap size: {min_overlap})"
                else:
                    plot_title = f"{separation_}:{category}\nsig. GO:terms overlap"
                plot_enriched_GO_overlap(sig_GO_terms, plot_filename= f"{out_path_figs}/upsetplot_GO_terms_{separation}_{category}.png", plot_title = plot_title, min_overlap = min_overlap)

            min_overlap = 5
            if min_overlap>0:
                plot_title = f"{separation_}:all categories\nsig. GO:terms overlap\n(min. overlap size: {min_overlap})"
            else:
                plot_title = f"{separation_}:all categories\nsig. GO:terms overlap"
            plot_enriched_GO_overlap(sep_GO_terms, plot_filename= f"{out_path_figs}/upsetplot_GO_terms_{separation}_all.png", plot_title = plot_title, min_overlap = min_overlap)
            