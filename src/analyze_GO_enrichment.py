"""
Analyze the significantly enriched GO terms in all the contrasts
"""

import pandas as pd
import upsetplot
import matplotlib.pyplot as plt
import go3


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

    


def get_interesting_GO_overlap_lists():
    """
    list of contrasts whose intersection of enriched GO terms is interesting
    """    
    out_dict = {
        "sex_separated" : {
            "females" :{
                "day14" : ["SL1_14 - SL3_14"],
                "day16" : ["SL1_16 - SL3_16"],
                "day18" : ["SL1_18 - SL3_18"],
                "SL3_day18-14" : ["SL3_18 - SL3_14"],
            },
            "males" : {
                "day14" : ["SL1_14 - SL3_14"],
                "day16" : ["SL1_16 - SL3_16"],
                "day18" : ["SL1_18 - SL3_18"],
                "SL1_day18-14" : ["SL1_18 - SL1_14"],
                "SL1_day18-16" : ["SL1_18 - SL1_16"],
                "SL1_day14-16" : ["SL1_14 - SL1_16"],
                "SL3_day18-14" : ["SL3_18 - SL3_14"],
                "SL3_day18-16" : ["SL3_18 - SL3_16"],
                "SL3_day14-16" : ["SL3_14 - SL3_16"],
            }
        },
        "line_separated" : {
            "SL1" : {
                "all" : ["F_14 - M_14","F_16 - M_16", "F_18 - M_18"],
                "early" : ["F_14 - M_14","F_16 - M_16"],
                "late" : ["F_16 - M_16", "F_18 - M_18"]
            },
            "SL3" : {
                "all" : ["F_14 - M_14","F_16 - M_16", "F_18 - M_18"],
                "early" : ["F_14 - M_14","F_16 - M_16"],
                "late" : ["F_16 - M_16", "F_18 - M_18"]
            }
        },
        "day_separated" : {
            "day14" : {
                "SL1" : ["F_1 - M_1"],
                "SL3" : ["F_3 - M_3"],
                "both" : ["F_1 - M_1","F_3 - M_3"]
            },
            "day16" : {
                "SL1" : ["F_1 - M_1"],
                "SL3" : ["F_3 - M_3"],
                "both" : ["F_1 - M_1","F_3 - M_3"]
            },
            "day18" : {
                "SL1" : ["F_1 - M_1"],
                "SL3" : ["F_3 - M_3"],
                "both" : ["F_1 - M_1","F_3 - M_3"]
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


def plot_enriched_GO_overlap(GO_Terms_dict:dict, plot_filename:str, plot_title  = "", min_overlap = 0, contrasts_of_interest = {}):
    """
    make an upseetplot of sig DE genes within each data separation
    """
    data = upsetplot.from_contents(GO_Terms_dict)

    if len(contrasts_of_interest) == 0:
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
    else:
        # make a table that gives a GO term list of all the contrasts of interest as they appear in the upsetplot
        table_filename = plot_filename.replace(".png", ".txt")
        membership_cols = data.index.names
        df = data.reset_index()
        with open(table_filename, "w") as table_out:
            for subset_name , contrasts_list in contrasts_of_interest.items():

                ## get intersection GO list
                intersection = df[df[contrasts_list].all(axis=1) & (df[membership_cols].sum(axis=1) == len(contrasts_list))] # all contrasts_list are true & sum of all is the same length as contrasts_list (therefore all others must be false)
                intersection = intersection["id"].tolist()

                ## write to table
                print(f"\t - {subset_name} : {len(intersection)} GO:terms")
                intersection_string = ",".join(intersection)
                contrasts_string = " AND ".join(contrasts_list)
                outfile_line = f"{subset_name}\t{contrasts_string}\t{intersection_string}\n"
                table_out.write(outfile_line)
        print(f"GO:terms of specified intersections saved as: {table_filename}")



if __name__ == "__main__":
    username = "miltr339"
    go_tables_paths = get_tables(username=username)
    contrasts_of_interest = get_interesting_GO_overlap_lists()
    out_path_figs = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures_python/GO_enrichment"
    
    ## don't plot only get lists of all the GO terms in each interaction
    plot = False

    ### TODO analyze go enrichment
    # either: REVIGO online tool: http://revigo.irb.hr/
    # or: GO3 library: https://go3.readthedocs.io/en/latest/examples.html
    #   * requires a GAF file https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/ unsure if i have that
    #   * unsure what to do now
    # if plot==False:
    #     go3.load_go_terms("")
    #     annots = go3.load_gaf(f"/Users/{username}/work/c_maculatus/C_mac_eggnog_diamond.emapper.annotations_geneIDs")

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
                
                if plot:
                    ## plot intersection sizes
                    plot_enriched_GO_overlap(sig_GO_terms, plot_filename= f"{out_path_figs}/upsetplot_GO_terms_{separation}_{category}.png", plot_title = plot_title, min_overlap = min_overlap, plot=plot)
                else:
                    ## do semantic clustering on contrasts of interest
                    contrasts_dict = contrasts_of_interest[separation][category]
                    plot_enriched_GO_overlap(sig_GO_terms, plot_filename= f"{out_path_figs}/upsetplot_GO_terms_{separation}_{category}.png", plot_title = plot_title, min_overlap = min_overlap, contrasts_of_interest=contrasts_dict)

            if plot:
                min_overlap = 5
                if min_overlap>0:
                    plot_title = f"{separation_}:all categories\nsig. GO:terms overlap\n(min. overlap size: {min_overlap})"
                else:
                    plot_title = f"{separation_}:all categories\nsig. GO:terms overlap"
                plot_enriched_GO_overlap(sep_GO_terms, plot_filename= f"{out_path_figs}/upsetplot_GO_terms_{separation}_all.png", plot_title = plot_title, min_overlap = min_overlap, plot=plot)
        
            