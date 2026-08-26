"""
Analyze the significantly enriched GO terms in all the contrasts
"""

import pandas as pd
import upsetplot
import matplotlib.pyplot as plt
from goatools.obo_parser import GODag
import warnings

def get_tables(username="miltr339"):
    """
    Tables from GO enrichment analysis in PhD_chapter4/src/GO_enrichment.Rmd
    """
    tables_dir = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/"

    out_dict = {
        "no_separation" : {
            "line_random" : {
                "(day14_F - day14_M) - (day16_F - day16_M) - (day18_F - day18_M)" : f"{tables_dir}full_dataset_line_ignored_day_sex_interaction_GO_enrichment.txt",
                "day14_F - day14_M" : f"{tables_dir}full_dataset_line_ignore_GO_enrichmentd_day14_F-M.csv",
                "day16_F - day16_M" : f"{tables_dir}full_dataset_line_ignore_GO_enrichmentd_day16_F-M.csv",
                "day18_F - day18_M" : f"{tables_dir}full_dataset_line_ignore_GO_enrichmentd_day18_F-M.csv",
            }
        },
        "sex_separated" : {
            "females" :{
                "SL1_14 - SL3_14" : f"{tables_dir}GO_enrichment_F_1-3_day14.txt",
                "SL1_16 - SL3_16" : f"{tables_dir}GO_enrichment_F_1-3_day16.txt",
                "SL1_18 - SL3_18" : f"{tables_dir}GO_enrichment_F_1-3_day18.txt",
                "SL3_18 - SL3_14" : f"{tables_dir}GO_enrichment_F_SL3_18_14.txt",
                # "day14 - day16" : f"{tables_dir}only_F_no_line_GO_enrichment_day_14-16.txt",
                "day14 - day18" : f"{tables_dir}only_F_no_line_GO_enrichment_day_14-18.txt",
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
                "day14 - day16" : f"{tables_dir}only_M_no_line_GO_enrichment_day_14-16.txt",
                "day16 - day18" : f"{tables_dir}only_M_no_line_GO_enrichment_day_16-18.txt",
                "day14 - day18" : f"{tables_dir}only_M_no_line_GO_enrichment_day_14-18.txt",
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
    Plot one upsetplot for each category within a separation, labels are the dict keys
    """    
    out_dict = {
        "no_separation" : {
            "line_random" : {
                "all days" : ["day14_F - day14_M","day16_F - day16_M","day18_F - day18_M"]
            }
        },
        "sex_separated" : {
            "females" :{
                "day14" : ["SL1_14 - SL3_14"],
                "day16" : ["SL1_16 - SL3_16"],
                "day18" : ["SL1_18 - SL3_18"],
                "SL3_day18-14" : ["SL3_18 - SL3_14"],
                "day18-14" : ["day18 - day14"],
                "day14-18" : ['day14 - day18'],
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
                "day18-14" : ["day18 - day14"],
                "day18-16" : ["day18 - day16"],
                "day14-16" : ["day14 - day16"],
                "day16-18" : ["day16 - day18"],
                "day14-18" : ["day14 - day18"],
                "all_days" : ["males: day14 - day16","males: day16 - day18","males: day14 - day18","females: day14 - day18"]
            }
        },
        "line_separated" : {
            "SL1" : {
                "day14" : ["F_14 - M_14"],
                "day16" : ["F_16 - M_16"],
                "day18" : ["F_18 - M_18"],
                "all" : ["F_14 - M_14","F_16 - M_16", "F_18 - M_18"],
                "early" : ["F_14 - M_14","F_16 - M_16"],
                "late" : ["F_16 - M_16", "F_18 - M_18"]
            },
            "SL3" : {
                "day14" : ["F_14 - M_14"],
                "day16" : ["F_16 - M_16"],
                "day18" : ["F_18 - M_18"],
                "all" : ["F_14 - M_14","F_16 - M_16", "F_18 - M_18"],
                "early" : ["F_14 - M_14","F_16 - M_16"],
                "late" : ["F_16 - M_16", "F_18 - M_18"]
            }
        },
        "day_separated" : {
            "day14" : {
                "males" : ["M_1 - M_3"],
                "females" : ["F_1 - F_3"],
                "SL1" : ["F_1 - M_1"],
                "SL3" : ["F_3 - M_3"],
                "both" : ["F_1 - M_1","F_3 - M_3"]
            },
            "day16" : {
                "males" : ["M_1 - M_3"],
                "females" : ["F_1 - F_3"],
                "SL1" : ["F_1 - M_1"],
                "SL3" : ["F_3 - M_3"],
                "both" : ["F_1 - M_1","F_3 - M_3"]
            },
            "day18" : {
                "males" : ["M_1 - M_3"],
                "females" : ["F_1 - F_3"],
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
        if len(plot_title)>0:
            plt.title(plot_title)
        plt.tight_layout()
        plt.savefig(plot_filename, dpi = 300, transparent = True)
        print(f"plot saved in current working directory as: {plot_filename}")
        plt.clf()
        plt.cla()
        plt.close()
    else:
        # make a table that gives a GO term list of all the contrasts of interest as they appear in the upsetplot
        print(f"make no plot only GO lists")
        membership_cols = data.index.names
        df = data.reset_index()
        table_filename = plot_filename.replace(".png", f".txt")
        with open(table_filename, "w") as table_out:

            for subset_name , contrasts_list in contrasts_of_interest.items():
                print(f"{subset_name} : [{contrasts_list}]")
                ## get intersection GO list
                try:
                    intersection = df[df[contrasts_list].all(axis=1) & (df[membership_cols].sum(axis=1) == len(contrasts_list))] # all contrasts_list are true & sum of all is the same length as contrasts_list (therefore all others must be false)
                except:
                    print(f"data columns: {df.columns}")
                    print(f"contrasts: {contrasts_list}")
                    print(f"member cols: {membership_cols}")
                    raise RuntimeError(f"some contrasts not found in columns")
                intersection = intersection["id"].tolist()

                ## write to table with functions
                table_filename_func = plot_filename.replace(".png", f"_{subset_name}.txt")
                table_filename_func = table_filename_func.replace(" ", "")
                GO_list = []
                for GO_term in intersection:
                    try:
                        go_name = go[GO_term].name 
                    except:
                        go_name = "OBSOLETE" # obsolete go terms are not read into 'go' and therefore should be skipped
                    GO_list.append(f"{GO_term}\t{go_name}\n")
                with open(table_filename_func, "w") as table_out_func:
                    for GO_line in GO_list:
                        table_out_func.write(GO_line)

                ## write just lists to summary outfile
                print(f"\t - {subset_name} : {len(intersection)} GO:terms ({table_filename_func})")
                contrasts_string = " AND ".join(contrasts_list)
                intersection_string = ",".join(intersection)
                outfile_line = f"{subset_name}\t{contrasts_string}\t{intersection_string}\n"
                table_out.write(outfile_line)

        
        print(f"GO:terms of specified intersections saved as: {table_filename}")


def make_table_sig_GO_terms(infile_path):
    outfile_path = infile_path.replace(".csv", "_sig_GO.txt")
    GO_file = pd.read_csv(infile_path)
    count_sig = 0
    with open(outfile_path, "w") as outfile:
        outfile.write(outfile_path.split("PhD_chapter4/")[-1]+"\n\n")
        for GO_term,pval in zip(GO_file["GO.ID"], GO_file["classicFisher"]):
            if pval>0.05:
                continue
            try:
                go_name = go[GO_term].name 
            except:
                go_name = "OBSOLETE" # obsolete go terms are not read into 'go' and therefore should be skipped
            GO_line = f"{GO_term}\t{go_name}\n"
            outfile.write(GO_line)
            count_sig+=1
    
    print(f"file with {count_sig} sig GO terms written to \n{outfile_path}")



def add_functional_information_to_geneIDs(infile_path, annotation_file, IDs_list = []):

    annotation_header = ["query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"]
    annotation_df = pd.read_csv(annotation_file, sep="\t", comment="#", names=annotation_header)

    if IDs_list == []:
        outfile_path = infile_path.replace(".txt", "_functional_annotation.txt")
        with open(infile_path, "r") as infile, open(outfile_path, "w") as outfile:
            infile_lines = infile.readlines()
            header = infile_lines[0].strip()
            header_ = f"{header}\tDescription\n"
            # outfile.write(header_)
            for i, line_full in enumerate(infile_lines[1:]):
                line = line_full.strip()
                geneID = line.split("\t")[-1]
                try:
                    desc = annotation_df.loc[annotation_df["query"] == geneID, "Description"].iloc[0]
                except:
                    desc = "---not_annotated---"
                if desc == "-":
                    desc = "no_description"
                if desc == "---not_annotated---":
                    name = "na"
                else:
                    try:
                        name = annotation_df.loc[annotation_df["query"] == geneID, "Preferred_name"].iloc[0]
                    except:
                        name = "-"

                print(f"{i+1}\t{geneID}\t{name}\t{desc}")
                outfile.write(f"{geneID}\t{name}\t{desc}\n")
    else:
        outfile_path = infile_path
        with open(outfile_path, "w") as outfile:
            # header_ = f"geneID\tDescription\n"
            # outfile.write(header_)
            for i, geneID in enumerate(IDs_list):
                try:
                    desc = annotation_df.loc[annotation_df["query"] == geneID, "Description"].iloc[0]
                except:
                    desc = "---not_annotated---"
                if desc == "-":
                    desc = "no_description"
                if desc == "---not_annotated---":
                    name = "na"
                else:
                    try:
                        name = annotation_df.loc[annotation_df["query"] == geneID, "Preferred_name"].iloc[0]
                    except:
                        name = "-"

                print(f"{i+1}\t{geneID}\t{name}\t{desc}")
                outfile.write(f"{geneID}\t{name}\t{desc}\n")

    print(f"file written to: \n{outfile_path}")


def get_genes_with_GO(infile_path, annotation_file, GOs_list = [], plot_file = "plot.png", plot_title = ""):

    annotation_header = ["query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"]
    annotation_df = pd.read_csv(annotation_file, sep="\t", comment="#", names=annotation_header)

    sigIDs_df = pd.read_csv(infile_path, sep=",")
    sigIDs = sigIDs_df.loc[sigIDs_df["sig_DE"]==1]["geneID"].tolist()

    GO_set = set(GOs_list)
    GO_geneIDs = {GO_term : [] for GO_term in GO_set}

    for i, geneID in enumerate(sigIDs):
        try:
            GO_terms = annotation_df.loc[annotation_df["query"] == geneID, "GOs"].iloc[0]
            if GO_terms != "-":
                gene_terms = set(GO_terms.split(","))
                shared = gene_terms & GO_set
                if len(shared)>0:
                    print(f"{i+1}\t{geneID}\t{shared}") 
                    for shared_GO in shared:
                        GO_geneIDs[shared_GO].append(geneID)       
        except:
            continue
            print(f"{i+1}\t{geneID}\tparsing failed!")

    from Y_expression_quantification import get_counts_paths,samples_group,plot_counts_sum_sets
    count_files = get_counts_paths(username=username)
    samples_group_dict = samples_group()
    print(f"{plot_title}")
    plot_counts_sum_sets(counts_table=count_files["no_log"], geneIDs_lists_dict = GO_geneIDs, 
                        outfile_name = plot_file, y_label= "normalized counts", errorbars=True, samples_group_dict = samples_group_dict, plot_title=plot_title)




if __name__ == "__main__":

    warnings.filterwarnings("ignore")
    username = "miltr339"
    go_tables_paths = get_tables(username=username)
    contrasts_of_interest = get_interesting_GO_overlap_lists()
    out_path_figs = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures_python/GO_enrichment"
    
    ########## if False: don't plot only get lists of all the GO terms in each interaction
    plot = False
    ##########

    # get functional information about GO terms
    if True:
        print(f"read GO-terms definitions...")
        go = GODag(f"/Users/{username}/work/go-basic.obo") # wget https://current.geneontology.org/ontology/go-basic.obo
        print(f"...done!")
        

    if False:
        for separation, seps_dict in go_tables_paths.items():
            if separation!="no_separation":
                print(f"ignore {separation}")
                continue
            ## if no special stuff below then just do an empty string
            suffix = "_line_ign" # for making only the model with line fixed effect and day contrasts plots, 

            print(f"\n=========================== {separation} ===========================")

            contrast_keys_list = []
            contrast_sep_list = []
            sep_GO_terms = {}
            for category, paths_dict in seps_dict.items():
                print(f"\n ------------------- {category} -------------------")
                
                sig_GO_terms = {}
                for contrast, table_path in paths_dict.items():
                    
                    print(f"{table_path}")
                    ### only plot the day ones
                    if len(suffix)==0 or "day" not in contrast:
                        continue
                    if "(" not in contrast:
                        contrast_keys_list.append(contrast)
                        GO_list = read_R_GO_table(table_filepath=table_path)
                        sig_GO_terms[contrast] = GO_list
                        sep_GO_terms[f"{category}: {contrast}"] = GO_list
                    if "day" in contrast:
                        contrast_sep_list.append(f"{category}: {contrast}")

                min_overlap = 5
                separation_  = separation.replace("_", "-")
                if len(suffix)>0:
                    min_overlap = 0
                    if category=="line_random":
                        plot_title = f"GO terms"
                    else:    
                        plot_title = f"{category}\nGO terms"
                elif min_overlap>0:
                    plot_title = f"{separation_}:{category}\nsig. GO:terms overlap\n(min. overlap size: {min_overlap})"
                else:
                    plot_title = f"{separation_}:{category}\nsig. GO:terms overlap"
                
                if plot==True and len(sig_GO_terms)>1:
                    ## plot intersection sizes
                    plot_enriched_GO_overlap(GO_Terms_dict = sig_GO_terms, plot_filename= f"{out_path_figs}/upsetplot_GO_terms_{separation}_{category}{suffix}.png", plot_title = plot_title, min_overlap = min_overlap)
                elif plot==False:
                    contrasts_dict = {name : contrast for name, contrast in contrasts_of_interest[separation][category].items() if any([True for c in contrast_keys_list if c in contrast])}
                    if len(contrasts_dict)==0:
                        print(f"no GO-enrichment lists for no contrasts of interest in {category}")
                        continue
                    print(f"CONTRASTS: {contrasts_dict}")
                    plot_enriched_GO_overlap(GO_Terms_dict = sig_GO_terms, plot_filename= f"{out_path_figs}/upsetplot_GO_terms_{separation}_{category}{suffix}.png", plot_title = plot_title, min_overlap = min_overlap, contrasts_of_interest=contrasts_dict)
                    contrast_keys_list = []
                else:
                    print(f"no GO-enrichment overlap for only one test: {sig_GO_terms.keys()}")

            if separation == "no_separation":
                continue
            print(f"\n ------------------- all categories in '{separation}' -------------------")

            if plot:
                min_overlap = 5
                if len(suffix)>0:
                    min_overlap = 0
                    plot_title = f"GO terms"
                elif min_overlap>0:
                    plot_title = f"{separation_}:all categories\nsig. GO:terms overlap\n(min. overlap size: {min_overlap})"
                else:
                    plot_title = f"{separation_}:all categories\nsig. GO:terms overlap"
                plot_enriched_GO_overlap(GO_Terms_dict = sep_GO_terms, plot_filename= f"{out_path_figs}/upsetplot_GO_terms_{separation}_all{suffix}.png", plot_title = plot_title, min_overlap = min_overlap)

                if separation=="day_separated":
                    print(sep_GO_terms.keys())
                    # sep_GO_terms_ = {sep : sep_GO_terms[sep] for sep in ['males: SL1_14 - SL3_14', 'males: SL1_16 - SL3_16', 'males: SL1_18 - SL3_18']}
                    sep_GO_terms_ = {sep : l_ for sep,l_ in sep_GO_terms.items() if "M_1 - M_3" in sep}
                    min_overlap_ = 0
                    plot_title_ = f""
                    plot_enriched_GO_overlap(GO_Terms_dict = sep_GO_terms_, plot_filename= f"{out_path_figs}/upsetplot_GO_terms_{separation}_male_line_bias.png", plot_title = plot_title_, min_overlap = min_overlap_)
            
            
            elif separation=="sex_separated" and plot==False:
                # Merge the individual day contrasts with the list of intersections from the relevant contrasts list
                contrasts_dict_ind = {c : [c] for c in contrast_sep_list}
                contrasts_dict_int = {name : contrast for name, contrast in contrasts_of_interest[separation]["males"].items() if any([True for c in contrast if "day" in c]) and len(contrast)>1}
                contrasts_dict = contrasts_dict_ind | contrasts_dict_int # merge dicts
                print(f"CONTRASTS: {contrasts_dict}")
                # contrasts_dict = {name : contrast for name, contrast in contrasts_of_interest[separation][category].items() if any([True for c in contrast_sep_list if c in contrast])}
                plot_enriched_GO_overlap(GO_Terms_dict = sep_GO_terms, plot_filename= f"{out_path_figs}/upsetplot_GO_terms_{separation}_all{suffix}.png", plot_title = plot_title, min_overlap = min_overlap, contrasts_of_interest=contrasts_dict)


    if False:
        ## make an output list based on the R csv table with full function names and only p<0.05
        all_days=f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/day_separated_all_lines_sex_bias_overlap_all_days_sigIDs.csv"
        make_table_sig_GO_terms(infile_path=all_days)
        late_days=f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/day_separated_all_lines_sex_bias_overlap_day1618_sigIDs.csv"
        make_table_sig_GO_terms(infile_path=late_days)
        early_days=f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/day_separated_all_lines_sex_bias_overlap_day14_sigIDs.csv"
        make_table_sig_GO_terms(infile_path=early_days)

    if False:
        ## add gene annotation information to a previously generated list of sig. geneIDs
        day_separated_line_bias_overlap = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/day_separated_male_line_bias_overlap_sigIDs.txt"
        annotation_path = f"/Users/{username}/work/c_maculatus/C_mac_eggnog_diamond.emapper.annotations_geneIDs"

        add_functional_information_to_geneIDs(infile_path = day_separated_line_bias_overlap, annotation_file=annotation_path)
        
    if True:
        SB_LB_genes = {
            "day14" : {
                "SL1" : ['gene-225158', 'gene-372264'],
                "SL3" : ['gene-371957', 'gene-39545', 'gene-130096', 'gene-166391', 'gene-130081', 'gene-158197', 'gene-414353', 'gene-115312', 'gene-122220', 'gene-204669', 'gene-143368', 'gene-398993', 'gene-149137', 'gene-100036', 'gene-227137', 'gene-246615', 'gene-406603', 'gene-277078', 'gene-372264'],
            },
            "day16" : {
                "SL1" : ['gene-372264', 'gene-240602', 'gene-390687', 'gene-24278', 'gene-13404', 'gene-63245', 'gene-55869', 'gene-23834', 'gene-23884', 'gene-122692', 'gene-279676', 'gene-104371', 'gene-329410', 'gene-431362', 'gene-411056', 'gene-127707', 'gene-80062', 'gene-351334', 'gene-48598', 'gene-30595', 'gene-90918'],
                "SL3" : ['gene-372264', 'gene-24185', 'gene-24290', 'gene-23840', 'gene-23597', 'gene-24203', 'gene-24120', 'gene-24132', 'gene-23538', 'gene-423321', 'gene-24088', 'gene-24221', 'gene-24167', 'gene-23514', 'gene-23365', 'gene-15763', 'gene-23689', 'gene-23413', 'gene-7220', 'gene-24052', 'gene-24079', 'gene-27466', 'gene-410366', 'gene-12075', 'gene-349163'],
            }
        }
        annotation_path = f"/Users/{username}/work/c_maculatus/C_mac_eggnog_diamond.emapper.annotations_geneIDs"
        for day, line_dict in SB_LB_genes.items():
            for line, IDs_list in line_dict.items():
                print(f"\n ------------------- {day}:{line} ({len(IDs_list)}) -------------------")

                day_separated_line_bias_overlap = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/day_separated_{day}_sex_{line}_and_line_M_biased_genes_functions.txt"
                add_functional_information_to_geneIDs(infile_path = day_separated_line_bias_overlap, annotation_file=annotation_path, IDs_list=IDs_list)

    if False:
        sig_IDs_list = {
            ## geneIDs that are significant in the day separated line-by-sex interaction
            ## from PhD_chapter4/data/sig_DE_genes_lists/sig_DE_list_day14_F-M_by_1-3.txt and other days
            "day_separated" : {
                "day14" : ["gene-237342","gene-181562","gene-58400"],
                "day18" : ["gene-426041","gene-48523","gene-73742","gene-2355","gene-237494","gene-279406","gene-237342","gene-276797","gene-367071"],
            }
        }
        annotation_path = f"/Users/{username}/work/c_maculatus/C_mac_eggnog_diamond.emapper.annotations_geneIDs"
        for sep, day_dict in sig_IDs_list.items():
            for day, IDs_list in day_dict.items():
                print(f"\n ------------------- {day} ({len(IDs_list)}) -------------------")

                day_separated_line_bias_overlap = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/day_separated_{day}_sex_by_line_interaction_sig_genes_functions.txt"
                add_functional_information_to_geneIDs(infile_path = day_separated_line_bias_overlap, annotation_file=annotation_path, IDs_list=IDs_list)

    if False:
        # check for smoothened signalling pathway in full data line ignored early and late day contrasts
        full_GO_path = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/full_dataset_line_ignore_GO_enrichmentd_M_16_18.csv"
        sig_genes_path = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/full_dataset_line_ignored_M_16_18.txt"
        annotation_path = f"/Users/{username}/work/c_maculatus/C_mac_eggnog_diamond.emapper.annotations_geneIDs"
        plot_dir = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures_python/counts_time_series"
        
        ### the function plots separate lines for each GO terms, but if a gene is annotated with more than one GO term, it is used for the median/SEM calculation in both!
        if False:
            smoothened_GO  = ["GO:0007224","GO:0008589","GO:0045879"]
            get_genes_with_GO(infile_path = sig_genes_path, annotation_file=annotation_path, GOs_list=smoothened_GO, 
                plot_file=f"{plot_dir}/full_dataset_line_ignored_M_16_18_smoothened_signalling_GO_terms.png", plot_title="GO-terms related to smoothened pathway")
        if True:
            aging_GO = ["GO:0007568","GO:0010259"]
            get_genes_with_GO(infile_path = sig_genes_path, annotation_file=annotation_path, GOs_list=aging_GO, 
                plot_file=f"{plot_dir}/full_dataset_line_ignored_M_16_18_aging_GO_terms.png", plot_title="GO-terms related to aging")

    if False:
        ## check ecdysone GO-terms
        GO_terms_ecdysone = {
            "14" : ["GO:0008205","GO:0006697","GO:0035072","GO:0035075"],
            "16" : ["GO:0008205","GO:0006697","GO:0035072"],
            "18" : ["GO:0008205"]
        }
        GO_terms_juvenile_hormone = ["GO:0006719","GO:0006716"] # only gene-540

        annotation_path = f"/Users/{username}/work/c_maculatus/C_mac_eggnog_diamond.emapper.annotations_geneIDs"
        plot_dir = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures_python/counts_time_series"

        for day in ["14","16","18"]:
            GO_terms_day_path = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/full_dataset_line_ignore_GO_enrichmentd_day{day}_F-M.csv"
            sig_genes_path = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/sig_DE_genes_lists/full_dataset_line_ignored_day{day}_F-M.txt"

            if True:
                GO_list = GO_terms_ecdysone[day]
                get_genes_with_GO(infile_path = sig_genes_path, annotation_file=annotation_path, GOs_list=GO_list, 
                    plot_file=f"{plot_dir}/full_dataset_line_ignored_ecdysone_day{day}_F-M.png", plot_title=f"GO-terms related to ecdysone on day {day}")
            
            # GO_list = GO_terms_juvenile_hormone
            # get_genes_with_GO(infile_path = sig_genes_path, annotation_file=annotation_path, GOs_list=GO_list, 
            #    plot_file=f"{plot_dir}/full_dataset_line_ignored_jh_day{day}_F-M.png", plot_title=f"GO-terms related to juvenile hormone on day {day}")

    if False:
        # plot the one gene for juvenile hormone
        from Y_expression_quantification import get_counts_paths,samples_group,plot_counts_sum_sets
        count_files = get_counts_paths(username=username)
        samples_group_dict = samples_group()

        GO_terms_juvenile_hormone = ["GO:0006719","GO:0006716"] # only gene-540, annotated as 'COesterase, Belongs to the type-B carboxylesterase lipase family'

        plot_counts_sum_sets(counts_table=count_files["no_log"], geneIDs_lists_dict = {"gene-540" : ["gene-540"]}, outfile_name = f"{plot_dir}/full_dataset_line_ignored_jh_all_days.png", 
            y_label= "normalized counts", errorbars=True, samples_group_dict = samples_group_dict, plot_title=f"gene annotated with juvenile hormone-related function")