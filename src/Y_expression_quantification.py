"""
Plot gene expression and compare Y and A
"""

from basic_plotting import plot_counts

def get_counts_paths(username = "milena"):
    counts_paths = {
        "no_log" : f"/Users/{username}/work/PhD_code/PhD_chapter4/data/gene_counts_normalized_nolog.tsv",
        "yes_log": f"/Users/{username}/work/PhD_code/PhD_chapter4/data/gene_counts_normalized.tsv",
    }
    return counts_paths

def get_y_information(username = "milena"):
    sex_chromosomes_superscaffolded = { 
        "X" : ['scaffold_10','scaffold_14','scaffold_23','scaffold_31','scaffold_34','scaffold_83'],
        "Y" : ['scaffold_26','scaffold_48','scaffold_103','scaffold_112','scaffold_164']}
    
    # grep "scaffold_26\t" C_mac_superscaffolded_yTor_LomeRNA_eggnog.gff | grep "\tgene" | cut -f9
    # no genes on 112 and 164
    y_contigs = {
        'scaffold_26' : ["gene-371557","gene-371566","gene-371572","gene-371578","gene-371593","gene-371599","gene-371605","gene-371611","gene-371617","gene-371627","gene-371633","gene-371651","gene-371660","gene-371666","gene-371676","gene-371686","gene-371698","gene-371704","gene-371719","gene-371725","gene-371738","gene-371744","gene-371750","gene-371763","gene-371769","gene-371775","gene-371781","gene-371787","gene-371793","gene-371799","gene-371805","yTor-A","yTor-B","gene-371844","yTor-C","gene-371883","gene-371889","gene-371898","gene-371907","gene-371913","gene-371922","gene-371936","gene-371957","gene-371984","gene-372026","gene-372053","gene-372068","gene-372207","gene-372216","gene-372237","gene-372243","gene-372264","gene-372302"],
        'scaffold_48' : ["gene-402962","gene-402968","gene-402989"],
        'scaffold_103' : ["gene-424042","gene-424048","gene-424063","gene-424069","gene-424075","gene-424085","gene-424095"],
        'scaffold_112' : [],
        'scaffold_164' : [],
        "all" : ["gene-402962","gene-402968","gene-402989","gene-424042","gene-424048","gene-424063","gene-424069","gene-424075","gene-424085","gene-424095","gene-371557","gene-371566","gene-371572","gene-371578","gene-371593","gene-371599","gene-371605","gene-371611","gene-371617","gene-371627","gene-371633","gene-371651","gene-371660","gene-371666","gene-371676","gene-371686","gene-371698","gene-371704","gene-371719","gene-371725","gene-371738","gene-371744","gene-371750","gene-371763","gene-371769","gene-371775","gene-371781","gene-371787","gene-371793","gene-371799","gene-371805","yTor-A","yTor-B","gene-371844","yTor-C","gene-371883","gene-371889","gene-371898","gene-371907","gene-371913","gene-371922","gene-371936","gene-371957","gene-371984","gene-372026","gene-372053","gene-372068","gene-372207","gene-372216","gene-372237","gene-372243","gene-372264","gene-372302"],
        }
    return sex_chromosomes_superscaffolded,y_contigs



if __name__ == "__main__":
    
    username = "milena"
    count_files = get_counts_paths(username=username)
    ex_chromosomes,y_contigs = get_y_information(username=username)
    out_path = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis"

    plot_counts(counts_table=count_files["no_log"], geneIDs_list = y_contigs["all"], 
                outfile_name = f"{out_path}/y_genes_mean_expression.png", y_label= "normalized counts" ,mean_per_sample=True)