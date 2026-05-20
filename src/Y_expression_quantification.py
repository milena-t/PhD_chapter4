"""
Plot gene expression and compare Y and A
"""

from basic_plotting import plot_counts

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats


def plot_counts_sum_sets(counts_table:str, geneIDs_lists_dict:dict, outfile_name:str, y_label="normalized counts", errorbars = False, samples_group_dict = {}):
    """
    plot the counts of the geneIDs in all samples
    except samples where all samples have zero counts
    gene_IDs_lists_dict is a dict of {'label' : [geneID_list]}, 
    and the mean counts with error bars of each list are plotted
    """
    assert len(geneIDs_lists_dict)>0
    counts_df = pd.read_csv(counts_table, sep="\t", index_col=0)
        
    # check if there are nonexpressed genes in the input list
    headers = counts_df.columns.tolist()

    # subsetting the df for each gene set to plot as separate lines
    dfs_dict = { geneID_label : counts_df.filter(items=geneIDs_list, axis="index") for geneID_label, geneIDs_list in geneIDs_lists_dict.items()}

    # check if there are samples that show no expression for any gene in any input list
    nonzero_samples2 = []
    for lab,df in dfs_dict.items():
        sums = df.sum()
        nonzero_samples2.extend(sums[sums>0].index.to_list())
    nonzero_samples = list(set(nonzero_samples2))
    print(f"out of {len(headers)} there are {len(nonzero_samples)} samples that have at least one count in one gene!")# \n{nonzero_samples}")
    
    dfs_dict = { geneID_label : df.filter(items=nonzero_samples, axis="columns") for geneID_label, df in dfs_dict.items()}

    ### plotting
    fig, ax = plt.subplots(1,1, figsize=(20, 10)) # for more than three rows

    fs = 25
    ps = fs*0.8# * 1/len(geneIDs_lists_dict) # point size
    lw=2
    linest = ":"

    # sort nonzero samples so that F/M and lines are separated and the plot is easier to interpret 
    females = sorted([sample for sample in nonzero_samples if "-F_" in sample])
    males = sorted([sample for sample in nonzero_samples if "-M_" in sample])
    nonzero_samples_sorted = males + females
    assert len(nonzero_samples) == len(nonzero_samples_sorted)

    colors_list = [
"#d2416e",
"#57b158",
"#c35dbd",
"#a6b541",
"#7965ca",
"#cd9445",
"#678ccc",
"#c85c3f",
"#4caf9a",
"#bd6b8f",
"#717935"]
    c = 0

    for list_label, gene_counts in dfs_dict.items():

        outfile_table = outfile_name.replace(".png", f"_{list_label}.txt")
        with open(outfile_table, "w" ) as out_table:

            print(f"--------------------- {list_label} ---------------------")
            ## make medians and standard errors
            if samples_group_dict == {}:
            # plot all samples individually, no groups
                medians_list = []
                errors_list = []
                tick_labels = [sample for sample  in nonzero_samples_sorted] # make fully new array, with just tick_labels=nonzero_samples_sorted i just kind of get a pointer and it messes up subsequent loops
                tick_pos = [i for i in range(len(nonzero_samples_sorted))]

                for i in range(len(nonzero_samples_sorted)):
                    sample = nonzero_samples_sorted[i]
                    try:
                        curr_counts = gene_counts[sample].to_list()
                        curr_med = np.median(curr_counts)
                        curr_sem = stats.sem(curr_counts)
                        count_expressed = len([count for count in curr_counts if count >0])
                    except:
                        curr_counts = []
                        curr_med = np.nan
                        curr_sem = np.nan
                        count_expressed = []

                    medians_list.append(curr_med)
                    errors_list.append(curr_sem)
                        
                    sample_ = sample.replace("WJ-3841-","SL").split("_")[0]
                    if len(geneIDs_lists_dict) ==1:
                        tick_labels[i] = f"{sample_} ({count_expressed})"
                    else:
                        tick_labels[i] = f"{sample_}"
                    
                    sample_outstr = sample.replace("WJ-3841-","SL")
                    out_str = f"{sample_outstr} ({len(curr_counts)} genes, {count_expressed} are expressed) :\t median: {curr_med:.3f}; std.err: {curr_sem:.3f}"
                    out_table.write(f"{out_str}\n")
                    # print(out_str)
                
                tick_cols = ["#000000" if "M" in sample else "#8A8A8A" for sample in nonzero_samples_sorted ]
                # for s,m in zip(tick_labels,medians_list):
                #     print(f"{s} : {m}")

            else:
                # use sample groups as defined in sample_groups_dict
                sample_group_labels = list(samples_group_dict.keys())
                medians_list = []
                errors_list = []
                tick_labels = [sample for sample  in sample_group_labels]
                tick_pos = [i for i in range(len(sample_group_labels))]

                for i in range(len(sample_group_labels)):
                    sample_group = sample_group_labels[i]
                    sample_group_list = samples_group_dict[sample_group]
                    curr_counts = []
                    count_samples = 0
                    count_expressed = 0
                    for i in range(len(sample_group_list)):
                        sample = sample_group_list[i]
                        try:
                            curr_counts.extend(gene_counts[sample].to_list())
                            count_samples +=1
                        except:
                            pass
                        
                    if len(curr_counts)>0:
                        curr_med = np.median(curr_counts)
                        curr_sem = stats.sem(curr_counts)
                        count_expressed = len([count for count in curr_counts if count >0])
                    else:
                        curr_med = np.nan
                        curr_sem = np.nan

                    medians_list.append(curr_med)
                    errors_list.append(curr_sem)

                    out_str =f"{sample_group} ({count_expressed} expr. genes from {count_samples} samples) :\t median: {curr_med:.3f}; std.err: {curr_sem:.3f}"
                    # print(out_str)
                    out_table.write(f"{out_str}\n")

                tick_cols = ["#000000" if "M" in sample else "#8A8A8A" for sample in tick_labels ]

                # for median, sem, sample in zip(medians_list, errors_list, tick_labels):
                #     print(f"{sample} : {median:.3f} +- {sem:.3f}")
            
            print(f"table outfile written to: {outfile_table}")

            if errorbars:
                ax.errorbar(tick_pos, medians_list, xerr = 0, yerr = errors_list, color=colors_list[c], linewidth =lw,
                            marker = ".", markersize=ps, linestyle = linest, label=list_label)
            else:
                ax.plot(tick_pos, medians_list,color=colors_list[c], linewidth =lw, 
                        marker = ".", markersize=ps, linestyle = linest, label=list_label)
            c +=1

    if "log" not in y_label:
        ymin, ymax = ax.get_ylim()
        ax.set_ylim([-0.1,ymax])

    
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels)
    for tick_label, color in zip(ax.get_xticklabels(), tick_cols):
        tick_label.set_color(color)
    
    ax.tick_params(axis='x', labelsize=fs*0.75,labelrotation=90)#, colors)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_ylabel(f"{y_label}", fontsize = fs)
    plt.tight_layout()
    if len(geneIDs_lists_dict)>1:
        plt.legend(fontsize = fs)
    plt.savefig(outfile_name, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {outfile_name}")


def samples_group():
    
    samples = {
        "SL1_14_M" : ["WJ-3841-1-14-10-M_S377","WJ-3841-1-14-14-M_S399","WJ-3841-1-14-15-M_S400","WJ-3841-1-14-19-M_S425","WJ-3841-1-14-20-M_S426","WJ-3841-1-14-23-M_S427","WJ-3841-1-14-24-M_S428","WJ-3841-1-14-7-M_S374_","WJ-3841-1-14-8-M_S375_","WJ-3841-1-14-9-M_S376_"],
        "SL1_16_M" : ["WJ-3841-1-16-10-M_S401","WJ-3841-1-16-12-M_S402","WJ-3841-1-16-14-M_S403","WJ-3841-1-16-18-M_S404","WJ-3841-1-16-20-M_S431","WJ-3841-1-16-24-M_S433","WJ-3841-1-16-25-M_S434","WJ-3841-1-16-6-M_S379_","WJ-3841-1-16-7-M_S380_","WJ-3841-1-16-8-M_S381_","WJ-3841-1-16-9-M_S382_"],
        "SL1_18_M" : ["WJ-3841-1-18-10-M_S383","WJ-3841-1-18-11-M_S406","WJ-3841-1-18-14-M_S407","WJ-3841-1-18-17-M_S408","WJ-3841-1-18-19-M_S437","WJ-3841-1-18-21-M_S438","WJ-3841-1-18-23-M_S439"],
        "SL3_14_M" : ["WJ-3841-3-14-10-M_S361","WJ-3841-3-14-11-M_S385","WJ-3841-3-14-13-M_S386","WJ-3841-3-14-14-M_S387","WJ-3841-3-14-17-M_S411","WJ-3841-3-14-22-M_S412","WJ-3841-3-14-23-M_S413","WJ-3841-3-14-5-M_S359_","WJ-3841-3-14-9-M_S360_"],
        "SL3_16_M" : ["WJ-3841-3-16-11-M_S365","WJ-3841-3-16-13-M_S390","WJ-3841-3-16-14-M_S391","WJ-3841-3-16-16-M_S392","WJ-3841-3-16-21-M_S416","WJ-3841-3-16-26-M_S417","WJ-3841-3-16-7-M_S364_"],
        "SL3_18_M" : ["WJ-3841-3-18-11-M_S395","WJ-3841-3-18-15-M_S396","WJ-3841-3-18-16-M_S397","WJ-3841-3-18-17-M_S398","WJ-3841-3-18-20-M_S421","WJ-3841-3-18-21-M_S422","WJ-3841-3-18-22-M_S423","WJ-3841-3-18-4-M_S367_","WJ-3841-3-18-6-M_S368_","WJ-3841-3-18-7-M_S369_","WJ-3841-3-18-8-M_S370_","WJ-3841-3-18-9-M_S371_"],
        "SL1_14_F" : ["WJ-3841-1-14-13-F_S429","WJ-3841-1-14-16-F_S430","WJ-3841-1-14-4-F_S378_"],
        "SL1_16_F" : ["WJ-3841-1-16-17-F_S405","WJ-3841-1-16-19-F_S435","WJ-3841-1-16-22-F_S436"],
        "SL1_18_F" : ["WJ-3841-1-18-12-F_S440","WJ-3841-1-18-13-F_S420","WJ-3841-1-18-5-F_S384_","WJ-3841-1-18-8-F_S409_","WJ-3841-1-18-9-F_S410_"],
        "SL3_14_F" : ["WJ-3841-3-14-12-F_S389","WJ-3841-3-14-16-F_S414","WJ-3841-3-14-18-F_S415","WJ-3841-3-14-4-F_S362_","WJ-3841-3-14-6-F_S363_","WJ-3841-3-14-7-F_S388_"],
        "SL3_16_F" : ["WJ-3841-3-16-10-F_S393","WJ-3841-3-16-12-F_S394","WJ-3841-3-16-15-F_S418","WJ-3841-3-16-17-F_S419","WJ-3841-3-16-18-F_S432","WJ-3841-3-16-9-F_S366_"],
        "SL3_18_F" : ["WJ-3841-3-18-10-F_S373","WJ-3841-3-18-24-F_S424","WJ-3841-3-18-5-F_S372_"],
    }
    return samples

def get_counts_paths(username = "milena"):
    counts_paths = {
        "no_log" : f"/Users/{username}/work/PhD_code/PhD_chapter4/data/gene_counts_normalized_nolog.tsv",
        "yes_log": f"/Users/{username}/work/PhD_code/PhD_chapter4/data/gene_counts_normalized.tsv",
        "raw" : f"/Users/{username}/work/PhD_code/PhD_chapter4/data/gene_counts_standard.txt",
    }
    return counts_paths

def get_y_information():
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
        "expressed" : ["gene-371805","yTor-A","gene-371844","yTor-C","gene-371913","gene-371922","gene-371957","gene-372053","gene-372068","gene-372216","gene-372264"]
        }
    x_contigs = {
        'scaffold_10' : ["gene-306335","gene-306353","gene-306359","gene-306484","gene-306496","gene-306523","gene-306550","gene-306577","gene-306589","gene-306648","gene-306675","gene-306713","gene-306736","gene-306977","gene-306995","gene-307001","gene-307013","gene-307098","gene-307178","gene-307190","gene-307264","gene-307421","gene-307427","gene-307433","gene-307442","gene-307448","gene-307460","gene-307469","gene-307478","gene-307484","gene-307496","gene-307508","gene-307537","gene-307555","gene-307679","gene-307691","gene-307700","gene-307736","gene-307746","gene-308117","gene-308126","gene-308132","gene-308150","gene-308177","gene-308198","gene-308216","gene-308231","gene-308282","gene-308300","gene-308321","gene-308424","gene-308454","gene-308487","gene-308508","gene-308549","gene-308564","gene-308597","gene-308603","gene-308697","gene-308708","gene-308738","gene-308795","gene-308810","gene-308822","gene-308840","gene-308850","gene-308862","gene-308874","gene-308880","gene-308887","gene-308897","gene-308921","gene-308930","gene-308945","gene-308981","gene-308990","gene-309049","gene-309070","gene-309076","gene-309252","gene-309273","gene-309288","gene-309315","gene-309327","gene-309339","gene-309387","gene-309399","gene-309408","gene-309414","gene-309420","gene-309471","gene-309516","gene-309522","gene-309556","gene-309583","gene-309642","gene-309648","gene-309660","gene-309708","gene-309714","gene-309744","gene-309768","gene-309828","gene-309843","gene-309870","gene-309876","gene-309885","gene-309915","gene-309966","gene-309993","gene-309999","gene-310012","gene-310065","gene-310115","gene-310121","gene-310152","gene-310158","gene-310176","gene-310194","gene-310230","gene-310251","gene-310287","gene-310338","gene-310383","gene-310422","gene-310482","gene-310604","gene-310639","gene-310648","gene-310660","gene-310714","gene-310729","gene-310850","gene-310868","gene-310889","gene-310951","gene-310975","gene-311116","gene-311131","gene-311152","gene-311193","gene-311211","gene-311273","gene-311303","gene-311343","gene-311367","gene-311396","gene-311426","gene-311444","gene-311450","gene-311462","gene-311477","gene-311521","gene-311551","gene-311563","gene-311581","gene-311593","gene-311620","gene-311647","gene-311674","gene-311725","gene-311734","gene-311749","gene-311758","gene-311767","gene-311796","gene-311808","gene-311817","gene-311904","gene-311910","gene-311934","gene-311952","gene-311970","gene-312110","gene-312140","gene-312229","gene-312264","gene-312270","gene-312288","gene-312294","gene-312314","gene-312326","gene-312362","gene-312409","gene-312416","gene-312440","gene-312452","gene-312479","gene-312544","gene-312583","gene-312604","gene-312715","gene-312748","gene-312890","gene-312908","gene-312923","gene-312938","gene-312968","gene-312983","gene-313040","gene-313064","gene-313139","gene-313154","gene-313196","gene-313303","gene-313324","gene-313330","gene-313455","gene-313478","gene-313493","gene-313589","gene-313618","gene-314377","gene-314386","gene-314478","gene-314490","gene-314511","gene-314576","gene-314588","gene-314599","gene-314611","gene-314715","gene-314807","gene-314825","gene-314834","gene-314864","gene-314873","gene-314900","gene-314921","gene-314945","gene-314951","gene-314980","gene-315004","gene-315034","gene-315114","gene-315138","gene-315170","gene-315215","gene-315254","gene-315287","gene-315379","gene-315455","gene-315467","gene-315482","gene-315515","gene-315539","gene-315551","gene-315569","gene-315578","gene-315593","gene-315603","gene-315630","gene-315640","gene-315808","gene-315814","gene-315826","gene-315967","gene-315991","gene-316014","gene-316126","gene-316255","gene-316282","gene-316300","gene-316413","gene-316434","gene-316459","gene-316474","gene-316544","gene-316589","gene-316686","gene-316713","gene-316731","gene-316755","gene-316779","gene-316827","gene-316836","gene-316880","gene-316886","gene-316951","gene-317002","gene-317103","gene-317127","gene-317184","gene-317205","gene-317217","gene-317235","gene-317256","gene-317283","gene-317289","gene-317336","gene-317372","gene-317405","gene-317564","gene-317573","gene-317951","gene-317969","gene-318475","gene-320235","gene-320352","gene-320375","gene-320404","gene-320470","gene-320479","gene-320494","gene-320541","gene-320568","gene-320594","gene-320647","gene-320770","gene-320782","gene-320803","gene-320812","gene-320821","gene-320833","gene-320854","gene-320878","gene-320978","gene-320990","gene-321011","gene-321017","gene-321078","gene-321107","gene-321177","gene-321219","gene-321346","gene-321367","gene-321450","gene-321510","gene-321546","gene-321587","gene-321629","gene-321638","gene-321653","gene-321691","gene-321706","gene-321733","gene-321745","gene-321775","gene-321790","gene-321814","gene-321841","gene-321850","gene-321859","gene-321868","gene-321913","gene-321925","gene-321991","gene-322044","gene-322092","gene-322136","gene-322154","gene-322160"],
        'scaffold_14' : ["gene-342132","gene-342138","gene-342144","gene-342150","gene-342165","gene-342177","gene-342183","gene-342189","gene-342198","gene-342204","gene-342210","gene-342216","gene-342222","gene-342228","gene-342240","gene-342246","gene-342252","gene-342261","gene-342270","gene-342279","gene-342285","gene-342291","gene-342297","gene-342303","gene-342309","gene-342315","gene-342321","gene-342330","gene-342339","gene-342348","gene-342354","gene-342360","gene-342366","gene-342372","gene-342378","gene-342384","gene-342390","gene-342399","gene-342408","gene-342417","gene-342423","gene-342429","gene-342438","gene-342447","gene-342456","gene-342465","gene-342471","gene-342477","gene-342483","gene-342489","gene-342495","gene-342501","gene-342510","gene-342519","gene-342525","gene-342531","gene-342540","gene-342546","gene-342552","gene-342561","gene-342567","gene-342573","gene-342579","gene-342594","gene-342639","gene-342660","gene-342667","gene-342682","gene-342688","gene-342694","gene-342718","gene-342950","gene-342962","gene-342971","gene-342977","gene-342983","gene-343084","gene-343093","gene-343099","gene-343117","gene-343123","gene-343165","gene-343194","gene-343203","gene-343215","gene-343242","gene-343369","gene-343376","gene-343453","gene-343500","gene-343547","gene-343580","gene-343648","gene-343660","gene-343684","gene-343696","gene-343708","gene-343717","gene-343726","gene-343747","gene-343838","gene-343924","gene-343933","gene-343951","gene-343971","gene-343989","gene-344019","gene-344052","gene-344096","gene-344111","gene-344117","gene-344147","gene-344156","gene-344189","gene-344201","gene-344210","gene-344537","gene-344549","gene-344555","gene-344564","gene-344945","gene-344978","gene-345019","gene-345031","gene-345061","gene-345067","gene-345073","gene-345094","gene-345100","gene-345115","gene-345129","gene-345135","gene-345162","gene-345180","gene-345195","gene-345222","gene-345246","gene-345258","gene-345273","gene-345288","gene-345306","gene-345372","gene-345392","gene-345419","gene-345515","gene-345524","gene-345545","gene-345557"],
        'scaffold_23' : ["gene-367203","gene-367224","gene-367230","gene-367236","gene-367251","gene-367267","gene-367277","gene-367351","gene-367360","gene-367366","gene-367381","gene-367393","gene-367399","gene-367408","gene-367414","gene-367420","gene-367427","gene-367433","gene-367439","gene-367448","gene-367454","gene-367460","gene-367469","gene-367478","gene-367484","gene-367490","gene-367502","gene-367517"],
        'scaffold_31' : ["gene-381161","gene-381182","gene-381188","gene-381194","gene-381200","gene-381209","gene-381218","gene-381237","gene-381246","gene-381259","gene-381280","gene-381292","gene-381304","gene-381310","gene-381322","gene-381331","gene-381337","gene-381349","gene-381355","gene-381361","gene-381370","gene-381376","gene-381394","gene-381436","gene-381442","gene-381451","gene-381491","gene-381512","gene-381547","gene-381574","gene-381597","gene-381624","gene-381651","gene-381698","gene-381757","gene-381766","gene-381772","gene-381799","gene-381805","gene-381988","gene-381997","gene-382021","gene-382055","gene-382088","gene-382095","gene-382158","gene-382176","gene-382294","gene-382300","gene-382309","gene-382321","gene-382327","gene-382333","gene-382345","gene-382351","gene-382357","gene-382363","gene-382369","gene-382375","gene-382384","gene-382396","gene-382402","gene-382408","gene-382414","gene-382432","gene-382438","gene-382444","gene-382455","gene-382464","gene-382473","gene-382479","gene-382485","gene-382491","gene-382500","gene-382506","gene-382512","gene-382518","gene-382524","gene-382530","gene-382536","gene-382547","gene-382559","gene-382571","gene-382577","gene-382583","gene-382589","gene-382601","gene-382607","gene-382613","gene-382619","gene-382625"],
        'scaffold_34' : ["gene-386046","gene-386052","gene-386058","gene-386118","gene-386151","gene-386178","gene-386234","gene-386270","gene-386276","gene-386303","gene-386309","gene-386321","gene-386336","gene-386375","gene-386470","gene-386717","gene-386821","gene-386830","gene-386836","gene-386842","gene-386848","gene-386933","gene-386954","gene-386960","gene-386969","gene-386975","gene-386981","gene-386987","gene-386993","gene-387002","gene-387081","gene-387251","gene-387272","gene-387443","gene-387497","gene-387749","gene-387992","gene-388013","gene-388076","gene-388126","gene-388138","gene-388159","gene-388171","gene-388183","gene-388204","gene-388216","gene-388234","gene-388249","gene-388261","gene-388296","gene-388322","gene-388340","gene-388416","gene-388422","gene-388434","gene-388479","gene-388488","gene-388556","gene-388577","gene-388598"],
        'scaffold_83' : ["gene-419623","gene-419632"],
        'all' : ["gene-306335","gene-306353","gene-306359","gene-306484","gene-306496","gene-306523","gene-306550","gene-306577","gene-306589","gene-306648","gene-306675","gene-306713","gene-306736","gene-306977","gene-306995","gene-307001","gene-307013","gene-307098","gene-307178","gene-307190","gene-307264","gene-307421","gene-307427","gene-307433","gene-307442","gene-307448","gene-307460","gene-307469","gene-307478","gene-307484","gene-307496","gene-307508","gene-307537","gene-307555","gene-307679","gene-307691","gene-307700","gene-307736","gene-307746","gene-308117","gene-308126","gene-308132","gene-308150","gene-308177","gene-308198","gene-308216","gene-308231","gene-308282","gene-308300","gene-308321","gene-308424","gene-308454","gene-308487","gene-308508","gene-308549","gene-308564","gene-308597","gene-308603","gene-308697","gene-308708","gene-308738","gene-308795","gene-308810","gene-308822","gene-308840","gene-308850","gene-308862","gene-308874","gene-308880","gene-308887","gene-308897","gene-308921","gene-308930","gene-308945","gene-308981","gene-308990","gene-309049","gene-309070","gene-309076","gene-309252","gene-309273","gene-309288","gene-309315","gene-309327","gene-309339","gene-309387","gene-309399","gene-309408","gene-309414","gene-309420","gene-309471","gene-309516","gene-309522","gene-309556","gene-309583","gene-309642","gene-309648","gene-309660","gene-309708","gene-309714","gene-309744","gene-309768","gene-309828","gene-309843","gene-309870","gene-309876","gene-309885","gene-309915","gene-309966","gene-309993","gene-309999","gene-310012","gene-310065","gene-310115","gene-310121","gene-310152","gene-310158","gene-310176","gene-310194","gene-310230","gene-310251","gene-310287","gene-310338","gene-310383","gene-310422","gene-310482","gene-310604","gene-310639","gene-310648","gene-310660","gene-310714","gene-310729","gene-310850","gene-310868","gene-310889","gene-310951","gene-310975","gene-311116","gene-311131","gene-311152","gene-311193","gene-311211","gene-311273","gene-311303","gene-311343","gene-311367","gene-311396","gene-311426","gene-311444","gene-311450","gene-311462","gene-311477","gene-311521","gene-311551","gene-311563","gene-311581","gene-311593","gene-311620","gene-311647","gene-311674","gene-311725","gene-311734","gene-311749","gene-311758","gene-311767","gene-311796","gene-311808","gene-311817","gene-311904","gene-311910","gene-311934","gene-311952","gene-311970","gene-312110","gene-312140","gene-312229","gene-312264","gene-312270","gene-312288","gene-312294","gene-312314","gene-312326","gene-312362","gene-312409","gene-312416","gene-312440","gene-312452","gene-312479","gene-312544","gene-312583","gene-312604","gene-312715","gene-312748","gene-312890","gene-312908","gene-312923","gene-312938","gene-312968","gene-312983","gene-313040","gene-313064","gene-313139","gene-313154","gene-313196","gene-313303","gene-313324","gene-313330","gene-313455","gene-313478","gene-313493","gene-313589","gene-313618","gene-314377","gene-314386","gene-314478","gene-314490","gene-314511","gene-314576","gene-314588","gene-314599","gene-314611","gene-314715","gene-314807","gene-314825","gene-314834","gene-314864","gene-314873","gene-314900","gene-314921","gene-314945","gene-314951","gene-314980","gene-315004","gene-315034","gene-315114","gene-315138","gene-315170","gene-315215","gene-315254","gene-315287","gene-315379","gene-315455","gene-315467","gene-315482","gene-315515","gene-315539","gene-315551","gene-315569","gene-315578","gene-315593","gene-315603","gene-315630","gene-315640","gene-315808","gene-315814","gene-315826","gene-315967","gene-315991","gene-316014","gene-316126","gene-316255","gene-316282","gene-316300","gene-316413","gene-316434","gene-316459","gene-316474","gene-316544","gene-316589","gene-316686","gene-316713","gene-316731","gene-316755","gene-316779","gene-316827","gene-316836","gene-316880","gene-316886","gene-316951","gene-317002","gene-317103","gene-317127","gene-317184","gene-317205","gene-317217","gene-317235","gene-317256","gene-317283","gene-317289","gene-317336","gene-317372","gene-317405","gene-317564","gene-317573","gene-317951","gene-317969","gene-318475","gene-320235","gene-320352","gene-320375","gene-320404","gene-320470","gene-320479","gene-320494","gene-320541","gene-320568","gene-320594","gene-320647","gene-320770","gene-320782","gene-320803","gene-320812","gene-320821","gene-320833","gene-320854","gene-320878","gene-320978","gene-320990","gene-321011","gene-321017","gene-321078","gene-321107","gene-321177","gene-321219","gene-321346","gene-321367","gene-321450","gene-321510","gene-321546","gene-321587","gene-321629","gene-321638","gene-321653","gene-321691","gene-321706","gene-321733","gene-321745","gene-321775","gene-321790","gene-321814","gene-321841","gene-321850","gene-321859","gene-321868","gene-321913","gene-321925","gene-321991","gene-322044","gene-322092","gene-322136","gene-322154","gene-322160","gene-342132","gene-342138","gene-342144","gene-342150","gene-342165","gene-342177","gene-342183","gene-342189","gene-342198","gene-342204","gene-342210","gene-342216","gene-342222","gene-342228","gene-342240","gene-342246","gene-342252","gene-342261","gene-342270","gene-342279","gene-342285","gene-342291","gene-342297","gene-342303","gene-342309","gene-342315","gene-342321","gene-342330","gene-342339","gene-342348","gene-342354","gene-342360","gene-342366","gene-342372","gene-342378","gene-342384","gene-342390","gene-342399","gene-342408","gene-342417","gene-342423","gene-342429","gene-342438","gene-342447","gene-342456","gene-342465","gene-342471","gene-342477","gene-342483","gene-342489","gene-342495","gene-342501","gene-342510","gene-342519","gene-342525","gene-342531","gene-342540","gene-342546","gene-342552","gene-342561","gene-342567","gene-342573","gene-342579","gene-342594","gene-342639","gene-342660","gene-342667","gene-342682","gene-342688","gene-342694","gene-342718","gene-342950","gene-342962","gene-342971","gene-342977","gene-342983","gene-343084","gene-343093","gene-343099","gene-343117","gene-343123","gene-343165","gene-343194","gene-343203","gene-343215","gene-343242","gene-343369","gene-343376","gene-343453","gene-343500","gene-343547","gene-343580","gene-343648","gene-343660","gene-343684","gene-343696","gene-343708","gene-343717","gene-343726","gene-343747","gene-343838","gene-343924","gene-343933","gene-343951","gene-343971","gene-343989","gene-344019","gene-344052","gene-344096","gene-344111","gene-344117","gene-344147","gene-344156","gene-344189","gene-344201","gene-344210","gene-344537","gene-344549","gene-344555","gene-344564","gene-344945","gene-344978","gene-345019","gene-345031","gene-345061","gene-345067","gene-345073","gene-345094","gene-345100","gene-345115","gene-345129","gene-345135","gene-345162","gene-345180","gene-345195","gene-345222","gene-345246","gene-345258","gene-345273","gene-345288","gene-345306","gene-345372","gene-345392","gene-345419","gene-345515","gene-345524","gene-345545","gene-345557","gene-367203","gene-367224","gene-367230","gene-367236","gene-367251","gene-367267","gene-367277","gene-367351","gene-367360","gene-367366","gene-367381","gene-367393","gene-367399","gene-367408","gene-367414","gene-367420","gene-367427","gene-367433","gene-367439","gene-367448","gene-367454","gene-367460","gene-367469","gene-367478","gene-367484","gene-367490","gene-367502","gene-367517","gene-381161","gene-381182","gene-381188","gene-381194","gene-381200","gene-381209","gene-381218","gene-381237","gene-381246","gene-381259","gene-381280","gene-381292","gene-381304","gene-381310","gene-381322","gene-381331","gene-381337","gene-381349","gene-381355","gene-381361","gene-381370","gene-381376","gene-381394","gene-381436","gene-381442","gene-381451","gene-381491","gene-381512","gene-381547","gene-381574","gene-381597","gene-381624","gene-381651","gene-381698","gene-381757","gene-381766","gene-381772","gene-381799","gene-381805","gene-381988","gene-381997","gene-382021","gene-382055","gene-382088","gene-382095","gene-382158","gene-382176","gene-382294","gene-382300","gene-382309","gene-382321","gene-382327","gene-382333","gene-382345","gene-382351","gene-382357","gene-382363","gene-382369","gene-382375","gene-382384","gene-382396","gene-382402","gene-382408","gene-382414","gene-382432","gene-382438","gene-382444","gene-382455","gene-382464","gene-382473","gene-382479","gene-382485","gene-382491","gene-382500","gene-382506","gene-382512","gene-382518","gene-382524","gene-382530","gene-382536","gene-382547","gene-382559","gene-382571","gene-382577","gene-382583","gene-382589","gene-382601","gene-382607","gene-382613","gene-382619","gene-382625","gene-386046","gene-386052","gene-386058","gene-386118","gene-386151","gene-386178","gene-386234","gene-386270","gene-386276","gene-386303","gene-386309","gene-386321","gene-386336","gene-386375","gene-386470","gene-386717","gene-386821","gene-386830","gene-386836","gene-386842","gene-386848","gene-386933","gene-386954","gene-386960","gene-386969","gene-386975","gene-386981","gene-386987","gene-386993","gene-387002","gene-387081","gene-387251","gene-387272","gene-387443","gene-387497","gene-387749","gene-387992","gene-388013","gene-388076","gene-388126","gene-388138","gene-388159","gene-388171","gene-388183","gene-388204","gene-388216","gene-388234","gene-388249","gene-388261","gene-388296","gene-388322","gene-388340","gene-388416","gene-388422","gene-388434","gene-388479","gene-388488","gene-388556","gene-388577","gene-388598","gene-419623","gene-419632"]
    }
    tor_related = {
        "FOXO" : ["gene-378716"],
        "SREBP" : ["gene-420407","gene-294541","gene-349478","gene-429040"],
        "Atg1" : ["gene-165615"],
        "Y-Tor" : ["yTor-all"],
        "A-Tor" : ["gene-30110"],

    }
    return sex_chromosomes_superscaffolded,y_contigs,x_contigs,tor_related


def filter_counts_file(counts_path, out_path, IDs_list):
    """
    filter a raw gene counts file to only include rows that are genes in IDs_list
    """
    line_counts = 0
    line_counts_all = 0
    with open(counts_path, "r") as counts_file, open(out_path, "w") as filtered_counts:
        for counts_line in counts_file.readlines():
            line_counts_all += 1
            geneID = counts_line.split()[0]
            if geneID in IDs_list:
                filtered_counts.write(counts_line)
                line_counts += 1
            elif geneID == "Geneid":
                filtered_counts.write(counts_line)
                line_counts += 1 
            elif counts_line[0] == "#":
                filtered_counts.write(counts_line)
                line_counts += 1
        print(f"filtering for {len(IDs_list)} genes : {line_counts} out of {line_counts_all} lines written to new file {out_path}")
            
                


if __name__ == "__main__":
    
    username = "miltr339"
    count_files = get_counts_paths(username=username)
    ex_chromosomes,y_contigs,x_contigs,tor_related = get_y_information()
    out_path = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/yTor_analysis"
    samples_group_dict = samples_group()
    samples_group_dict_rev = {item: key for key, values in samples_group_dict.items() for item in values}


    if True:
        print(f"\n")
        plot_counts_sum_sets(counts_table=count_files["no_log"], geneIDs_lists_dict = {"Y" : y_contigs["all"]}, 
                    outfile_name = f"{out_path}/y_genes_mean_expression.png", y_label= "normalized counts", errorbars=True)
        print(f"\n")
        plot_counts_sum_sets(counts_table=count_files["no_log"], geneIDs_lists_dict = {"Y" : y_contigs["all"], "X" : x_contigs["all"]}, 
                    outfile_name = f"{out_path}/y_x_genes_mean_expression.png", y_label= "normalized counts", errorbars=False)
    if True:
        print(f"\n")
        plot_counts_sum_sets(counts_table=count_files["no_log"], geneIDs_lists_dict = {"Y" : y_contigs["all"]}, 
                    outfile_name = f"{out_path}/y_genes_mean_expression_sample_groups.png", y_label= "normalized counts", errorbars=True, samples_group_dict = samples_group_dict)
        print(f"\n")
        plot_counts_sum_sets(counts_table=count_files["no_log"], geneIDs_lists_dict = {"Y" : y_contigs["all"], "X" : x_contigs["all"]}, 
                    outfile_name = f"{out_path}/y_x_genes_mean_expression_sample_groups.png", y_label= "normalized counts", errorbars=False, samples_group_dict = samples_group_dict)
    if True:
        print(f"\n")
        plot_counts_sum_sets(counts_table=count_files["no_log"], geneIDs_lists_dict = tor_related, 
                    outfile_name = f"{out_path}/tor_pathway_related_expression.png", y_label= "normalized counts", errorbars=True)
        print(f"\n")
        plot_counts_sum_sets(counts_table=count_files["no_log"], geneIDs_lists_dict = tor_related, 
                    outfile_name = f"{out_path}/tor_pathway_related_expression_sample_groups.png", y_label= "normalized counts", errorbars=True, samples_group_dict = samples_group_dict)
    
    ### test stuff with edgeR on downsampled counts files, did not work! only here for posterity just in case
    # filter_counts_file(counts_path=count_files["raw"], out_path = count_files["raw"].replace(".txt", "_only_Y.txt"), IDs_list=y_contigs["all"])
    # filter_counts_file(counts_path=count_files["raw"], out_path = count_files["raw"].replace(".txt", "_only_X.txt"), IDs_list=x_contigs["all"])