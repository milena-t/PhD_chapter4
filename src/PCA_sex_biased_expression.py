""" 
plot for the differential expression analysis in combination with DE_analysis_edgeR.Rmd
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


def counts_paths(username="miltr339"):
    counts_path = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/gene_counts_standard.txt"
    VST_path = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/gene_counts_normalized.tsv"
    return counts_path,VST_path

def metadata_paths(username="miltr339"):
    metadata = f"/Users/{username}/work/PhD_code/PhD_chapter4/data/larvae_metadata.csv"
    return metadata


colors_dict = {
    "sex" : {
        "M" : "#4570B0",
        "F" : "#C32B09",
    },
    "line" : {
        1 : "#8D94BA", # lavender grey
        3 : "#54457F" , # dusty grape
    },
    "day" : {
        14 : "#003B36", # evergreen
        16 : "#7EA172", # sage green
        18 : "#DE871C", # amber earth
    }
}

points_dict = {
    "sex" : {
        "M" : "v", # down triangle
        "F" : "o", # default circle
    },
    "line" : {
        1 : "X", # bold X
        3 : "D" , # diamond shape
    },
    "day" : {
        14 : "1", # tri_up
        16 : "2" , # tri_down
        18 : "3" , # tri_left
    }
}

def plot_PCA_vst_counts(counts_path:str, metadata_path:str, plot_path:str="", colors_dict=colors_dict, points_dict=points_dict, condition="line", excl_list = []):
    """
    make a PCA of the DE counts after they are vst-normalized
    """
    if plot_path=="":
        plot_path=counts_path.replace(".txt", "_PCA.png")
    
    norm_counts = pd.read_csv(counts_path, sep="\t", comment="#", index_col=0)
    if len(excl_list)>0:
        print(f"{norm_counts.shape}")
        norm_counts = norm_counts.drop(index=excl_list)
        print(f"{norm_counts.shape}")
        excl_str = "excl_female_line-bias"
    else:
        excl_str = ""
    # read metadata and sort/order according to counts
    metadata = pd.read_csv(metadata_path)
    metadata = metadata.set_index('ID')
    metadata = metadata.reindex(norm_counts.columns)
    print(metadata)
    # categories = list(metadata.columns)

    # create transpose
    vst_t = norm_counts.T
    vst_t.columns = vst_t.columns.astype(str) # convert all column names to string
    # pca = PCA(n_components=2)
    pca = PCA()
    pca_scores = pca.fit_transform(vst_t)
    # print(pca_scores)

    pca_df = metadata.copy()
    pca_df['PC1'] = pca_scores[:,0]
    pca_df['PC2'] = pca_scores[:,1]
    # print(pca_df)

    # Explained variance
    pc1_var = pca.explained_variance_ratio_[0]
    pc2_var = pca.explained_variance_ratio_[1]

    ### plotting
    fig, ax = plt.subplots(1,1, figsize=(15, 10)) # for more than three rows

    fs = 25
    ps = fs*15 # point size

    for i, row in pca_df.iterrows():
        ax.scatter(row.loc["PC1"], row.loc["PC2"], marker=points_dict["sex"][row.loc['sex']], s=ps, color=colors_dict[condition][row.loc[condition]])
    # ax.scatter(pca_df["PC1"], pca_df["PC2"], color=color_by_sex_vec, s=100, marker=marker_by_organ_vec)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)", fontsize = fs)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)", fontsize = fs)
    ax.tick_params(axis='x', labelsize=fs)
    ax.tick_params(axis='y', labelsize=fs)
    ylim = ax.get_ylim()
    ax.set_ylim(ylim)
    xlim = ax.get_xlim()
    ax.set_xlim(xlim)

    ## make legend
    sex_labels = { sex : points_dict['sex'][sex] for sex in pca_df['sex']}
    org_labels = { org : colors_dict[condition][org] for org in pca_df[condition]}
    yleg = ylim[0]-1e6
    xleg = xlim[0]-1e6

    for key,value in sex_labels.items():
        ax.scatter(xleg,yleg,marker=value,s=ps,label=key, color='black')
    for key,value in org_labels.items():
        if key>10:
            label_ = f"Day {key}"
        elif key ==1:
            label_ = "small-Y"
        elif key==3:
            label_ = "large-Y"
        ax.scatter(xleg,yleg,color=value,s=ps,label=label_, marker="s")
    ax.legend(fontsize=fs)
    if excl_str != "":
        excl_str_ = excl_str.replace("_", " ")
        plt.suptitle(f"Differential expression \nbased on {len(metadata)} samples\n{excl_str_}", fontsize=fs*1.25)
        plot_path = plot_path.replace(".png", f"_{excl_str}.png")
    else:
        plt.suptitle(f"Differential expression \nbased on {len(metadata)} samples", fontsize=fs*1.25)
    
    # Adjust layout to prevent overlap
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    plt.savefig(plot_path, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {plot_path}")



def plot_PCA_separation(counts_path:str, metadata_path:str, plot_path:str="", colors_dict=colors_dict, points_dict=points_dict, sex="", day = "", line = ""):
    """
    make a PCA of the DE counts for only one sex by day and line
    """
    if plot_path=="":
        plot_path=counts_path.replace(".txt", "_PCA.png")
    
    if sex !="":
        cat_name = "male"
        if sex == "F":
            cat_name = "female"
        elif sex != "M":
            raise RuntimeError(f"invalid sex specifier, you picked '{sex}', pick 'M' or 'F'!")
    elif day != "":
        if day not in ["14","16","18"]:
            raise RuntimeError(f"invalid day specifier, you picked '{day}', pick '14', '16' or '18'!")
        cat_name = f"day {day}"
    elif line != "":
        if line not in ["1","3"]:
            raise RuntimeError(f"invalid day specifier, you picked '{line}', pick '1' or '3'!")
        if line == "1":
            cat_name = f"small-Y"
        if line == "3":
            cat_name = f"large-Y"
    
    norm_counts = pd.read_csv(counts_path, sep="\t", comment="#", index_col=0)
    # read metadata and sort/order according to counts
    metadata = pd.read_csv(metadata_path)
    metadata = metadata.set_index('ID')
    metadata = metadata.reindex(norm_counts.columns)
    ## filter only {sex} samples

    if sex != "":
        categories = [sample for sample in list(norm_counts.columns) if f"-{sex}_" in sample]
    elif day != "":
        line1 = [sample for sample in list(norm_counts.columns) if f"-1-{day}-" in sample] #otherwise you also get other days where the sample number happens to be 14/16/18
        line3 = [sample for sample in list(norm_counts.columns) if f"-3-{day}-" in sample]
        categories = line1+line3
    elif line != "":
        categories = [sample for sample in list(norm_counts.columns) if f"WJ-3841-{line}-" in sample]

    print(f"{len(list(norm_counts.columns))} categories before filtering, {len(categories)} after filtering")
    norm_counts = norm_counts[categories]
    metadata = metadata.filter(categories, axis="index")

    # create transpose
    vst_t = norm_counts.T
    vst_t.columns = vst_t.columns.astype(str) # convert all column names to string
    # pca = PCA(n_components=2)
    pca = PCA()
    pca_scores = pca.fit_transform(vst_t)
    # print(pca_scores)

    pca_df = metadata.copy()
    pca_df['PC1'] = pca_scores[:,0]
    pca_df['PC2'] = pca_scores[:,1]
    # print(pca_df)

    # Explained variance
    pc1_var = pca.explained_variance_ratio_[0]
    pc2_var = pca.explained_variance_ratio_[1]

    ### plotting
    fig, ax = plt.subplots(1,1, figsize=(15, 10)) # for more than three rows

    fs = 25
    ps = fs*15 # point size

    if sex != "":
        for i, row in pca_df.iterrows():
            ax.scatter(row.loc["PC1"], row.loc["PC2"], marker=points_dict["line"][row.loc['line']], s=ps, color=colors_dict["day"][row.loc["day"]])
        # ax.scatter(pca_df["PC1"], pca_df["PC2"], color=color_by_sex_vec, s=100, marker=marker_by_organ_vec)
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)", fontsize = fs)
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)", fontsize = fs)
        ax.tick_params(axis='x', labelsize=fs)
        ax.tick_params(axis='y', labelsize=fs)
        ylim = ax.get_ylim()
        ax.set_ylim(ylim)
        xlim = ax.get_xlim()
        ax.set_xlim(xlim)

        ## make legend
        sex_labels = { sex : points_dict['line'][sex] for sex in pca_df['line']}
        org_labels = { org : colors_dict["day"][org] for org in pca_df["day"]}
    
        yleg = ylim[0]-1e6
        xleg = xlim[0]-1e6

        for key,value in sex_labels.items():
            key = int(key)
            if key>10:
                label_ = f"Day {key}"
            elif key ==1:
                label_ = "small-Y"
            elif key==3:
                label_ = "large-Y"
            ax.scatter(xleg,yleg,marker=value,s=ps,label=label_, color='black')
        for key,value in org_labels.items():
            key = int(key)
            if key>10:
                label_ = f"Day {key}"
            elif key ==1:
                label_ = "small-Y"
            elif key==3:
                label_ = "large-Y"
            ax.scatter(xleg,yleg,color=value,s=ps,label=label_, marker="s")

    elif line != "":
        for i, row in pca_df.iterrows():
            ax.scatter(row.loc["PC1"], row.loc["PC2"], marker=points_dict["sex"][row.loc['sex']], s=ps, color=colors_dict["day"][row.loc["day"]])
        # ax.scatter(pca_df["PC1"], pca_df["PC2"], color=color_by_sex_vec, s=100, marker=marker_by_organ_vec)
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)", fontsize = fs)
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)", fontsize = fs)
        ax.tick_params(axis='x', labelsize=fs)
        ax.tick_params(axis='y', labelsize=fs)
        ylim = ax.get_ylim()
        ax.set_ylim(ylim)
        xlim = ax.get_xlim()
        ax.set_xlim(xlim)

        ## make legend
        sex_labels = { sex : points_dict['sex'][sex] for sex in pca_df['sex']}
        org_labels = { org : colors_dict["day"][org] for org in pca_df["day"]}
    
        yleg = ylim[0]-1e6
        xleg = xlim[0]-1e6

        for key,value in sex_labels.items():
            ax.scatter(xleg,yleg,marker=value,s=ps,label=key, color='black')
        for key,value in org_labels.items():
            if key>10:
                label_ = f"Day {key}"
            elif key ==1:
                label_ = "small-Y"
            elif key==3:
                label_ = "large-Y"
            # print(f"label: {key}, {type(key)} --> {label_}")
            ax.scatter(xleg,yleg,color=value,s=ps,label=label_, marker="s")

    elif day != "":
        for i, row in pca_df.iterrows():
            ax.scatter(row.loc["PC1"], row.loc["PC2"], marker=points_dict["sex"][row.loc['sex']], s=ps, color=colors_dict["line"][row.loc["line"]])
        # ax.scatter(pca_df["PC1"], pca_df["PC2"], color=color_by_sex_vec, s=100, marker=marker_by_organ_vec)
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)", fontsize = fs)
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)", fontsize = fs)
        ax.tick_params(axis='x', labelsize=fs)
        ax.tick_params(axis='y', labelsize=fs)
        ylim = ax.get_ylim()
        ax.set_ylim(ylim)
        xlim = ax.get_xlim()
        ax.set_xlim(xlim)

        ## make legend
        sex_labels = { sex : points_dict['sex'][sex] for sex in pca_df['sex']}
        org_labels = { org : colors_dict["line"][org] for org in pca_df["line"]}
    
        yleg = ylim[0]-1e6
        xleg = xlim[0]-1e6

        for key,value in sex_labels.items():
            try:
                key = int(key)
                if key>10:
                    label_ = f"Day {key}"
                elif key ==1:
                    label_ = "small-Y"
                elif key==3:
                    label_ = "large-Y"
            except:
                if key == "F" or key=="M":
                    label_ = key
            ax.scatter(xleg,yleg,marker=value,s=ps,label=label_, color='black')
        for key,value in org_labels.items():
            try:
                key = int(key)
                if key>10:
                    label_ = f"Day {key}"
                elif key ==1:
                    label_ = "small-Y"
                elif key==3:
                    label_ = "large-Y"
            except:
                if key == "F" or key=="M":
                    label_ = key
            ax.scatter(xleg,yleg,color=value,s=ps,label=label_, marker="s")
    ax.legend(fontsize=fs)
    plt.suptitle(f"Differential expression \nbased on {len(metadata)} {cat_name} samples", fontsize=fs*1.25)
    
    # Adjust layout to prevent overlap
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    plt.savefig(plot_path, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {plot_path}")



if __name__ == "__main__":

    username = "miltr339"


    ### plot PCA
    if True:
        counts_path, vst_path = counts_paths(username=username)
        metadata_path = metadata_paths(username=username)
        if True:
            # genes that are sig. line-biased in males and females which we think are batch effects, 
            # try to exclude them from the PCA to see if we can reduce the line-difference in females.
            excl_list = list(set(['gene-390678', 'gene-326825', 'gene-428756', 'gene-392159', 'gene-237318', 'gene-80484', 'gene-222486', 'gene-326909', 'gene-391198', 'gene-220544', 'gene-225236', 'gene-323148', 'gene-90157', 'gene-403706', 'gene-428738', 'gene-328764', 'gene-241055', 'gene-390956', 'gene-241108', 'gene-403583', 'gene-323803', 'gene-80466', 'gene-89234', 'gene-237881', 'gene-224357', 'gene-226245', 'gene-395080', 'gene-220249', 'gene-430044', 'gene-222365', 'gene-240929', 'gene-392186', 'gene-430032', 'gene-224782', 'gene-239553', 'gene-240623', 'gene-224743', 'gene-403700', 'gene-117712', 'gene-224860', 'gene-407253', 'gene-225635', 'gene-403902', 'gene-225709', 'gene-225030', 'gene-406468', 'gene-395143', 'gene-322927', 'gene-222600', 'gene-240871', 'gene-84949', 'gene-81551', 'gene-224845', 'gene-431701', 'gene-240833', 'gene-222383', 'gene-225738', 'gene-224227', 'gene-243308', 'gene-224277', 'gene-241126', 'gene-240910', 'gene-240935', 'gene-392248', 'gene-283443', 'gene-403809', 'gene-222344', 'gene-218086', 'gene-390637', 'gene-392224', 'gene-225629', 'gene-224697', 'gene-224682', 'gene-406796', 'gene-214979', 'gene-222531', 'gene-231925', 'gene-430068', 'gene-225173', 'gene-80359', 'gene-241001', 'gene-81640', 'gene-81599', 'gene-225140', 'gene-391222', 'gene-225720', 'gene-219019', 'gene-222332', 'gene-430263', 'gene-390616', 'gene-120952', 'gene-222555', 'gene-224968', 'gene-224614', 'gene-224250', 'gene-428747', 'gene-224307', 'gene-222501', 'gene-81427', 'gene-222430', 'gene-221953', 'gene-81572', 'gene-225107', 'gene-430314', 'gene-327441', 'gene-222350', 'gene-240983', 'gene-240691', 'gene-84970', 'gene-322912', 'gene-224201', 'gene-224956', 'gene-428765', 'gene-224593', 'gene-392290', 'gene-326873', 'gene-224896', 'gene-326849', 'gene-403818', 'gene-326810', 'gene-403851', 'gene-88715', 'gene-234575', 'gene-220028', 'gene-260693', 'gene-225325', 'gene-224875', 'gene-222519', 'gene-430080', 'gene-403652', 'gene-390678', 'gene-431030', 'gene-231854', 'gene-428756', 'gene-282886', 'gene-282746', 'gene-222486', 'gene-222159', 'gene-220544', 'gene-391198', 'gene-326909', 'gene-225236', 'gene-223491', 'gene-323148', 'gene-90157', 'gene-400402', 'gene-428738', 'gene-224890', 'gene-241055', 'gene-390956', 'gene-282347', 'gene-241108', 'gene-238407', 'gene-323803', 'gene-89234', 'gene-282491', 'gene-237881', 'gene-224357', 'gene-226245', 'gene-282853', 'gene-395080', 'gene-220249', 'gene-380466', 'gene-282398', 'gene-399475', 'gene-430044', 'gene-286545', 'gene-222365', 'gene-240929', 'gene-392186', 'gene-430032', 'gene-224782', 'gene-428113', 'gene-241682', 'gene-240623', 'gene-224743', 'gene-222519', 'gene-224860', 'gene-242512', 'gene-225635', 'gene-225709', 'gene-225030', 'g14784', 'gene-223773', 'gene-395143', 'gene-322927', 'gene-222600', 'gene-227308', 'gene-240871', 'gene-84949', 'gene-224845', 'gene-282590', 'gene-431701', 'gene-240833', 'gene-222383', 'gene-225738', 'gene-224227', 'gene-282551', 'gene-399317', 'gene-224277', 'gene-241126', 'gene-240638', 'gene-240910', 'gene-225158', 'gene-403809', 'gene-222344', 'gene-390637', 'gene-225629', 'gene-399424', 'gene-224697', 'gene-400393', 'gene-224682', 'gene-243299', 'gene-406796', 'gene-241262', 'gene-282620', 'gene-222531', 'gene-430068', 'gene-225173', 'gene-400384', 'gene-80359', 'gene-227370', 'gene-241001', 'gene-81640', 'gene-225140', 'gene-391222', 'gene-428104', 'gene-225720', 'gene-222332', 'gene-430263', 'gene-282458', 'gene-390616', 'gene-120952', 'gene-399484', 'gene-222555', 'gene-224968', 'gene-282362', 'gene-428747', 'gene-224307', 'gene-222501', 'gene-81427', 'gene-239506', 'gene-222430', 'gene-221953', 'gene-282524', 'gene-81572', 'gene-225107', 'gene-399270', 'gene-430314', 'gene-222350', 'gene-240983', 'gene-240691', 'gene-84970', 'gene-238849', 'gene-322912', 'gene-224201', 'gene-224956', 'gene-428765', 'gene-224593', 'gene-224896', 'gene-240860', 'gene-326849', 'gene-282784', 'gene-403818', 'gene-282665', 'gene-220028', 'gene-260693', 'gene-225325', 'gene-224875', 'gene-392290', 'gene-430080', 'gene-282701']))
            for condition in ["line", "day"]:
                plot_PCA_vst_counts(
                    counts_path=vst_path, 
                    metadata_path=metadata_path, 
                    plot_path=f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures/PCA_sex_{condition}_all_counts.png",
                    # excl_list = excl_list,
                    condition=condition)
        if False:
            for sex in ["M", "F"] :
                plot_PCA_separation(
                    counts_path=vst_path, 
                    metadata_path=metadata_path, 
                    plot_path=f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures/PCA_{sex}_day_line.png",
                    sex=sex)
            for day in ["14","16","18"]:
                plot_PCA_separation(
                    counts_path=vst_path, 
                    metadata_path=metadata_path, 
                    plot_path=f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures/PCA_sex_day{day}_line.png",
                    day=day)
            for line in ["1","3"]:
                plot_PCA_separation(
                    counts_path=vst_path, 
                    metadata_path=metadata_path, 
                    plot_path=f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures/PCA_sex_day_SL{line}.png",
                    line=line)