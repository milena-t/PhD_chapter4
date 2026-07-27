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

def plot_PCA_vst_counts(counts_path:str, metadata_path:str, plot_path:str="", colors_dict=colors_dict, points_dict=points_dict, condition="line"):
    """
    make a PCA of the DE counts after they are vst-normalized
    """
    if plot_path=="":
        plot_path=counts_path.replace(".txt", "_PCA.png")
    
    norm_counts = pd.read_csv(counts_path, sep="\t", comment="#", index_col=0)
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
        if False:
            for condition in ["line", "day"]:
                plot_PCA_vst_counts(
                    counts_path=vst_path, 
                    metadata_path=metadata_path, 
                    plot_path=f"/Users/{username}/work/PhD_code/PhD_chapter4/data/DE_figures/PCA_sex_{condition}_all_counts.png",
                    condition=condition)
        if True:
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