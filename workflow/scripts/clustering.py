import os
import numpy as np
import pandas as pd
import yaml
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial import distance
from sklearn.metrics import silhouette_score
#import rpy2.robjects as ro
#from rpy2.robjects import pandas2ri

import pickle
from utils import pickler

from plots import save_empty

def get_dendrogram(df,
                    linkage,
                    thresh = 0,
                    Plot = False, 
                    title = 'Hierarchical Clustering Dendrogram',
                    figsize = (50, 20)):

    if Plot:
        plt.figure(figsize=figsize)

    dn = hierarchy.dendrogram(
        linkage,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
        color_threshold=thresh,
        labels = df.columns,
        get_leaves = True,
        no_plot = not Plot
    )

    if Plot:
        plt.title(title)
        plt.xlabel('sample index')
        plt.ylabel('distance')
        plt.axhline(y=thresh, c='k', linestyle='--')
        #plt.savefig(f'{outpath}/{dataname}-dendrogram.png', facecolor='w', transparent=False)
        plt.show()

    return dn

def hierarchical_clustering(heatmap_data, distance_threshold):

    linkage_matrix = hierarchy.linkage(distance.squareform(1 - heatmap_data.values), method='average')

    # Clusters
    clusters = hierarchy.fcluster(linkage_matrix, t=distance_threshold, criterion='distance')
    cluster_df = pd.DataFrame(index = heatmap_data.index, data=clusters, columns=["GO_Cluster"])

    # Heatmap

    unique_clusters = np.unique(clusters)
    cluster_colors = sns.color_palette('tab10', len(unique_clusters))
    import matplotlib.colors as mcolors
    cluster_color_map = {cluster: mcolors.to_hex(cluster_colors[i]) for i, cluster in enumerate(unique_clusters)}
    row_colors_cl = pd.Series(clusters, index=heatmap_data.index).map(cluster_color_map)
    g = sns.clustermap(heatmap_data,row_linkage=linkage_matrix,col_linkage=linkage_matrix,robust=True,row_colors=row_colors_cl,col_colors=row_colors_cl) #, vmin=vmin2, vmax=vmax2, method=distmethod)
    g.ax_row_dendrogram.axvline(distance_threshold, c='magenta', linestyle='--',lw=1.5, alpha=0.5)

    return g, cluster_df


def find_optimal_threshold(heatmap_data, steps = 50):

    n_samples = len(heatmap_data)

    # Create linkage matrix
    linkage_matrix = hierarchy.linkage(distance.squareform(1 - heatmap_data.values), method='average')

    # Range of distance_threshold values to try
    thresholds = np.linspace(0.1, 1.0, steps)  # Adjust this range depending on your data

    # Store silhouette scores for each threshold
    silhouette_scores = []

    for threshold in thresholds:
        # Get clusters at the current distance threshold
        clusters = hierarchy.fcluster(linkage_matrix, t=threshold, criterion='distance')
        # Calculate silhouette score (requires more than 1 cluster)
        if n_samples > len(set(clusters)) > 1:
            score = silhouette_score(heatmap_data.values, clusters, metric='euclidean')
            silhouette_scores.append(score)
        else:
            silhouette_scores.append(-1)  # In case there's only 1 cluster

    # Find the threshold with the best silhouette score
    optimal_threshold = thresholds[np.argmax(silhouette_scores)]

    print(f"Optimal distance threshold: {optimal_threshold}")

    return optimal_threshold

# needs different environment!
# def get_clusters_from_depth_df(depth_df, ont, org, qval = 0.05, figpath="", max_thresh=0):

#     pandas2ri.activate()
#     ro.r['source']('workflow/scripts/gosemsim.R')
#     semsim_r = ro.r['semsim']
#     depth_df_r = pandas2ri.py2rpy(depth_df)
#     df_out_r = semsim_r(depth_df_r, ont, org, qval)
#     sim_matrix = pandas2ri.rpy2py(df_out_r)

#     thresh = find_optimal_threshold(sim_matrix)
#     thresh = max(max_thresh, thresh)
#     heat, cluster_df = hierarchical_clustering(sim_matrix, thresh)

#     if figpath != "":
#         heat.figure.suptitle(f"Threshold: {thresh:.2f}")
#         heat.figure.tight_layout()
#         heat.savefig(figpath)

#     return cluster_df

def get_clusters_from_sim_matrix(sim_matrix_path, figpath="", max_thresh=0):

    sim_matrix = pd.read_csv(sim_matrix_path, index_col=0)

    if len(sim_matrix) < 1:
        print("No terms found for GO Clusters, saving empty:", figpath)
        save_empty(figpath)
        return
    
    thresh = find_optimal_threshold(sim_matrix)
    thresh = min(max_thresh, thresh)
    heat, cluster_df = hierarchical_clustering(sim_matrix, thresh)

    if figpath != "":
        heat.figure.suptitle(f"Threshold: {thresh:.2f}")
        heat.figure.tight_layout()
        heat.savefig(figpath)

    return cluster_df

def append_GO_clusters_to_depth_df(sim_matrix_cache_folder, depth_df, figpath, lib, project_name, max_thresh=0):
        d = depth_df
        d["Top_GO_Cluster"] = ""
        d["GO_Cluster"] = ""
        for subont in ["BP","CC","MF"]:
            print(subont)
            sim_matrix_path = os.path.join(sim_matrix_cache_folder, f"sim_matrix_{subont}.csv")
            figfile = os.path.join(figpath,f"go.heat.{subont}.{lib}.{project_name}.png")
            clusters = get_clusters_from_sim_matrix(sim_matrix_path, figpath=figfile, max_thresh=max_thresh)
            if clusters is not None:
                d.loc[clusters.index, "GO_Cluster"] = clusters["GO_Cluster"].astype(str) + f"_{subont}"

        # Hacky way of removing ties
        eps = np.random.uniform(0,1,len(d))
        d["epsilon"] = (1+eps/1e5) * d["Combined FDR"]
        d['Top_GO_Cluster'] = d.groupby('GO_Cluster')['epsilon'].transform(lambda x: x == x.min())
        del d["epsilon"]

        #d['Top_GO_Cluster'] = d.groupby('GO_Cluster')['Combined FDR'].transform(lambda x: x == x.min())
        # For cases where there are multiple True values (ties), keep only the first occurrence # doesn't work?
        #d['Top_GO_Cluster'] = d.groupby('GO_Cluster').apply(lambda x: x.assign(Top_GO_Cluster=(x['Top_GO_Cluster'].cumsum() == 1))).reset_index(drop=True)['Top_GO_Cluster']

        lastcols = ["Configurations","Genes"]
        cols = [col for col in d.columns if col not in lastcols] + lastcols
        d = d[cols]
        return d

if __name__ == "__main__":

    config_file_path = os.path.join("config", "config.yaml")
    with open(config_file_path, 'r') as stream:
        config = yaml.safe_load(stream)

    libs = config.get('libraries', [])
    lib_names = {lib.split(".gmt")[0]: lib for lib in libs}
    project_name = config.get('project_name')
    qval = config.get('qval')
    go_sem_sim_max_distance = 100#config.get("go_sem_sim_max_distance")
    organismKEGG = config.get("organismKEGG")
    match(organismKEGG):
        case "mmu": 
            org = "Org.Mm.eg.db"
        case "hsa": 
            org = "Org.Hs.eg.db"
        case _: 
            raise Exception(f"Organism not supported: {organismKEGG}")
        
    savepath = os.path.join("results",project_name,"combined")
    sim_matrix_cache_folder = os.path.join("results",project_name,".cache")
    figpath = os.path.join("results",project_name,"figures")

    for lib in lib_names:
        if not lib.startswith("GO"): 
            continue
        depth_df_file = os.path.join(savepath, f"syn.depth.{lib}.{project_name}.csv")
        depth_df = pd.read_csv(depth_df_file, index_col=0)
        depth_df = append_GO_clusters_to_depth_df(sim_matrix_cache_folder, depth_df, figpath, lib, project_name, max_thresh=go_sem_sim_max_distance)
        depth_df.to_csv(depth_df_file)

        summary_dict_file = os.path.join(savepath, f"syn.summary_dict.{project_name}.txt")
        with open(summary_dict_file, "rb") as f:
            summary_dict = pickle.load(f)
        summary_dict[lib]["depth_df"] = depth_df
        pickler(summary_dict, summary_dict_file)