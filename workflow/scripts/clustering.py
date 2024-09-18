import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial import distance

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

def sort_cluster_colors(dn):
    """
    Sort leaves_color_list by index so that they can be clustered again by
    sns.clustermap and used as row colors
    """
    zipped = zip(dn["leaves_color_list"], dn["leaves"])
    zipped = sorted(zipped, key = lambda zipped: zipped[1])
    return [c[0] for c in zipped]

def hierarchical_clustering(heatmap_data, distance_threshold):

  linkage_matrix = hierarchy.linkage(distance.squareform(1 - heatmap_data.values), method='average')

  # Clusters
  clusters = hierarchy.fcluster(linkage_matrix, t=distance_threshold, criterion='distance')
  cluster_df = pd.DataFrame(index = heatmap_data.index, data=clusters, columns=["GO_Cluster"])

  # Heatmap
  dn = get_dendrogram(heatmap_data,linkage_matrix,thresh=distance_threshold,Plot=False, title='Hierarchical Clustering Dendrogram')
  row_colors_cl = sort_cluster_colors(dn)
  g = sns.clustermap(heatmap_data,row_linkage=linkage_matrix,col_linkage=linkage_matrix,robust=True,row_colors=row_colors_cl,col_colors=row_colors_cl) #, vmin=vmin2, vmax=vmax2, method=distmethod)
  g.ax_row_dendrogram.axvline(distance_threshold, c='magenta', linestyle='--',lw=1.5, alpha=0.5)

  return g, cluster_df