import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import matplotlib.patches as mpatches
from .cdhit import cd_hit


def clustermap(df_feature):
    g = sns.clustermap(df_feature, method="ward", yticklabels=[], xticklabels=[],)
    return g


def labeled_clustermap(
    df_feature,
    annotation: pd.Series,
    xlabels = [],
    ylabels = [], 
    colors: list = ["cyan", "magenta", "yellow", "green", "orange"],
    legend_loc: str = "upper right",
    legend_bbox: tuple = (1.05, 1.25),
    legend_frame: bool = True,
    legend_fontsize: int = 10,
):
    assert len(df_feature.index) == len(annotation)

    color_map = dict(zip(sorted(annotation.unique()), colors))

    g = sns.clustermap(
        df_feature,
        method="ward",
        xticklabels=xlabels,
        yticklabels=ylabels,
        row_colors=annotation.map(color_map),
    )
    legend = g.ax_heatmap.legend(
        loc=legend_loc,
        bbox_to_anchor=legend_bbox,
        handles=[mpatches.Patch(color=c, label=l) for l, c in color_map.items()],
        frameon=legend_frame,
    )
    legend.set_title(title=annotation.name, prop={"size": legend_fontsize})

    return color_map, g


def get_clusters(df_features, n_clusters=2):

    cluster = AgglomerativeClustering(n_clusters=n_clusters, linkage="ward")

    cluster.fit(df_features)

    cluster_list = []
    for label in np.unique(cluster.labels_):
        cluster_list.append(df_features.index[cluster.labels_ == label].tolist())

    return cluster_list


def cluster_samples_plot(sequences: pd.Series, labels: pd.Series):
    records = []
    for clustering_threshold in [40, 50, 60, 70, 80, 90, 100]:
        cluster_proteins = cd_hit(sequences, identity_threshold=clustering_threshold)
        labels_clustered = labels.loc[cluster_proteins]
        sample_counts = labels_clustered.value_counts()
        for label in sample_counts.index.unique():
            records.append([clustering_threshold, label, sample_counts.loc[label]])

    df_records = pd.DataFrame.from_records(
        records, columns=["Threshold", "Substrate", "Proteins"]
    )
    sns.lineplot(data=df_records, x="Threshold", y="Proteins", hue="Substrate")
