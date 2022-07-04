import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from main.cdhit import cd_hit


def labeled_clustermap(df_feature, annotation: pd.Series):
    color_map = dict(
        zip(
            sorted(annotation.unique()),
            ["cyan", "magenta", "yellow", "green", "orange"],
        )
    )
    sns.clustermap(
        df_feature,
        method="ward",
        yticklabels=[],
        xticklabels=[],
        row_colors=annotation.map(color_map),
    )
    return color_map


def get_clusters(df_features, n_clusters):

    cluster = AgglomerativeClustering(n_clusters=2, linkage="ward")

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
