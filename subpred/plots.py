import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import f_classif, chi2
from sklearn.pipeline import make_pipeline
from .cdhit import cd_hit

# plot long form table returned from eval.full_test
def plot_full_test(df_scores):
    g = sns.barplot(data=df_scores, x="label", y="F1 score", hue="dataset")
    g.set(ylim=(0, 1))
    return g


def get_feature_score(df_feature, labels: pd.Series, method: str = "f_classif"):
    func = None
    if method == "f_classif":
        func = f_classif
    elif method == "chi2":
        func = chi2
    else:
        raise ValueError(f"Inalid method: {method}")

    df_score = pd.DataFrame(
        {
            "Feature": df_feature.columns.tolist(),
            "Normalized score": func(df_feature, labels)[0],
            "Measure": [f"Feature importance ({method})"] * df_feature.shape[1],
        }
    )
    df_score["Normalized score"] = df_score["Normalized score"] / sum(
        df_score["Normalized score"]
    )
    return df_score


def feature_importance_plot(
    df_feature,
    labels: pd.Series,
    feature_name: str = "AA Frequency",
    method: str = "f_classif",
    figsize=(10, 7),
):
    long_dfs = []
    for label in labels.unique():
        df_feature_label = df_feature.loc[labels[labels == label].index.tolist()]
        df_feature_label_long = df_feature_label.melt(
            var_name="Feature", value_name="Normalized score"
        ).assign(Measure=f"{feature_name} ({label})")
        long_dfs.append(df_feature_label_long)

    df_score = get_feature_score(df_feature, labels, method=method)
    long_dfs.insert(0, df_score)

    df_plot = pd.concat(long_dfs).reset_index(drop=True)
    plt.figure(figsize=figsize)
    return sns.barplot(
        data=df_plot,
        x="Feature",
        y="Normalized score",
        hue="Measure",
        order=df_score.sort_values("Normalized score", ascending=False).Feature,
    )


def perform_pca(df_feature, labels: pd.Series, n_components: int = 2):
    pipe = make_pipeline(StandardScaler(), PCA(n_components=n_components))
    df_pca = pd.DataFrame(
        data=pipe.fit_transform(df_feature),
        columns=[f"PC{pc}" for pc in range(1, n_components + 1)],
        index=df_feature.index,
    )
    df_pca[labels.name] = labels
    return df_pca


def pca_plot_2d(df_feature, labels: pd.Series, figsize=(10, 6)):
    df_pca2 = perform_pca(df_feature, labels, n_components=2)
    plt.figure(figsize=figsize)
    sns.set_style("darkgrid")
    return sns.scatterplot(data=df_pca2, x="PC1", y="PC2", hue=labels.name)


def pca_plot_3d(df_feature, labels, figsize=(10, 10)):
    df_pca3 = perform_pca(df_feature, labels, n_components=3)
    fig = plt.figure(figsize=figsize)
    axes = plt.axes(projection="3d")
    axes.set_xlabel("PC1")
    axes.set_ylabel("PC2")
    axes.set_zlabel("PC3")
    for label in df_pca3[labels.name].unique():
        axes.scatter3D(
            df_pca3.PC1[df_pca3[labels.name] == label],
            df_pca3.PC2[df_pca3[labels.name] == label],
            df_pca3.PC3[df_pca3[labels.name] == label],
            label=label,
        )
    axes.legend()
    return axes


def corr_heatmap(df_feature):
    return sns.heatmap(df_feature.corr(), cmap="YlGnBu", vmin=-1, vmax=1,)


def clustermap(df_feature):
    return sns.clustermap(df_feature, method="ward", yticklabels=[], xticklabels=[],)


def labeled_clustermap(
    df_feature,
    annotation: pd.Series,
    xlabels=[],
    ylabels=[],
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
    # seems to be deterministic
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
