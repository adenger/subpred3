import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.svm import SVC
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.feature_selection import f_classif, chi2
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score, f1_score, make_scorer
from imblearn.under_sampling import RandomUnderSampler
from joblib import Parallel, delayed

from .cdhit import cd_hit
from .util import perform_pca, get_feature_score, get_clusters

# plot long form table returned from eval.full_test
def plot_full_test(df_scores):
    g = sns.barplot(data=df_scores, x="label", y="F1 score", hue="dataset")
    g.set(ylim=(0, 1))
    return g


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


def pca_plot_2d(
    df_feature,
    labels: pd.Series,
    labels2: pd.Series = None,
    figsize: tuple = (10, 6),
    marker_size: int = None,
    marker_opacity:float = None
):
    # Important: labels parameter must have a name.
    df_pca2 = perform_pca(df_feature, labels, n_components=2)
    plt.figure(figsize=figsize)
    sns.set_style("darkgrid")
    kwargs = {}
    if marker_size:
        kwargs.update({"s": marker_size})
    if marker_opacity:
        assert marker_opacity >= 0 and marker_opacity <= 1
        kwargs.update({"alpha": marker_opacity})
    return sns.scatterplot(
        data=df_pca2, x="PC1", y="PC2", hue=labels.name, style=labels2, **kwargs
    )


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


def corr_heatmap(df_feature, figsize=(15, 10)):
    plt.figure(figsize=figsize)
    return sns.heatmap(
        df_feature.corr(), cmap="YlGnBu", vmin=-1, vmax=1, annot=True, fmt=".2f"
    )


def clustermap(df_feature):
    return sns.clustermap(
        df_feature,
        method="ward",
        yticklabels=[],
        xticklabels=[],
    )


def labeled_clustermap(
    df_feature,
    annotation: pd.Series,
    xlabels=[],
    ylabels=[],
    colors: list = ["cyan", "magenta", "yellow", "green", "orange", "red"],
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


# TODO refactor with other function?
def downsampling_plot(
    df_features,
    labels: pd.Series,
    min_features=16,
    random_seeds=list(range(50)),
    figsize=(10, 7),
):
    """
    Perform downsampling and evaluation, until number of samples reaches min_features.
    """

    def test_case(df_features, labels, n_samples, random_seed):
        df_features_sampled, _, labels_sampled, _ = train_test_split(
            df_features,
            labels,
            stratify=labels,
            random_state=random_seed,
            train_size=n_samples,
            shuffle=True,
        )
        label_encoder = LabelEncoder()
        X = df_features_sampled.to_numpy()
        y = label_encoder.fit_transform(labels_sampled)

        estimator = make_pipeline(StandardScaler(), SVC(class_weight="balanced"))
        scorers_dict = dict()
        labels_unique = labels.unique()
        labels_unique_numerical = label_encoder.transform(labels_unique)
        for label, label_numerical in zip(labels_unique, labels_unique_numerical):
            scorers_dict[f"F1 {label}"] = make_scorer(
                f1_score, pos_label=label_numerical
            )

        scorers_dict["F1 macro"] = make_scorer(f1_score, average="macro")

        records = []
        for scorer_name, scorer in sorted(scorers_dict.items()):
            cv_results = cross_val_score(estimator, X, y, n_jobs=1, scoring=scorer)
            for cv_result in cv_results:
                records.append([scorer_name, n_samples, cv_result])
        return records

    n_samples_list = list(range(min_features, df_features.shape[0] - 1))
    records_list = Parallel(n_jobs=-1)(
        delayed(test_case)(df_features, labels, n_samples, random_seed)
        for n_samples in n_samples_list
        for random_seed in random_seeds
    )
    records = []
    for sl in records_list:
        for ssl in sl:
            records.append(ssl)
    results_df = pd.DataFrame.from_records(
        records,
        columns=["Score name", "Total Samples", "Score"],
    )
    plt.figure(figsize=figsize)
    g = sns.lineplot(
        data=results_df,
        x="Total Samples",
        y="Score",
        hue="Score name",
    )
    return g


def downsample_majority_class_plot(
    df_feature_unclustered,
    labels: pd.Series,
    min_class_sample_fractions=[x / 100 for x in range(40, 101, 5)],
    random_seeds=list(range(50)),
    figsize=(10, 7),
    n_jobs=4,
    include_macro: bool = True,
):
    """
    Binary classification plot.
    Performs downsampling of majority class, until it contains less samples than minority class.
    A simple evaluation is performed at every fraction, and results are saved in plot.
    The process is repeated for multiple random seeds, to rule out bias introduced by undersampling
    """

    def test_case(df_features, labels, min_class_sample_fraction, random_seed):
        rus = RandomUnderSampler(
            random_state=random_seed, sampling_strategy=min_class_sample_fraction
        )
        df_features_sampled, labels_sampled = rus.fit_resample(df_features, labels)

        label_encoder = LabelEncoder()
        X = df_features_sampled.to_numpy()
        y = label_encoder.fit_transform(labels_sampled)

        estimator = make_pipeline(StandardScaler(), SVC(class_weight="balanced"))

        labels_unique = sorted(labels.unique())
        labels_numerical = label_encoder.transform(labels_unique)

        scorers_dict = dict()
        for label, label_numerical in zip(labels_unique, labels_numerical):
            scorer_name = f"F1 {label}"
            scorer = make_scorer(f1_score, pos_label=label_numerical)
            scorers_dict[scorer_name] = scorer

        if include_macro:
            scorers_dict["F1 macro"] = make_scorer(f1_score, average="macro")

        records = []
        for scorer_name, scorer in sorted(scorers_dict.items()):
            cv_results = cross_val_score(estimator, X, y, n_jobs=1, scoring=scorer)
            for cv_result in cv_results:
                records.append([scorer_name, min_class_sample_fraction, cv_result])
        return records

    records_list = Parallel(n_jobs=n_jobs)(
        delayed(test_case)(
            df_feature_unclustered, labels, min_class_sample_fraction, random_seed
        )
        for min_class_sample_fraction in min_class_sample_fractions
        for random_seed in random_seeds
    )
    records = []
    for sl in records_list:
        for ssl in sl:
            records.append(ssl)

    labels_unique = sorted(labels.unique())
    x_lab = f"|{labels_unique[0]}|/|{labels_unique[1]}|"
    results_df = pd.DataFrame.from_records(
        records,
        columns=["Score name", x_lab, "Score"],
    )
    plt.figure(figsize=figsize)
    g = sns.lineplot(
        data=results_df,
        x=x_lab,
        y="Score",
        hue="Score name",
    )
    return g
