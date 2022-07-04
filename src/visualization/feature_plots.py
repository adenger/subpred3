import argparse
import pandas as pd
import logging
import seaborn as sns
import matplotlib.pyplot as plt
from yellowbrick.features import pca_decomposition, radviz
from sklearn.preprocessing import LabelEncoder

log = logging.getLogger("FEATUREPLOT")



def __get_label_colors(df_labels):
    labels_unique = df_labels.labels.unique()
    if len(labels_unique) > 6:
        log.warning("not enough colors for all samples, add more.")
    colors = ["red", "green", "blue", "yellow", "magenta", "cyan"]
    colors_dict = {labels_unique[i]: colors[i] for i in range(len(labels_unique))}
    colors_labels = [colors_dict[label] for label in df_labels.labels]
    return colors_labels, colors_dict


def __pandas_to_sklearn(df_features, df_labels):
    X = df_features.to_numpy()
    label_enc = LabelEncoder()
    y = label_enc.fit_transform(df_labels.values.ravel())
    feature_names = df_features.columns.values
    return X, y, feature_names, label_enc


def create_plot(
    input_features: str, input_labels: str, type, log_file=None, verbose=True
):

    df_features = pd.read_table(input_features, index_col=0)
    df_labels = pd.read_table(input_labels, index_col=0)

    logging.basicConfig(
        format="[%(levelname)s] %(filename)s: %(message)s",
        level=logging.DEBUG if verbose else logging.INFO,
        filename=log_file,
    )

    log.debug(df_labels.value_counts())

    if not df_features.index.identical(df_labels.index):
        log.warning("indices don't match.")
        common_index = df_features.index[
            df_features.index.isin(df_labels.index)
        ].tolist()
        df_features = df_features.loc[common_index]
        df_labels = df_labels.loc[common_index]

    label_colors, colors_dict = __get_label_colors(df_labels)
    X, y, feature_names, label_enc = __pandas_to_sklearn(df_features, df_labels)

    print(colors_dict)

    if type == "sample":
        sns.set(font_scale=0.5)
        sns.clustermap(
            data=df_features,
            xticklabels=df_features.columns,
            # yticklabels=None,
            row_colors=label_colors,
            method="ward",
        )
        sns.set(font_scale=1)
    elif type == "corr_samples":
        sns.set(font_scale=0.5)
        sns.clustermap(
            data=df_features.transpose().corr(),
            xticklabels=df_labels.index,
            yticklabels=df_labels.index,
            col_colors=label_colors,
            row_colors=label_colors,
            method="ward",
        )
        sns.set(font_scale=1)
    elif type == "corr_features":
        sns.set(font_scale=0.5)
        sns.clustermap(
            data=df_features.corr(),
            # xticklabels=df_features.columns,
            # yticklabels=df_features.columns,
            method="ward",
        )
        sns.set(font_scale=1)
    elif type == "pca":
        pca_decomposition(
            X,
            y,
            features=feature_names,
            classes=label_enc.classes_,
        )
    elif type == "radviz":
        radviz(
            X,
            y,
            features=feature_names,
            classes=label_enc.classes_,
        )
    else:
        raise ValueError(f"invalid plot type: {type}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create plot from feature tsv")
    parser.add_argument(
        "--input-features",
        required=True,
    )
    parser.add_argument(
        "--input-labels",
        required=True,
    )
    parser.add_argument(
        "--type",
        choices={"sample", "corr_samples", "pca", "radviz", "corr_features"},
        required=True,
    )
    parser.add_argument(
        "--log-file",
    )
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    create_plot(
        args.input_features, args.input_labels, args.type, args.log_file, args.verbose
    )
