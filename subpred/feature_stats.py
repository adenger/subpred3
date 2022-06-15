import pandas as pd

def get_feature_stats(
    df_features, df_labels_, labels=["Amino-acid transport", "Sugar transport"]
):
    df_stats = pd.concat(
        {
            "corr": df_features.corrwith(
                df_labels_.labels.transform(lambda x: 1.0 if x == labels[1] else 0.0)
            ),
            "mean": df_features.mean(),
            "std": df_features.std(),
        },
        axis=1,
    )

    df_stats["sum"] = df_stats.sum(axis=1)
    df_stats["corr_abs"] = df_stats["corr"].abs()

    df_stats["mean0"] = df_features.loc[
        df_labels_[df_labels_.labels == labels[0]].index
    ].mean()
    df_stats["mean1"] = df_features.loc[
        df_labels_[df_labels_.labels == labels[1]].index
    ].mean()

    df_stats["median0"] = df_features.loc[
        df_labels_[df_labels_.labels == labels[0]].index
    ].median()
    df_stats["median1"] = df_features.loc[
        df_labels_[df_labels_.labels == labels[1]].index
    ].median()

    df_stats["mediandiff"] = (df_stats["median0"] - df_stats["median1"]).abs()
    df_stats = df_stats.sort_values("mediandiff", ascending=False)
    return df_stats


def get_pca_contributions(gsearch, feature_names):
    return (
        pd.DataFrame.from_records(
            [
                feature_names,
                gsearch.best_estimator_["pca"].components_[0],
                gsearch.best_estimator_["pca"].components_[1],
                gsearch.best_estimator_["pca"].components_[2],
            ],
            index=["Feature", "Contrib. PC1", "Contrib. PC2", "Contrib. PC3"],
        )
        .T.set_index("Feature")
        .sort_values("Contrib. PC1", ascending=False)
    )


def get_removed_features(gsearch, feature_names):
    return feature_names[~gsearch.best_estimator_["selectkbest"].get_support()]
