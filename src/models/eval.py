import argparse
import pandas as pd
import numpy as np
import logging
from pathlib import Path
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.svm import SVC
from sklearn.feature_selection import SelectKBest, VarianceThreshold
from sklearn.model_selection import cross_validate, GridSearchCV
from sklearn.decomposition import PCA

log = logging.getLogger("EVAL")


def nested_crossval(
    X: np.array,
    y: np.array,
    params: dict,
    pipe: Pipeline,
    scoring: str,
    cv_splits: int = 5,
) -> float:
    search = GridSearchCV(
        pipe, params, cv=cv_splits, scoring=scoring, n_jobs=-1
    )  
    log.debug(f"Pipeline: \n{pipe}")
    log.debug(f"Parameters: \n{params}")
    cv_dict = cross_validate(
        search,
        X,
        y,
        cv=cv_splits,
        return_estimator=True,
        return_train_score=True,
        n_jobs=1,
        scoring=scoring,
    )

    scores_df = pd.DataFrame.from_dict(
        {a: b for a, b in cv_dict.items() if a in {"test_score", "train_score"}}
    ).transpose()
    scores_df_mean = scores_df.mean(axis=1)
    scores_df_std = scores_df.std(axis=1)
    scores_df["mean"] = scores_df_mean
    scores_df["std"] = scores_df_std
    scores_df = scores_df.round(3)
    log.debug(f"Scores: \n{scores_df}")

    params_df = pd.DataFrame.from_records(
        [cv_dict["estimator"][i].best_params_ for i in range(cv_splits)]
    ).transpose()
    log.debug(f"Best Parameters: \n{params_df}")

    return pd.concat([scores_df, params_df]).fillna("")


# def evaluate_batch(
#     input_files: list,
#     input_labels_file: str,
#     min_samples_per_class: int,
#     output_file: str,
# ):

#     log.debug("#" * 60)
#     log.debug("Machine Learning Pipeline")
#     results_df = pd.DataFrame(columns=list(range(5)) + ["mean", "std"])

#     for input_file_name in input_files:
#         log.debug("=" * 60)
#         log.debug(f"Performing evaluation for {input_file_name}:")

#         # read data
#         df_features = pd.read_table(input_file_name, index_col=0)
#         df_labels = pd.read_table(input_labels_file, index_col=0)

#         if not df_features.index.identical(df_labels.index):
#             log.debug("WARNING: Indices of labels and features did not match")
#             common_index = df_features.index[
#                 df_features.index.isin(df_labels.index)
#             ].tolist()
#             df_features = df_features.loc[common_index]
#             df_labels = df_labels.loc[common_index]

#         feature_names = df_features.columns
#         sample_names = df_features.index

#         X = df_features.to_numpy()
#         log.debug(f"Number of features: {df_features.shape[0]}")
#         log.debug(f"Labels in dataset: \n{df_labels.value_counts().to_string()}")
#         if min(df_labels.value_counts()) < min_samples_per_class:
#             raise ValueError(
#                 "Not enough samples for at least one class, skipping file."
#             )
#         label_encoder = LabelEncoder()
#         y = label_encoder.fit_transform(df_labels.to_numpy().ravel())

#         scores_df, params_df = evaluate(
#             X, y, sample_names=sample_names, feature_names=feature_names
#         )
#         scores_df["file"] = Path(input_file_name).stem
#         params_df["file"] = Path(input_file_name).stem
#         results_df = results_df.append(scores_df)
#         results_df = results_df.append(params_df)
#     log.debug("=" * 60)
#     log.debug(f"Overall results:\n{results_df.to_string(na_rep='-')}")

#     results_df.to_csv(output_file, sep="\t")


# if __name__ == "__main__":

#     pd.set_option("display.max_rows", 500)
#     pd.set_option("display.max_columns", 500)
#     pd.set_option("display.width", 1000)

#     parser = argparse.ArgumentParser(
#         description="train ML algorithm for every file in folder and write table with CV resuls"
#     )
#     parser.add_argument(
#         "-i",
#         "--input-files",
#         nargs="+",
#         required=True,
#     )
#     parser.add_argument(
#         "--input-labels",
#         required=True,
#     )
#     parser.add_argument(
#         "-a",
#         "--ml-algorithm",
#         default="svm",
#     )

#     parser.add_argument(
#         "--output-table",
#         required=True,
#     )
#     parser.add_argument(
#         "--output-log",
#         type=str,
#         required=True,
#     )
#     parser.add_argument("--verbose", action="store_true")
#     parser.add_argument("--min-samples-per-class", default=20, type=int)

#     args = parser.parse_args()

#     logging.basicConfig(
#         format="[%(levelname)s] %(filename)s: %(message)s",
#         level=logging.DEBUG if args.verbose else logging.INFO,
#         filename=args.output_log,
#     )

#     evaluate_batch(
#         input_files=args.input_files,
#         input_labels_file=args.input_labels,
#         min_samples_per_class=args.min_samples_per_class,
#         output_file=args.output_table,
#     )

# TODO GridsearchCV in cross validate?
# TODO XGBoost? Multilabel Binarizer?
# TODO most important features
# TODO logging
# TODO more scores (for binary, multi)
# TODO optimize param grid
# TODO Constant features throw errors with number of features
