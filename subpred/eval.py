from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, VarianceThreshold
from sklearn.linear_model import SGDClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.svm import SVC, LinearSVC
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    balanced_accuracy_score,
    f1_score,
)
from sklearn.model_selection import (
    LeaveOneOut,
    cross_val_score,
    cross_val_predict,
    train_test_split,
    GridSearchCV,
)
from joblib import Parallel, delayed
import numpy as np
import pandas as pd

from .custom_transformers import PSSMSelector


def __encode_labels(labels: pd.Series) -> np.ndarray:
    labels_unique = labels.unique()
    labels_unique = np.sort(labels_unique)
    labels_dict = {label: pos for pos, label in enumerate(labels_unique)}
    return labels.map(labels_dict).values


def preprocess_pandas(df_features: pd.DataFrame, labels: pd.Series, return_names=False):
    X = df_features.to_numpy()
    y = __encode_labels(labels)
    ret = (X, y)
    if return_names:
        feature_names = df_features.columns.to_numpy(dtype=str)
        sample_names = df_features.index.to_numpy(dtype=str)
        ret += (feature_names, sample_names)
    return ret


def get_independent_test_set(X, y, sample_names=None, random_state=42, test_size=0.2):
    if sample_names is not None:
        return train_test_split(
            X,
            y,
            sample_names,
            stratify=y,
            random_state=random_state,
            shuffle=True,
            test_size=test_size,
        )
    else:
        return train_test_split(
            X,
            y,
            stratify=y,
            random_state=random_state,
            shuffle=True,
            test_size=test_size,
        )


def __get_cv_method(cross_val_method: str):
    if cross_val_method == "5fold":
        return 5
    elif cross_val_method.upper() == "LOOCV":
        return LeaveOneOut()
    else:
        raise ValueError(f"Unsupported CV method: {cross_val_method}")


def optimize_hyperparams(
    X_train,
    y_train,
    feature_transformer=None,
    feature_names=None,
    kernel="rbf",
    C=[1, 0.1, 10],
    gamma=["scale", 0.01, 0.1, 1],
    class_weight=["balanced", None],
    dim_reduction=None,
    decision_function_shape=None,
    verbose=True,
    remove_zero_var=False,
    select_k_steps=1,
    max_k_to_select: int = None,
    cross_val_method: str = "5fold",
    n_jobs: int = -1,
):
    pipe_list = list()
    param_grid = dict()
    if feature_transformer == "pssm":
        if feature_names is None:
            raise ValueError("feature names need to be provided for PSSM filtering")
        pipe_list.append(PSSMSelector(feature_names=feature_names))
        param_grid.update(
            {
                "pssmselector__uniref_threshold": [50, 90, "all"],
                "pssmselector__iterations": [1, 3, "all"],
            }
        )
    if remove_zero_var:
        pipe_list.append(VarianceThreshold(threshold=0))
    pipe_list.append(StandardScaler())

    if dim_reduction == "pca":
        pipe_list.append(PCA())
        param_grid.update({"pca__n_components": np.linspace(0.8, 0.99, 20)})
        pipe_list.append(StandardScaler())
    elif dim_reduction == "kbest":
        pipe_list.append(SelectKBest())
        param_grid.update(
            {
                "selectkbest__k": range(
                    1,
                    max_k_to_select if max_k_to_select else X_train.shape[1],
                    select_k_steps,
                )
            }
        )
        pipe_list.append(StandardScaler())
    if kernel == "rbf":
        pipe_list.append(SVC())
        param_grid.update(
            {"svc__class_weight": class_weight, "svc__C": C, "svc__gamma": gamma,}
        )
        if len(np.unique(y_train)) > 2:
            param_grid.update(
                {
                    "svc__decision_function_shape": ["ovo", "ovr"]
                    if not decision_function_shape
                    else decision_function_shape,
                }
            )
    elif kernel == "linear":
        pipe_list.append(LinearSVC())
        param_grid.update(
            {
                "linearsvc__class_weight": class_weight,
                "linearsvc__C": C,
                "linearsvc__dual": [True, False],
                "linearsvc__max_iter": [1e8],
            }
        )
        if len(np.unique(y_train)) > 2:
            param_grid.update(
                {
                    "linearsvc__multi_class": ["ovr", "crammer_singer"]
                    if not decision_function_shape
                    else decision_function_shape,
                }
            )

    pipe = make_pipeline(*pipe_list)
    gsearch = GridSearchCV(
        estimator=pipe,
        param_grid=param_grid,
        cv=__get_cv_method(cross_val_method),
        scoring="f1_macro",
        n_jobs=n_jobs,
        return_train_score=True,
        refit=True,
    )
    gsearch.fit(X_train, y_train)
    if verbose:
        print(gsearch.best_params_)
        print(gsearch.best_score_.round(3))
    return gsearch


def get_confusion_matrix(X_test, y_test, clf, labels: pd.Series = None):
    y_pred = clf.predict(X_test)

    df_confusion_matrix = pd.DataFrame(confusion_matrix(y_test, y_pred))

    if labels is not None:
        labels_unique = labels.unique()
        labels_unique = np.sort(labels_unique)
        df_confusion_matrix = df_confusion_matrix.set_index(labels_unique)
        df_confusion_matrix.columns = labels_unique
        df_confusion_matrix.index = df_confusion_matrix.index.rename("observed")
        df_confusion_matrix.columns = df_confusion_matrix.columns.rename("predicted")
    return df_confusion_matrix


def get_classification_report(
    X_test, y_test, clf, labels: pd.Series = None, add_balanced_accuracy: bool = False
):
    y_pred = clf.predict(X_test)

    report_dict = classification_report(y_true=y_test, y_pred=y_pred, output_dict=True,)

    df_report = pd.DataFrame.from_dict(report_dict).T
    df_report = df_report.astype({"support": "int"})
    df_report = df_report.round(3)
    df_report = df_report.drop("accuracy")

    if labels is not None:
        # Since numerical labels are assigned alphabetically by __encode_labels
        labels_unique = labels.unique()
        labels_unique = np.sort(labels_unique)
        labels_unique_dict = {
            str(pos): label for pos, label in enumerate(labels_unique)
        }
        df_report = df_report.rename(index=labels_unique_dict)

    if add_balanced_accuracy:
        balanced_accuracy = balanced_accuracy_score(y_true=y_test, y_pred=y_pred).round(
            3
        )
        # Accuracy is calculated on all samples in the dataset
        support = len(y_test)
        df_report.loc["balanced accuracy"] = ["", "", balanced_accuracy, support]

    return df_report


def quick_test(df_features, labels: pd.Series):
    X, y, feature_names, sample_names = preprocess_pandas(
        df_features, labels, return_names=True
    )
    gsearch = optimize_hyperparams(X, y, dim_reduction=None)
    print("Default ", gsearch.best_score_.round(3))
    gsearch = optimize_hyperparams(X, y, dim_reduction="pca")
    print("PCA", gsearch.best_score_.round(3))
    gsearch = optimize_hyperparams(X, y, dim_reduction="kbest")
    print("Kbest", gsearch.best_score_.round(3))


# kwargs are passed to gridsearch in inner loop
def nested_loocv(
    df_features: pd.DataFrame,
    labels: pd.Series,
    verbose: bool = False,
    decimal_places: int = 3,
    n_jobs_outer_loop: int = -1,
    n_jobs_inner_loop: int = 1,
    **kwargs,
):
    X, y, feature_names, sample_names = preprocess_pandas(
        df_features, labels, return_names=True
    )

    
    def outer_loop(
        X: np.ndarray,
        y: np.ndarray,
        feature_names: np.ndarray,
        train_index: np.ndarray,
        test_index: np.ndarray,
    ):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        # inner loop
        gsearch = optimize_hyperparams(
            X_train,
            y_train,
            verbose=verbose,
            feature_names=feature_names,
            cross_val_method="LOOCV",
            n_jobs=n_jobs_inner_loop,
            **kwargs,
        )
        best_estimator = gsearch.best_estimator_

        train_score = gsearch.best_score_

        y_pred = best_estimator.predict(X_test)

        # y_pred and y_test are both numpy arrays, this does not work for multi-label clf!
        return [train_score, y_pred[0], y_test[0]]

    splits = [(train_index, test_index) for train_index, test_index in LeaveOneOut().split(X)]

    # setting all numpy arrays to read-only, to avoid race conditions
    results = Parallel(n_jobs=n_jobs_outer_loop, mmap_mode="r")(
        delayed(outer_loop)(X, y, feature_names, train_index, test_index)
        for train_index, test_index in splits
    )

    df_results = pd.DataFrame.from_records(results, columns=["train_scores", "y_pred", "y_true"])

    df_results_mean = pd.DataFrame(columns=["train", "test"])
    f1_train = np.mean(df_results.train_scores)
    f1_test = f1_score(y_true=df_results.y_true, y_pred=df_results.y_pred, average="macro")
    df_results_mean.loc["F1 (macro)"] = [f1_train, f1_test]

    return df_results_mean.round(decimal_places)


# kwargs are passed to hyperparam optimization
def full_test(
    df_features: pd.DataFrame,
    labels: pd.Series,
    repetitions: int = 10,
    cross_val_method: str = "5fold",
    verbose: bool = False,
    **kwargs,
):
    X, y, feature_names, sample_names = preprocess_pandas(
        df_features, labels, return_names=True
    )
    df_params = pd.DataFrame()
    scores = []
    for random_seed in range(repetitions):
        (
            X_train,
            X_test,
            y_train,
            y_test,
            sample_names_train,
            sample_names_test,
        ) = get_independent_test_set(
            X, y, sample_names, test_size=0.2, random_state=random_seed
        )

        gsearch = optimize_hyperparams(
            X_train,
            y_train,
            verbose=verbose,
            feature_names=feature_names,
            cross_val_method=cross_val_method,
            **kwargs,
        )
        df_params[random_seed] = pd.Series(gsearch.best_params_, dtype="object")
        best_estimator = gsearch.best_estimator_

        y_pred_train = cross_val_predict(
            best_estimator,
            X_train,
            y_train,
            cv=__get_cv_method(cross_val_method),
            n_jobs=-1,
        )
        df_report_train = get_classification_report(
            X_train, y_pred_train, best_estimator, labels
        )
        # print(get_confusion_matrix(X_test, y_test, best_estimator, labels))
        df_report_test = get_classification_report(
            X_test, y_test, best_estimator, labels
        )
        for label in labels.unique():
            scores.append([label, df_report_train.loc[label, "f1-score"], "train"])
            scores.append([label, df_report_test.loc[label, "f1-score"], "test"])
        # print(df_report)

    df_scores = pd.DataFrame(scores, columns=["label", "F1 score", "dataset"])
    return df_scores, df_params


# def print_validation_results(y_true_, y_pred_, labels_unique):
#     report_dict = classification_report(
#         y_true=y_true_, y_pred=y_pred_, output_dict=True
#     )
#     labels_unique = sorted(labels_unique)
#     tmp_dict = {
#         labels_unique[i]: report_dict[str(i)] for i in range(len(labels_unique))
#     }
#     tmp_dict.update(
#         {"Macro": report_dict["macro avg"], "Weighted": report_dict["weighted avg"]}
#     )
#     # TODO multiclass
#     report_df = pd.DataFrame.from_dict(tmp_dict)
#     confusion_matrix_df = pd.DataFrame(
#         confusion_matrix(y_true_, y_pred_), columns=labels_unique, index=labels_unique,
#     )
#     return report_df, confusion_matrix_df


# def validate_test_set(X_train, y_train, X_test, y_test, gsearch):
#     best_estimator = gsearch.best_estimator_
#     best_scores = cross_val_score(
#         estimator=clone(best_estimator), X=X_train, y=y_train, scoring="f1_macro"
#     )
#     print(f"Train scores: {best_scores.mean().round(3)}+-{best_scores.std().round(3)}")

#     y_pred = best_estimator.predict(X_test)
#     y_true = y_test.copy()

#     report_df, confusion_matrix_df = print_validation_results(
#         y_true, y_pred, labels_unique=["Amino", "Sugar"]
#     )
#     return report_df.round(3), confusion_matrix_df


def models_quick_compare(X_train, y_train):
    records = []
    for estimator in [
        LinearSVC(random_state=0, max_iter=1e6),
        LinearSVC(random_state=0, max_iter=1e6, class_weight="balanced"),
        SVC(random_state=0,),
        SVC(random_state=0, class_weight="balanced"),
        GaussianNB(),
        KNeighborsClassifier(),
        RandomForestClassifier(random_state=0,),
        RandomForestClassifier(random_state=0, class_weight="balanced"),
        SGDClassifier(random_state=0),
    ]:
        pipeline = make_pipeline(StandardScaler(), estimator)
        scores = cross_val_score(pipeline, X_train, y_train, scoring="f1_macro", cv=5)
        records.append([str(estimator)] + scores.tolist())

    df_cv = pd.DataFrame.from_records(
        records, columns=["est"] + list(range(5))
    ).set_index("est")
    mean = df_cv.mean(axis=1)
    std = df_cv.std(axis=1)
    df_cv = df_cv.assign(mean=mean)
    df_cv = df_cv.assign(std=std)
    return df_cv.round(3).sort_index()

