import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
import re

class PSSMSelector(BaseEstimator, TransformerMixin):
    def __init__(self, feature_names, uniref_threshold="all", iterations="all"):
        self.feature_names = feature_names
        self.uniref_threshold = uniref_threshold
        self.iterations = iterations
        self.column_regex = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY]_\d\d_\d$")

    def fit(self, X, y=None):
        if self.uniref_threshold in {50, 90}:
            has_uniref = (
                np.char.find(self.feature_names, str(self.uniref_threshold)) >= 0
            )
        elif self.uniref_threshold == "all":
            has_uniref = np.array([True] * len(self.feature_names))
        else:
            raise ValueError(f"Incorrect uniref threshold {self.uniref_threshold}")

        if self.iterations in {1, 3}:
            has_iterations = np.char.find(self.feature_names, str(self.iterations)) >= 0
        elif self.iterations == "all":
            has_iterations = np.array([True] * len(self.feature_names))
        else:
            raise ValueError(f"Incorrect iteration count: {self.iterations}")

        is_no_pssm_feature = [not bool(self.column_regex.match(col_name)) for col_name in self.feature_names]
        self.mask = np.bitwise_or(
            np.bitwise_and(has_uniref, has_iterations), is_no_pssm_feature
        )
        # self.mask = np.bitwise_and(has_uniref, has_iterations)
        return self

    def transform(self, X, y=None):
        X = np.array(X)
        X = X[:, self.mask]
        return X


class CoexpParameterSelector(BaseEstimator, TransformerMixin):
    def __init__(self, feature_names, normalized_expression=True, neighbors=(5, 3)):
        self.feature_names = feature_names
        self.neighbors = neighbors
        self.normalized_expression = normalized_expression

    def fit(self, X, y=None):
        self.mask = np.char.endswith(
            self.feature_names,
            "{}_{}_{}".format(
                self.neighbors[0],
                self.neighbors[1],
                "norm" if self.normalized_expression else "notnorm",
            ),
        )
        return self

    def transform(self, X, y=None):
        X = np.array(X)
        X = X[:, self.mask]
        return X
