# %%
import sys

import pandas as pd

sys.path.append("../../src")
from main.transporter_dataset import create_dataset
from main.eval import full_test
from main.compositions import calculate_aac, calculate_paac
from main.pssm import calculate_pssms_notebook
import matplotlib.pyplot as plt
import seaborn as sns

LOG_FILE = "../../logs/athaliana_amino_sugar_potassium.log"


# %%
sodium = ["Q2UVJ5", "O04034", "Q9FKP1", "Q9LJI2", "Q9SYG9", "Q84TI7"]
gdu = ["O81775", "Q9SW07", "Q9FHH5", "Q8S8A0", "Q3E965", "Q3EAV6", "Q3E8L0"]
df_all = create_dataset(
    keywords_substrate_filter=[
        "Amino-acid transport",
        "Electron transport",
        "Sugar transport",
        "Potassium transport",
    ],
    keywords_component_filter=["Transmembrane"],
    keywords_transport_filter=["Transport"],
    input_file="../../data/raw/swissprot/uniprot-reviewed_yes.tab.gz",
    multi_substrate="integrate",
    outliers=gdu + sodium,
    verbose=True,
    tax_ids_filter=[3702],
    output_log=LOG_FILE,
    sequence_clustering=70,
)


# %% [markdown]
# # Feature generation

# %%
df_aac = calculate_aac(df_all.sequence)
df_paac = calculate_paac(df_all.sequence)
df_pssm = calculate_pssms_notebook(df_all.sequence)
df_combined = pd.concat([df_aac, df_paac, df_pssm], axis=1)
labels = df_all.keywords_transport
labels.value_counts()

# %% [markdown]
# # Function

# %%
def test_case(df_features, labels: pd.Series, test_name: str, **kwargs):
    for dim_reduction in [None, "pca", "kbest"]:
        print("#" * 60)
        print("Feature reduction: ", dim_reduction)
        print("#" * 60)
        df_scores, df_params = full_test(
            df_features, labels, kernel="rbf", dim_reduction=dim_reduction, **kwargs
        )
        df_scores.label = df_scores.label.map(lambda x: x.replace(" transport", ""))

        df_params.to_csv(
            f"results/{test_name}_{dim_reduction}_params.tsv", sep="\t"
        )        
        df_scores_piv = df_scores.groupby(["label", "dataset"]).mean().unstack(1)
        df_scores_piv.loc["mean"] = df_scores_piv.mean()
        # display(df_scores_piv.round(3))
        df_scores_piv.round(3).to_csv(
            f"results/{test_name}_{dim_reduction}_scores.tsv", sep="\t"
        )

        g = sns.barplot(data=df_scores, y="F1 score", x="label", hue="dataset")
        g.set_ylim((0, 1))
        g.set_xlabel("")
        plt.savefig(f"results/{test_name}_{dim_reduction}_barplot.png")
        plt.close()
        # plt.show()

        g = sns.boxplot(data=df_scores, y="F1 score", x="label", hue="dataset")
        g.set_xlabel("")
        g.set_ylim((0, 1))
        plt.savefig(f"results/{test_name}_{dim_reduction}_boxplot.png")
        plt.close()
        # plt.show()


# %% [markdown]
# Eval AAC

# %%
test_case(
    df_aac,
    labels,
    test_name="AAC",
    decision_function_shape=["ovo"],
    gamma=["scale"],
    class_weight=["balanced"],
    C=[1, 10, 100],
)

# %% [markdown]
# # Eval PAAC

# %%
test_case(
    df_paac,
    labels,
    test_name="PAAC",
    decision_function_shape=["ovo"],
    gamma=["scale"],
    C=[1, 10, 100],
    class_weight=["balanced"],
)

# %% [markdown]
# # PSSM

# %%
test_case(
    df_pssm,
    labels,
    test_name="PSSM",
    decision_function_shape=["ovo"],
    gamma=["scale"],
    C=[1, 10, 100],
    class_weight=["balanced"],
    feature_filter="pssm",
    select_k_steps=10,
)


# %% [markdown]
# # Combined

# %%
test_case(
    pd.concat([df_aac, df_paac, df_pssm], axis=1),
    labels,
    test_name="Combined",
    decision_function_shape=["ovo"],
    gamma=["scale"],
    C=[1, 10, 100],
    class_weight=["balanced"],
    feature_filter="pssm",
    select_k_steps=10,
)
