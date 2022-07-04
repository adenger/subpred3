import argparse
import pandas as pd
from soft import read_soft
from qnorm import quantile_normalize


def __read_file(soft_file_name: str) -> pd.DataFrame:
    soft_data = read_soft(soft_file_name)
    soft_df = (
        pd.DataFrame(soft_data[1:], columns=soft_data[0])
        .replace("null", 0.0)
        .drop("ID_REF", axis=1)
        .rename(columns={"IDENTIFIER": "Gene_Name"})
    )
    return soft_df


def __annotate(soft_df: pd.DataFrame, mapping_file_name: str) -> pd.DataFrame:
    mapping_df = pd.read_csv(
        mapping_file_name, sep="\t", names=["Uniprot", "Type", "Value"]
    )
    gene_name_df = (
        mapping_df[mapping_df.Type == "Gene_Name"]
        .drop("Type", axis=1)
        .rename(columns={"Value": "Gene_Name"})
    )
    soft_df_uniprot = soft_df.merge(gene_name_df, on="Gene_Name", how="left").drop(
        "Gene_Name", axis=1
    )
    soft_df_uniprot = soft_df_uniprot[
        soft_df_uniprot.columns[~soft_df_uniprot.columns.str.startswith("GSM")].tolist()
        + soft_df_uniprot.columns[
            soft_df_uniprot.columns.str.startswith("GSM")
        ].tolist()
    ]
    soft_df_uniprot = soft_df_uniprot.astype(
        {
            col_name: "float"
            for col_name in soft_df_uniprot.columns
            if col_name != "Uniprot"
        }
    )
    return soft_df_uniprot


def __preprocess_values(soft_df, remove_zero_var: bool, exponent: int):
    soft_df_id = soft_df[soft_df.columns[~soft_df.columns.str.startswith("GSM")]]
    soft_df_gsm = soft_df[soft_df.columns[soft_df.columns.str.startswith("GSM")]]

    if remove_zero_var:
        zero_var_mask = soft_df_gsm.apply(
            lambda row: row.mean() == 0.0 and row.var() == 0.0, axis=1
        )
        soft_df_id = soft_df_id[~zero_var_mask]
        soft_df_gsm = soft_df_gsm[~zero_var_mask]

    if exponent > 1:
        soft_df_gsm = soft_df_gsm ** exponent

    return pd.concat([soft_df_id, soft_df_gsm], axis=1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create dataset")
    parser.add_argument("-i", "--soft-file", help="Input Soft file", required=True)
    parser.add_argument(
        "-m",
        "--mapping-file",
        help="Uniprot id mapping dat file for organism",
        required=True,
    )
    parser.add_argument("-o", "--output-file", help="Output tsv file", required=True)
    parser.add_argument(
        "-e", "--exp", help="exponent for data (revert log)", default=1, type=int
    )
    parser.add_argument(
        "--quantile-normalize",
        help="normalize data before writing",
        action="store_true",
    )
    parser.add_argument(
        "--remove-zero-var",
        help="remove genes with zero variance",
        action="store_true",
    )

    args = parser.parse_args()

    soft_file_name = args.soft_file
    mapping_file_name = args.mapping_file
    out_file_name = args.output_file
    normalize_flag = args.quantile_normalize
    exponent = int(args.exp)
    remove_zero_var = args.remove_zero_var

    soft_df = __read_file(soft_file_name)

    soft_df = __annotate(soft_df=soft_df, mapping_file_name=mapping_file_name)

    soft_df = __preprocess_values(
        soft_df=soft_df, remove_zero_var=remove_zero_var, exponent=exponent
    )

    soft_df = soft_df.groupby("Uniprot").aggregate("median")

    if normalize_flag:
        soft_df = quantile_normalize(soft_df, axis=1)

    soft_df.to_csv(out_file_name, sep="\t")
