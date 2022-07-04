import pandas as pd
import argparse
from qnorm import quantile_normalize
from joblib import Parallel, delayed
import os


def read_gtex_file(input_folder: str):
    exp_df = pd.read_csv(
        f"{input_folder}/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz",
        sep="\t",
        skiprows=2,
        # low_memory=False,
        index_col=0,
    )
    exp_df.drop("Description", axis=1, inplace=True)

    exp_df = exp_df[
        ~exp_df.apply(lambda row: row.mean() == 0.0 and row.var() == 0.0, axis=1)
    ]
    exp_df["Ensembl"] = exp_df.index
    exp_df["Ensembl"] = exp_df.Ensembl.apply(lambda e: e.split(".")[0])
    exp_df.set_index("Ensembl", inplace=True, verify_integrity=True)
    return exp_df


def get_hgnc_df(input_folder):
    hgnc_df = pd.read_csv(f"{input_folder}/hgnc.tsv", sep="\t")

    hgnc_df_filtered = hgnc_df[["Ensembl", "Uniprot"]]
    hgnc_df_filtered = hgnc_df_filtered[
        ~hgnc_df_filtered.Ensembl.isnull() & ~hgnc_df_filtered.Uniprot.isnull()
    ]
    hgnc_df_filtered = hgnc_df_filtered.set_index("Ensembl", verify_integrity=True)

    return hgnc_df_filtered


def write_tissue_df(tissue_name, sample_annotation_df, exp_df):
    samples_tissue = sample_annotation_df[
        sample_annotation_df.Tissue == tissue_name
    ].index.to_list()
    print(tissue_name)
    if len(samples_tissue) >= 20:
        exp_df_tissue = exp_df[
            exp_df.columns[exp_df.columns.isin(samples_tissue)].to_list()
        ]
        exp_df_tissue.to_csv(
            f"{output_folder}/gtex_{tissue_name.lower().replace(' ', '_')}.tsv",
            sep="\t",
        )
        exp_df_tissue_norm = quantile_normalize(exp_df_tissue)
        exp_df_tissue_norm.to_csv(
            f"{output_folder}/gtex_{tissue_name.lower().replace(' ', '_')}_norm.tsv",
            sep="\t",
        )
    else:
        print(f"skipped {tissue_name} ({len(samples_tissue)} samples)")


def write_output(exp_df: pd.DataFrame, output_folder):
    sample_annotation_df = pd.read_csv(
        f"{input_folder}/GTEx_v7_Annotations_SampleAttributesDS.txt",
        sep="\t",
        index_col=0,
        usecols=["SAMPID", "SMTS"],
    )
    sample_annotation_df.index.rename("Sample", inplace=True)
    sample_annotation_df.rename(columns={"SMTS": "Tissue"}, inplace=True)

    exp_df = exp_df[
        exp_df.columns[exp_df.columns.isin(sample_annotation_df.index)].to_list()
    ]
    sample_annotation_df = sample_annotation_df[
        sample_annotation_df.index.isin(exp_df.columns)
    ]

    print("starting parallel processing")
    Parallel(n_jobs=-1)(
        delayed(write_tissue_df)(tissue_name, sample_annotation_df, exp_df)
        for tissue_name in sample_annotation_df.Tissue.unique()
    )
    # for tissue_name in sample_annotation_df.Tissue.unique():
    #     samples_tissue = sample_annotation_df[
    #         sample_annotation_df.Tissue == tissue_name
    #     ].index.to_list()
    #     print(tissue_name)
    #     if len(samples_tissue) >= 20:
    #         exp_df_tissue = exp_df[
    #             exp_df.columns[exp_df.columns.isin(samples_tissue)].to_list()
    #         ]
    #         exp_df_tissue.to_csv(
    #             f"{output_folder}/gtex_{tissue_name.lower().replace(' ', '_')}.tsv.gz",
    #             sep="\t",
    #         )
    #         exp_df_tissue_norm = quantile_normalize(exp_df_tissue)
    #         exp_df_tissue_norm.to_csv(
    #             f"{output_folder}/gtex_{tissue_name.lower().replace(' ', '_')}_norm.tsv.gz",
    #             sep="\t",
    #         )
    #     else:
    #         print(f"..skipped ({len(samples_tissue)} samples)")

    # exp_df.to_pickle(f"{output_folder}/gtex.tsv.pickle")
    # sample_annotation_df.to_csv(f"{output_folder}/gtex_annotation.tsv.gz", sep="\t")
    os.makedirs(f"{output_folder}/annotation", exist_ok=True)
    sample_annotation_df.to_csv(
        f"{output_folder}/annotation/gtex_annotation.tsv", sep="\t"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create dataset")
    parser.add_argument("-i", "--input_folder", help="Input Folder", required=True)
    parser.add_argument("-o", "--output_folder", help="Output Folder", required=True)

    args = parser.parse_args()
    input_folder = args.input_folder
    output_folder = args.output_folder

    print("reading expression data...")
    exp_df = read_gtex_file(input_folder=input_folder)

    print("reading hgnc mapping data...")
    hgnc_df = get_hgnc_df(input_folder=input_folder)

    print("processing data...")
    exp_df = exp_df.join(hgnc_df, how="left")
    exp_df = exp_df[~exp_df.Uniprot.isnull()]
    exp_df = exp_df.groupby("Uniprot").aggregate("median")

    print("writing output files...")
    write_output(exp_df=exp_df, output_folder=output_folder)
