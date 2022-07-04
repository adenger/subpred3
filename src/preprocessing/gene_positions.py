import pandas as pd
import argparse
from gff3 import read_gff3


def process_genome(
    gff_df: pd.DataFrame, mapping_file_name: str, mapping_column: str
) -> pd.DataFrame:
    gff_df = gff_df[~gff_df.gene_id.isnull()]
    # gff_df = gff_df[gff_df.type == "gene"]
    # gff_df = gff_df[gff_df.biotype == "protein_coding"]
    gff_df = gff_df[["seqid", "start", "end", "strand", "gene_id"]]

    uniprot_df = pd.read_table(
        mapping_file_name, header=None, names=["Uniprot", "DB", "ID"]
    )
    uniprot_df = (
        uniprot_df[uniprot_df.DB == mapping_column]
        .drop("DB", axis=1)
        .rename(columns={"ID": "gene_id"})
    )
    gff_df = gff_df.merge(uniprot_df, how="left", on="gene_id")
    gff_df = gff_df[~gff_df.Uniprot.isnull()]

    gff_df = gff_df.rename(columns={"seqid": "chr"})
    gff_df = gff_df[["Uniprot", "chr", "start", "end", "strand"]]

    gff_df = gff_df.drop_duplicates()

    return gff_df


# def process_genome_human(gff_df: pd.DataFrame):
#     gff_df = gff_df[~gff_df.gene_id.isnull()]
#     # gff_df = gff_df[gff_df.type == "gene"]
#     # gff_df = gff_df[gff_df.biotype == "protein_coding"]
#     gff_df = gff_df[["seqid", "start", "end", "strand", "gene_id"]]

#     uniprot_df = pd.read_table(
#         mapping_file_name, header=None, names=["Uniprot", "DB", "ID"]
#     )
#     uniprot_df = (
#         uniprot_df[uniprot_df.DB == "Ensembl"]  # TODO
#         .drop("DB", axis=1)
#         .rename(columns={"ID": "gene_id"})
#     )
#     gff_df = gff_df.merge(uniprot_df, how="left", on="gene_id")
#     gff_df = gff_df[~gff_df.Uniprot.isnull()]

#     gff_df = gff_df.rename(columns={"seqid": "chr"})
#     gff_df = gff_df[["Uniprot", "chr", "start", "end", "strand"]]

#     return gff_df


# def process_genome_athaliana(gff_df: pd.DataFrame):
#     gff_df = gff_df[~gff_df.gene_id.isnull()]
#     # gff_df = gff_df[gff_df.type == "gene"]
#     # gff_df = gff_df[gff_df.biotype == "protein_coding"]
#     gff_df = gff_df[["seqid", "start", "end", "strand", "gene_id"]]

#     uniprot_df = pd.read_table(
#         mapping_file_name,
#         header=None,
#         names=["Uniprot", "DB", "ID"],
#     )
#     uniprot_df = (
#         uniprot_df[uniprot_df.DB == "EnsemblGenome"]
#         .drop("DB", axis=1)
#         .rename(columns={"ID": "gene_id"})
#     )

#     gff_df = gff_df.merge(uniprot_df, how="left", on="gene_id")
#     gff_df = gff_df[~gff_df.Uniprot.isnull()]

#     gff_df = gff_df.rename(columns={"seqid": "chr"})
#     gff_df = gff_df[["Uniprot", "chr", "start", "end", "strand"]]

#     return gff_df


# def process_genome_ecoli(gff_df: pd.DataFrame):
#     gff_df = gff_df[~gff_df.gene_id.isnull()]
#     # gff_df = gff_df[gff_df.type == "gene"]
#     # gff_df = gff_df[gff_df.biotype == "protein_coding"]
#     gff_df = gff_df[["seqid", "start", "end", "strand", "gene_id"]]

#     uniprot_df = pd.read_table(
#         mapping_file_name,
#         header=None,
#         names=["Uniprot", "DB", "ID"],
#     )
#     uniprot_df = (
#         uniprot_df[uniprot_df.DB == "EnsemblGenome"]
#         .drop("DB", axis=1)
#         .rename(columns={"ID": "gene_id"})
#     )

#     gff_df = gff_df.merge(uniprot_df, how="left", on="gene_id")
#     gff_df = gff_df[~gff_df.Uniprot.isnull()]

#     gff_df = gff_df.rename(columns={"seqid": "chr"})
#     gff_df = gff_df[["Uniprot", "chr", "start", "end", "strand"]]

#     return gff_df


# def process_genome_yeast(gff_df: pd.DataFrame):
#     gff_df = gff_df[~gff_df.gene_id.isnull()]
#     # gff_df = gff_df[gff_df.type == "gene"]
#     # gff_df = gff_df[gff_df.biotype == "protein_coding"]
#     gff_df = gff_df[["seqid", "start", "end", "strand", "gene_id"]]

#     uniprot_df = pd.read_table(
#         mapping_file_name,
#         header=None,
#         names=["Uniprot", "DB", "ID"],
#     )
#     uniprot_df = (
#         uniprot_df[uniprot_df.DB == "EnsemblGenome"]
#         .drop("DB", axis=1)
#         .rename(columns={"ID": "gene_id"})
#     )

#     gff_df = gff_df.merge(uniprot_df, how="left", on="gene_id")
#     gff_df = gff_df[~gff_df.Uniprot.isnull()]

#     gff_df = gff_df.rename(columns={"seqid": "chr"})
#     gff_df = gff_df[["Uniprot", "chr", "start", "end", "strand"]]

#     return gff_df


parser = argparse.ArgumentParser(description="Create dataset")
parser.add_argument("-i", "--input_file", help="Input gff file", required=True)
parser.add_argument("-m", "--mapping_file", help="Uniprot mapping file", required=True)
parser.add_argument(
    "-c",
    "--mapping_column",
    help="Column in uniprot mapping file that corresponds to id in genome file",
    required=True,
)
parser.add_argument("-o", "--output_file", help="Output file", required=True)

args = parser.parse_args()

input_file_name = args.input_file
output_file_name = args.output_file
mapping_file_name = args.mapping_file
mapping_column_name = args.mapping_column


gff_df = read_gff3(input_file_name)

gff_df = process_genome(gff_df, mapping_file_name, mapping_column_name)

gff_df.to_csv(output_file_name, index=None, sep="\t")
