import pandas as pd
from soft import read_soft

# from qnorm import quantile_normalize
from pathlib import Path
import argparse
from expression import aggregate_maxmean

parser = argparse.ArgumentParser(
    description="Create gene expression tsv from arrayexpress file."
)
parser.add_argument(
    "-e",
    "--expression-file",
    help="Input Expression file in E-TABM format",
    required=True,
)
parser.add_argument(
    "-g",
    "--gpl-file",
    help="GPL Soft file describing microarray platform",
    required=True,
)
parser.add_argument(
    "-m",
    "--mapping-file",
    help="Uniprot id mapping dat file for organism",
    required=True,
)
parser.add_argument(
    "-a", "--annotation-file", help="SDRF file with annotations", required=True,
)
parser.add_argument(
    "--ecotypes", help="ecotypes to filter for", type=str, nargs="*", default=None
)
parser.add_argument(
    "--organism-parts",
    help="organism parts to filter for",
    type=str,
    nargs="*",
    default=None,
)
parser.add_argument(
    "-o",
    "--output-file",
    help="Output tsv file for gene expression data",
    required=True,
)
parser.add_argument("--quantile-normalize", default=False, action="store_true")
parser.add_argument("--aggregate", default=False, action="store_true")

args = parser.parse_args()

###################################
# Reading expression data         #
###################################

df_exp = pd.read_table(args.expression_file, low_memory=False, index_col=0)
df_exp = df_exp[
    df_exp.columns[
        df_exp.loc["Composite Element REF"] == "RMA:mean_normalized_log_ratio"
    ]
].drop("Composite Element REF")
df_exp = df_exp.reset_index()

###################################
# Identifier conversion           #
###################################

# Affy probe id to Entrez gene id
gpl_data = read_soft(args.gpl_file)
gpl_df = pd.DataFrame.from_records(gpl_data[1:], columns=gpl_data[0])
gpl_df_entrez = gpl_df[["ID", "Gene ID"]]

df_exp = df_exp.merge(gpl_df_entrez, left_on="Scan REF", right_on="ID").drop(
    ["ID", "Scan REF"], axis=1
)
df_exp = df_exp[df_exp["Gene ID"] != ""]

df_exp["Gene ID"] = df_exp["Gene ID"].str.split("///")
df_exp = df_exp.explode("Gene ID")
df_exp = df_exp.astype({"Gene ID": int})

df_exp = df_exp.rename(columns={"Gene ID": "ID"})

# Entrez to Uniprot
df_idmapping = pd.read_table(
    args.mapping_file, header=None, names=["Uniprot", "DB", "ID"], low_memory=False,
)

df_idmapping_entrez = df_idmapping[df_idmapping.DB == "GeneID"].drop("DB", axis=1)
df_idmapping_entrez = df_idmapping_entrez.astype({"ID": int})

df_exp = df_exp.merge(df_idmapping_entrez, on="ID", how="left").drop("ID", axis=1)
df_exp = df_exp[~df_exp.Uniprot.isnull()].reset_index(drop=True)


###################################
# Rename columns                  #
###################################

df_exp.columns = df_exp.columns.str.replace(" ", "")
df_exp = df_exp.rename(columns=lambda c: "_".join(c.split("_")[:2]))
df_exp = df_exp.reindex(df_exp.columns.sort_values(), axis=1)

###################################
# Aggregation                     #
###################################

df_exp = df_exp.astype(
    {colname: float for colname in df_exp.columns if colname != "Uniprot"}
)
if args.aggregate:
    df_exp = aggregate_maxmean(df_exp, col_name="Uniprot")

###################################
# Filtering annotations           #
###################################

# Reading annotation file
df_annotation_raw = pd.read_table(args.annotation_file, index_col=0)
df_annotation = pd.DataFrame()

df_annotation["age_days"] = df_annotation_raw["Factor Value [time]"]
df_annotation["ecotype"] = df_annotation_raw["Factor Value [ecotype]"]
df_annotation["organism_part"] = df_annotation_raw["Factor Value [organism part]"]
df_annotation["developmental_stage"] = df_annotation_raw[
    "Characteristics [developmental stage]"
]
df_annotation["sample_name"] = df_annotation_raw["Scan Name"]

df_annotation.sample_name = df_annotation.sample_name.apply(
    lambda x: "_".join(x.split("_")[:-1]).replace(" ", "")
)

df_annotation = (
    df_annotation.drop_duplicates().set_index("sample_name", drop=True).sort_index()
)

df_annotation = df_annotation.loc[df_exp.columns]

# TODO filter
# Four datasets: columbia_all, columbia_four, columbia_


if args.ecotypes:
    df_annotation = df_annotation[df_annotation.ecotype.isin(args.ecotypes)]

if args.organism_parts:
    df_annotation = df_annotation[df_annotation.organism_part.isin(args.organism_parts)]

df_exp = df_exp[df_annotation.index]

###################################
# Writing output                  #
###################################

output_path = Path(args.output_file)
if not output_path.parent.exists():
    output_path.parent.mkdir(parents=True, exist_ok=True)

# if args.quantile_normalize:
#     df_exp = quantile_normalize(df_exp)
print(df_exp.shape)
# if median is not calculated, then Uniprot never becomes index, so don't write index.
df_exp.to_csv(output_path, sep="\t", index=args.aggregate)
