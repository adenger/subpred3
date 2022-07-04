import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Create dataset")
parser.add_argument("-i", "--input_file", help="Input file", required=True)
parser.add_argument("-o", "--output_file", help="Output file", required=True)

args = parser.parse_args()

input_file_name = args.input_file
output_file_name = args.output_file

# http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/
goa_df = pd.read_table(
    input_file_name,
    comment="!",
    header=None,
    low_memory=True,
    engine="c",
    dtype="str",
    names=[
        "DB",
        "DB_Object_ID",
        "DB_Object_Symbol",
        "Qualifier",
        "GO_ID",
        "DB:Reference",
        "Evidence_Code",
        "With_or_From",
        "Aspect",
        "DB_Object_Name",
        "DB_Object_Synonym",
        "DB_Object_Type",
        "Taxon_and_Interacting_taxon",
        "Date",
        "Assigned_By",
        "Annotation_Extension",
        "Gene_Product_Form_ID",
    ],
    usecols=[
        # "DB",
        "DB_Object_ID",
        # "DB_Object_Symbol",
        "Qualifier",
        "GO_ID",
        # "DB:Reference",
        "Evidence_Code",
        # "With_or_From",
        "Aspect",
        # "DB_Object_Name",
        # "DB_Object_Synonym",
        # "DB_Object_Type",
        # "Taxon_and_Interacting_taxon",
        # "Date",
        # "Assigned_By",
        # "Annotation_Extension",
        # "Gene_Product_Form_ID",
    ],
)
# goa_df = goa_df[goa_df.db == "UniProtKB"]

# goa_df = goa_df[~goa_df.taxon.str.contains("\|")]
# goa_df = goa_df.assign(
#     tax_id = goa_df.taxon.str.split(":", expand=True)[1]
# )

goa_df = goa_df.rename(
    columns={
        "DB_Object_ID": "Uniprot",
        "Qualifier": "qualifier",
        "GO_ID": "go_id",
        "Evidence_Code": "evidence_code",
        "Aspect": "ontology"
    }
)

goa_df = goa_df[~goa_df.qualifier.str.contains("NOT")]
# goa_df = goa_df[~goa_df.qualifier.str.contains("upstream")]
goa_df = goa_df[~goa_df.ontology.isnull()]

goa_df = goa_df.drop_duplicates()

goa_df.to_csv(output_file_name, sep="\t", index=False)
