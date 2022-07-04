import argparse
import pandas as pd
from util.fasta import read_fasta


def fasta_to_labels(input_fasta: str, output_tsv: str, tax_id: int = None):

    fasta_data = read_fasta(input_fasta)

    df = pd.DataFrame.from_records(
        [
            (header.split("|")[1], int(header.split("|")[3]), header.split("|")[5])
            for header, _ in fasta_data
        ],
        columns=["Uniprot", "tax_id", "labels"],
    ).set_index("Uniprot", verify_integrity=True)

    if tax_id:
        df = df[df.tax_id == tax_id]
    df = df.drop("tax_id", axis=1)

    df.to_csv(output_tsv, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="extract labels from fasta files")

    parser.add_argument(
        "-i",
        "--input-fasta",
        required=True,
    )

    parser.add_argument(
        "-o",
        "--output-tsv",
        required=True,
    )

    parser.add_argument("-t", "--tax-id", type=int)

    args = parser.parse_args()

    fasta_to_labels(args.input_fasta, args.output_tsv, args.tax_id)
