from collections import OrderedDict
import numpy as np
import pandas as pd


def __get_amino_acids():
    return list("ACDEFGHIKLMNPQRSTVWY")


def __amino_acid_comp(sequence: str, alphabet: str = "ACDEFGHIKLMNPQRSTVWY"):
    counter = OrderedDict({aa: 0 for aa in alphabet})
    for c in sequence:
        counter[c] += 1
    return np.divide(np.array(list(counter.values())), len(sequence))


def calculate_aac(sequences: pd.Series):

    return pd.DataFrame(
        data=sequences.apply(__amino_acid_comp).tolist(),
        index=sequences.index,
        columns=__get_amino_acids(),
    )


def __dipeptide_comp(sequence: str, alphabet: str = "ACDEFGHIKLMNPQRSTVWY"):
    counter = OrderedDict({aa1 + aa2: 0 for aa1 in alphabet for aa2 in alphabet})
    for i in range(len(sequence) - 1):
        peptide = sequence[i : i + 2]
        counter[peptide] += 1
    return np.array([b / (len(sequence) - 1) for _, b in sorted(counter.items())])


def calculate_paac(sequences: pd.Series):
    return pd.DataFrame(
        data=sequences.apply(__dipeptide_comp).tolist(),
        index=sequences.index,
        columns=[
            aa1 + aa2 for aa1 in __get_amino_acids() for aa2 in __get_amino_acids()
        ],
    )


# def calculate_composition_feature(input_fasta: str, output_tsv: str, feature_type: str):
#     fasta_data = read_fasta(input_fasta)
#     fasta_data_records = []
#     for header, sequence in fasta_data:
#         fasta_data_records.append(
#             {"Uniprot": header.split("|")[1], "sequence": sequence,}
#         )

#     sequence_df = pd.DataFrame.from_records(fasta_data_records).set_index("Uniprot")

#     if feature_type == "aac":
#         sequence_df = sequence_df.sequence.apply(__amino_acid_comp)
#         sequence_df = pd.DataFrame(
#             sequence_df.to_list(), index=sequence_df.index, columns=__get_amino_acids()
#         )
#     elif feature_type == "paac":
#         sequence_df = sequence_df.sequence.apply(__dipeptide_comp)
#         sequence_df = pd.DataFrame(
#             sequence_df.to_list(),
#             index=sequence_df.index,
#             columns=[
#                 aa1 + aa2 for aa1 in __get_amino_acids() for aa2 in __get_amino_acids()
#             ],
#         )
#     else:
#         raise ValueError(f"Invalid feature type: {feature_type}")

#     sequence_df.to_csv(output_tsv, sep="\t")


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(
#         description="Create amino acid composition feature dataframe from protein sequences in fasta file"
#     )

#     parser.add_argument(
#         "-i",
#         "--input-file",
#         help="FASTA file with Uniprot accession and Sequence",
#         required=True,
#     )

#     parser.add_argument(
#         "-o", "--output-file", required=True,
#     )

#     parser.add_argument("-t", "--type", choices=["aac", "paac"], default="aac")

#     args = parser.parse_args()

#     calculate_composition_feature(args.input_file, args.output_file, args.type)
