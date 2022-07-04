import pandas as pd
import argparse
from util.fasta import read_fasta, write_fasta
import os
from sklearn.preprocessing import minmax_scale
import subprocess
import platform


def __process_pssm_file(pssm_file_name, tmp_folder_name):
    # print(pssm_file_name, tmp_folder_name)
    with open(pssm_file_name) as pssm_file:
        next(pssm_file)
        next(pssm_file)

        amino_acids = pssm_file.readline().strip().split()[:20]
        amino_acid_to_sum_vector = {
            amino_acid: [0.0] * 20 for amino_acid in amino_acids
        }

        for line in pssm_file:
            if line == "\n":  # end of file, before overall scores
                break

            values = line.strip().split()
            amino_acid = values[1]
            if amino_acid not in amino_acids:
                raise ValueError(
                    f"unexpected amino acid in pssm file {pssm_file_name}: {amino_acid}"
                )

            scores = [float(score) for score in values[2:22]]
            sum_vector = amino_acid_to_sum_vector.get(amino_acid)

            if len(scores) != 20:
                raise ValueError(
                    f"incomplete PSSM file: {pssm_file_name}. Delete from folder {tmp_folder_name} and recompute"
                )

            for pos in range(20):
                sum_vector[pos] += scores[pos]
            amino_acid_to_sum_vector[amino_acid] = sum_vector

        pssm = []
        for _, sum_vector in sorted(amino_acid_to_sum_vector.items()):
            pssm.extend(sum_vector)

        pssm = minmax_scale(pssm).tolist()  # scale to [0,1]

        return pssm


def __create_pssm_file(
    psiblast_location: str,
    fasta_file_name: str,
    pssm_file_name: str,
    blastdb_location: str,
    iterations: int,
    evalue: float,
    threads: int,
) -> None:
    log_file_name = f"{pssm_file_name}.log"
    subprocess.run(
        "{} -query {} -db {} -num_iterations {} -inclusion_ethresh {} -num_threads {} -save_pssm_after_last_round\
             -out_ascii_pssm {} -out {} -comp_based_stats {}".format(
            psiblast_location,
            fasta_file_name,
            blastdb_location,
            iterations,
            evalue,
            threads,
            pssm_file_name,
            log_file_name,
            2
            if iterations == 1
            else 1,  # default is 2, but not supported when matrix is PSSM instead of BLOSUM
        ),
        check=True,
        shell=False if platform.system() == "Windows" else True,
        # stdout=subprocess.DEVNULL,
        # stderr=subprocess.DEVNULL,
    )


def __get_pssm_feature(
    accession: str,
    sequence: str,
    blastdb_fasta_file: str,
    pssm_folder_path: str,
    psiblast_location: str,
    iterations: int = 1,
    evalue: float = 0.002,
    threads: int = 1,
    verbose: bool = False,
) -> list:

    if not os.path.exists(pssm_folder_path):
        os.makedirs(pssm_folder_path)

    fasta_file_name = f"{pssm_folder_path}/{accession}.fasta"
    pssm_file_name = f"{pssm_folder_path}/{accession}.pssm"

    pssm = []
    if os.path.isfile(pssm_file_name):
        pssm = __process_pssm_file(pssm_file_name, pssm_folder_path)
        if verbose:
            print(f"PSSM for accession {accession} was found in tmp folder")
    else:
        if verbose:
            print(
                f"PSSM for accession {accession} was not found in tmp folder, calling psiblast"
            )
        write_fasta(
            fasta_file_name=fasta_file_name, fasta_data=[(">" + accession, sequence)]
        )
        __create_pssm_file(
            psiblast_location=psiblast_location,
            fasta_file_name=fasta_file_name,
            pssm_file_name=pssm_file_name,
            blastdb_location=blastdb_fasta_file,
            iterations=iterations,
            evalue=evalue,
            threads=threads,
        )
        pssm = __process_pssm_file(pssm_file_name, pssm_folder_path)
        if verbose:
            print(f"PSSM for accession {accession} was generated")

    return pssm


def calculate_pssm_feature(
    input_fasta: str,
    output_tsv: str,
    tmp_folder: str,
    blast_db: str,
    iterations: int,
    psiblast_executable: str,
    psiblast_threads: int = 4,
    verbose: bool = False,
):

    fasta_data = read_fasta(input_fasta)
    fasta_data_records = []
    for header, sequence in fasta_data:
        fasta_data_records.append(
            {
                "Uniprot": header.split("|")[1],
                "sequence": sequence,
            }
        )
    sequence_df = pd.DataFrame.from_records(fasta_data_records).set_index("Uniprot")

    features = [
        __get_pssm_feature(
            accession=sequence_df.index[i],
            sequence=sequence_df.sequence[i],
            blastdb_fasta_file=blast_db,
            pssm_folder_path=tmp_folder,
            iterations=iterations,
            psiblast_location=psiblast_executable,
            threads=psiblast_threads,
            verbose=verbose,
        )
        for i in range(len(sequence_df))
    ]

    pssm_aa_order = "ARNDCQEGHILKMFPSTWYV"
    pssm_aa_substitutions = [
        aa1 + aa2 for aa1 in pssm_aa_order for aa2 in pssm_aa_order
    ]

    pssm_df = pd.DataFrame(
        data=features, index=sequence_df.index, columns=pssm_aa_substitutions
    )

    pssm_df.to_csv(output_tsv, sep="\t")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Create pssm composition feature dataframe from protein sequences in fasta file"
    )
    parser.add_argument(
        "-i",
        "--input-file",
        help="FASTA file with Uniprot accession and Sequence",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-file",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--tmp-folder",
        required=True,
    )
    parser.add_argument(
        "-db",
        "--blast-db",
        help="BLAST db fasta file",
        required=True,
    )
    parser.add_argument(
        "-it",
        "--iterations",
        required=True,
    )
    parser.add_argument(
        "-e",
        "--psiblast-executable",
        required=True,
    )
    parser.add_argument("--psiblast-threads", default=4)
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    calculate_pssm_feature(
        input_fasta=args.input_file,
        output_tsv=args.output_file,
        tmp_folder=args.tmp_folder,
        blast_db=args.blast_db,
        iterations=args.iterations,
        psiblast_executable=args.psiblast_executable,
        psiblast_threads=args.psiblast_threads,
        verbose=args.verbose,
    )
