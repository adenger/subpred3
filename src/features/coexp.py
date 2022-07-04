import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from util.fasta import read_fasta
import logging
from joblib.parallel import Parallel, delayed

log = logging.getLogger(__name__)


class ExpressionFeature:
    def __init__(
        self,
        gene_expression_file_path: str,
        gene_pos_file_path: str,
        go_file_path: str,
        genome_neighbor_count: int,
        genome_selected_neighbor_count: int,
        go_percentage_threshold: float,
        # go_curated_terms_only: bool = False,
        feature_type: str,
        go_ontology: set,
        verbose: bool = False,
    ):
        """
        Input data columns (all organism-specific)
            GO file
                Index:None
                Uniprot
                qualifier
                go_id
                evidence_code
                ontology
            Gene positions file
                Index:None
                Uniprot
                chr
                start
                end
                strand
            gene expression file (with and without norm)
                Index: Uniprot; Accessions need to be unique
                <Sample names>
        """
        if feature_type not in ["binary", "percentage"]:
            raise ValueError(f"Invalid feature type: {feature_type}")

        self.genome_neighbor_count = genome_neighbor_count
        self.genome_selected_neighbor_count = genome_selected_neighbor_count
        self.go_percentage_threshold = go_percentage_threshold

        self.feature_type = feature_type

        """
        Read expression data
        """
        log.debug(
            f"reading gene expression dataset {Path(gene_expression_file_path).stem}"
        )
        self.expression = pd.read_table(gene_expression_file_path, index_col=0)
        log.debug(f"found {len(self.expression.columns)} samples.")
        if len(self.expression.index.unique()) != len(self.expression.index):
            log.warning("Warning: Expression data index is not unique!")

        self.tissue_name = Path(gene_expression_file_path).stem

        """
        GO Data
        """

        go_data = pd.read_table(go_file_path)
        go_data = go_data[go_data.ontology.isin(go_ontology)]
        go_data = go_data[["Uniprot", "go_id"]]
        self.accession_to_go_terms = (
            go_data.groupby("Uniprot")["go_id"].apply(set).to_dict()
        )

        """
        Read Gene Positions
        """
        gene_positions = pd.read_table(gene_pos_file_path)
        gene_positions = gene_positions[gene_positions.Uniprot.isin(go_data.Uniprot)]
        gene_positions = gene_positions[
            gene_positions.Uniprot.isin(self.expression.index)
        ]
        gene_positions = gene_positions[~gene_positions.Uniprot.duplicated()]

        def chr_to_neighbors(df):
            accession_to_neighbors = {}
            min_index = 0
            max_index = len(df.index) - 1
            for index in range(len(df.index)):
                index_shifted = max(index, min_index + genome_neighbor_count)
                index_shifted = min(index_shifted, max_index - genome_neighbor_count)
                start_index = index_shifted - genome_neighbor_count
                end_index = index_shifted + genome_neighbor_count
                accession = df.iloc[index].Uniprot
                neighbors = df.Uniprot[start_index : (end_index + 1)].to_list()
                neighbors = list(filter(lambda x: x != accession, neighbors))
                accession_to_neighbors[accession] = neighbors
            return accession_to_neighbors

        chr_to_accession_to_neighbors = gene_positions.groupby("chr").apply(
            chr_to_neighbors
        )

        # merge dicts from all chromosomes. each gene only has one chromosome
        self.accession_to_neighbors = dict()
        for accession_to_neighbors_chr in chr_to_accession_to_neighbors:
            self.accession_to_neighbors = {
                **self.accession_to_neighbors,
                **accession_to_neighbors_chr,
            }

        """
        Check which accessions have data
        """
        self.accessions_with_data = set(self.accession_to_neighbors.keys()) & set(
            self.accession_to_go_terms.keys() & set(self.expression.index.tolist())
        )
        log.debug(f"Feature is available for {len(self.accessions_with_data)} proteins")

        """
        Remove proteins for which not all data points are available
        """
        self.accession_to_neighbors = {
            accession: neighbors
            for accession, neighbors in self.accession_to_neighbors.items()
            if accession in self.accessions_with_data
        }
        self.accession_to_go_terms = {
            accession: go_terms
            for accession, go_terms in self.accession_to_go_terms.items()
            if accession in self.accessions_with_data
        }
        self.expression = self.expression[
            self.expression.index.isin(self.accessions_with_data)
        ]

        self.is_fitted = False

    # Returns neighbor_count neighbors to each side of the accession
    # or the next 2*gene_count if gene is at the end of the chromosome
    def __get_genome_neighbors(
        self,
        accession: str,
    ) -> list:
        return self.accession_to_neighbors[accession]

    def __get_go_terms(self, accessions: list) -> list:
        # single accession
        if type(accessions) is str:
            accessions = [accessions]
        go_terms = set()
        for accession in accessions:
            go_terms_accession = self.accession_to_go_terms[accession]
            go_terms |= go_terms_accession
        return go_terms

    def train_calculate(self, training_df: pd.DataFrame) -> pd.Series:
        """Calculates the feature vectors for proteins in training_df by
        iteratively removing a protein from the dataset, training ExpressionFeature with the
        remaining proteins, and calculating the feature vector. Finally, train() is called
        with all proteins in training_df.

        calculate() can then be called on the proteins in the test set, without any dependency.

        Args:
            training_df (pd.DataFrame): Dataframe containing Uniprot accessions (index)
                and transported substrates (Substrates column)
            tissue (str): The GTEx tissue for which to calculate the feature

        Returns:
            pd.Series(np.array): Series with Uniprot accession as index,
                feature vectors (np.array) as elements
        """
        # calculates feature for every accession in training_df iteratively by calling train and calculate
        # then trains class on entire training_df
        training_df = training_df[training_df.index.isin(self.accessions_with_data)]
        log.debug(
            f"Feature is available for {training_df.shape[0]} proteins from training dataset"
        )
        log.debug(training_df.value_counts())

        feature_vector_series_list = list()
        for accession in training_df.index:
            training_df_filtered = training_df[training_df.index != accession]
            self.train(training_df=training_df_filtered)
            feature_vector_series = self.calculate(accession=accession)
            feature_vector_series_list.append(feature_vector_series)

        feature_vectors_df = pd.DataFrame(
            data=feature_vector_series_list, index=training_df.index
        )
        # self.train(training_df=training_df)
        return feature_vectors_df

    def train(self, training_df: pd.DataFrame):
        """
        Calculates the GO term profiles for the substrate classes in training_df.
        This is done independently of the calculation, to avoid sharing information
        between training and test sets.

        Args:
            training_df (pd.DataFrame): The dataframe used for training.
                Index: Uniprot Accessions of transporters
                Substrate: The transported substrates
        """
        training_df = training_df[training_df.index.isin(self.accessions_with_data)]

        training_df = training_df.assign(
            neighbors=training_df.index.map(
                lambda accession: self.__get_genome_neighbors(accession),
            )
        )

        training_df = training_df.assign(
            GO=training_df.index.map(lambda accession: self.__get_go_terms(accession))
        )

        go_profiles_transporters = dict()
        go_profiles_neighbors = dict()
        for substrate in training_df.Substrate.unique():
            go_terms_substrate = set()
            for go_terms_set in training_df[training_df.Substrate == substrate].GO:
                go_terms_substrate |= go_terms_set
            go_profiles_transporters[substrate] = go_terms_substrate

            neighbors_substrate = set()
            for neighbors in training_df[training_df.Substrate == substrate].neighbors:
                neighbors_substrate |= set(neighbors)
            go_terms_neighbors_substrate = self.__get_go_terms(neighbors_substrate)
            go_profiles_neighbors[substrate] = go_terms_neighbors_substrate

        self.is_fitted = True
        self.go_profiles_transporters = go_profiles_transporters
        self.go_profiles_neighbors = go_profiles_neighbors
        self.training_data_accessions = training_df.index.tolist()

    def calculate(self, accession: str) -> np.array:
        """Calculates the COEXP feature for a protein, using the expression data of "tissue".
        The train() function has to be called first.

        Args:
            accession (str): The transport protein
            tissue (str): The name of a tissue in the GTEx gene expression data.
                Available options can be found in self.tissue_names.
            feature_type (string): Feature processing. Options: "binary" or "percentage".
                Can be original approach from Tran2018 (convert to binary with cutoff)
                or percentages.

        Raises:
            RuntimeError: train() was not called first
            ValueError: Accession was part of the training data (i.e. substrate is known)
            ValueError: Invalid name for tissue
            ValueError: Not enough data found to calculate feature for protein
            ValueError: Feature type is not "binary" or "percentage"

        Returns:
            np.array(float): Feature vector for accession
        """
        if not self.is_fitted:
            raise RuntimeError("Call train function first.")

        if accession in self.training_data_accessions:
            raise ValueError(f"Accession {accession} was in training dataset!")

        # if not tissue in self.tissue_names:
        #     raise ValueError(f"Tissue {tissue} was not found in the expression data")

        if accession not in self.accessions_with_data:
            # Accessions with data are stored in self.accessions_with_data. filter beforehand
            raise ValueError(f"No data for accession {accession}")

        # get neighbors of current accession
        accession_neighbors = self.__get_genome_neighbors(accession)

        # get expression data of neighbors and accession.
        expression_data = self.expression.loc[accession_neighbors + [accession]]

        # calculate pairwise correlation between neighbors and accession, sort by highest
        corr = (
            expression_data.loc[accession_neighbors]
            .corrwith(expression_data.loc[accession], axis=1, method="pearson")
            .sort_values(ascending=False)
        )

        # get neighbors with highest correlation
        selected_neighbors = corr.iloc[
            0 : self.genome_selected_neighbor_count
        ].index.tolist()

        # get go terms of selected neighbors
        selected_neighbors_go_terms = {
            selected_neighbor: self.__get_go_terms(selected_neighbor)
            for selected_neighbor in selected_neighbors
        }

        feature_vector = list()
        feature_names = list()
        # sorting to make it deterministic
        for substrate in sorted(self.go_profiles_transporters.keys()):
            go_profile_transporters = self.go_profiles_transporters[substrate]
            go_profile_neighbors = self.go_profiles_neighbors[substrate]
            for selected_neighbor_pos, selected_neighbor in enumerate(
                sorted(selected_neighbors_go_terms.keys())
            ):
                selected_neighbor_go_terms = selected_neighbors_go_terms[
                    selected_neighbor
                ]
                percentage_transporters = len(
                    selected_neighbor_go_terms & go_profile_transporters
                ) / len(selected_neighbor_go_terms)

                if self.feature_type == "binary":
                    feature_vector.append(
                        1
                        if percentage_transporters >= self.go_percentage_threshold
                        else 0
                    )
                else:  # percentage
                    feature_vector.append(percentage_transporters)

                feature_names.append(
                    f"neighbor{selected_neighbor_pos}_{substrate}_transporters"
                )

                percentage_transporter_neighbors = len(
                    selected_neighbor_go_terms & go_profile_neighbors
                ) / len(selected_neighbor_go_terms)

                if self.feature_type == "binary":
                    feature_vector.append(
                        1
                        if percentage_transporter_neighbors
                        >= self.go_percentage_threshold
                        else 0
                    )
                else:  # percentage
                    feature_vector.append(percentage_transporter_neighbors)

                feature_names.append(
                    f"neighbor{selected_neighbor_pos}_{substrate}_neighbors"
                )

        return pd.Series(data=feature_vector, index=feature_names)


def __write_coexp_df(
    gene_expression_file_path: Path,
    gene_pos_file: str,
    go_file: str,
    output_folder: str,
    neighbor_count: int,
    selected_neighbor_count: int,
    go_percentage_threshold: float,
    feature_type: str,
    ontologies: list,
    substrate_df: pd.DataFrame,
):
    log.debug("#" * 80)

    log.debug(f"gene expression file {gene_expression_file_path}")

    coexp_calculator = ExpressionFeature(
        gene_expression_file_path=gene_expression_file_path,
        gene_pos_file_path=gene_pos_file,
        go_file_path=go_file,
        genome_neighbor_count=neighbor_count,
        genome_selected_neighbor_count=selected_neighbor_count,
        go_percentage_threshold=go_percentage_threshold,
        feature_type=feature_type,
        go_ontology=set(ontologies),
    )
    results_df = coexp_calculator.train_calculate(training_df=substrate_df)

    results_filename = "{}_n{}_s{}_p{}_{}_{}".format(
        str(gene_expression_file_path.name).replace(".tsv", ""),
        neighbor_count,
        selected_neighbor_count,
        go_percentage_threshold,
        feature_type,
        "".join(sorted(list(set(ontologies)))),
    )

    output_path = Path(output_folder)

    log.debug(f"writing results to {output_path}")

    if not output_path.exists():
        log.warning(f"Output folder {output_path} did not exist, creating folder.")
        output_path.mkdir(parents=True, exist_ok=True)

    if feature_type == "binary":
        results_df = results_df.astype(int)
    elif feature_type == "percentage":
        results_df = results_df.astype(float)

    results_df.to_csv(f"{output_folder}/{results_filename}.tsv", sep="\t")


# TODO tissue sample counts beforehand?
# TODO one chromosome per gene (gene positions)
# TODO option to train on whole DF
# TODO take and write individual files


def calculate_coexp_feature(
    gene_expression_folder: str,
    fasta_file_training: str,
    log_file: str,
    gene_pos_file: str,
    go_file: str,
    output_folder: str,
    neighbor_count: int,
    selected_neighbor_count: int,
    go_percentage_threshold: float,
    feature_type: str,
    ontologies: list = ["P", "F", "C"],
    verbose: bool = False,
    n_threads: int = 4,
):
    logging.basicConfig(
        format="[%(levelname)s] %(filename)s: %(message)s",
        level=logging.DEBUG if verbose else logging.INFO,
        filename=log_file,
    )

    fasta_data_list = read_fasta(fasta_file_training)
    fasta_substrate_records = []
    for header, _ in fasta_data_list:
        header_elements = header.split("|")
        fasta_substrate_records.append(
            {
                "Uniprot": header_elements[1],
                "Substrate": header_elements[5],
            }
        )
    substrate_df = pd.DataFrame.from_records(fasta_substrate_records)
    substrate_df = substrate_df.set_index("Uniprot", verify_integrity=True)
    log.debug(f"Read {substrate_df.shape[0]} proteins from fasta file")
    log.debug(substrate_df.value_counts())
    substrate_df = substrate_df[~substrate_df.Substrate.str.contains(";")]
    log.debug(substrate_df.value_counts())

    Parallel(n_jobs=n_threads)(
        delayed(__write_coexp_df)(
            gene_expression_file_path=gene_expression_file_path,
            gene_pos_file=gene_pos_file,
            go_file=go_file,
            output_folder=output_folder,
            neighbor_count=neighbor_count,
            selected_neighbor_count=selected_neighbor_count,
            go_percentage_threshold=go_percentage_threshold,
            feature_type=feature_type,
            ontologies=ontologies,
            substrate_df=substrate_df,
        )
        for gene_expression_file_path in Path(gene_expression_folder).glob("*.tsv")
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create coexp feature")
    """
    Files
    """
    parser.add_argument(
        "-e",
        "--gene-expression-folder",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--fasta-file-training",
        help="FASTA file with protein training data",
        required=True,
    )
    parser.add_argument("--log-file", help="log file to write output to", default=None)
    parser.add_argument(
        "-p",
        "--gene-pos-file",
        required=True,
    )
    parser.add_argument("-g", "--go-file", required=True)
    parser.add_argument("-o", "--output-folder", help="Output Folder", required=True)
    """
    Coexp parameters
    """
    parser.add_argument("--neighbor-count", default=5, type=int)
    parser.add_argument("--selected-neighbor-count", default=3, type=int)
    parser.add_argument("--go-percentage-threshold", default=0.8, type=float)
    parser.add_argument(
        "--feature-type", choices=["percentage", "binary"], default="binary", type=str
    )
    parser.add_argument(
        "--ontologies", choices=["P", "F", "C"], nargs="+", default=["P", "F", "C"]
    )
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--n-threads", type=int, default=-1)

    args = parser.parse_args()

    calculate_coexp_feature(
        gene_expression_folder=args.gene_expression_folder,
        fasta_file_training=args.fasta_file_training,
        log_file=args.log_file,
        gene_pos_file=args.gene_pos_file,
        go_file=args.go_file,
        output_folder=args.output_folder,
        neighbor_count=args.neighbor_count,
        selected_neighbor_count=args.selected_neighbor_count,
        go_percentage_threshold=args.go_percentage_threshold,
        feature_type=args.feature_type,
        ontologies=args.ontologies,
        verbose=args.verbose,
        n_threads=args.n_threads,
    )
