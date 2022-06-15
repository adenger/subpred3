from scipy.stats import pearsonr
import pandas as pd
import numpy as np


class Coexp:
    #################
    # Preprocessing #
    #################

    def __init__(
        self,
        expression_tsv: str,
        go_tsv: str,
        gene_pos_tsv: str,
        swissprot_tsv: str,
        chromosomes: list,
        tax_id: int,
        n_neighbors: int,
    ):
        df_exp, df_go_long, df_gene_pos_long = self.__read_tables(
            expression_tsv, go_tsv, gene_pos_tsv, chromosomes
        )

        proteins_swissprot = self.__get_swissprot_proteins(
            swissprot_tsv=swissprot_tsv, tax_id=tax_id
        )

        (
            df_exp,
            df_go_long,
            df_gene_pos_long,
            proteins_enough_data,
        ) = self.__filter_datasets(
            df_exp, df_go_long, df_gene_pos_long, proteins_swissprot
        )
        df_neighbors_long = self.__get_neighbors(
            df_gene_pos_long=df_gene_pos_long, n_neighbors=n_neighbors
        )

        self.__proteins_whitelist = proteins_enough_data
        self.__df_exp = df_exp
        self.__df_go_long = df_go_long
        self.__df_neighbors = df_neighbors_long

    def __read_tables(
        self, expression_tsv: str, go_tsv: str, gene_pos_tsv: str, chromosomes: list
    ):
        df_exp = pd.read_table(expression_tsv, index_col=0)
        df_go_long = pd.read_table(go_tsv)
        df_gene_pos_long = pd.read_table(gene_pos_tsv).sort_values(
            ["chr", "start", "strand"]
        )
        df_gene_pos_long = df_gene_pos_long[df_gene_pos_long.chr.isin(chromosomes)]
        return df_exp, df_go_long, df_gene_pos_long

    def __get_swissprot_proteins(self, swissprot_tsv: str, tax_id: int):
        df_swissprot = pd.read_table(
            swissprot_tsv,
            index_col=0,
            usecols=["Entry", "Organism ID", "Protein existence", "Fragment"],
        )
        df_swissprot_filtered = df_swissprot.copy()
        df_swissprot_filtered = df_swissprot_filtered[
            (df_swissprot_filtered["Organism ID"] == tax_id)
        ]
        df_swissprot_filtered = df_swissprot_filtered[
            ~df_swissprot_filtered["Protein existence"].isin(["Predicted", "Uncertain"])
        ]
        df_swissprot_filtered = df_swissprot_filtered[
            df_swissprot_filtered.Fragment.isnull()
        ]
        return df_swissprot_filtered.index.values

    def __filter_datasets(
        self,
        df_exp: pd.DataFrame,
        df_go_long: pd.DataFrame,
        df_gene_pos_long: pd.DataFrame,
        proteins_swissprot,
    ):
        proteins_enough_data = (
            set(proteins_swissprot)
            & set(df_exp.index)
            & set(df_go_long.Uniprot)
            & set(df_gene_pos_long.Uniprot)
        )
        df_exp = df_exp[df_exp.index.isin(proteins_enough_data)]
        df_gene_pos_long = df_gene_pos_long[
            df_gene_pos_long.Uniprot.isin(proteins_enough_data)
        ].reset_index(drop=True)
        df_go_long = df_go_long[
            df_go_long.Uniprot.isin(proteins_enough_data)
        ].reset_index(drop=True)
        return df_exp, df_go_long, df_gene_pos_long, proteins_enough_data

    def __get_neighbors(
        self, df_gene_pos_long: pd.DataFrame, n_neighbors: int
    ) -> pd.DataFrame:
        records_neighbors = []

        for chromosome in df_gene_pos_long.chr.unique():
            df_chr = df_gene_pos_long[df_gene_pos_long.chr == chromosome].reset_index(
                drop=True
            )
            for idx in df_chr.index:
                accession = df_chr.at[idx, "Uniprot"]
                idx_shifted = max(idx, n_neighbors)
                idx_shifted = min(idx_shifted, df_chr.shape[0] - 1 - n_neighbors)
                neighbors = df_chr.iloc[
                    idx_shifted - n_neighbors : idx_shifted + 1 + n_neighbors
                ].Uniprot
                neighbors = neighbors[neighbors != accession]
                for neighbor in neighbors:
                    records_neighbors.append([accession, neighbor, chromosome])
        df_neighbors_long = pd.DataFrame(
            records_neighbors, columns=["Uniprot", "Neighbor", "Chromosome"]
        )
        return df_neighbors_long

    #################
    # Calculation   #
    #################

    def __get_go_profiles(self, df_training: pd.DataFrame) -> dict:
        go_profiles = {}
        for name, df_substrate in df_training.groupby("Label"):
            label = name.replace(" ", "_")
            transporters = df_substrate.Uniprot.unique()
            neighbors = df_substrate.Neighbor.unique()

            go_profiles[label + "_transporters"] = self.__df_go_long[
                self.__df_go_long.Uniprot.isin(transporters)
            ].go_id.unique()
            go_profiles[label + "_neighbors"] = self.__df_go_long[
                self.__df_go_long.Uniprot.isin(neighbors)
            ].go_id.unique()
        return go_profiles

    def get_selected_neighbors(self, accession: str, n_selected: int) -> pd.DataFrame:
        accession_selected_neighbors = self.__df_neighbors[
            self.__df_neighbors.Uniprot == accession
        ]
        accession_selected_neighbors = accession_selected_neighbors.drop_duplicates()
        accession_selected_neighbors = accession_selected_neighbors.assign(
            coexp=accession_selected_neighbors.apply(
                lambda row: pearsonr(
                    self.__df_exp.loc[row.Uniprot], self.__df_exp.loc[row.Neighbor]
                )[0],
                axis=1,
            )
        )

        accession_selected_neighbors = accession_selected_neighbors.sort_values(
            "coexp", ascending=False
        ).iloc[:n_selected]

        accession_selected_neighbors_go = (
            accession_selected_neighbors[["Neighbor", "coexp"]]
            .merge(
                self.__df_go_long, how="left", left_on="Neighbor", right_on="Uniprot"
            )
            .drop(["Uniprot", "coexp"], axis=1)
            .drop_duplicates()
        )

        selected_neighbors_wide = (
            accession_selected_neighbors_go.groupby("Neighbor").apply(
                lambda df_neighbor: df_neighbor.go_id.unique()
            )
            # .sort_index()
        )
        return selected_neighbors_wide

    def __get_percentages(
        self, go_profiles: dict, selected_neighbors_go: pd.DataFrame
    ) -> pd.DataFrame:
        records = []
        for label, profile in sorted(go_profiles.items()):
            i = 0
            for neighbor, neighbor_go in selected_neighbors_go.iteritems():
                percentage = (
                    neighbor_go[np.isin(neighbor_go, profile)].size / neighbor_go.size
                )
                records.append([f"{label}_{i}", neighbor, percentage])
                i += 1
        df_percentages = pd.DataFrame.from_records(
            records, columns=["profile", "neighbor", "percentage"]
        )
        return df_percentages

    def __check_inputs(
        self, accession, training_accessions, training_labels, n_selected
    ):
        assert training_accessions.size == training_labels.size
        assert accession not in training_accessions
        assert isinstance(accession, str)
        assert isinstance(training_accessions, np.ndarray)
        assert isinstance(training_labels, np.ndarray)
        assert isinstance(n_selected, int)

    def __percentage_to_aggregate(self, df_percentages: pd.DataFrame) -> pd.Series:
        df_percentages.profile = (
            df_percentages.profile.str.split("_")
            .str[:3]
            .transform(lambda x: "_".join(x))
        )
        df_percentages = df_percentages.groupby("profile").percentage.max()
        df_percentages = df_percentages.reset_index()
        return df_percentages

    def get_feature(
        self,
        accession: str,
        training_accessions: np.ndarray,
        training_labels: np.ndarray,
        n_selected: int,
        threshold: float,
        aggregate: bool = False,
        binary: bool = True,
    ):
        self.__check_inputs(accession, training_accessions, training_labels, n_selected)
        if accession not in self.__proteins_whitelist:
            return None
        df_training = pd.DataFrame(training_accessions, columns=["Uniprot"])
        df_training["Label"] = training_labels

        df_training = df_training.merge(self.__df_neighbors, on="Uniprot", how="left")
        df_training = df_training[~df_training.Neighbor.isnull()].reset_index(drop=True)

        go_profiles = self.__get_go_profiles(df_training)

        selected_neighbors_go = self.get_selected_neighbors(accession, n_selected)
        df_percentages = self.__get_percentages(
            go_profiles=go_profiles, selected_neighbors_go=selected_neighbors_go
        )

        if aggregate:
            df_percentages = self.__percentage_to_aggregate(df_percentages)
        if binary:
            df_percentages["percentage"] = np.where(
                df_percentages.percentage > threshold, 1, 0
            )

        df_percentages = df_percentages.set_index("profile", verify_integrity=True)
        result = df_percentages.percentage.rename(accession)

        return result

    def get_features(
        self,
        training_accessions: np.ndarray,
        training_labels: np.ndarray,
        n_selected: int,
        threshold: float,
        aggregate: bool = False,
        binary: bool = True,
    ):
        res = []
        for accession in training_accessions:
            mask = training_accessions != accession
            accessions_subset = training_accessions[mask]
            labels_subset = training_labels[mask]

            feature = self.get_feature(
                accession=accession,
                training_accessions=accessions_subset,
                training_labels=labels_subset,
                n_selected=n_selected,
                threshold=threshold,
                aggregate=aggregate,
                binary=binary,
            )
            if isinstance(feature, pd.Series):
                res.append(feature)
        res = pd.concat(res, axis=1).transpose()
        return res
