import os
import pandas as pd
import argparse
import logging
from typing import List

from main.fasta import write_fasta
from main.cdhit import cd_hit


def create_dataset(
    keywords_substrate_filter: List[str],
    keywords_component_filter: List[str],
    keywords_transport_filter: List[str],
    input_file: str,
    multi_substrate: str = "keep",
    outliers: List[str] = None,
    verbose: bool = False,
    tax_ids_filter: List[int] = None,
    output_tsv: str = None,
    output_fasta: str = None,
    output_log: str = None,
    sequence_clustering: int = None,
):

    # pd.set_option("display.max_rows", 500)
    # pd.set_option("display.max_columns", 500)
    # pd.set_option("display.width", 1000)

    #########################
    # Logging               #
    #########################
    if os.path.exists(output_log):
        with open(output_log, "w"):
            pass

    log = logging.getLogger("DATASET")
    logging.basicConfig(
        format="[%(levelname)s] %(filename)s: %(message)s",
        level=logging.DEBUG if verbose else logging.INFO,
        filename=output_log,
    )

    #########################
    # Reading               #
    #########################

    log.debug("#" * 60)

    log.debug(
        "Parameters:\n{}".format(
            pd.DataFrame.from_dict(
                {k: str(v) for k, v in locals().items() if k != "log"}, orient="index"
            )
        )
    )

    log.debug(f"Reading file {input_file}")

    df = pd.read_table(input_file, index_col=0)

    log.debug(f"Found {df.shape[0]} rows and {df.shape[1]} columns")

    #########################
    # Perform basic cleanup #
    #########################

    log.debug("=" * 60)
    df = df.rename(
        columns={
            "Keyword ID": "keyword_ids",
            "Gene ontology IDs": "go_ids",
            "Gene ontology (GO)": "go_terms",
            "Cross-reference (TCDB)": "tcdb_id",
            "Taxonomic lineage IDs": "tax_id",
        },
    )
    df.columns = df.columns.map(lambda c: c.lower().replace(" ", "_"))
    df.index = df.index.rename("Uniprot")
    df.tcdb_id = df.tcdb_id.str.replace(";", "").str.strip()
    df.sequence = df.sequence.str.replace("X", "").str.replace("U", "")

    log.debug(f"Removing {df[df.keywords.isnull()].shape[0]} rows without any keywords")
    df = df[~df.keywords.isnull()]

    # Mostly peptides, apparently. Like Pollen
    log.debug(f"Removing {df.gene_names.isnull().sum()} proteins without gene symbols")
    df = df[~df.isnull()]

    #######################################
    # Removing proteins invalid sequences #
    #######################################

    log.debug("=" * 60)
    log.debug("Filtering invalid protein sequences:")
    log.debug("-" * 60)
    log.debug(
        f"Protein existence value counts: \n{df.protein_existence.value_counts(dropna=False).to_string()}"
    )
    log.debug(
        f"Sequence fragment value counts: \n{df.fragment.value_counts(dropna=False).to_string()}"
    )
    df = df[
        df.protein_existence.isin(
            {"Evidence at protein level", "Evidence at transcript level"}
        )
    ]
    df = df[df.fragment.isnull()]
    log.debug(f"Removed proteins with invalid sequences, {df.shape[0]} left")
    df = df.drop(["protein_existence", "fragment"], axis=1)

    ######################
    # Taxonomy filtering #
    ######################

    log.debug("=" * 60)
    if tax_ids_filter:
        log.debug(
            f"Removing {df.organism_id.isnull().sum()} proteins without taxonomy ids"
        )
        df = df[~df.organism_id.isnull()]
        tax_ids_keep = set(tax_ids_filter)
        log.debug(f"filtering for tax_ids {tax_ids_keep}")
        log.debug(
            f"tax_id value counts: \n{df.organism_id.value_counts()[list(tax_ids_keep)].to_string()}"
        )
        df = df[df.organism_id.isin(tax_ids_keep)]
        for tax_id in tax_ids_keep:
            if tax_id not in df.organism_id.value_counts().index:
                raise RuntimeError(f"No proteins found for tax id {tax_id}")
        log.debug(f"Proteins left in the dataset: {df.shape[0]}")
    else:
        log.debug("Not filtering for taxonomy ids.")

    ######################
    # Keyword annotation #
    ######################

    keywords_transport = {
        "Ion transport",
        "Anion exchange",
        "Protein transport",
        "Sodium/potassium transport",
        "Polysaccharide transport",
        "Bacteriocin transport",
        "Peptide transport",
        "Translocation",
        "Bacterial flagellum protein export",
        "Amino-acid transport",
        "Electron transport",
        "Lipid transport",
        "mRNA transport",
        "Neurotransmitter transport",
        "Oxygen transport",
        "Phosphate transport",
        "Ammonia transport",
        "Phosphonate transport",
        "Viral movement protein",
        "Sulfate transport",
        "Sugar transport",
        "Calcium transport",
        "Cobalt transport",
        "Copper transport",
        "Hydrogen ion transport",
        "Iron transport",
        "Zinc transport",
        "Nickel transport",
        "Potassium transport",
        "Sodium transport",
        "Chloride",
    }
    df["keywords_transport"] = df.keywords.str.split(";").apply(
        lambda keywords: ";".join(
            [keyword for keyword in keywords if keyword in keywords_transport]
        )
    )
    keywords_transport_related = {
        "Antibiotic resistance",
        "Transport",
        "Symport",
        "Antiport",
        "ER-Golgi transport",
        "Ion channel",
        "Calcium channel",
        "Potassium channel",
        "Chloride channel",
        "Sodium channel",
        "Viral ion channel",
        "Voltage-gated channel",
        "Ligand-gated ion channel",
        "Porin",
        "Nuclear pore complex",
        "Respiratory chain",
    }
    df["keywords_transport_related"] = df.keywords.str.split(";").apply(
        lambda keywords: ";".join(
            [keyword for keyword in keywords if keyword in keywords_transport_related]
        )
    )
    keywords_location = {
        "Membrane",
        "Cell membrane",
        "Cell inner membrane",
        "Cell outer membrane",
        "Transmembrane",
        "Nucleus",
        "Mitochondrion",
        "Endoplasmic reticulum",
        "Plastid outer membrane",
        "Plastid inner membrane",
        "Mitochondrion outer membrane",
        "Mitochondrion inner membrane",
        "Postsynaptic cell membrane",
    }
    df["keywords_location"] = df.keywords.str.split(";").apply(
        lambda keywords: ";".join(
            [keyword for keyword in keywords if keyword in keywords_location]
        )
    )

    ################################################
    # Removing outliers                            #
    ################################################

    log.debug("=" * 60)
    if outliers:
        log.debug(
            "removing outliers: \n{}".format(
                df.loc[
                    df.index.isin(outliers),
                    [
                        "gene_names",
                        "protein_names",
                        "keywords_transport",
                        "keywords_transport_related",
                        "keywords_location",
                    ],
                ]
            )
        )

        df = df[~df.index.isin(outliers)]
    else:
        log.debug("No outliers were removed.")

    ################################################
    # Printing stats before filtering              #
    ################################################

    def get_value_counts_keywords(kw_series, split=False):
        if split:
            keywords_list = []
            for keyword_str in kw_series.values:
                keywords_split = keyword_str.split(";")
                for keyword_split in keywords_split:
                    keywords_list.append(keyword_split)
            return (
                pd.Series(keywords_list)
                .value_counts(dropna=False)
                .rename(index={"": "NONE"})
                .to_string()
            )
        else:
            return (
                kw_series.value_counts(dropna=False)
                .rename(index={"": "NONE"})
                .to_string()
            )

    log.debug("=" * 60)
    log.debug("Printing dataset stats before filtering:")
    log.debug("-" * 60)
    log.debug(
        f"Substrate keywords:\n{get_value_counts_keywords(df.keywords_transport)}"
    )
    log.debug("-" * 60)
    log.debug(
        f"Substrate keywords individual:\n{get_value_counts_keywords(df.keywords_transport, split=True)}"
    )
    log.debug("-" * 60)
    log.debug(f"Location keywords:\n{get_value_counts_keywords(df.keywords_location)}")
    log.debug("-" * 60)
    log.debug(
        f"Location keywords individual:\n{get_value_counts_keywords(df.keywords_location, split=True)}"
    )
    log.debug("-" * 60)
    log.debug(
        f"Transport-related keywords:\n{get_value_counts_keywords(df.keywords_transport_related)}"
    )
    log.debug("-" * 60)
    log.debug(
        "Transport-related keywords individual:\n{}".format(
            get_value_counts_keywords(df.keywords_transport_related, split=True)
        )
    )

    ################################################
    # Removing proteins without necessary keywords #
    ################################################
    log.debug("=" * 60)
    log.debug("Filtering...")
    log.debug("-" * 60)
    log.debug(
        f"Removing {df[df.keywords_transport == ''].shape[0]} rows without any substrate keywords"
    )
    df = df[df.keywords_transport != ""]
    log.debug(f"Proteins left after filtering: {df.shape[0]}")

    ###############################
    # Transport keyword filtering #
    ###############################
    log.debug("-" * 60)
    log.debug(f"Filtering for transport-related keywords {keywords_transport_filter}")
    substrate_keywords = {kw.strip() for kw in keywords_transport_filter}
    df = df[
        df.keywords_transport_related.str.split(";").apply(
            lambda l: len(set(l) & substrate_keywords) > 0
        )
    ]
    log.debug(f"Proteins left after filtering: {df.shape[0]}")

    #################################
    # Compartment keyword filtering #
    #################################
    log.debug("-" * 60)
    log.debug(f"Filtering for location keywords {keywords_component_filter}")
    df = df[
        df.keywords_location.str.split(";").apply(
            lambda l: len(set(l) & {kw.strip() for kw in keywords_component_filter}) > 0
        )
    ]
    log.debug(f"Proteins left after filtering: {df.shape[0]}")

    ###############################
    # Substrate keyword filtering #
    ###############################

    log.debug("-" * 60)
    log.debug(f"Filtering for substrates {keywords_substrate_filter}...")
    log.debug(
        f"Applying strategy '{keywords_substrate_filter}' for multiple substrates"
    )
    keywords_filter_set = {kw.strip() for kw in keywords_substrate_filter}

    if multi_substrate == "keep":
        # Keep protein if at least one of the substrates is in parameter
        df = df[
            df.keywords_transport.str.split(";").apply(
                lambda l: len(set(l) & keywords_filter_set) > 0
            )
        ]
    elif multi_substrate in {"remove", "integrate"}:
        if multi_substrate == "integrate":
            # Remove all other keywords. Only those proteins where exactly one desired keywords is left are kept
            df.keywords_transport = df.keywords_transport.str.split(";").apply(
                lambda kw_list: ";".join(
                    [kw for kw in kw_list if kw in keywords_filter_set]
                )
            )
        # Only keep protein if it is annotated with one substrate, and that substrate is in parameter
        df = df[df.keywords_transport.apply(lambda s: s.strip() in keywords_filter_set)]
    else:
        # Should not happen, handled by argparse
        raise ValueError("Invalid parameter for multi_substrate")

    log.debug(f"Proteins left after filtering: {df.shape[0]}")

    ##################
    # Printing stats #
    ##################

    log.debug("=" * 60)
    log.debug("Dataset stats:")
    log.debug("-" * 60)
    log.debug(
        f"Final substrate counts: \n{get_value_counts_keywords(df.keywords_transport)}"
    )
    log.debug("-" * 60)

    #     # counts for tcdb classes

    def get_tcdb_stats(df):
        df_tcdb = df[~df.tcdb_id.isnull()]
        log.debug(f"Transporters without tcdb entry: {df.shape[0] - df_tcdb.shape[0]}")
        df_tcdb = pd.DataFrame.from_dict(
            {
                "tcdb_class": range(1, 10),
                "count": [
                    df_tcdb[df_tcdb.tcdb_id.str.startswith(str(i))].shape[0]
                    for i in range(1, 10)
                ],
                "class_name": [
                    "Channels/Pores",
                    "Electrochemical Potential-driven Transporters",
                    "Primary Active Transporters",
                    "Group Translocators",
                    "Transmembrane Electron Carriers",
                    "",
                    "",
                    "Accessory Factors Involved in Transport",
                    "Incompletely Characterized Transport Systems",
                ],
            }
        )
        df_tcdb = df_tcdb.set_index("tcdb_class")
        return df_tcdb

    log.debug("TCDB Stats all keywords:\n{}".format(get_tcdb_stats(df)))
    for keyword in df.keywords_transport.unique():
        log.debug("-" * 60)
        log.debug(
            "TCDB Stats for keyword {}:\n{}".format(
                keyword, get_tcdb_stats(df[df.keywords_transport == keyword])
            )
        )

    ########################
    # Sequence clustering  #
    ########################

    if sequence_clustering:
        log.debug(
            f"performing {sequence_clustering}% sequence clustering. proteins before {len(df.index)}"
        )
        cluster_repr = cd_hit(
            df.sequence, identity_threshold=sequence_clustering, verbose=verbose
        )
        df = df.loc[cluster_repr]
        log.debug(f"finished sequence clustering. proteins after {len(df.index)}")

    ########################
    # TCDB class field     #
    ########################

    df = df.assign(tcdb_class=df.tcdb_id.fillna("0.0").apply(lambda x: x[:3]))

    ########################
    # Writing Output Files #
    ########################

    df_fasta = df[
        [
            "keywords_transport",
            "gene_names",
            "protein_names",
            "tcdb_id",
            "organism_id",
            "sequence",
        ]
    ]

    if output_fasta:
        df_fasta = df[
            [
                "keywords_transport",
                "gene_names",
                "protein_names",
                "tcdb_id",
                "organism_id",
                "sequence",
            ]
        ]
        fasta_data = (
            df_fasta.reset_index()
            .apply(
                lambda row: (
                    (">sp" + "|{}" * 6).format(
                        row.Uniprot,
                        ";".join(row.gene_names.split()),
                        row.organism_id,
                        row.tcdb_id,
                        row.keywords_transport,
                        row.protein_names,
                    ),
                    row.sequence,
                ),
                axis=1,
                result_type="reduce",
            )
            .tolist()
        )
        log.debug(f"writing output fasta file to {output_fasta}...")
        write_fasta(fasta_file_name=output_fasta, fasta_data=fasta_data)

    df_tsv = df[
        [
            "keywords_transport",
            "keywords_location",
            "keywords_transport_related",
            "gene_names",
            "protein_names",
            "tcdb_id",
            "tcdb_class",
            "organism_id",
            "sequence",
        ]
    ]
    if output_tsv:
        df_tsv.to_csv(output_tsv, sep="\t")
    return df_tsv


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--keywords-substrate",
        type=str,
        required=True,
        help="List of Uniprot keywords, in Quotes and separated by semicolons",
    )
    parser.add_argument(
        "--keywords-component",
        type=str,
        required=True,
        help="List of Uniprot keywords, in Quotes and separated by semicolons",
    )
    parser.add_argument(
        "--keywords-transport",
        type=str,
        required=True,
        help="List of Uniprot keywords, in Quotes and separated by semicolons",
    )
    parser.add_argument("--input-file", type=str, required=True)
    parser.add_argument("--output-tsv", type=str)
    parser.add_argument("--output-fasta", type=str)
    parser.add_argument("--output-log", type=str)
    parser.add_argument(
        "--tax-ids", type=int, nargs="+", help="tax id(s) to filter for", default=None,
    )
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument(
        "--multi-substrate",
        choices=["keep", "remove", "integrate"],
        default="keep",
        help="How to treat proteins with multiple substrates. \
            Remove them, keep them, or integrate them into \
            single-substrate classes if the other substrates are not in the dataset",
    )
    args = parser.parse_args()
    create_dataset(
        keywords_substrate_filter=args.keywords_substrate.split(";"),
        keywords_component_filter=args.keywords_component.split(";"),
        keywords_transport_filter=args.keywords_transport.split(";"),
        input_file=args.input_file,
        multi_substrate=args.multi_substrate,
        verbose=args.verbose,
        tax_ids_filter=args.tax_ids,
        output_tsv=args.output_tsv,
        output_fasta=args.output_fasta,
        output_log=args.output_log,
    )
