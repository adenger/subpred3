import pandas as pd


def __read_df(filepath: str) -> pd.DataFrame:
    gff_df = pd.read_table(
        filepath,
        comment="#",
        low_memory=False,
        header=None,
        names=[
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attribute",
        ],
        usecols=[
            "seqid",
            "source",
            "type",
            "start",
            "end",
            # "score",
            "strand",
            # "phase",
            "attribute",
        ],
    )
    return gff_df


def __split_gff_attributes(attribute_col):
    attribute_list = attribute_col.str.split(";").to_list()
    attribute_records = []
    for attribute_sublist in attribute_list:
        attribute_sublist_dict = {}
        for attribute_subslist_element in attribute_sublist:
            attribute_subslist_element_items = attribute_subslist_element.split("=")
            attribute_sublist_dict[
                attribute_subslist_element_items[0]
            ] = attribute_subslist_element_items[1]
        attribute_records.append(attribute_sublist_dict)

    attribute_df = pd.DataFrame.from_records(attribute_records)
    return attribute_df


def read_gff3(filepath: str):

    gff_df = __read_df(filepath=filepath)
    attribute_df = __split_gff_attributes(gff_df.attribute)

    gff_df = gff_df.join(attribute_df)
    gff_df = gff_df.drop("attribute", axis=1)

    return gff_df
