import pandas as pd


def aggregate_maxmean(df, col_name):
    index_elements_duplicated = df.loc[df[col_name].duplicated(), col_name].unique()

    df_dupl = df[df[col_name].isin(index_elements_duplicated)]
    df_dupl = (
        df_dupl.groupby(col_name)
        .apply(
            lambda df_gr: df_gr.loc[df_gr.drop(col_name, axis=1).mean(axis=1).idxmax()]
        )
        .drop(col_name, axis=1)
    )

    df_single = df[~df[col_name].isin(index_elements_duplicated)].set_index(
        col_name, verify_integrity=True
    )
    df_maxmean = pd.concat([df_single, df_dupl])
    df_maxmean = df_maxmean.sort_index()
    return df_maxmean
