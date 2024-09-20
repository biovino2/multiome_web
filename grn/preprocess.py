"""Preprocess cell oracle data for plotting GRNs.

08/10/24    Ben Iovino    CZ-Biohub
"""

import os
import celloracle as co
import pandas as pd
import seaborn as sns
import sys

sys.setrecursionlimit(10000)


def load_links(path: str, files: 'list[str]') -> 'dict[str: co.Links]':
    """Returns a dictionary of cell oracle links objects.

    Args:
        path: Path to the links.
        files: List of file names to load.

    Returns:
        dict_links: Dictionary of cell oracle links objects.
    """

    dict_links = {}
    for dataset in files:
        if dataset.endswith(".txt"):
            continue
        file_name = f"{path}/links/{dataset}_celltype_GRNs.celloracle.links"
        dict_links[dataset] = co.load_hdf5(file_name)

    return dict_links


def get_dicts(links: 'dict[str: co.Links]') -> 'tuple[dict, dict]':
    """Returns two dictionaries, one of network scores and one of filtered GRNs, with keys for each
    timepoint.

    Args:
        dict_links: Dictionary of cell oracle links objects.

    Returns:
        dict_merged_score (dict[str: pd.Dataframe]: Dictionary of network scores.
            ex. {'TDR126': pd.Dataframe}
        dict_filtered GRNs (dict[str: dict[str: pd.Dataframe]]): Dictionary of filtered GRNs.
            ex. {'TDR126': {'PSM': pd.Dataframe}}
    """

    merged_score = {}
    filtered_GRNs = {}
    for dataset in links.keys():
        filtered_GRNs[dataset] = links[dataset].filtered_links

    return merged_score, filtered_GRNs


def process_celltypes(path: str, timepoints: 'list[str]', celltypes: 'list[str]', filtered_GRNs: 'dict[str: dict[str: pd.Dataframe]]'):
    """Saves clustered heatmaps to csv files by cell type for each time point.

    Args:
        path (str): Path to save the csv files.
        timepoints (list[str]): List of timepoints.
        filtered_GRNs (dict[str: dict[str: pd.Dataframe]]): Dictionary of filtered GRNs.
    """

    # For each celltype, get counts for each timepoint and linkages for all timepoints
    vmax, vmin = 0.1, -0.1
    for ct in celltypes:

        # Step 1. collect all sources and targets across all timepoints
        all_sources = set()
        all_targets = set()
        for timepoint in timepoints:
            df = filtered_GRNs[timepoint][ct]
            all_sources.update(df['source'].unique())
            all_targets.update(df['target'].unique())

        # Step 2: Recreate each df_counts DataFrame
        df_counts_union = {}
        for timepoint in timepoints:
            df = filtered_GRNs[timepoint][ct]
    
            # Pivot the original DataFrame
            df_pivot = df.pivot(index='target',
                                columns='source',
                                values='coef_mean').reindex(index=all_targets,
                                                            columns=all_sources).fillna(0)
            df_counts_union[timepoint] = df_pivot

        # compute the linkages from the first and the last timepoints, by augmenting the "time" components
        df_counts1 = df_counts_union["TDR126"]
        df_counts2 = df_counts_union["TDR127"]
        df_counts3 = df_counts_union["TDR128"]
        df_counts4 = df_counts_union["TDR118"]
        df_counts5 = df_counts_union["TDR125"]
        df_counts6 = df_counts_union["TDR124"]

        # concatenate over the rows/cols
        df_counts_rows = pd.concat([df_counts1, df_counts2, df_counts3,
                            df_counts4, df_counts5, df_counts6], axis=1)
        df_counts_cols = pd.concat([df_counts1, df_counts2, df_counts3,
                            df_counts4, df_counts5, df_counts6], axis=0)
        
        # create a clustered heatmap for the "rows"
        g1 = sns.clustermap(df_counts_rows, method='ward', metric='euclidean', 
                    cmap='coolwarm', standard_scale=None,
                    row_cluster=True, col_cluster=True, 
                    xticklabels=False, yticklabels=False, 
                    vmax=vmax, vmin=vmin)

        # create a clustered heatmap for the "cols"
        g2 = sns.clustermap(df_counts_cols, method='ward', metric='euclidean', 
                    cmap='coolwarm', standard_scale=None,
                    row_cluster=True, col_cluster=True, 
                    xticklabels=2, yticklabels=2, 
                    vmax=vmax, vmin=vmin)

        # Reorder the dataframes
        row_order = g1.dendrogram_row.reordered_ind
        col_order = g2.dendrogram_col.reordered_ind
        for timepoint in timepoints:
            df_counts_union[timepoint] = df_counts_union[timepoint].iloc[row_order, col_order]

        # Save reordered dataframes
        if not os.path.exists(f'{path}/{ct}'):
            os.makedirs(f'{path}/{ct}')
        if not os.path.exists(f'{path}/{ct}/tp'):
            os.makedirs(f'{path}/{ct}/tp')
        for timepoint in timepoints:
            df_counts_union[timepoint].to_csv(f"{path}/tp/{ct}/{ct}_{timepoint}.csv")


def process_timepoints(path, timepoints, celltypes, filtered_GRNs):
    """
    """

    # For each timepoint, get counts for each celltype and linkages for all celltypes
    vmax, vmin = 0.1, -0.1
    for tp in timepoints:

        all_sources = set()
        all_targets = set()

        for celltype in celltypes:
            df = filtered_GRNs[tp][celltype]
            all_sources.update(df['source'].unique())
            all_targets.update(df['target'].unique())

        # Step 2: Recreate each df_counts DataFrame
        df_counts_union = {}

        for celltype in celltypes:
            df = filtered_GRNs[tp][celltype]
            # Pivot the original DataFrame
            df_pivot = df.pivot(index='target', columns='source', values='coef_mean').reindex(index=all_targets, columns=all_sources).fillna(0)
            df_counts_union[celltype] = df_pivot

        # compute the linkages from the first and the last timepoints, by augmenting the "time" components
        df_counts1 = df_counts_union["neural_posterior"]
        df_counts2 = df_counts_union["spinal_cord"]
        df_counts3 = df_counts_union["NMPs"]
        df_counts4 = df_counts_union["tail_bud"]
        df_counts5 = df_counts_union["PSM"]
        df_counts6 = df_counts_union["somites"]

        # concatenate over the columns
        df_counts_rows = pd.concat([df_counts1, df_counts2, df_counts3,
                            df_counts4, df_counts5, df_counts6], axis=1)

        # concatenate over the rows
        df_counts_cols = pd.concat([df_counts1, df_counts2, df_counts3,
                            df_counts4, df_counts5, df_counts6], axis=0)

        # create a clustered heatmap for the "rows"
        g1 = sns.clustermap(df_counts_rows, method='ward', metric='euclidean', 
                    cmap='coolwarm', standard_scale=None,
                    row_cluster=True, col_cluster=True, 
                    xticklabels=False, yticklabels=False, 
                    vmax=vmax, vmin=vmin)

        # create a clustered heatmap for the "cols"
        g2 = sns.clustermap(df_counts_cols, method='ward', metric='euclidean', 
                    cmap='coolwarm', standard_scale=None,
                    row_cluster=True, col_cluster=True, 
                    xticklabels=False, yticklabels=False, 
                    vmax=vmax, vmin=vmin)
        
        # Reorder the dataframes
        row_order = g1.dendrogram_row.reordered_ind
        col_order = g2.dendrogram_col.reordered_ind
        for celltype in celltypes:
            df_counts_union[celltype] = df_counts_union[celltype].iloc[row_order, col_order]

        # Save reordered dataframes
        if not os.path.exists(f'{path}/{tp}'):
            os.makedirs(f'{path}/{tp}')
        if not os.path.exists(f'{path}/{tp}/ct'):
            os.makedirs(f'{path}/{tp}/ct')
        for celltype in celltypes:
            df_counts_union[celltype].to_csv(f"{path}/ct/{celltype}/{celltype}_{tp}.csv")


def main():
    """
    """

    path = 'grn/data'
    if not os.path.exists(path):
        os.makedirs(path)
    timepoints = ['TDR126', 'TDR127', 'TDR128', 'TDR118', 'TDR125', 'TDR124']
    celltypes = ['NMPs', 'PSM', 'neural_posterior', 'somites', 'spinal_cord', 'tail_bud']
    
    #  Load the pruned links and get network scores + GRNs
    links = load_links(path, timepoints)
    _, filtered_GRNs = get_dicts(links)
    process_celltypes(path, timepoints, celltypes, filtered_GRNs)
    process_timepoints(path, timepoints, celltypes, filtered_GRNs)


if __name__ == '__main__':
    main()

