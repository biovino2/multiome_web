"""Preprocess cell oracle data for plotting GRNs.

08/10/24    Ben Iovino    CZ-Biohub
"""

import os
import celloracle as co
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys

sys.setrecursionlimit(10000)


def load_links(files: 'list[str]') -> 'dict[str: co.Links]':
    """Returns a dictionary of cell oracle links objects.

    Args:
        files: List of file names to load.

    Returns:
        dict_links: Dictionary of cell oracle links objects.
    """

    dict_links = {}
    for dataset in files:
        if dataset.endswith(".txt"):
            continue
        file_name = f"data/links/{dataset}_celltype_GRNs.celloracle.links"
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


def main():
    """
    """

    if not os.path.exists('data/'):
        os.makedirs('data')

    #  Load the pruned links and get network scores + GRNs
    timepoints = ['TDR126', 'TDR127', 'TDR128', 'TDR118', 'TDR125', 'TDR124']
    
    links = load_links(timepoints)
    merged_score, filtered_GRNs = get_dicts(links)

    # For each celltype, get counts for each timepoint and linkages for all timepoints
    celltypes = ['NMPs', 'PSM', 'fast_muscle', 'neural_posterior', 'somites', 'spinal_cord', 'tail_bud']
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

            # Save to celltype_timepoint to csv
            if not os.path.exists(f'data/{ct}'):
                os.makedirs(f'data/{ct}')
            df_pivot.to_csv(f"data/{ct}/{ct}_{timepoint}.csv")


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

        # extract the row/col indices
        row_linkage = g1.dendrogram_row.linkage
        col_linkage = g2.dendrogram_col.linkage

        # Save linkages to npz
        np.savez(f'data/{ct}/{ct}_row_linkage.npz', row_linkage)
        np.savez(f'data/{ct}/{ct}_col_linkage.npz', col_linkage)


if __name__ == '__main__':
    main()
