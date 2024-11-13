"""Preprocess cell oracle data for plotting GRNs.

08/10/24    Ben Iovino    CZ-Biohub
"""

import os
import celloracle as co
import pandas as pd
import seaborn as sns
import sys
import zipfile
from util import get_timepoints_abbr

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
    path = '/hpc/projects/data.science/yangjoon.kim/zebrahub_multiome/data/processed_data/09_NMPs_subsetted_v2'
    for dataset in files:
        if dataset.endswith(".txt"):
            continue
        file_name = f"{path}/{dataset}/08_{dataset}_celltype_GRNs.celloracle.links"
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


def cluster_counts(df_grn: pd.DataFrame) -> pd.DataFrame:
    """Returns a clustered dataframe of counts. Rows are clustered by TF family and clustered
    within each family.

    Args:
        df_grn: Dataframe of counts.

    Returns:
        final_data: Clustered dataframe
    """

    df_grn.fillna(0, inplace=True)
    df_grn = df_grn.loc[(df_grn == 0).sum(axis=1) < 5]

    # Take family of TF with regex
    reg = r'^(.*?)(?=\d)'
    df_grn['TF'] = df_grn.index.str.split('_').str[0]
    df_grn['TF_family'] = df_grn['TF'].str.extract(reg)

    # If NaN in TF_family, just take name of TF
    df_grn['TF_family'] = df_grn['TF_family'].fillna(df_grn['TF'])
    df_grn['TF_family'] = df_grn['TF_family'].str[:3]

    # Order TF's by frequency
    tf_counts = df_grn['TF_family'].value_counts()
    ordered_tfs = tf_counts.index

    # Cluster each TF family separately
    clustered_subgroups = {}
    for tf_name in ordered_tfs:
        tf_data = df_grn[df_grn['TF_family'] == tf_name].drop(columns=['TF_family', 'TF'])
        if tf_data.shape[0] == 1:  # Skip if row has only one entry, can't cluster
            continue
        clustermap = sns.clustermap(tf_data, row_cluster=True, col_cluster=False)
    
        # Retrieve the ordered rows based on clustering
        ordered_indices = clustermap.dendrogram_row.reordered_ind
        clustered_tf_data = tf_data.iloc[ordered_indices]
        clustered_tf_data = clustered_tf_data.iloc[::-1]
        clustered_subgroups[tf_name] = clustered_tf_data

    # Concatenate clustered subgroups in the order of TF frequency
    final_data = pd.concat(clustered_subgroups.values())
    final_data = final_data.iloc[::-1]

    return final_data


def cluster_timepoints(filtered_GRNs: dict, path: str, timepoints: 'list[str]', celltypes: 'list[str]'):
    """
    """

    if not os.path.exists(f'{path}/tp'):
        os.makedirs(f'{path}/tp')
    for tp in timepoints:
        df_grn = pd.DataFrame()
        for ct in celltypes:
            df = filtered_GRNs[tp][ct]

            # Refactor dataframe for source_target and coef_mean
            df['TF-gene'] = df['source'] + '_' + df['target']
            df_filt = df[['TF-gene', 'coef_mean']].set_index('TF-gene')
            df_filt.rename(columns={'coef_mean': ct}, inplace=True)
            df_grn = pd.concat([df_grn, df_filt], axis=1)

        # Cluster and save to csv
        df_grn_clustered = cluster_counts(df_grn)
        df_grn_clustered.to_csv(f'{path}/tp/{tp}.csv')


def cluster_celltypes(filtered_GRNs: dict, path: str, timepoints: 'list[str]', celltypes: 'list[str]'):
    """
    """

    if not os.path.exists(f'{path}/ct'):
        os.makedirs(f'{path}/ct')
    for ct in celltypes:
        df_grn = pd.DataFrame()
        for tp in timepoints:
            df = filtered_GRNs[tp][ct]

            # Refactor dataframe for source_target and coef_mean
            df['TF-gene'] = df['source'] + '_' + df['target']
            df_filt = df[['TF-gene', 'coef_mean']].set_index('TF-gene')
            df_filt.rename(columns={'coef_mean': tp}, inplace=True)
            df_grn = pd.concat([df_grn, df_filt], axis=1)

        # Cluster and save to csv
        df_grn_clustered = cluster_counts(df_grn)
        df_grn_clustered.to_csv(f'{path}/ct/{ct}.csv')


def get_scores(path: str, metric: str, top_n: int):
    """Saves top n scores for each cell type at each time point.

    Args:
        path (str): Path to save the csv files.
        metric (str): Metric to sort by.
        top_n (int): Number of top scores to save.
    """

    if not os.path.exists(f'{path}/scores'):
        os.makedirs(f'{path}/scores')

    timepoints = get_timepoints_abbr()
    celltypes = ['NMPs', 'PSM', 'neural_posterior', 'somites', 'spinal_cord', 'tail_bud']
    for dataset in list(timepoints.keys()):

        links = co.load_hdf5(f'/hpc/projects/data.science/benjamin.iovino/multiome_web/src/grn/data/links/08_{dataset}_celltype_GRNs.celloracle.links')
        scores = links.merged_score[['degree_centrality_all', 'cluster']]

        # For each cell type, sort by degree_centrality_all and get top 30
        top_scores = {}
        for ct in celltypes:
            df = scores[scores['cluster'] == ct]
            df = df.sort_values(by=metric, ascending=False)
            top_scores[ct] = df.head(top_n)

            if not os.path.exists(f'{path}/scores/{timepoints[dataset]}'):
                os.makedirs(f'{path}/scores/{timepoints[dataset]}')

            top_scores[ct].to_csv(f'{path}/scores/{timepoints[dataset]}/{ct}.csv')


def main():
    """
    """

    path = os.path.dirname(os.path.abspath(__file__))+'/data'
    if not os.path.exists(path):
        os.makedirs(path)
    timepoints = ['TDR126', 'TDR127', 'TDR128', 'TDR118', 'TDR125', 'TDR124']
    celltypes = ['neural_posterior', 'spinal_cord', 'NMPs', 'tail_bud', 'PSM', 'somites']
    
    #  Load the pruned links and get network scores + GRNs
    links = load_links(path, timepoints)
    _, filtered_GRNs = get_dicts(links)
    cluster_timepoints(filtered_GRNs, path, timepoints, celltypes)
    cluster_celltypes(filtered_GRNs, path, timepoints, celltypes)
    

    '''
    # Save data as zipfile for download
    with zipfile.ZipFile(f'{path}/data.zip', 'w') as z:
        for root, _, files in os.walk(path):
            if root.endswith('links'):
                continue
            for file in files:
                z.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), path))
    '''


if __name__ == '__main__':
    main()
