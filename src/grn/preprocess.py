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


def load_links(files: 'list[str]') -> 'dict[str: co.Links]':
    """Returns a dictionary of cell oracle links objects.

    Args:
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


def cluster_counts(df_grn: pd.DataFrame, drop: int) -> pd.DataFrame:
    """Returns a clustered dataframe. Rows are ordered by TF family and clustering is performed
    within each family.

    Args:
        df_grn: Dataframe with rows as TF-gene pairs and cols as cell types/time points.
        drop: Minimum number of non-zero entries in a row (TF-gene) to keep it in the dataframe.

    Returns:
        df_all: Clustered dataframe
    """

    df_grn.fillna(0, inplace=True)
    if drop:
        df_grn = df_grn.loc[(df_grn == 0).sum(axis=1) < 6-drop]  # works for 6 columns..must fix

    # Take family of TF and order by frequency
    r'^[A-Za-z]{0,3}|^[A-Za-z]*?\d'
    df_grn['TF'] = df_grn.index.str.split('_').str[0]
    df_grn['TF_family'] = df_grn['TF'].str.extract(r'(^[A-Za-z]{0,3}|^[A-Za-z]*?\d)')
    tf_counts = df_grn['TF_family'].value_counts()
    ordered_tfs = tf_counts.index

    # Cluster each TF family separately
    clustered_subgroups = {}
    for tf_name in ordered_tfs:
        df_tf = df_grn[df_grn['TF_family'] == tf_name].drop(columns=['TF_family', 'TF'])
        if df_tf.shape[0] == 1:  # Skip if row has only one entry, can't cluster
            continue
        clustermap = sns.clustermap(df_tf, row_cluster=True, col_cluster=False)
    
        # Retrieve the ordered rows based on clustering
        ordered_indices = clustermap.dendrogram_row.reordered_ind
        clustered_tf_data = df_tf.iloc[ordered_indices]
        clustered_tf_data = clustered_tf_data.iloc[::-1]
        clustered_subgroups[tf_name] = clustered_tf_data

    # Concatenate clustered subgroups in the order of TF frequency
    df_all = pd.concat(clustered_subgroups.values())
    df_all = df_all.iloc[::-1]
    df_all['TF'] = df_all.index.str.split('_').str[0]
    df_all['TF_family'] = df_all['TF'].str.extract(r'(^[A-Za-z]{0,3}|^[A-Za-z]*?\d)')

    return df_all


def order_tfs(df_grn: pd.DataFrame, category: str) -> pd.DataFrame:
    """Returns a dataframe similar to the input, but with TFs ordered by time point/cell type
    within each family.

    Args:
        df_grn: Dataframe with rows as TF-gene pairs and cols as cell types/time points.
        category: 'timepoint' or 'celltype'.

    Returns:
        df_grn_max: Dataframe with TFs ordered by time point/cell type within each family.
    """

    # Change column names to integers for ordering (necessary for cell types)
    if category == 'timepoint':
        mapping1 = {'10hpf': 1, '12hpf': 2, '14hpf': 3, '16hpf': 4, '19hpf': 5, '24hpf': 6}
        mapping2 = {1: '10hpf', 2: '12hpf', 3: '14hpf', 4: '16hpf', 5: '19hpf', 6: '24hpf'}
    if category == 'celltype':
        mapping1 = {'neural_posterior': 1, 'spinal_cord': 2, 'NMPs': 3, 'tail_bud': 4, 'PSM': 5, 'somites': 6}
        mapping2 = {1: 'neural_posterior', 2: 'spinal_cord', 3: 'NMPs', 4: 'tail_bud', 5: 'PSM', 6: 'somites'}
    df_grn.rename(columns=mapping1, inplace=True)

    # Sum up values for each TF (count any non-zero as 1) and find the max for ordering
    df_sum = pd.DataFrame()
    df_sum = df_grn.groupby('TF').apply(lambda x: (x.iloc[:, :] != 0).sum()).drop(columns=['TF', 'TF_family'])
    df_sum['max'] = df_sum.idxmax(axis=1)

    # Map max value to original dataframe
    df_grn_max = df_grn.copy()
    df_grn_max['max'] = df_grn['TF'].map(df_sum['max'])

    # Order by max value within each TF family
    df_grn_ordered = pd.DataFrame()
    for tf_fam in df_grn_max['TF_family'].unique():
    
        # Copy rows with same TF family and order within family by max value
        df_tf_fam = df_grn_max[df_grn_max['TF_family'] == tf_fam]
        df_tf_fam = df_tf_fam.sort_values(by='max')[::-1]
        df_grn_ordered = pd.concat([df_grn_ordered, df_tf_fam])  # Maintain family structure in new df

    # Reverse mappings
    df_grn_ordered.rename(columns=mapping2, inplace=True)

    return df_grn_ordered


def cluster_timepoints(filtered_GRNs: dict, path: str, timepoints: 'list[str]', celltypes: 'list[str]'):
    """Saves a pandas dataframe of the GRN for each time point across all cell types.

    Args:
        filtered_GRNs (dict): Dictionary of filtered GRNs.
        path (str): Path to save the csv files.
        timepoints (list[str]): List of time points.
        celltypes (list[str]): List of cell types.
    """

    tp_dict = get_timepoints_abbr()
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
        df_grn_clustered = cluster_counts(df_grn, None)
        df_grn_ordered = order_tfs(df_grn_clustered, 'timepoint')
        df_grn_ordered.to_csv(f'{path}/tp/{tp_dict[tp]}.csv')


def cluster_celltypes(filtered_GRNs: dict, path: str, timepoints: 'list[str]', celltypes: 'list[str]'):
    """Saves a pandas dataframe of the GRN for each cell type across all time points.

    Args:
        filtered_GRNs (dict): Dictionary of filtered GRNs.
        path (str): Path to save the csv files.
        timepoints (list[str]): List of time points.
        celltypes (list[str]): List of cell types.
    """

    tp_dict = get_timepoints_abbr()
    if not os.path.exists(f'{path}/ct'):
        os.makedirs(f'{path}/ct')
    for ct in celltypes:
        df_grn = pd.DataFrame()
        for tp in timepoints:
            df = filtered_GRNs[tp][ct]

            # Refactor dataframe for source_target and coef_mean
            df['TF-gene'] = df['source'] + '_' + df['target']
            df_filt = df[['TF-gene', 'coef_mean']].set_index('TF-gene')
            df_filt.rename(columns={'coef_mean': tp_dict[tp]}, inplace=True)
            df_grn = pd.concat([df_grn, df_filt], axis=1)

        # Cluster and save to csv
        df_grn_clustered = cluster_counts(df_grn, None)
        df_grn_ordered = order_tfs(df_grn_clustered, 'celltype')
        df_grn_ordered.to_csv(f'{path}/ct/{ct}.csv')


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
    links = load_links(timepoints)
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
