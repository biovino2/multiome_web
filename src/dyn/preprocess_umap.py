"""Preprocess adata to get dataframe of UMAP coordinates and cluster assignments.

12/04/24    Ben Iovino    CZ-Biohub
"""

import scanpy as sc
import pandas as pd

def get_umap_df(adata: 'sc.AnnData') -> pd.DataFrame:
    """Returns a dataframe of UMAP coordinates and cluster assignments.

    Args:
        adata (sc.AnnData): The annotated data object.

    Returns:
        pd.DataFrame: The dataframe of UMAP coordinates and cluster assignments.
    """

    df_umap = pd.DataFrame(index=adata.obs.index, data=adata.obsm['X_umap'], columns=['umap1', 'umap2'])
    df_umap['Leiden'] = adata.obs['leiden_0.7']

    # Map colors to cluster
    cluster_colors = adata.uns['leiden_0.7_colors']
    color_dict = {i: color for i, color in enumerate(cluster_colors)}
    df_umap['color'] = df_umap['Leiden'].map(color_dict)

    return df_umap


def main():
    """
    """

    adata_file_path = "/hpc/projects/data.science/yangjoon.kim/zebrahub_multiome/data/processed_data/annotations/genes_by_ct_tp_top_5068genes.h5ad"
    adata = sc.read(adata_file_path)
    df_umap = get_umap_df(adata)
    df_umap.to_csv('src/dyn/data/umap.csv')


if __name__ == '__main__':
    main()
