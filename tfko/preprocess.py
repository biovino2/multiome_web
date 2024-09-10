"""Preprocess cell oracle data for plotting knock outs.

09/05/24    Ben Iovino    CZ-Biohub
"""

import anndata as ad
import celloracle as co
import numpy as np
import pandas as pd
import os


def get_adata(path: str, timepoints: 'list[str]'):
    """Loads celloracle .oracle objects and writes them to anndata objects. Celloracle refers to
    the active regulatory units as genes, but I believe they are actually transcription factors.

    Args:
        path (str): Path to the celloracle .oracle objects.
        timepoints (list): List of timepoints to load.
    """

    genes: 'dict[str, list]' = {}
    for timepoint in timepoints:
        filename = f'{path}/15_{timepoint}_in_silico_KO_trans_probs_added.celloracle.oracle'
        co_object = co.load_hdf5(filename)

        # Update list of genes
        genes[timepoint] = co_object.active_regulatory_genes

        # Subset for adata, few obs columns and all obsm
        adata = co_object.adata
        adata_subset = ad.AnnData(
            X=adata.X,
            obs=adata.obs[['manual_annotation', 'SEACell']],
            obsm=adata.obsm
        )
        adata_subset.write(f'data/{timepoint}_KO.h5ad')

    # Find intersection of genes
    common_tfs = set(genes[timepoints[0]])
    for timepoint in timepoints[1:]:
        common_tfs = common_tfs.intersection(genes[timepoint])
    with open(f'{path}/common_tfs.txt', 'w') as f:
        for gene in common_tfs:
            f.write(f'{gene}\n')


def average_metacells(path: str, timepoints: 'list[str]', common_tfs: list[str]):
    """Writes the average UMAP position and transition vector of metacells based off of each
    metacell's population. This is performed for each timepoint, only for the knockouts that
    are shared across all timepoints.

    Args:
        path (str): Path to the adata
        timepoints (list): List of timepoints.
        common_tfs (list): List of common transcription factors.
    """

    # Write to npz
    if not os.path.exists(f'{path}/metacells'):
        os.makedirs(f'{path}/metacells')

    for timepoint in timepoints:

        # Open adata
        adata = ad.read(f'{path}/{timepoint}_KO.h5ad')
        adata.obs['SEACell'] = adata.obs['SEACell'].astype('object')
        adata.obs['manual_annotation'] = adata.obs['manual_annotation'].astype('object')
        X_umap = adata.obsm['X_umap_aligned']

        # Convert metacell column to categorical if it's not already
        if not pd.api.types.is_categorical_dtype(adata.obs['SEACell']):
            metacells = pd.Categorical(adata.obs['SEACell'])
        else:
            metacells = adata.obs['SEACell']

        # 2D transition vectors for each metacell
        for ko in common_tfs:
            if ko == 'WT_global_nmps':
                V_cell = adata.obsm[f'{ko}_umap_aligned']
            else:
                V_cell = adata.obsm[f'{ko}_KO_umap_aligned']

            # X_metacell is the average UMAP position of the metacells
            # V_metacell is the average transition vector of the metacells
            n_metacells = len(metacells.categories)
            X_metacell = np.zeros((n_metacells, 2))
            V_metacell = np.zeros((n_metacells, 2))

            for i, category in enumerate(metacells.categories):
                mask = metacells == category
                X_metacell[i] = X_umap[mask].mean(axis =0)
                V_metacell[i] = V_cell[mask].mean(axis=0)

            np.savez(f'{path}/metacells/{timepoint}_{ko}_metacells.npz',
                     X_metacell=X_metacell,
                     V_metacell=V_metacell)


def main():
    """
    """

    path = 'tfko/data'
    if not os.path.exists(path):
        os.makedirs(path)

    # Subset large celloracle objects to small adata objects
    timepoints = ['TDR126', 'TDR127', 'TDR128', 'TDR118', 'TDR125', 'TDR124']
    get_adata(path, timepoints)

    # Read common tfs and average metacells
    with open('tfko/data/common_tfs.txt', 'r') as file:  # faster than rerunning get_adata()
        common_tfs = file.read().splitlines()
    common_tfs.append('WT_global_nmps')
    average_metacells(path, timepoints, common_tfs)


if __name__ == "__main__":
    main()

