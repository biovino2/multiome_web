"""Preprocess cell oracle data for plotting knock outs.

09/05/24    Ben Iovino    CZ-Biohub
"""

import anndata as ad
import celloracle as co
import os


def get_adata(path: str,timepoints: 'list[str]'):
    """Loads a celloracle .oracle objects and writes them to anndata objects.

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
    common_genes = set(genes[timepoints[0]])
    for timepoint in timepoints[1:]:
        common_genes = common_genes.intersection(genes[timepoint])
    with open(f'{path}/common_genes.txt', 'w') as f:
        for gene in common_genes:
            f.write(f'{gene}\n')


def main():
    """
    """

    path = 'data'
    if not os.path.exists("data"):
        os.makedirs("data")

    timepoints = ['TDR126', 'TDR127', 'TDR128', 'TDR118', 'TDR125', 'TDR124']
    get_adata(path, timepoints)


if __name__ == "__main__":
    main()

