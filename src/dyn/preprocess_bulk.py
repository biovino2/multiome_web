"""Preprocess adata for each timepoint for plotting pseudobulk correlations.

10/17/24    Ben Iovino    CZ-Biohub
"""

import numpy as np
import pandas as pd
import scanpy as sc


def replace_periods_with_underscores(df: pd.DataFrame) -> pd.DataFrame:
    """Replace periods with underscores in column names.

    Args:
        df (pd.DataFrame): The dataframe.

    Returns:
        pd.DataFrame: The updated dataframe.
    """

    df.columns = df.columns.str.replace('.', '_', regex=False)

    return df


def get_objects(path: str, objects: str) -> tuple:
    """Return two adata objects, one for RNA and one for ATAC.

    Args:
        path (str): The path to the data.
        objects (str): The list of object names to read.

    Returns:
        tuple (sc.AnnData, sc.AnnData): The adata objects.
    """

    adata_RNA = sc.read_h5ad(f'{path}/{objects[0]}')
    adata_ATAC = sc.read_h5ad(f'{path}/{objects[1]}')

    # filter out the "low_quality_cells" - either using the adata_RNA cell_ids, or the latest annotation csv file
    adata_ATAC = adata_ATAC[adata_ATAC.obs_names.isin(adata_RNA.obs_names)]

    # remove unnecessary fields from adata_ATAC 
    columns_to_drop = adata_ATAC.obs.columns[adata_ATAC.obs.columns.str.startswith("prediction")]
    adata_ATAC.obs = adata_ATAC.obs.drop(columns=columns_to_drop)

    # replace the "." within the strings in adata.obs.columns to "_" (this is for exCellxgene)
    adata_ATAC.obs = replace_periods_with_underscores(adata_ATAC.obs)
    del adata_ATAC.raw

    return adata_RNA, adata_ATAC


def process_objects(adata: sc.AnnData, common_genes: list[str], annotation: pd.Series) -> tuple:
    """Returns two pandas dataframes, one containing the mean values for each gene in the adata
    objects and the second containing the standard error of values for each gene.

    Args:
        adata (sc.AnnData): The adata object with all genes.
        common_genes (list[str]): The list of common genes between all adata objects.
        annotation (pd.Series): The manual annotations for each cell.

    Returns:
        tuple (pd.DataFrame, pd.DataFrame): The pandas dataframes.
    """

    # subset the adata_RNA for the list of genes (this is to reduce the computing resource/time)
    adata_sub = adata[:, adata.var_names.isin(common_genes)]

    # create a dataframe of cells-by-gene.activity (Signac)
    count_matrix = pd.DataFrame(adata_sub.X.todense(),
                             index=adata_sub.obs.index,
                             columns=adata_sub.var_names)

    # transfer the "dataset" category labels to count_matrice df
    count_matrix["dataset"] = adata_sub.obs["dataset"]
    count_matrix["annotation_ML_coarse"] = annotation

    # add the "timepoints" category for the column
    dict_timepoints = {"TDR126": "10hpf",
                   "TDR127": "12hpf",
                   "TDR128": "14hpf",
                   "TDR118": "16hpf",
                   "TDR119": "16hpf",
                   "TDR125": "19hpf",
                   "TDR124": "24hpf"}
    count_matrix["timepoints"] = count_matrix["dataset"].map(dict_timepoints)

    # Compute mean values across all cells per "timepoints" (TDR118 and TDR119 will be merged here)
    numeric_columns = count_matrix.select_dtypes(include=[np.number]).columns.tolist()
    timepoints_by_genes = count_matrix.groupby('timepoints')[numeric_columns].mean()
    timepoints_by_genes_sem = count_matrix.groupby('timepoints')[numeric_columns].sem()

    # numeric_timepoints to re-order the rows (timepoints_by_genes)
    timepoints_by_genes['numeric_timepoints'] = timepoints_by_genes.index.str.extract('(\d+)').astype(int).values
    timepoints_by_genes_sem['numeric_timepoints'] = timepoints_by_genes_sem.index.str.extract('(\d+)').astype(int).values

    # Sort by the numeric timepoints to ensure correct order in plot
    timepoints_by_genes_sorted = timepoints_by_genes.sort_values('numeric_timepoints')
    timepoints_by_genes_sem_sorted = timepoints_by_genes_sem.sort_values('numeric_timepoints')

    return timepoints_by_genes_sorted, timepoints_by_genes_sem_sorted


def main():
    """
    """

    path = '/hpc/projects/data.science/yangjoon.kim/zebrahub_multiome/data/processed_data/01_Signac_processed'
    objects = ['integrated_RNA_ATAC_counts_RNA_master_filtered.h5ad',
                'integrated_RNA_ATAC_counts_gene_activity.h5ad']
    adata_RNA, adata_ATAC = get_objects(path, objects)
    annotation = adata_RNA.obs['annotation_ML_coarse']

    # Process objects
    common_genes = np.intersect1d(adata_RNA.var_names, adata_ATAC.var_names)
    gex_values, gex_errors = process_objects(adata_RNA, common_genes, annotation)
    atac_values, atac_errors = process_objects(adata_ATAC, common_genes, annotation)

    # Save each to csv
    gex_values.to_csv('src/dyn/data/bulk/gex_values.csv')
    gex_errors.to_csv('src/dyn/data/bulk/gex_errors.csv')
    atac_values.to_csv('src/dyn/data/bulk/atac_values.csv')
    atac_errors.to_csv('src/dyn/data/bulk/atac_errors.csv')

if __name__ == '__main__':
    main()
