"""Preprocess adata for each timepoint for plotting correlation scatterplots.

08/16/24    Ben Iovino    CZ-Biohub
"""

import polars as pl
import scanpy as sc
import numpy as np
import os
import pandas as pd
import pickle


def define_color_dict() -> 'dict[str:str]':
    """Return a dictionary of colors for each cell type.

    Returns:
        dict[str, str]: The dictionary of colors.
    """

    cell_type_color_dict = {
        'NMPs': '#8dd3c7',
        'PSM': '#008080',
        'differentiating_neurons': '#bebada',
        'endocrine_pancreas': '#fb8072',
        'endoderm': '#80b1d3',
        'enteric_neurons': '#fdb462',
        'epidermis': '#b3de69',
        'fast_muscle': '#df4b9b',
        'floor_plate': '#d9d9d9',
        'hatching_gland': '#bc80bd',
        'heart_myocardium': '#ccebc5',
        'hemangioblasts': '#ffed6f',
        'hematopoietic_vasculature': '#e41a1c',
        'hindbrain': '#377eb8',
        'lateral_plate_mesoderm': '#4daf4a',
        'midbrain_hindbrain_boundary': '#984ea3',
        'muscle': '#ff7f00',
        'neural': '#e6ab02',
        'neural_crest': '#a65628',
        'neural_floor_plate': '#66a61e',
        'neural_optic': '#999999',
        'neural_posterior': '#393b7f',
        'neural_telencephalon': '#fdcdac',
        'neurons': '#cbd5e8',
        'notochord': '#f4cae4',
        'optic_cup': '#c0c000',
        'pharyngeal_arches': '#fff2ae',
        'primordial_germ_cells': '#f1e2cc',
        'pronephros': '#cccccc',
        'somites': '#1b9e77',
        'spinal_cord': '#d95f02',
        'tail_bud': '#7570b3'
    }

    return cell_type_color_dict


def get_datasets(path: str, datasets: 'list[str]') -> 'dict[str:sc.AnnData]':
    """Returns a dictionary of adata objects for each dataset.

    Args:
        path (str): The path to the data.
        datasets (list[str]): The list of dataset names.

    Returns:
        dict[str, sc.AnnData]: The dictionary of adata objects.
    """

    # Import all datasets (adata) and save in a dictionary
    dict_meta_ad = {}
    for data_id in datasets:

        # define the sample_id (withou "reseq" in the handle)
        sample_id = data_id.replace("reseq","")
    
        # import the single-cell adata's annotation (with adata.obs["SEACell"] annotation)
        df_seacells = pd.read_csv(f'{path}/{data_id}_seacells_obs_annotation_ML_coarse.csv', index_col=0)
    
        # Group by SEACell and find the most prevalent annotation for each SEACell
        most_prevalent = df_seacells.groupby('SEACell')['annotation_ML_coarse'].agg(lambda x: x.value_counts().idxmax())
        prevalent_count = df_seacells.groupby('SEACell')['annotation_ML_coarse'].agg(lambda x: x.value_counts().max())
        total_counts = df_seacells.groupby('SEACell')['annotation_ML_coarse'].count()
        prevalent_freq = prevalent_count/total_counts

        # define a dataframe
        df_prevalent = pd.DataFrame({
            'SEACell': most_prevalent.index,
            'most_prevalent_annotation': most_prevalent.values,
            'prevalent_frequency': prevalent_freq.values
        })
        
        df_prevalent.set_index("SEACell", inplace=True)

        # Convert the result to a dictionary
        seacell_to_annotation = df_prevalent["most_prevalent_annotation"].to_dict()

        # import adata (aggregated over seacells)
        rna_meta_ad = sc.read_h5ad(f'{path}/{data_id}/{sample_id}_RNA_seacells_aggre.h5ad')
        atac_meta_ad = sc.read_h5ad(f'{path}/{data_id}/{sample_id}_ATAC_seacells_aggre.h5ad')
    
        # First, we'll subset the features that are shared between RNA and ATAC (gene names)
        shared_genes = np.intersect1d(rna_meta_ad.var_names, atac_meta_ad.var_names)

        # subset the RNA and ATAC objects for the shared genes
        rna_meta_ad = rna_meta_ad[:, shared_genes]
        atac_meta_ad = atac_meta_ad[:, shared_genes]
    
        # transfer the "most prevalent celltype" annotation to the aggregated objects
        rna_meta_ad.obs["celltype"] = rna_meta_ad.obs_names.map(seacell_to_annotation)
        atac_meta_ad.obs["celltype"] = atac_meta_ad.obs_names.map(seacell_to_annotation)
    
        # save in a dictionary
        dict_meta_ad[f"{sample_id}_rna"] = rna_meta_ad
        dict_meta_ad[f"{sample_id}_atac"] = atac_meta_ad

    return dict_meta_ad


def get_gene_dict(dict_meta: 'dict[str:sc.AnnData]') -> 'dict[str:dict[tuple[np.ndarray, np.ndarray, list[str]]]]':
    """Returns a dictionary of genes and their corresponding expression/availability values.
    
    Args:
        dict_meta (dict[str:sc.AnnData]): The dictionary of adata objects.

    Returns:
        dict[str, dict[tuple[np.ndarray, np.ndarray, list[str]]]]: The dictionary of genes and their values.
            ex. gene_dict = {'meox1: {'TDR126': (rna_expr, atac_expr, colors),
                                    'TDR127': (rna_expr, atac_expr, colors), ...
                                    }
                            'hbbe3': {'TDR126': (rna_expr, atac_expr, colors),
                                    'TDR127': (rna_expr, atac_expr, colors), ...
                                }
                            }
    """

    datasets = set([key.split("_")[0] for key in dict_meta.keys()])
    gene_dict = {}
    for dataset in datasets:

        # Get shared genes between RNA and ATAC adata objects
        rna_ad = dict_meta[f"{dataset}_rna"]
        atac_ad = dict_meta[f"{dataset}_atac"]
        shared_genes = np.intersect1d(rna_ad.var_names, atac_ad.var_names)

        # Append gene expression values to gene dictionary
        for gene in shared_genes:
            rna_expr = rna_ad[:, gene].X.toarray().flatten()
            atac_expr = atac_ad[:, gene].X.toarray().flatten()
            if gene not in gene_dict.keys():
                gene_dict[gene] = {}
            
            # Append cell type colors to gene dictionary as well
            celltypes = rna_ad.obs["celltype"]
            colors = celltypes.map(define_color_dict()).to_list()  # Map colors to SEACell annnotations
            gene_dict[gene][dataset] = (rna_expr, atac_expr, colors)

    return gene_dict


def write_gene_csv(gene_dict: 'dict[str:dict[tuple[np.ndarray, np.ndarray, list[str]]]]', path: str):
    """Writes each gene to a csv file.

    Args:
        gene_dict (dict[str, dict[tuple[np.ndarray, np.ndarray, list[str]]]): The dictionary of genes and their values.
        path (str): The path to save the csv files.
    """

    if not os.path.exists('src/dyn/data/corr'):
        os.makedirs('src/dyn/data/corr')

    timepoints = ['TDR126', 'TDR127', 'TDR128', 'TDR118', 'TDR125', 'TDR124']
    for gene in list(gene_dict.keys()):
        gene = gene.replace('/', '-')  # a few genes have '/' in their name, who would do such a thing?

        # Extract values from gene dictionary
        rna_full, atac_full, colors_full, tp_full = [], [], [], []
        for tp in timepoints:
            try:  # Skip if gene is not present in a timepoint
                rna, atac, colors = gene_dict[gene][tp]
                colors = colors[:len(rna)]
            except KeyError:
                continue

            # Add each value in rna and atac array to rna_full and atac_full lists
            rna_full.extend(rna)
            atac_full.extend(atac)
            colors_full.extend(colors)
            tp_full.extend([tp] * len(rna))

        # Create a polars dataframe and save to file
        df = pl.DataFrame({
            'RNA': rna_full,
            'ATAC': atac_full,
            'Color': colors_full,
            'Timepoint': tp_full
            })
        df.write_csv(f'src/dyn/data/corr/{gene}.csv')


def main():
    """
    """

    # Get meta dictionary containing all datasets
    path = '/hpc/projects/data.science/yangjoon.kim/zebrahub_multiome/data/processed_data/05_SEACells_processed'
    datasets = ['TDR126', 'TDR127', 'TDR128', 'TDR118reseq', 'TDR125reseq', 'TDR124reseq']
    dict_meta = get_datasets(path, datasets)

    # Get gene dictionary
    gene_dict = get_gene_dict(dict_meta)
    pickle.dump(gene_dict, open("src/dyn/data/gene_dict.pkl", "wb"))

    # Write each gene to a csv file
    write_gene_csv(gene_dict, 'src/dyn/data/corr')

    # Write gene names to file
    with open('src/dyn/data/corr/gene_names.txt', 'w') as f:
        for gene in list(gene_dict.keys()):
            f.write(f"{gene}\n")


if __name__ == "__main__":
    main()
