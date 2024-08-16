"""Preprocess adata for each timepoint for plotting correlation scatterplots.

08/16/24    Ben Iovino    CZ-Biohub
"""

import scanpy as sc
import numpy as np
import pandas as pd
import pickle


def get_datasets(datasets: 'list[str]') -> 'dict[str:sc.AnnData]':
    """Returns a dictionary of adata objects for each dataset.

    Returns:
        dict[str, sc.AnnData]: The dictionary of adata objects.
    """

    # Import all datasets (adata) and save in a dictionary
    dict_meta_ad = {}
    path = "corr/data/"
    for data_id in datasets:

        # define the sample_id (withou "reseq" in the handle)
        sample_id = data_id.replace("reseq","")
    
        # import the single-cell adata's annotation (with adata.obs["SEACell"] annotation)
        df_seacells = pd.read_csv(f'{path}/{data_id}/{data_id}_seacells_obs_annotation_ML_coarse.csv', index_col=0)
    
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
        rna_meta_ad = sc.read_h5ad(path + f"{data_id}/{sample_id}_RNA_seacells_aggre.h5ad")
        atac_meta_ad = sc.read_h5ad(path + f"{data_id}/{sample_id}_ATAC_seacells_aggre.h5ad")
    
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


def main():
    """
    """

    # Get meta dictionary containing all datasets
    datasets = ['TDR126', 'TDR127', 'TDR128', 'TDR118reseq', 'TDR125reseq', 'TDR124reseq']
    dict_meta_ad = get_datasets(datasets)

    # Save to pickle
    with open('corr/data/dict_meta.pkl', 'wb') as f:
        pickle.dump(dict_meta_ad, f)


if __name__ == "__main__":
    main()
