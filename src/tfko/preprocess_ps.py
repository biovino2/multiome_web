"""Preprocess cell oracle data for calculating perturbation scores.

09/23/24    Ben Iovino    CZ-Biohub
"""

from anndata import AnnData
import celloracle as co
import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine
from tqdm import tqdm


def compute_metacell_transitions(adata: AnnData, trans_key="T_fwd_WT", metacell_key="SEACell"):
    """Returns dataframe of metacell-to-metacell transition probabilities.

    Args:
        adata (AnnData): Annotated data
        trans_key (str): Key for transition matrix
        metacell_key (str): Key for metacell assignment
    """

    # Get the cell-cell transition matrix
    T_cell = adata.obsp[trans_key]
    
    # Get metacell labels
    metacells = adata.obs[metacell_key]
    
    # Get unique metacells
    unique_metacells = metacells.unique()
    
    # Initialize the metacell transition matrix
    n_metacells = len(unique_metacells)
    T_metacell = np.zeros((n_metacells, n_metacells))
    
    # Create a mapping of metacell to cell indices
    metacell_to_indices = {mc: np.where(metacells == mc)[0] for mc in unique_metacells}
    
    # Compute metacell transitions
    for i, source_metacell in enumerate(unique_metacells):
        source_indices = metacell_to_indices[source_metacell]
        for j, target_metacell in enumerate(unique_metacells):
            target_indices = metacell_to_indices[target_metacell]
            
            # Extract the submatrix of transitions from source to target metacell
            submatrix = T_cell[source_indices][:, target_indices]
            
            # Sum all transitions and normalize by the number of source cells
            T_metacell[i, j] = submatrix.sum() / len(source_indices)
    
    # Create a DataFrame for easier interpretation
    T_metacell_df = pd.DataFrame(T_metacell, index=unique_metacells, columns=unique_metacells)
    
    # Normalize rows to sum to 1
    T_metacell_df = T_metacell_df.div(T_metacell_df.sum(axis=1), axis=0)
    
    return T_metacell_df


def compute_row_cosine_similarities(df_wt: pd.DataFrame, df_ko: pd.DataFrame) -> pd.Series:
    """
    Compute cosine similarities between corresponding rows of two dataframes.
    
    Args:
        df_wt (pd.DataFrame): Transition probability matrix for WT
        df_ko (pd.DataFrame): Transition probability matrix for KO
    
    Returns:
        pd.Series: Cosine similarities for each row
    """

    assert df_wt.index.equals(df_ko.index), "Dataframes must have the same index"
    assert df_wt.columns.equals(df_ko.columns), "Dataframes must have the same columns"
    
    similarities = {}
    for idx in df_wt.index:
        wt_row = df_wt.loc[idx].values
        ko_row = df_ko.loc[idx].values
        
        # Compute cosine similarity (1 - cosine distance)
        similarity = 1 - cosine(wt_row, ko_row)
        similarities[idx] = similarity
    
    return pd.Series(similarities, name="cos_sim")


def calc_cos_sims(path: str, list_datasets: 'list[str]'):
    """Calculates cosine similarity scores between wildtype and knockout transition probabilities.

    Args:
        path (str): Path to save the dataframes.
        list_datasets (list): List of datasets to process.
    """

    dict_cos_sims = {}
    oracle_path = "/hpc/projects/data.science/yangjoon.kim/zebrahub_multiome/data/processed_data/09_NMPs_subsetted_v2/"
    
    for data_id in tqdm(list_datasets, desc="Processing Datasets"):
    
        # load the oracle object and metacell info
        oracle = co.load_hdf5(f'{oracle_path}/14_{data_id}_in_silico_KO_trans_probs_added.celloracle.oracle')
        metacell = pd.read_csv(f'{oracle_path}/metacells/{data_id}_seacells_obs_manual_annotation_30cells.csv', index_col=0)
    
        # add the metacell information to the oracle.adata
        metacell_dict = metacell.SEACell.to_dict()
        oracle.adata.obs["SEACell"] = oracle.adata.obs_names.map(metacell_dict)
        oracle.adata.obs.head()
    
        # Calculate most prevalent cell type for each metacell
        most_prevalent = oracle.adata.obs.groupby("SEACell")["manual_annotation"].agg(lambda x: x.value_counts().idxmax())
        most_prevalent

        # average the 2D embedding and 2D transition vectors across "metacells"
        trans_probs_metacell_WT = compute_metacell_transitions(oracle.adata, 
                                                        trans_key="T_fwd_WT_global_nmps", 
                                                        metacell_key="SEACell")

        # Initialize an empty DataFrame with celltypes as the index
        metacells = trans_probs_metacell_WT.index
        cosine_sim_df = pd.DataFrame(index=metacells)
    
        # Compute cosine similarities for each gene knockout
        for gene in tqdm(oracle.active_regulatory_genes, desc=f"Processing Genes for {data_id}"):

            # Compute transition probabilities for the current gene knockout
            trans_key = f"T_fwd_{gene}_KO"
            trans_probs_metacell_KO = compute_metacell_transitions(oracle.adata, trans_key=trans_key, 
                                                                metacell_key="SEACell")
        
            # Compute cosine similarities between wildtype and knockout
            cosine_similarities = compute_row_cosine_similarities(trans_probs_metacell_WT, trans_probs_metacell_KO)
        
            # Add the cosine similarities as a new column to the DataFrame
            cosine_sim_df[gene] = cosine_similarities

        # Display the resulting DataFrame (metacell-by-genes)
        cosine_sim_df["celltype"] = cosine_sim_df.index.map(most_prevalent)

        # average the cosine similarities across cell types
        cosine_sim_df_avg = cosine_sim_df.groupby("celltype").mean()
        df_averaged = cosine_sim_df_avg.reset_index()
    
        # save the dataframes
        df_averaged.to_csv(f'{path}/{data_id}_sim.csv')
        dict_cos_sims[data_id] = df_averaged
        dict_cos_sims[f"{data_id}_metacells"] = cosine_sim_df
    
        print(f"{data_id} is completed")


def main():
    """
    """

    save_path = 'src/tfko/data'
    list_datasets = ['TDR118', 'TDR124','TDR125', 'TDR126', 'TDR127', 'TDR128']
    calc_cos_sims(save_path, list_datasets)


if __name__ == "__main__":
    main()
