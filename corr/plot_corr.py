"""Plot correlation scatterplots between ATAC and RNA data.

Ben Iovino  08/08/24    CZ-Biohub
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr
import streamlit as st


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
    

    rna_meta_ad_filtered = rna_meta_ad[:,rna_meta_ad.var_names.isin(atac_meta_ad.var_names)]
    atac_meta_ad_filtered = atac_meta_ad[:,atac_meta_ad.var_names.isin(rna_meta_ad.var_names)]

    return dict_meta_ad, rna_meta_ad_filtered, atac_meta_ad_filtered


def plot_genes(genes: 'list[str]', dict_meta: 'dict[str:sc.AnnData]', rna_meta, atac_meta) -> plt.Figure:
    """Returns one figure containing scatter plots for each gene.

    Args:
        genes (list[str]): The list of genes to plot.
        dict_meta (dict[str:sc.AnnData]): The dictionary of adata objects.

    Returns:
        plt.Figure: The figure object.
    """

    list_timepoints = ['TDR126', 'TDR127', 'TDR128', 'TDR118', 'TDR125', 'TDR124']

    # for each gene, generate the time-resolved scatter plots for the RNA/ATAC correlation
    for gene in genes:
        fig, axs = plt.subplots(1, 6, figsize=(30,5))

        # Loop over all timepoints
        for index, sample_id in enumerate(list_timepoints):
            
            # extract the rna_ad and atac_ad
            rna_ad = dict_meta[f"{sample_id}_rna"]
            atac_ad = dict_meta[f"{sample_id}_atac"]
        
            # compute the correlation between RNA and ATAC
            expr_rna = rna_ad[:, gene].X.toarray().flatten()
            expr_atac = atac_ad[:, gene].X.toarray().flatten()
        
            # if either exp1 or expr2 are "constant", then
            if np.all(expr_rna == expr_rna[0]) or np.all(expr_atac== expr_atac[0]):
                correlation = np.nan
            else:
                correlation, _ = pearsonr(expr_atac, expr_rna)
        
            # Get the cell type for each cell
            cell_types = rna_ad.obs["celltype"]
        
            # Map the cell type to colors
            colors = cell_types.map(define_color_dict())
        
            # Plot the scatter plot for the current timepoint
            axs[index].scatter(expr_atac, expr_rna, color=colors)
            axs[index].set_title(f"{sample_id} (r={correlation:.2f})")
            axs[index].set_xlabel('ATAC expression')
            axs[index].set_ylabel('RNA expression')
            axs[index].grid(False)

        plt.suptitle(f"Scatter plots for {gene}", fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to fit the main title

    return fig


def add_selectbox():
    st.session_state.selectboxes.append(len(st.session_state.selectboxes))


def main():
    """
    """

    # define the list of datasets (data_id)
    list_datasets = ['TDR126', 'TDR127', 'TDR128',
                'TDR118reseq', 'TDR125reseq', 'TDR124reseq']
    dict_meta, rna_meta, atac_meta = get_datasets(list_datasets)
    gene_names = rna_meta.var_names

    # Page setup
    st.set_page_config(layout="wide")
    st.sidebar.markdown('# Settings')
    st.title('Meta-cell ATAC-RNA correlation')
    st.write('For each gene, we plot plot a time-resolved scatter plot of ATAC and RNA expression.')
    if "selectboxes" not in st.session_state:
        st.session_state.selectboxes = [0]

    if st.sidebar.button("Add gene"):
        st.session_state.selectboxes.append(len(st.session_state.selectboxes))

    selected_genes = []
    for i, key in enumerate(st.session_state.selectboxes):
        st.sidebar.selectbox(
            'Select a gene to plot',
            gene_names,
            key=key
        )
        selected_genes.append(st.session_state[key])
    
    # Sidebar with selectboxes
    for gene in selected_genes:
        fig = plot_genes([gene], dict_meta, rna_meta, atac_meta)
        st.pyplot(fig)


if __name__ == "__main__":
    main()
