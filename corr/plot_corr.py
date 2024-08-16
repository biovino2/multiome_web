"""Plot correlation scatterplots between ATAC and RNA data.

Ben Iovino  08/08/24    CZ-Biohub
"""

import matplotlib.pyplot as plt
import pickle as pkl
import numpy as np
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


def plot_genes(gene: str, dict_meta: 'dict[str:sc.AnnData]') -> plt.Figure:
    """Returns one figure containing scatter plots for each gene.

    Args:
        genes (str): Gene to plot
        dict_meta (dict[str:sc.AnnData]): The dictionary of adata objects.

    Returns:
        plt.Figure: The figure object.
    """

    timepoints = {'TDR126': '0 hours post fertilization',
                    'TDR127': '5 hours post fertilization',
                    'TDR128': '10 hours post fertilization',
                    'TDR118': '15 hours post fertilization',
                    'TDR125': '20 hours post fertilization',
                    'TDR124': '30 hours fertilization'}
    fig, axs = plt.subplots(1, 6, figsize=(30,5))

    # Loop over all timepoints
    for index, sample_id in enumerate(list(timepoints.keys())):
            
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
        axs[index].set_title(f"{timepoints[sample_id]} (r={correlation:.2f})")
        axs[index].set_xlabel('ATAC expression')
        axs[index].set_ylabel('RNA expression')
        axs[index].grid(False)

    plt.suptitle(f"Scatter plots for {gene}", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to fit the main title

    return fig


def main():
    """
    """

    dict_meta = pkl.load(open('corr/data/dict_meta.pkl', 'rb'))
    
    # Get all gene names
    gene_names = set()
    for key in dict_meta.keys():
        gene_names.update(dict_meta[key].var_names)

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
    for key in st.session_state.selectboxes:
        st.sidebar.selectbox(
            'Select a gene to plot',
            gene_names,
            key=key
        )
        selected_genes.append(st.session_state[key])
    
    # Sidebar with selectboxes
    for gene in selected_genes:
        fig = plot_genes(gene, dict_meta)
        st.pyplot(fig)

if __name__ == "__main__":
    main()
