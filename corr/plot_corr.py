"""Plot correlation scatterplots between ATAC and RNA data.

Ben Iovino  08/08/24    CZ-Biohub
"""

import matplotlib.pyplot as plt
import pickle as pkl
import numpy as np
import scanpy as sc
from scipy.stats import pearsonr
import streamlit as st


def plot_genes(gene: str, gene_dict: 'dict[str:sc.AnnData]') -> plt.Figure:
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

        expr_rna = gene_dict[gene][sample_id][0]
        expr_atac = gene_dict[gene][sample_id][1]
        colors = gene_dict[gene][sample_id][2]
        
        # if either exp1 or expr2 are "constant", then
        if np.all(expr_rna == expr_rna[0]) or np.all(expr_atac== expr_atac[0]):
            correlation = np.nan
        else:
            correlation, _ = pearsonr(expr_atac, expr_rna)

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

    gene_dict = pkl.load(open("corr/data/gene_dict.pkl", "rb"))
    gene_names = list(gene_dict.keys())

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
        fig = plot_genes(gene, gene_dict)
        st.pyplot(fig)

if __name__ == "__main__":
    main()
