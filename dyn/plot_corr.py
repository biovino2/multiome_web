"""Plot correlation scatterplots between ATAC and RNA data.

Ben Iovino  08/08/24    CZ-Biohub
"""

import matplotlib.pyplot as plt
import os
import pickle as pkl
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import scanpy as sc
from scipy.stats import pearsonr
import streamlit as st
from preprocess_corr import define_color_dict


def st_setup() -> str:
    """Initializes streamlit session.

    Returns:
        str: The selected cell type.
    """

    st.set_page_config(layout="wide")
    st.sidebar.markdown('# Settings')
    st.title('Gene Activity vs. Expression')
    st.write('For each time point, we plot plot time-resolved scatter plots of gene activity, which \
              represents the signal from the peaks within the gene body, and gene expression, which \
             represents the normalized transcript counts. Each point represents a metacell (SEACell), \
              colored by cell type. The correlation between gene activity and expression is also \
              calculated.')
    
    # Remove extra space at top of the page
    st.markdown(
    """
    <style>
        /* Remove padding from main block */
        .block-container {
            padding-top: 2rem;
        }
    </style>
    """,
    unsafe_allow_html=True
    )

    # Initialize drop down boxes
    if "selectboxes" not in st.session_state:
        st.session_state.selectboxes = [0]

    # Set up buttons for adding/removing genes
    with st.sidebar:
        add, reset = st.columns([1, 1])
        with add:
            if st.button('Add Gene'):
                st.session_state.selectboxes.append(len(st.session_state.selectboxes))
        with reset:
            if st.button('Reset'):
                st.session_state.selectboxes = [0]

    # Create selectbox for cell type selection
    celltypes = list(define_color_dict().keys())
    celltypes.insert(0, 'All')
    celltype = st.sidebar.selectbox(
        'Isolate a cell type',
        celltypes,
        index=celltypes.index('All')
    )

    return celltype


def save_config() -> dict:
    """Returns a config to save plotly figure as SVG.

    Returns:
        config (dict): The configuration.
    """

    config = {
    'toImageButtonOptions': {
        'format': 'svg', # one of png, svg, jpeg, webp
        'filename': 'correlation_scatterplot',
        'height': None,
        'width': None,
        'scale': 1 # Multiply title/legend/axis/canvas sizes by this factor
    },
    'displayModeBar': True
    }

    return config


def subset_data(gene_dict: dict, gene: str, sample_id: str, celltype: str) -> 'tuple[np.ndarray, np.ndarray, list[str], list[str]]':
    """Returns gene expression and activity values as arrays, as well as colors for each cell type.

    Args:
        gene_dict (dict): The dictionary of genes and their values for each time point.
        gene (str): The gene to plot.
        sample_id (str): The timepoint.
        celltype (str): The selected cell type.

    Returns:
        tuple (np.ndarray, np.ndarray, list[str], list[str]): The gene expression values, activity values, colors, and cell type.
    """

    # Get data for timepoint
    expr_rna = gene_dict[gene][sample_id][0]
    expr_atac = gene_dict[gene][sample_id][1]
    colors = gene_dict[gene][sample_id][2]

    # Map colors to cell types
    color_dict = define_color_dict()
    inv_color_dict = {v: k for k, v in color_dict.items()}
    cell_colors = [inv_color_dict[color] for color in colors]

    # Isolate points related to celltype
    if celltype != 'All':
        indices = [i for i, cell in enumerate(cell_colors) if cell == celltype]
        expr_rna = expr_rna[indices]
        expr_atac = expr_atac[indices]
        colors = [colors[i] for i in indices]
        cell_colors = [cell_colors[i] for i in indices]

    return expr_rna, expr_atac, colors, cell_colors


def calc_corr(expr_rna: np.ndarray, expr_atac: np.ndarray) -> float:
    """Returns the correlation between two arrays.

    Args:
        expr_rna (np.ndarray): RNA expression values.
        expr_atac (np.ndarray): ATAC expression values.

    Returns:
        float: The correlation value.
    """

    if len(expr_rna) == 0 or len(expr_atac) == 0:
        return np.nan

    # Calculate correlation
    if np.all(expr_rna == expr_rna[0]) or np.all(expr_atac== expr_atac[0]):
        return np.nan
    else:
        correlation, _ = pearsonr(expr_rna, expr_atac)
        return correlation
    

def get_margins(expr_rna: np.ndarray, expr_atac: np.ndarray, rna_max: float, atac_max: float) -> 'tuple[float, float]':
    """Returns maximum values for RNA and ATAC expression.

    Args:
        expr_rna (np.ndarray): RNA expression values.
        expr_atac (np.ndarray): ATAC expression values.
        rna_max (float): The current maximum RNA value.
        atac_max (float): The current maximum ATAC value.

    Returns:
        tuple (float, float): The updated maximum values.
    """

    # No values if no cells are present at timepoint
    if len(expr_rna) == 0 or len(expr_atac) == 0:
        return rna_max, atac_max
    
    # Update margins of the plot
    if np.max(expr_rna) > rna_max:
        rna_max = np.max(expr_rna)
    if np.max(expr_atac) > atac_max:
        atac_max = np.max(expr_atac)

    return rna_max, atac_max


def plot_genes(celltype: str, gene: str, gene_dict: 'dict[str:sc.AnnData]') -> plt.Figure:
    """Returns one figure containing scatter plots for each gene.

    Args:
        celltype (str): The selected cell type.
        genes (str): Gene to plot
        dict_meta (dict[str:sc.AnnData]): The dictionary of adata objects.

    Returns:
        plt.Figure: The figure object.
    """

    timepoints = {'TDR126': '10 hours post fertilization',
                    'TDR127': '12 hours post fertilization',
                    'TDR128': '14 hours post fertilization',
                    'TDR118': '16 hours post fertilization',
                    'TDR125': '19 hours post fertilization',
                    'TDR124': '24 hours post fertilization'}
    fig = make_subplots(rows=2, cols=3, subplot_titles=list(timepoints.values()))

    # Loop over all timepoints
    corr_scores = []
    rna_max, atac_max = 0, 0
    for index, sample_id in enumerate(list(timepoints.keys())):

        # Plot scatter plot with plotly
        expr_rna, expr_atac, colors, type_colors = subset_data(gene_dict, gene, sample_id, celltype)
        scatter = go.Scatter(
                        x=expr_atac,
                        y=expr_rna,
                        mode='markers',
                        marker=dict(color=colors),
                        hoverinfo='text',
                        text=type_colors,
                    )
        
        # Calculate correlation and update margins with max rna/atac values
        corr_scores.append(calc_corr(expr_rna, expr_atac))
        rna_max, atac_max = get_margins(expr_rna, expr_atac, rna_max, atac_max)

        # Add go.Scatter object to figure
        fig.add_trace(scatter, row=(index//3)+1, col=(index%3)+1)

    # Update layout
    fig.update_layout(
        margin=dict(l=10, r=10, t=70, b=0),
        showlegend=False,
        xaxis=dict(title_standoff=5)
    )

    # Update titles
    for i, annotation in enumerate(fig.layout.annotations):
        annotation.text = f"{list(timepoints.values())[i]}<br>Correlation: {corr_scores[i]:.2f}"
        annotation.font.size = 14

    # Update axes labels and titles
    for i in range (0, len(corr_scores)):
        fig.update_xaxes(title='ATAC', row=i//3+1, col=i%3+1, title_font=dict(size=10), range=[0, atac_max], title_standoff=4)
        fig.update_yaxes(title='RNA', row=i//3+1, col=i%3+1, title_font=dict(size=10), range=[0, rna_max], title_standoff=8)

    return fig


# Main

# Load data and set up streamlit
path = os.path.dirname(os.path.abspath(__file__))+'/data'
gene_dict = pkl.load(open(f"{path}/gene_dict.pkl", "rb"))
gene_names = list(gene_dict.keys())
celltype = st_setup()

# Create selectboxes
selected_genes = []
default = gene_names.index('slc4a1a')
for i, key in enumerate(st.session_state.selectboxes):
    st.sidebar.selectbox(
        f'Select gene {i+1}',
        gene_names,
        key=key,
        index=default
    )
selected_genes.append(st.session_state[key])

# Plot each figure
for gene in selected_genes:
    fig = plot_genes(celltype, gene, gene_dict)
    st.markdown(f'### {gene}')
    st.plotly_chart(fig, config=save_config())
