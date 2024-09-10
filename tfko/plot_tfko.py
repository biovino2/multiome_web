"""Plot UMAPS at single and meta-cell level for a given time point and transcription factor knockout.

Ben Iovino  09/05/24    CZ-Biohub
"""

import anndata as ad
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import scanpy as sc
import streamlit as st
import numpy as np


def average_metacells(adata: ad.AnnData, knockout='WT') -> 'tuple[np.ndarray, np.ndarray]':
    """Returns the average UMAP position and transition vector of metacells based off each cell
    belonging to a metacell.

    Args:
        adata (ad.AnnData): The AnnData object.
        knockout (str): The transcription factor knockout (default is wildtype (WT)).

    Returns:
        tuple: A tuple containing the average UMAP position and transition vector of metacells.
    """

    if knockout != 'WT_global_nmps':
        knockout += '_KO'

    # Cell-level UMAP and 2D transition vectors
    X_umap = adata.obsm['X_umap_aligned']
    V_cell = adata.obsm[f'{knockout}_umap_aligned'] 
    
    # Convert metacell column to categorical if it's not already
    if not pd.api.types.is_categorical_dtype(adata.obs['SEACell']):
        metacells = pd.Categorical(adata.obs['SEACell'])
    else:
        metacells = adata.obs['SEACell']
    
    # X_metacell is the average UMAP position of the metacells
    # V_metacell is the average transition vector of the metacells
    n_metacells = len(metacells.categories)
    X_metacell = np.zeros((n_metacells, 2))
    V_metacell = np.zeros((n_metacells, 2))
    
    for i, category in enumerate(metacells.categories):
        mask = metacells == category
        X_metacell[i] = X_umap[mask].mean(axis =0)
        V_metacell[i] = V_cell[mask].mean(axis=0)
    
    return X_metacell, V_metacell


def plot_single_cells(cell_colors: 'dict[str, str]', umap_data: pd.DataFrame) -> go.Figure:
    """Returns a scatter plot containing single cells mapped in UMAP space.

    Args:
        cell_colors (dict): The color codes for cell types.
        umap_dat (pd.DataFrame): The df containing UMAP coordinates, SEACell, and cell type info.

    Returns:
        go.Figure: A plotly figure.
    """

    fig = go.Figure()
    for cell_type, color in cell_colors.items():
        filtered_data = umap_data[umap_data['celltype'] == cell_type]
        fig.add_trace(go.Scatter(
            x=filtered_data[0],
            y=filtered_data[1],
            mode='markers',
            name=cell_type,
            marker=dict(color=color, size=6, opacity=0.7),
            showlegend=False
        ))

    return fig


def plot_meta_cells(fig: go.Figure, cell_colors: 'dict[str, str]', umap_data: pd.DataFrame) -> go.Figure:
    """Returns a scatter plot containing meta cells mapped in UMAP space.

    Args:
        fig (go.Figure): The plotly figure object
        cell_colors (dict): The color codes for cell types.
        umap_dat (pd.DataFrame): The df containing UMAP coordinates, SEACell, and cell type info.

    Returns:
        go.Figure: A plotly figure.
    """

    # Prepare metacell data
    most_prevalent = umap_data.groupby('SEACell')['celltype'].agg(lambda x: x.value_counts().index[0])
    celltypes = umap_data['celltype']
    umap_data = umap_data.drop(['celltype'], axis=1)
    mcs = umap_data.groupby('SEACell').mean().reset_index()
    mcs['celltype'] = most_prevalent.values
    umap_data['celltype'] = celltypes

    # Plot metacells
    for cell_type, color in cell_colors.items():
        filtered_mcs = mcs[mcs['celltype'] == cell_type]
        fig.add_trace(go.Scatter(
            x=filtered_mcs[0],
            y=filtered_mcs[1],
            mode='markers',
            name=cell_type,
            marker=dict(color=color, size=10, line=dict(color='black', width=1.25)),
            showlegend=False
        ))

    return fig


def plot_trans_vecs(fig: go.Figure, X_metacell: np.ndarray, V_metacell: np.ndarray, ref: int) -> go.Figure:
    """Returns a scatter plot containing arrows representing the transition probabilities of
    metacells in UMAP space.

    Args:
        fig (go.Figure): The plotly figure object.
        X_metacell (np.ndarray): The average UMAP position of the metacells.
        V_metacell (np.ndarray): The average transition vector of the metacells.
        ref (int): The reference point.

    Returns:
        go.Figure: A plotly figure.
    """

    # Arrow lines
    for i in range(X_metacell.shape[0]):
        start_x = X_metacell[i, 0]
        start_y = X_metacell[i, 1]
    
        # Normalize the velocity vector similar to how quiver would do
        magnitude = np.sqrt(V_metacell[i, 0]**2 + V_metacell[i, 1]**2)
        scale_factor = 1 / 2
        end_x = start_x + (V_metacell[i, 0] * scale_factor / magnitude)
        end_y = start_y + (V_metacell[i, 1] * scale_factor / magnitude)
        
        # Add arrow heads
        fig.add_annotation(
            x=end_x, 
            y=end_y, 
            ax=start_x, 
            ay=start_y, 
            xref=f'x{ref}', 
            yref=f'y{ref}', 
            axref=f'x{ref}', 
            ayref=f'y{ref}',
            showarrow=True, 
            arrowhead=2,  # Adjust the type of arrowhead
            arrowsize=0.55,  # Adjust the size of the arrowhead
            arrowwidth=2,
            arrowcolor="black"
        )

    return fig


def plot_cells(timepoint: str, tf_names) -> go.Figure:
    """Returns a plotly figure with single cells and metacells, along with metacell transition
    probabilities given the TF knockout, plotted on UMAP coordinates

    Args:
        timepoint (str): The time point.
        knockout (str): The TF knockout.
        timepoints (dict): The dictionary of time points.

    Returns:
        go.Figure: A plotly figure.
    """
    
    cell_colors: 'dict[str, str]' = {
        'NMPs': '#8dd3c7',
        'PSM': '#008080',
        'fast_muscle': '#df4b9b',
        'neural_posterior': '#393b7f',
        'somites': '#1b9e77',
        'spinal_cord': '#d95f02',
        'tail_bud': '#7570b3'
    }

    # Obtain selected knockouts
    selected_knockouts = []
    for i, key in enumerate(st.session_state.selectboxes):
        st.sidebar.selectbox(
            'Select a knockout to plot',
            tf_names,
            key=key
        )
        selected_knockouts.append(st.session_state[key])

    # Create figure (subplots)
    fig = make_subplots(
        rows=len(selected_knockouts),
        cols=1,
        subplot_titles=[f"{ko}" for ko in selected_knockouts],
        shared_xaxes=True,
        shared_yaxes=True,
        )
    
    # Generate plot for each knockout
    for i, knockout in enumerate(selected_knockouts):

        # Load data
        adata = sc.read_h5ad(f"tfko/data/{timepoint}_KO.h5ad")
        adata.obs['SEACell'] = adata.obs['SEACell'].astype('object')
        adata.obs['manual_annotation'] = adata.obs['manual_annotation'].astype('object')
        X_metacell, V_metacell = average_metacells(adata, knockout)
    
        # Prepare data for plotting
        umap_coords = pd.DataFrame(adata.obsm['X_umap_aligned'], columns=[0, 1], index=adata.obs_names)
        umap_data = umap_coords.join(adata.obs[['SEACell', 'manual_annotation']])
        umap_data = umap_data.rename(columns={'manual_annotation': 'celltype'})

        # Plot each component
        plot = plot_single_cells(cell_colors, umap_data)
        plot = plot_meta_cells(plot, cell_colors, umap_data)
        plot = plot_trans_vecs(plot, X_metacell, V_metacell, i+1)

        # Add plot to figure
        for trace in plot.data:
            fig.add_trace(trace, row=i+1, col=1)
            fig.update_xaxes(row=i+1, col=1, matches='x',
                            showticklabels=False, showgrid=False, zeroline=False)
            fig.update_yaxes(row=i+1, col=1, matches='y',
                            showticklabels=False, showgrid=False, zeroline=False)
        for annotation in plot.layout.annotations:
            fig.add_annotation(annotation, row=i+1, col=1)

    fig.update_layout(height=1000, width=2000)

    return fig


def st_setup(timepoints: list[str]) -> 'tuple[str, str]':
    """Initializes streamlit session.

    Args:
        tf_names (list): The list of transcription factors to knockout.
    """

    st.set_page_config(layout="wide")
    st.title('Transcription Factor Knockout')
    st.write('')
    st.sidebar.markdown('# Settings')

    # Initialize drop down boxes
    if "selectboxes" not in st.session_state:
        st.session_state.selectboxes = [0]

   # Set up buttons for adding/removing timepoints
    with st.sidebar:
        add, reset = st.columns([1, 1])
        with add:
            if st.button('Add Knockout'):
                st.session_state.selectboxes.append(len(st.session_state.selectboxes))
        with reset:
            if st.button('Reset'):
                st.session_state.selectboxes = [0]

    # Get timepoint to show
    timepoint = st.sidebar.selectbox(
        'Select a time point to plot',
        timepoints
    )
    
    return timepoint


def main():
    """
    """

    with open('tfko/data/common_genes.txt', 'r') as file:
        tf_names = file.read().splitlines()
    tf_names.append('WT_global_nmps')

    # Set up streamlit
    timepoint = 'TDR125'
    timepoints = {'10 hours post fertilization': 'TDR126',
                '12 hours post fertilization': 'TDR127',
                '14 hours post fertilization': 'TDR128',
                '16 hours post fertilization': 'TDR118',
                '19 hours post fertilization': 'TDR125',
                '24 hours post fertilization': 'TDR124'}
    timepoint = st_setup(list(timepoints.keys()))

    # Generate plot
    fig = plot_cells(timepoints[timepoint], tf_names)    
    st.plotly_chart(fig)


if __name__ == '__main__':
    main()
