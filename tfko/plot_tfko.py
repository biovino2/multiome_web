"""Plot UMAPS at single and meta-cell level for a given time point and transcription factor knockout.

Ben Iovino  09/05/24    CZ-Biohub
"""

import anndata as ad
import plotly.graph_objects as go
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


def plot_cells(adata: ad.AnnData, X_metacell: np.ndarray,
                V_metacell: np.ndarray, timepoint: str, knockout: str) -> go.Figure:
    """Returns a plotly figure with single cells and metacells, along with metacell transition
    probabilities given the TF knockout, plotted on UMAP coordinates

    Args:
        adata (ad.AnnData): The AnnData object.
        X_metacell (np.ndarray): The average UMAP position of the metacells.
        V_metacell (np.ndarray): The average transition vector of the metacells.
        timepoint (str): The time point.
        knockout (str): The TF knockout.

    Returns:
        go.Figure: A plotly figure.
    """
    
    cell_type_color_dict = {
        'NMPs': '#8dd3c7',
        'PSM': '#008080',
        'fast_muscle': '#df4b9b',
        'neural_posterior': '#393b7f',
        'somites': '#1b9e77',
        'spinal_cord': '#d95f02',
        'tail_bud': '#7570b3'
    }
    
    # Prepare data for plotting
    umap_coords = pd.DataFrame(adata.obsm['X_umap_aligned'], columns=[0, 1], index=adata.obs_names)
    umap_data = umap_coords.join(adata.obs[['SEACell', 'manual_annotation']])
    umap_data = umap_data.rename(columns={'manual_annotation': 'celltype'})

    # Compute the most prevalent cell type in each metacell
    most_prevalent = umap_data.groupby('SEACell')['celltype'].agg(lambda x: x.value_counts().index[0])

    # Prepare metacell data
    celltypes = umap_data['celltype']
    umap_data = umap_data.drop(['celltype'], axis=1)
    mcs = umap_data.groupby('SEACell').mean().reset_index()
    mcs['celltype'] = most_prevalent.values
    umap_data['celltype'] = celltypes

    # Create a Plotly figure
    fig = go.Figure()

    # Plot single cells
    for cell_type, color in cell_type_color_dict.items():
        filtered_data = umap_data[umap_data['celltype'] == cell_type]
        fig.add_trace(go.Scatter(
            x=filtered_data[0],
            y=filtered_data[1],
            mode='markers',
            name=cell_type,
            marker=dict(color=color, size=6, opacity=0.7),
            showlegend=False
        ))

    # Plot metacells
    for cell_type, color in cell_type_color_dict.items():
        filtered_mcs = mcs[mcs['celltype'] == cell_type]
        fig.add_trace(go.Scatter(
            x=filtered_mcs[0],
            y=filtered_mcs[1],
            mode='markers',
            name=cell_type,
            marker=dict(color=color, size=10, line=dict(color='black', width=1.25)),
            showlegend=False
        ))

    # Plot transition vectors (arrows with arrowheads)
    for i in range(X_metacell.shape[0]):
        start_x = X_metacell[i, 0]
        start_y = X_metacell[i, 1]
        end_x = start_x + V_metacell[i, 0] * 15
        end_y = start_y + V_metacell[i, 1] * 15
        
        fig.add_annotation(
            x=end_x, 
            y=end_y, 
            ax=start_x, 
            ay=start_y, 
            xref="x", 
            yref="y", 
            axref="x", 
            ayref="y",
            showarrow=True, 
            arrowhead=2,  # Adjust the type of arrowhead
            arrowsize=0.55,  # Adjust the size of the arrowhead
            arrowwidth=2,
            arrowcolor="black"
        )

    # Customize layout
    fig.update_layout(
        width=800, height=600,
        title=f"{knockout} KO at {timepoint}",
        xaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
        yaxis=dict(showgrid=False, showticklabels=False, zeroline=False)
    )

    return fig


def st_setup(timepoints: list[str], tf_names: list[str]) -> 'tuple[str, str]':
    """Initializes streamlit session.

    Args:
        tf_names (list): The list of transcription factors to knockout.
    """

    st.set_page_config(layout="wide")
    st.title('Transcription Factor Knockout')
    st.write('')
    st.sidebar.markdown('# Settings')

    # Get TF to knockout
    default = tf_names.index('WT_global_nmps')
    trans_fac = st.sidebar.selectbox(
        'Select a transcription factor to knockout',
        tf_names,
        index=default
    )

    # Get timepoint to show
    # Get celltype
    timepoint = st.sidebar.selectbox(
        'Select a time point to plot',
        timepoints
    )
    
    return timepoint, trans_fac


def main():
    """
    """

    with open('tfko/data/common_genes.txt', 'r') as file:
        tf_names = file.read().splitlines()
    tf_names.append('WT_global_nmps')

    timepoint = 'TDR125'
    timepoints = {'10 hours post fertilization': 'TDR126',
                '12 hours post fertilization': 'TDR127',
                '14 hours post fertilization': 'TDR128',
                '16 hours post fertilization': 'TDR118',
                '19 hours post fertilization': 'TDR125',
                '24 hours post fertilization': 'TDR124'}
    
    timepoint, trans_fac = st_setup(list(timepoints.keys()), tf_names)

    # Load data
    adata = sc.read_h5ad(f"tfko/data/{timepoints[timepoint]}_KO.h5ad")
    adata.obs['SEACell'] = adata.obs['SEACell'].astype('object')
    adata.obs['manual_annotation'] = adata.obs['manual_annotation'].astype('object')

    # Calculate metacell averages and generate plot
    X_metacell, V_metacell = average_metacells(adata, trans_fac)
    fig = plot_cells(adata, X_metacell, V_metacell, timepoint, trans_fac)    
    st.plotly_chart(fig)


if __name__ == '__main__':
    main()
