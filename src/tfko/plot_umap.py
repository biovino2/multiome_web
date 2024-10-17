"""Plot UMAPS at single and meta-cell level for a given time point and transcription factor knockout.

Ben Iovino  09/05/24    CZ-Biohub
"""

import os
import plotly.graph_objects as go
import pandas as pd
import scanpy as sc
import streamlit as st
import numpy as np
from util import get_datasets


def map_ko(name: str) -> str:
    """Changes name of control. Lazy solution to avoid changing all the names in preprocessing.

    Args:
        name (str): The name of the knockout.

    Returns:
        str: The new knockout name.
    """

    if name == 'WT_global_nmps':
        return 'Control (No Knockout)'
    if name == 'Control (No Knockout)':
        return 'WT_global_nmps'
    return name


def get_hover_text(cell_type: str, color: str, num: int) -> str:
    """Returns HTML hover text for a cell type.

    Args:
        cell_type (str): The cell type.
        color (str): The color of the cell type.
        num (int): The number of cells of that type.
    """

    hover_texts = [
        f'<span style="color:{color};">{cell_type}</span>' 
        for _ in range(num)
    ]

    return hover_texts


def plot_single_cells(fig: go.Figure, cell_colors: 'dict[str, str]', umap_data: pd.DataFrame) -> go.Figure:
    """Returns a scatter plot containing single cells mapped in UMAP space.

    Args:
        fig (go.Figure): The plotly figure object.
        cell_colors (dict): The color codes for cell types.
        umap_dat (pd.DataFrame): The df containing UMAP coordinates, SEACell, and cell type info.

    Returns:
        go.Figure: A plotly figure.
    """

    for cell_type, color in cell_colors.items():
        filtered_data = umap_data[umap_data['celltype'] == cell_type]
        hover_texts = get_hover_text(cell_type, color, len(filtered_data))

        # Add all cells of same type to plot
        fig.add_trace(go.Scatter(
            x=filtered_data[0],
            y=filtered_data[1],
            mode='markers',
            name=cell_type,
            marker=dict(color=color, size=6, opacity=0.7),
            showlegend=False,
            customdata=hover_texts,  # Pass the custom HTML hover text
            hovertemplate='%{customdata}<extra></extra>'

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
        hover_texts = get_hover_text(cell_type, color, len(filtered_mcs))

        # Add all metacells of same type to plot
        fig.add_trace(go.Scatter(
            x=filtered_mcs[0],
            y=filtered_mcs[1],
            mode='markers',
            name=cell_type,
            marker=dict(color=color, size=10, line=dict(color='black', width=1.25)),
            showlegend=False,
            customdata=hover_texts,  # Pass the custom HTML hover text
            hovertemplate='%{customdata}<extra></extra>'  # Display hover text with HTML styling
        ))

    return fig


def plot_trans_vecs(fig: go.Figure, X_metacell: np.ndarray, V_metacell: np.ndarray, color: str) -> go.Figure:
    """Returns a scatter plot containing arrows representing the transition probabilities of
    metacells in UMAP space.

    Args:
        fig (go.Figure): The plotly figure object.
        X_metacell (np.ndarray): The average UMAP position of the metacells.
        V_metacell (np.ndarray): The average transition vector of the metacells.
        color (str): The color of the arrows.

    Returns:
        go.Figure: A plotly figure.
    """

    lines_x = []
    lines_y = []
    arrowheads_x = []
    arrowheads_y = []

    # Arrow lines
    for i in range(X_metacell.shape[0]):
        start_x = X_metacell[i, 0]
        start_y = X_metacell[i, 1]
    
        # Calculate end point of the line
        magnitude = np.sqrt(V_metacell[i, 0]**2 + V_metacell[i, 1]**2)
        end_x = start_x + (V_metacell[i, 0] * 0.5 / magnitude)
        end_y = start_y + (V_metacell[i, 1] * 0.5 / magnitude)
    
        # Line between start and end, None to break the line
        lines_x.extend([start_x, end_x, None])
        lines_y.extend([start_y, end_y, None])
    
        # Arrowhead direction
        arrowhead_length = 0.05
        arrowhead_angle = np.pi / 5  # 30 degrees for arrowhead angle

        # Calculate angle of the line
        angle = np.arctan2(end_y - start_y, end_x - start_x)
    
        # Create two arrowhead points (left and right)
        left_x = end_x - arrowhead_length * np.cos(angle - arrowhead_angle)
        left_y = end_y - arrowhead_length * np.sin(angle - arrowhead_angle)
        right_x = end_x - arrowhead_length * np.cos(angle + arrowhead_angle)
        right_y = end_y - arrowhead_length * np.sin(angle + arrowhead_angle)
    
        # Use None to break line between each arrowhead
        arrowheads_x.extend([end_x, left_x, None, end_x, right_x, None])
        arrowheads_y.extend([end_y, left_y, None, end_y, right_y, None])

    # Add the lines as scattergl
    fig.add_trace(go.Scattergl(
        x=lines_x, 
        y=lines_y,
        mode='lines',
        line=dict(color=color, width=2),
        hoverinfo='none',
        showlegend=False
    ))
    
    # Add arrowheads as scattergl
    fig.add_trace(go.Scattergl(
        x=arrowheads_x,
        y=arrowheads_y,
        mode='lines',
        line=dict(color=color, width=2),
        hoverinfo='none',
        showlegend=False
    ))

    return fig


def plot_ct_legend(fig: go.Figure, cell_colors: 'dict[str, str]') -> go.Figure:
    """Returns figure with legend for cell type colors.

    Args:
        fig (go.Figure): The plotly figure object.
        cell_colors (dict): The color codes for cell types.
    """

    for cell_type, color in cell_colors.items():  # Legend for cell types
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(color=color, size=10),
            showlegend=True,
            legendgroup='group1',
            name=cell_type
        ))

    return fig


def plot_arr_legend(fig: go.Figure, arrow_colors: 'dict[str, str]', num_ko: int, knockout: str) -> go.Figure:
    """Returns figure with legend for arrow colors.

    Args:
        fig (go.Figure): The plotly figure object.
        arrow_colors (dict): The color codes for arrows.
        num_ko (int): The number of knockouts plotted.
        knockout (str): The knockout plotted.

    Returns:
        go.Figure: A plotly figure.
    """

    fig.add_trace(go.Scatter(  # Legend for arrow colors
        x=[None],
        y=[None],
        mode='lines+markers',
        marker=dict(color=arrow_colors[num_ko], size=10, symbol='triangle-right'),
        line=dict(width=2),
        showlegend=True,
        legendgroup='group2',
        name=knockout
    ))

    return fig


def plot_cells(path: str, knockouts: list[str], timepoint: str) -> go.Figure:
    """Returns a plotly figure with single cells and metacells, along with metacell transition
    probabilities given the TF knockout, plotted on UMAP coordinates

    Args:
        path (str): The path to the data.
        knockouts (list): The list of TF knockouts.
        timepoint (str): The time point.

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

    arrow_colors: 'dict[str, str]' = {
        0: '#000000',
        1: '#FF0000',
        2: '#FFFF00',
        3: '#00FFFF',
        4: '#FF00FF',
    }

    # Load cell/metacell data
    adata = sc.read_h5ad(f"{path}/{timepoint}_KO.h5ad")
    adata.obs['SEACell'] = adata.obs['SEACell'].astype('object')
    adata.obs['manual_annotation'] = adata.obs['manual_annotation'].astype('object')

    # Prepare data for plotting
    umap_coords = pd.DataFrame(adata.obsm['X_umap_aligned'], columns=[0, 1], index=adata.obs_names)
    umap_data = umap_coords.join(adata.obs[['SEACell', 'manual_annotation']])
    umap_data = umap_data.rename(columns={'manual_annotation': 'celltype'})

    # Plot cells and metacells
    fig = go.Figure()
    fig = plot_single_cells(fig, cell_colors, umap_data)
    fig = plot_meta_cells(fig, cell_colors, umap_data)
    fig = plot_ct_legend(fig, cell_colors)

    # Plot transition vectors for each knockout
    for i, ko in enumerate(knockouts):
        ko = map_ko(ko)

        # Load and plot transition vectors
        X_metacell, V_metacell = np.load(f"{path}/metacells/{timepoint}_{ko}_metacells.npz").values()
        fig = plot_trans_vecs(fig, X_metacell, V_metacell, arrow_colors[i])
        fig = plot_arr_legend(fig, arrow_colors, i, map_ko(ko))

    # Update layout
    fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False)
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False)
    fig.update_layout(height=550, width=1000, margin=dict(l=10, r=10, t=20, b=0))

    return fig


def st_setup(timepoints: list[str]) -> str:
    """Initializes streamlit session.

    Args:
        timepoints (list): The list of time points.

    Returns:
        str: The time point.
    """

    st.set_page_config(layout="wide")
    st.markdown("# *in silico* Genetic Perturbation of Transcription Factors")
    st.write('For each time point, we visualize the mesodermal and neuro-ectodermal cells in UMAP space, with each cell colored according to its cell type.\
            We also display metacells, computed using SEACells (Persad et al., 2023), as larger points, \
            colored by their most prevalent cell type. The arrows on metacells represent 2D-projected transition probabilities, \
            which are dependent on the transcription factor knockout being modeled.')
    st.sidebar.markdown('# Settings')

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

   # Set up buttons for adding/removing knockouts
    with st.sidebar:
        add, reset = st.columns([1, 1])
        with add:
            if st.button('Add Knockout'):
                st.session_state.selectboxes.append(len(st.session_state.selectboxes))
            if len(st.session_state.selectboxes) > 5:  # Max 5 knockouts
                st.session_state.selectboxes.pop()
        with reset:
            if st.button('Reset'):
                st.session_state.selectboxes = [0]

    # Get timepoint to show
    timepoint = st.sidebar.selectbox(
        'Select time point',
        timepoints
    )
    
    return timepoint


### Main

# Get names of knockouts
path = os.path.dirname(os.path.abspath(__file__))+'/data'
with open(f'{path}/common_tfs.txt', 'r') as file:
    tf_names = file.read().splitlines()
tf_names.insert(0, 'Control (No Knockout)')

# Set up streamlit
datasets = get_datasets()
timepoint = st_setup(list(datasets.keys()))
st.markdown(f'### {timepoint}')

# Generate plot based on selected knockouts
selected_knockouts = []
default = 'Control (No Knockout)'
for i, key in enumerate(st.session_state.selectboxes):
    st.sidebar.selectbox(
        f'Select knockout {i+1}',
        tf_names,
        index=tf_names.index(default),
        key=key
    )
    selected_knockouts.append(st.session_state[key])
fig = plot_cells(path, selected_knockouts, datasets[timepoint])
st.plotly_chart(fig, use_container_width=True)
