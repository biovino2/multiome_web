"""Plot UMAPS at single and meta-cell level for a given time point and transcription factor knockout.

Ben Iovino  09/05/24    CZ-Biohub
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import scanpy as sc
import streamlit as st
import numpy as np


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
        line=dict(color='black', width=2),
        showlegend=False
    ))
    
    # Add arrowheads as scattergl
    fig.add_trace(go.Scattergl(
        x=arrowheads_x,
        y=arrowheads_y,
        mode='lines',
        line=dict(color='black', width=2),
        showlegend=False
    ))

    return fig


def plot_cells(knockouts: list[str], timepoint: str) -> go.Figure:
    """Returns a plotly figure with single cells and metacells, along with metacell transition
    probabilities given the TF knockout, plotted on UMAP coordinates

    Args:
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

    # Create figure (subplots)
    fig = make_subplots(
        rows=len(knockouts),
        cols=1,
        subplot_titles=[f"{ko}" for ko in knockouts],
        shared_xaxes=True,
        shared_yaxes=True,
        vertical_spacing=0,
        )
    
    # Generate plot for each knockout
    for i, knockout in enumerate(knockouts):

        # Map control to WT_global_nmps (pretty hacky)
        if knockout == 'Control (No Knockout)':
            knockout = 'WT_global_nmps'

        # Load data
        adata = sc.read_h5ad(f"tfko/data/{timepoint}_KO.h5ad")
        adata.obs['SEACell'] = adata.obs['SEACell'].astype('object')
        adata.obs['manual_annotation'] = adata.obs['manual_annotation'].astype('object')
        X_metacell, V_metacell = np.load(f"tfko/data/metacells/{timepoint}_{knockout}_metacells.npz").values()
    
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

    # Scale height and width with number of knockouts
    height = 600 * len(knockouts)
    fig.update_layout(height=height, width=1000, margin=dict(l=10, r=10, t=20, b=0))

    return fig


def st_setup(timepoints: list[str]) -> str:
    """Initializes streamlit session.

    Args:
        timepoints (list): The list of time points.

    Returns:
        str: The time point.
    """

    st.set_page_config(layout="wide")
    st.title('In-Silico Transcription Factor Knockout')
    st.write('For any time point (hours post fertilization), we plot the mesodermal and neuro-ectodermal cells in UMAP space,' \
            ' colored by their cell type. We also plot the metacells (SEACells), depicted as larger points, colored by their most' \
            ' prevalent cell type. The arrows represent the transition probabilities of the metacells given the transcription factor' \
            ' knockout')
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
        with reset:
            if st.button('Reset'):
                st.session_state.selectboxes = [0]

    # Get timepoint to show
    timepoint = st.sidebar.selectbox(
        'Select a time point to plot',
        timepoints
    )
    
    return timepoint


# Main
with open('tfko/data/common_tfs.txt', 'r') as file:
    tf_names = file.read().splitlines()
tf_names.append('Control (No Knockout)')

# Set up streamlit
timepoint = 'TDR125'
timepoints = {'10 hours post fertilization': 'TDR126',
                '12 hours post fertilization': 'TDR127',
                '14 hours post fertilization': 'TDR128',
                '16 hours post fertilization': 'TDR118',
                '19 hours post fertilization': 'TDR125',
                '24 hours post fertilization': 'TDR124'}
timepoint = st_setup(list(timepoints.keys()))
st.markdown(f'### {timepoint}')

# Generate plot
selected_knockouts = []
default = 'Control (No Knockout)'
for key in st.session_state.selectboxes:
    st.sidebar.selectbox(
        'Select a knockout to plot',
        tf_names,
        index=tf_names.index(default),
        key=key
    )
    selected_knockouts.append(st.session_state[key])
fig = plot_cells(selected_knockouts, timepoints[timepoint])
st.plotly_chart(fig, use_container_width=True)
