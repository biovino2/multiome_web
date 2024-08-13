"""Plot GRN cluster maps.

Ben Iovino  08/09/24    CZ-Biohub
"""

import numpy as np
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import seaborn as sns
import streamlit as st


def load_data(celltype: str, timepoint: str) -> 'tuple[pd.DataFrame, np.array, np.array]':
    """Returns counts for celltype at a timepoint, as well as row and column linkages.

    Args:
        celltype: The celltype to load.
        timepoint: The timepoint to load.

    Returns:
        df_counts (pd.DataFrame): The counts for the celltype at the timepoint.
        row_linkage (np.array): The row dendrogram linkage.
        col_linkage (np.array): The column dendrogram linkage.
    """

    path = 'grn/data'
    df_counts = pd.read_csv(f"{path}/{celltype}/{celltype}_{timepoint}.csv", index_col=0)
    row_linkage = np.load(f"{path}/{celltype}/{celltype}_row_linkage.npz")['arr_0']
    col_linkage = np.load(f"{path}/{celltype}/{celltype}_col_linkage.npz")['arr_0']

    return df_counts, row_linkage, col_linkage


def st_setup(celltypes: list[str]):
    """Initializes streamlit session.

    Args:
        celltypes: The list of cell types.
    """

    st.set_page_config(layout="wide")
    st.title('Time-Resolved Gene Regulatory Networks (GRNs)')
    st.write('For each timepoint, we plot the CellOracle GRNs for each cell type.')
    st.sidebar.markdown('# Settings')

    # Initialize drop down boxes
    if "selectboxes" not in st.session_state:
        st.session_state.selectboxes = [0]

    # Set up buttons for adding/removing timepoints
    with st.sidebar:
        add, reset = st.columns([1, 1])
        with add:
            if st.button('Add Timepoint'):
                st.session_state.selectboxes.append(len(st.session_state.selectboxes))
        with reset:
            if st.button('Reset'):
                st.session_state.selectboxes = [0]

    # Get celltype
    celltype = st.sidebar.selectbox(
        'Select a cell type to plot',
        celltypes
    )

    return celltype


def custom_scale() -> 'list[list[float, str]]':
    """Returns custom color scale for plotly.

    Returns:
        colorscale (list): The custom color scale.
    """

    # Blue, grey, to red
    colorscale = [
        [0, 'rgb(0,0,255)'],
        [0.1, 'rgb(188,188,188)'],
        [0.2, 'rgb(255,0,0)']
        ]
    
    return colorscale


def plot_grn(celltype: str, timepoint: str) -> 'px.imshow':
    """Returns a plotly figure of the GRN clustermap.

    Args:
        celltype (str): The celltype to plot.
        timepoint (str): The timepoint to plot.

    Returns:
        fig (plotly.express.imshow): The plotly figure.
    """

    # Seaborn clustermap
    vmax, vmin = 0.1, -0.1
    df_counts, row_linkage, col_linkage = load_data(celltype, timepoint)

    # Create clustermap
    g = sns.clustermap(df_counts, method='ward', metric='euclidean', 
                       row_cluster=True, col_cluster=True, 
                       xticklabels=df_counts.columns.tolist(),
                       yticklabels=df_counts.index.tolist(), 
                       vmax=vmax, vmin=vmin, 
                       row_linkage=row_linkage, col_linkage=col_linkage)

    # Convert seaborn plot to plotly and plot on page
    fig = px.imshow(g.data2d,
                    labels=dict(x='Transcription Factor', y='Gene'))

    return fig


def main():
    """
    """
    
    timepoints = {'0 hours post fertilization': 'TDR126',
              '5 hours post fertilization': 'TDR127',
              '10 hours post fertilization': 'TDR128',
              '15 hours post fertilization': 'TDR118',
              '20 hours post fertilization': 'TDR125',
              '30 hours post fertilization': 'TDR124'}
    celltypes = ['fast_muscle', 'neural_posterior', 'NMPs', 'PSM', 'somites', 'spinal_cord', 'tail_bud']
    celltype = st_setup(celltypes)

    # Graph selected timepoints
    selected_timepoints = []
    for i, key in enumerate(st.session_state.selectboxes):
        st.sidebar.selectbox(
            'Select a timepoint to plot',
            list(timepoints.keys()),
            key=key
        )
        selected_timepoints.append(st.session_state[key])
        fig = make_subplots(rows=1,
                            cols=len(selected_timepoints),
                            horizontal_spacing=.1,
                            subplot_titles=[f"{tp}" for tp in selected_timepoints])

    # Create a single plot containing all GRNs
    for i, timepoint in enumerate(selected_timepoints):
        plot = plot_grn(celltype, timepoints[timepoint])
        for trace in plot['data']:
            fig.add_trace(trace, row=1, col=i+1)
        fig.update_xaxes(tickfont=dict(size=8), row=1, col=i+1)
        fig.update_yaxes(tickfont=dict(size=8), row=1, col=i+1)

    # Display complete figure
    fig.update_layout(height=800, width=3000, showlegend=False)
    fig.update_layout(coloraxis=dict(colorscale=custom_scale()))
    st.plotly_chart(fig)


if __name__ == '__main__':
    main()
