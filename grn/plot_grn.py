"""Plot GRN cluster maps.

Ben Iovino  08/09/24    CZ-Biohub
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st


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


def load_data(celltype: str, timepoint: str) -> 'tuple[pd.DataFrame, np.array, np.array]':
    """Returns counts for celltype at a timepoint.

    Args:
        celltype: The celltype to load.
        timepoint: The timepoint to load.

    Returns:
        df_counts (pd.DataFrame): The counts for the celltype at the timepoint.
    """

    path = 'grn/data'
    df_counts = pd.read_csv(f"{path}/{celltype}/{celltype}_{timepoint}.csv", index_col=0)

    return df_counts


def plot_grn(celltype: str, timepoint: str) -> 'px.imshow':
    """Returns a plotly figure of the GRN clustermap.

    Args:
        celltype (str): The celltype to plot.
        timepoint (str): The timepoint to plot.

    Returns:
        fig (plotly.express.imshow): The plotly figure.
    """

    # Create heatmap directly using graph_objects
    df_counts = load_data(celltype, timepoint)
    fig = go.Figure(data=go.Heatmap(
        z=df_counts.values,
        x=df_counts.columns,
        y=df_counts.index,
        colorscale='balance',
        zmin=-0.1, zmax=0.1
    ))

    # Customize layout
    fig.update_layout(
        xaxis_title="Transcription Factor",
        yaxis_title="Gene",
    )

    return fig


def make_figure(celltype: str, timepoints: dict[str:str]):
    """Returns a plotly figure containing a clustermap representing the GRN for each timepoint.

    Args:
        celltype (str): The celltype to plot.
        timepoints (dict[str:str]): The dictionary of timepoints.

    Returns:
        fig (plotly.graph_objs.Figure): The plotly figure.
    """

    # Obtain selected timepoints
    selected_timepoints = []
    for i, key in enumerate(st.session_state.selectboxes):
        st.sidebar.selectbox(
            'Select a timepoint to plot',
            list(timepoints.keys()),
            key=key
        )
        selected_timepoints.append(st.session_state[key])

    # Create figure (subplots)
    fig = make_subplots(
        rows=1,
        cols=len(selected_timepoints),
        horizontal_spacing=.1,
        subplot_titles=[f"{tp}" for tp in selected_timepoints]
        )

    # Generate clustermap for each subplot
    for i, timepoint in enumerate(selected_timepoints):
        plot = plot_grn(celltype, timepoints[timepoint])
        for trace in plot['data']:
            fig.add_trace(trace, row=1, col=i+1)
        fig.update_xaxes(tickfont=dict(size=12), row=1, col=i+1)
        fig.update_yaxes(tickfont=dict(size=12), row=1, col=i+1)

    # Figure layout
    fig.update_layout(height=800, width=3000, showlegend=False)
    fig.update_layout(coloraxis=dict(colorscale='RdBu_r'))

    return fig


def save_config() -> dict:
    """Returns a config to save plotly figure as SVG.

    Returns:
        config (dict): The configuration.
    """

    config = {
        'toImageButtonOptions': {
        'format': 'svg', # one of png, svg, jpeg, webp
        'filename': 'grn_clustermap',
        'height': None,
        'width': None,
        'scale': 1 # Multiply title/legend/axis/canvas sizes by this factor
    },
    'displayModeBar': True
    }

    return config


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
    fig = make_figure(celltype, timepoints)
    st.plotly_chart(fig, config=save_config())


if __name__ == '__main__':
    main()
