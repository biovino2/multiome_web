"""Plot GRN heat maps for each time point given a cell type.

Ben Iovino  08/22/24    CZ-Biohub
"""

import streamlit as st
from plot_grn import plot_grn, save_config
from plotly.subplots import make_subplots


def st_setup(celltypes: list[str]) -> str:
    """Initializes streamlit session.

    Args:
        celltypes: The list of cell types options.

    Returns:
        celltype (str): The selected cell type.
    """

    st.set_page_config(layout="wide")
    st.title('Time-Resolved Gene Regulatory Networks (GRNs)')
    st.write('For each time point (hours post fertilization), we plot the CellOracle predicted GRNs for each cell type. \
             The GRNs are represented as a heatmap, with the transcription factors on the x-axis and the genes on the y-axis. \
             The color of each cell represents the strength of the interaction between the transcription factor and the gene, \
             with red indicating activation and blue indicating repression, as predicted by CellOracle.')
    st.write('You can select the cell type and time points to display using the sidebar on the left.')
    st.sidebar.markdown('# Settings')

    # Initialize drop down boxes
    if "selectboxes" not in st.session_state:
        st.session_state.selectboxes = [0]

    # Set up buttons for adding/removing timepoints
    with st.sidebar:
        add, reset = st.columns([1, 1])
        with add:
            if st.button('Add Time Point'):
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


def make_figure(celltype: str, timepoints: dict[str:str]):
    """Returns a plotly figure containing a heatmap representing the GRN for each time point.

    Args:
        celltype (str): The celltype to plot.
        timepoints (dict[str:str]): The dictionary of time points.

    Returns:
        fig (plotly.graph_objs.Figure): The plotly figure.
    """

    # Obtain selected timepoints
    selected_timepoints = []
    for i, key in enumerate(st.session_state.selectboxes):
        st.sidebar.selectbox(
            'Select a time point to plot',
            list(timepoints.keys()),
            key=key
        )
        selected_timepoints.append(st.session_state[key])

    # Create figure (subplots)
    fig = make_subplots(
        rows=1,
        cols=len(selected_timepoints),
        subplot_titles=[f"{tp}" for tp in selected_timepoints],
        shared_yaxes=True,
        shared_xaxes=True
        )

    # Generate heatmap for each subplot
    for i, timepoint in enumerate(selected_timepoints):
        plot = plot_grn(celltype, timepoints[timepoint])
        fig.add_trace(plot, row=1, col=i+1)
        fig.update_xaxes(tickfont=dict(size=12), row=1, col=i+1, matches='x')
        fig.update_yaxes(tickfont=dict(size=12), row=1, col=i+1, matches='y')

    # Figure layout
    fig.update_layout(height=700, width=2000, showlegend=False)
    fig.update_layout(coloraxis=dict(colorscale='RdBu_r'))

    return fig

# Main
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
