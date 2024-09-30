"""Plot GRN heat maps for each time point given a cell type.

Ben Iovino  08/22/24    CZ-Biohub
"""

import os
import streamlit as st
from plot_grn import plot_grn, save_config
from plotly.subplots import make_subplots
from util import get_datasets


def st_setup(celltypes: list[str]) -> str:
    """Initializes streamlit session.

    Args:
        celltypes: The list of cell types options.

    Returns:
        celltype (str): The selected cell type.
    """

    st.set_page_config(layout="wide")
    st.title('Time-Resolved Gene Regulatory Networks (GRNs)')
    st.write('For each time point, we plot the CellOracle predicted GRNs for each cell type. The \
             GRNs are represented as a heatmap, with the transcription factors on the x-axis and \
             the genes on the y-axis. The color of each cell represents the strength of the \
             interaction between the transcription factor and the gene, with red indicating \
             activation and blue indicating repression.')
             
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
        'Select cell type',
        celltypes
    )

    return celltype


def make_figure(path: str, celltype: str, timepoints: dict[str:str]):
    """Returns a plotly figure containing a heatmap representing the GRN for each time point.

    Args:
        path (str): The path to the data.
        celltype (str): The celltype to plot.
        timepoints (dict[str:str]): The dictionary of time points.

    Returns:
        fig (plotly.graph_objs.Figure): The plotly figure.
    """

    # Obtain selected timepoints
    selected_timepoints = []
    for i, key in enumerate(st.session_state.selectboxes):
        st.sidebar.selectbox(
            f'Select time point {i+1}',
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
        shared_xaxes=True,
        )

    # Generate heatmap for each subplot
    for i, timepoint in enumerate(selected_timepoints):
        plot = plot_grn(path, celltype, timepoints[timepoint], 'timepoint')
        fig.add_trace(plot, row=1, col=i+1)

        # Update hover template and axis labels
        fig.data[-1].update(
            hovertemplate='TF: %{x}<br>Gene: %{y}<br>Strength: %{z}<extra></extra>',
        hoverlabel=dict(namelength=0))
        fig.update_xaxes(tickfont=dict(size=12), row=1, col=i+1, matches='x')
        fig.update_yaxes(tickfont=dict(size=12), row=1, col=i+1, matches='y')

    # Figure layout
    fig.update_layout(height=500, width=2000, showlegend=False, margin=dict(l=0, r=0, t=50, b=0))
    fig.update_layout(coloraxis=dict(colorscale='RdBu_r'))

    return fig

### Main

path = os.path.dirname(os.path.abspath(__file__))+'/data'

# Plot figure
datasets = get_datasets()
celltypes = ['neural_posterior', 'NMPs', 'PSM', 'somites', 'spinal_cord', 'tail_bud']
celltype = st_setup(celltypes)
st.markdown(f'### {celltype}')
fig = make_figure(path, celltype, datasets)
st.plotly_chart(fig, config=save_config())

# Download button
with open(path + '/data.zip', 'rb') as file:
    data = file.read()
st.sidebar.download_button('Download Data', data=data, file_name=path + '/data.zip')
