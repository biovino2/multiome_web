"""Plot GRN heat maps for each time point given a cell type.

Ben Iovino  08/22/24    CZ-Biohub
"""

import os
import streamlit as st
from plot_grn import plot_grn, save_config
from plotly.subplots import make_subplots
from util import get_datasets


def st_setup(lineages: list[str] = None) -> str:
    """Initializes streamlit session.

    Args:
        lineages: The list of lineage options.

    Returns:
        lineage (str): The selected lineage.
    """

    st.set_page_config(layout="wide")
    st.title('Time-Resolved Gene Regulatory Networks (GRNs)')
    st.write('For each time point, we present the Gene Regulatory Networks (GRNs), computed by CellOracle (Kamimoto et al., 2023), \
            for every cell-type. These GRNs are visualized as heatmaps, with time points on the x-axis and transcription factor-gene \
            pairs on the y-axis. The color intensity in each element represents the strength of the interaction between a transcription \
            factor and a gene, with red indicating activation and blue indicating repression.')
             
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

    # Get lineage
    # lineage = st.sidebar.selectbox(
    #    'Select lineage',
    #    lineages
    # )

    # return lineage


def make_figure(path: str, timepoints: dict[str:str]):
    """Returns a plotly figure containing a heatmap representing the GRN for each time point.

    Args:
        path (str): The path to the data.
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
        horizontal_spacing=0.075,
        cols=len(selected_timepoints),
        subplot_titles=[f"{tp}" for tp in selected_timepoints]
        )

    # Generate heatmap for each subplot
    for i, timepoint in enumerate(selected_timepoints):
        plot, fams = plot_grn(path, celltype=None, timepoint=timepoints[timepoint])
        fig.add_trace(plot, row=1, col=i+1)

        # Update hover template and axis labels
        fig.data[-1].update(
            hovertemplate='<b>Cell type:</b> %{x}<br>'
            '<b>TF-Gene:</b> %{y}<br>'
            '<b>Strength:</b> %{z}<extra></extra>',
        hoverlabel=dict(namelength=0))
        fig.update_xaxes(tickfont=dict(size=12), row=1, col=i+1)
        fig.update_yaxes(tickfont=dict(size=12),
                        tickvals=list(fams.values()),
                        ticktext=list(fams.keys()),
                        row=1, col=i+1)

    # Figure layout
    fig.update_layout(height=500, width=2000, showlegend=False, margin=dict(l=0, r=0, t=50, b=0))
    fig.update_layout(coloraxis=dict(colorscale='RdBu_r'))

    return fig, selected_timepoints


### Main

path = os.path.dirname(os.path.abspath(__file__))+'/data'

# Plot figure
st_setup()
datasets = get_datasets()
fig, selected_timepoints = make_figure(path, datasets)
st.plotly_chart(fig, config=save_config())

# Download button
#with open(path + '/data.zip', 'rb') as file:
#    data = file.read()
#st.sidebar.download_button('Download Data', data=data, file_name=path + '/data.zip')
