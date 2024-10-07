"""Plot GRN heat maps for each cell type given a time point.

Ben Iovino  08/22/24    CZ-Biohub
"""

import os
import streamlit as st
from plot_grn import plot_grn, save_config, plot_scores
from plotly.subplots import make_subplots
from util import get_datasets


def st_setup(timepoints: list[str]) -> str:
    """Initializes streamlit session.

    Args:
        timepoints: The list of time point options.

    Returns:
        timepoint (str): The selected time point.
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

    # Set up buttons for adding/removing time points
    with st.sidebar:
        add, reset = st.columns([1, 1])
        with add:
            if st.button('Add Cell Type'):
                st.session_state.selectboxes.append(len(st.session_state.selectboxes))
        with reset:
            if st.button('Reset'):
                st.session_state.selectboxes = [0]

    # Get time point
    timepoint = st.sidebar.selectbox(
        'Select time point',
        timepoints
    )

    return timepoint


def make_figure(path: str, timepoint: str, celltypes: 'list[str]'):
    """Returns a plotly figure containing a heatmap representing the GRN for each time point.

    Args:
        path (str): The path to the data.
        timepoint (str): The time point to plot.
        celltypes (list[str]): The list of cell types.

    Returns:
        fig (plotly.graph_objs.Figure): The plotly figure.
    """

    # Obtain selected time points
    selected_celltypes = []
    for i, key in enumerate(st.session_state.selectboxes):
        st.sidebar.selectbox(
            f'Select cell type {i+1}',
            celltypes,
            key=key
        )
        selected_celltypes.append(st.session_state[key])

    # Create figure (subplots)
    fig = make_subplots(
        rows=1,
        horizontal_spacing=0.075,
        cols=len(selected_celltypes),
        subplot_titles=[f"{(' '.join(ct.split('_')))}" for ct in selected_celltypes],
        shared_xaxes=True,
        shared_yaxes=True
        )

    # Generate heatmap for each subplot
    for i, celltype in enumerate(selected_celltypes):
        plot = plot_grn(path, celltype, timepoint, 'celltype')
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

    return fig, selected_celltypes

### Main

path = os.path.dirname(os.path.abspath(__file__))+'/data'

# Plot figure
datasets = get_datasets()
celltypes = ['neural_posterior', 'NMPs', 'PSM', 'somites', 'spinal_cord', 'tail_bud']
timepoint = st_setup(list(datasets.keys()))
st.markdown(f'### {timepoint}')
fig, selected_celltypes = make_figure(path, datasets[timepoint], celltypes)
scores = plot_scores(path, selected_celltypes, [timepoint])
st.plotly_chart(fig, config=save_config())

# Download button
with open(path + '/data.zip', 'rb') as file:
    data = file.read()
st.sidebar.download_button('Download Data', data=data, file_name=path + '/data.zip')

# Button to display scatter plots
if st.checkbox('Show degree centrality scores'):
    st.plotly_chart(scores, config=save_config())
