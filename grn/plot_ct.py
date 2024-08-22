"""Plot GRN heat maps for each cell type given a time point.

Ben Iovino  08/22/24    CZ-Biohub
"""

import streamlit as st
from plot_grn import plot_grn, save_config
from plotly.subplots import make_subplots


def st_setup(timepoints: list[str]) -> str:
    """Initializes streamlit session.

    Args:
        timepoints: The list of time point options.

    Returns:
        timepoint (str): The selected time point.
    """

    st.set_page_config(layout="wide")
    st.title('Time-Resolved Gene Regulatory Networks (GRNs)')
    st.write('For each cell type, we plot the CellOracle predicted GRNs for each time point (hours post fertilization). \
             The GRNs are represented as a heatmap, with the transcription factors on the x-axis and the genes on the y-axis. \
             The color of each cell represents the strength of the interaction between the transcription factor and the gene, \
             with red indicating activation and blue indicating repression, as predicted by CellOracle.')
    st.write('You can select the time point and cell types to display using the sidebar on the left.')
    st.sidebar.markdown('# Settings')

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
        'Select a time point to plot',
        timepoints
    )

    return timepoint


def make_figure(timepoint: str, celltypes: 'list[str]'):
    """Returns a plotly figure containing a heatmap representing the GRN for each time point.

    Args:
        timepoint (str): The time point to plot.
        celltypes (list[str]): The list of cell types.

    Returns:
        fig (plotly.graph_objs.Figure): The plotly figure.
    """

    # Obtain selected time points
    selected_celltypes = []
    for i, key in enumerate(st.session_state.selectboxes):
        st.sidebar.selectbox(
            'Select a cell type to plot',
            celltypes,
            key=key
        )
        selected_celltypes.append(st.session_state[key])

    # Create figure (subplots)
    fig = make_subplots(
        rows=1,
        horizontal_spacing=0.1,
        cols=len(selected_celltypes),
        subplot_titles=[f"{(' '.join(ct.split('_')))}" for ct in selected_celltypes]
        )

    # Generate heatmap for each subplot
    for i, celltype in enumerate(selected_celltypes):
        plot = plot_grn(celltype, timepoint)
        fig.add_trace(plot, row=1, col=i+1)
        fig.update_xaxes(tickfont=dict(size=12), row=1, col=i+1)
        fig.update_yaxes(tickfont=dict(size=12), row=1, col=i+1)

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
timepoint = st_setup(list(timepoints.keys()))
fig = make_figure(timepoints[timepoint], celltypes)
st.plotly_chart(fig, config=save_config())