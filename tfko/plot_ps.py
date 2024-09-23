"""Plot perturbation score scatter plots for a given time point.

Ben Iovino  09/16/24    CZ-Biohub
"""

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st
import sys


def st_setup():
    """Initializes streamlit session.
    """

    st.set_page_config(layout="wide")
    st.title('In-Silico Transcription Factor Knockout')
    st.write('For all time points, we plot the perturbation scores of each transcription factor knockout for mesodermal and neuro-ectodermal cells. ' \
            'The perturbation score is a measure of how much the knockout influences the direction of the cell type transition.')
    
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
    

def make_figure(path: str, timepoints: dict[str:str], inv_tp: dict[str:str]) -> go.Figure:
    """Returns a plotly figure containing scatter plots of perturbation scores for each time point.

    Args:
        path (str): The path to the data.
        timepoints (dict[str:str]): The dictionary of time points.
        inv_tp (dict[str:str]): The inverted dictionary of time points.
    Returns:
        fig (plotly.graph_objs.Figure): The plotly figure.
    """

    fig = make_subplots(
        rows=2,
        cols=3,
        subplot_titles=list(timepoints.keys()))

    # Keep track of PS scores for every time point for each gene
    meso_scores: 'dict[dict[str, float]]' = {}  # ex. {meox1: {TDR118: 0.42, ...}, ...}
    ne_scores: 'dict[dict[str, float]]' = {}

    # Plot scatter plots for each timepoint
    for i, tp in enumerate(timepoints.values()):
        cos_sim = pd.read_csv(f'{path}/{tp}_sim.csv', index_col=0)

        # Subset mesoderm cells and compute perturbation score
        df_meso = cos_sim[cos_sim.celltype.isin(["PSM","fast_muscle","somites"])].drop(columns="celltype")
        df_meso_avg = df_meso.median(axis=0)
        df_meso_avg = 1 - df_meso_avg

        # Subset neuro-ectoderm celltypes
        df_ne = cos_sim[cos_sim.celltype.isin(["neural_posterior","spinal_cord"])].drop(columns="celltype")
        df_ne_avg = df_ne.median(axis=0)
        df_ne_avg = 1 - df_ne_avg

        # Merge dataframes
        df_merged = df_meso_avg.to_frame(name="meso").join(df_ne_avg.to_frame(name="ne"))

        # Add genes to dictionary
        for gene in df_merged.index:
            meso_scores[inv_tp[tp]] = meso_scores.get(inv_tp[tp], {})
            ne_scores[inv_tp[tp]] = ne_scores.get(inv_tp[tp], {})
            meso_scores[inv_tp[tp]][gene] = df_merged.loc[gene, "meso"]
            ne_scores[inv_tp[tp]][gene] = df_merged.loc[gene, "ne"]

        # Scatter plot
        scatter = go.Scatter(
            x=df_merged['meso'],
            y=df_merged['ne'],
            mode='markers',
            marker=dict(size=5, color='dodgerblue', opacity=0.6),
            hoverinfo='text',
            text=df_merged.index,
            showlegend=False
        )

        # Add to figure, update hover info
        fig.add_trace(scatter, row=(i//3)+1, col=(i%3)+1)
        fig.data[-1].update(
            hovertemplate='<b>TF</b>: %{text}<br><b>Mesoderm</b>: %{x}<br><b>Neuro-ectoderm</b>: %{y}',
        hoverlabel=dict(namelength=0))

    # Update axes labels and titles
    for i in range (0, len(timepoints)):
        fig.update_xaxes(title='Mesoderm', row=i//3+1, col=i%3+1, range=[0.4, 0.80])
        fig.update_yaxes(title='Neuro-ectoderm', row=i//3+1, col=i%3+1, range=[0.4, 0.65])

    # Update figure size
    fig.update_layout(width=1200, height=550,
        margin=dict(l=10, r=10, t=70, b=0),
        showlegend=False )
    
    # Convert dictionary to dataframe
    meso_df = pd.DataFrame.from_dict(meso_scores)
    ne_df = pd.DataFrame.from_dict(ne_scores)

    return fig, meso_df, ne_df


# Main
arg: 'list[str]' = sys.argv[0].split('/')[:-1]
path = '/'.join(arg) + '/data'

timepoints = {'10 hours post fertilization': 'TDR126',
                '12 hours post fertilization': 'TDR127',
                '14 hours post fertilization': 'TDR128',
                '16 hours post fertilization': 'TDR118',
                '19 hours post fertilization': 'TDR125',
                '24 hours post fertilization': 'TDR124'}
inv_tp = {'TDR126': '10 hpf',
            'TDR127': '12 hpf',
            'TDR128': '14 hpf',
            'TDR118': '16 hpf',
            'TDR125': '19 hpf',
            'TDR124': '24 hpf'}
st_setup()
fig, meso_df, ne_df = make_figure(path, timepoints, inv_tp)
st.plotly_chart(fig)

# Add two checkboxes to display dataframes, one in each column
meso, ne = st.columns([1, 1])
with meso:
    if st.checkbox('Show mesoderm perturbation scores'):
        st.dataframe(meso_df)
with ne:
    if st.checkbox('Show neuro-ectoderm perturbation scores'):
        st.dataframe(ne_df)
