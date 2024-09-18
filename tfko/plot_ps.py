"""Plot perturbation score scatter plots for a given time point.

Ben Iovino  09/16/24    CZ-Biohub
"""

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st


def st_setup():
    """Initializes streamlit session.
    """

    st.set_page_config(layout="wide")
    st.title('In-Silico Transcription Factor Knockout')
    st.write('For all time points, we plot the perturbation scores of each transcription factor knockout for mesodermal and neuro-ectodermal cells. ' \
            'The perturbation score is a measure of how much the knockout influences the direction of the cell type transition.')
    

def make_figure(timepoints: dict[str:str]) -> go.Figure:
    """Returns a plotly figure containing scatter plots of perturbation scores for each time point.

    Args:
        timepoints (dict[str:str]): The dictionary of time points.

    Returns:
        fig (plotly.graph_objs.Figure): The plotly figure.
    """

    fig = make_subplots(rows=2, cols=3, subplot_titles=list(timepoints.keys()))

    # Plot scatter plots for each timepoint
    for i, tp in enumerate(timepoints.values()):

        cos_sim = pd.read_csv(f'tfko/data/{tp}_sim.csv', index_col=0)

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
        fig.add_trace(scatter, row=(i//3)+1, col=(i%3)+1)

    # Update layout
    fig.update_layout(
        margin=dict(l=10, r=10, t=70, b=0),
        showlegend=False
    )

    # Update axes labels and titles
    for i in range (0, len(timepoints)):
        fig.update_xaxes(title='Mesoderm', row=i//3+1, col=i%3+1, range=[0.4, 0.80])
        fig.update_yaxes(title='Neuro-ectoderm', row=i//3+1, col=i%3+1, range=[0.4, 0.65])
    fig.update_layout(width=1200, height=700)

    return fig


# Main
timepoints = {'10 hours post fertilization': 'TDR126',
                '12 hours post fertilization': 'TDR127',
                '14 hours post fertilization': 'TDR128',
                '16 hours post fertilization': 'TDR118',
                '19 hours post fertilization': 'TDR125',
                '24 hours post fertilization': 'TDR124'}
st_setup()
fig = make_figure(timepoints)
st.plotly_chart(fig)
