"""Master script for plotting GRNs.

Ben Iovino  08/09/24    CZ-Biohub
"""

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st
from util import get_timepoints_abbr


def load_data(path: str, celltype: str, timepoint: str, control: str) -> 'pd.DataFrame':
    """Returns counts for celltype at a timepoint.

    Args:
        path (str): The path to the data.
        celltype (str): The celltype to load.
        timepoint (str): The timepoint to load.
        control (str): The control variable (time point or cell type)

    Returns:
        df_counts (pd.DataFrame): The counts for the celltype at the timepoint.
    """

    abbr = get_timepoints_abbr()

    # Plotting multiple time points for each cell type
    if control == 'timepoint':
        path += '/ct'  # this is very confusing
        df_counts = pd.read_csv(f"{path}/{celltype}/{celltype}_{abbr[timepoint]}.csv", index_col=0)

    # Plotting multiple cell types for each time point
    if control == 'celltype':
        path += '/tp'  # still confusing
        df_counts = pd.read_csv(f"{path}/{abbr[timepoint]}/{abbr[timepoint]}_{celltype}.csv", index_col=0)
    df_counts = df_counts.transpose()  # genes on y-axis, TFs on x-axis

    return df_counts


def plot_grn(path: str, celltype: str, timepoint: str, control: str) -> 'go.Heatmap':
    """Returns a plotly Heatmap of the GRN.

    Args:
        path (str): The path to the data.
        celltype (str): The celltype to plot.
        timepoint (str): The timepoint to plot.
        control (str): The control variable (time point or cell type)

    Returns:
        fig (go.Heatmap): The plotly figure.
    """

    # Load data
    df_counts = load_data(path, celltype, timepoint, control)

    # Create heatmap figure
    fig = go.Heatmap(
        z=df_counts.values,
        x=df_counts.columns,
        y=df_counts.index,
        colorscale='balance',
        zmin=-0.1, zmax=0.1
    )

    return fig


def plot_scores(path: str, celltypes: 'list[str]', timepoints: 'list[str]') -> 'go.Heatmap':
    """Returns a scatter plot of the network scores.

    Args:
        path (str): The path to the data.
        celltypes (list[str]): The list of cell types.
        timepoint (list[str]): The list of time points.
    """

    # Create figure depending on size of lists
    fig = make_subplots(
        rows = 1,
        cols = max(len(celltypes), len(timepoints)),
    )

    # Change timepoints to abbreviations
    timepoints = [f'{tp.split()[0]}hpf' for tp in timepoints]

    for i, ct in enumerate(celltypes):
        for j, tp in enumerate(timepoints):
            df_scores = pd.read_csv(f"{path}/scores/{tp}/{ct}.csv", index_col=0)

            # Create scatter plot where x is score, y is gene name
            fig.add_trace(
                go.Scatter(
                    x=df_scores['degree_centrality_all'],
                    y=df_scores.index,
                    mode='markers',
                    marker=dict(
                        size=10,
                        color=df_scores['degree_centrality_all'],
                    ),
                    name=f"{ct} {tp}"
                ),
                row=1, col=max(i, j)+1
            )

            # Reverse y-axis
            fig.update_yaxes(autorange='reversed', row=1, col=max(i, j)+1)
            fig.update_xaxes(title='Degree Centrality', row=1, col=max(i, j)+1)
            fig.update_layout(height=500, width=2000,
                                showlegend=False,
                                margin=dict(l=0, r=0, t=40, b=0),
                                yaxis=dict(tickmode='linear'))

            # Update hoverinfo
            fig.data[-1].update(
                hovertemplate='TF: %{y}<br>Degree Centrality: %{x}<extra></extra>',
                hoverlabel=dict(namelength=0))

    return fig


def save_config() -> dict:
    """Returns a config to save plotly figure as SVG.

    Returns:
        config (dict): The save configuration.
    """

    config = {
        'toImageButtonOptions': {
        'format': 'svg', # one of png, svg, jpeg, webp
        'filename': 'grn_heatmap',
        'height': None,
        'width': None,
        'scale': 1 # Multiply title/legend/axis/canvas sizes by this factor
    },
    'displayModeBar': True
    }

    return config


def main():
    """Runs selected page.
    """
    
    pg = st.navigation([
        st.Page('plot_tp.py', title='Time Points'),
        st.Page('plot_ct.py', title='Cell Types'),
        ])
    pg.run()


if __name__ == '__main__':
    main()
