"""Master script for plotting GRNs.

Ben Iovino  08/09/24    CZ-Biohub
"""

import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from util import get_timepoints_abbr


def load_data(path: str, celltype: str, timepoint: str) -> 'pd.DataFrame':
    """Returns counts for celltype at a timepoint.

    Args:
        path (str): The path to the data.
        celltype (str): The celltype to load.
        timepoint (str): The timepoint to load.

    Returns:
        df_grn (pd.DataFrame): The GRN for the given celltype and timepoint.
    """

    abbr = get_timepoints_abbr()

    # Plotting multiple time points for each cell type
    if timepoint:
        df_grn = pd.read_csv(f"{path}/tp/{abbr[timepoint]}.csv", index_col=0)

    # Plotting multiple cell types for each time point
    if celltype:
        df_grn = pd.read_csv(f"{path}/ct/{celltype}.csv", index_col=0)

    return df_grn


def plot_grn(path: str, celltype: str, timepoint: str) -> 'go.Heatmap':
    """Returns a plotly Heatmap of the GRN.

    Args:
        path (str): The path to the data.
        celltype (str): The celltype to plot.
        timepoint (str): The timepoint to plot.

    Returns:
        fig (go.Heatmap): The plotly figure.
    """

    # Load data
    df_grn = load_data(path, celltype, timepoint)

    # Create heatmap figure
    fig = go.Heatmap(
        z=df_grn.values,
        x=df_grn.columns[:-3],  # last two columns are TF names
        y=df_grn.index,
        colorscale='balance',
        zmin=-0.1, zmax=0.1
    )

    # Change y-tick labels so that they occur when a new TF family starts
    fams = {}
    for i, index in enumerate(df_grn.index):
        fams[df_grn.loc[index, 'TF_family']] = i

    return fig, fams


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
