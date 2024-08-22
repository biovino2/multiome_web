"""Master script for plotting GRNs.

Ben Iovino  08/09/24    CZ-Biohub
"""

import pandas as pd
import plotly.graph_objects as go
import streamlit as st


def load_data(celltype: str, timepoint: str) -> 'pd.DataFrame':
    """Returns counts for celltype at a timepoint.

    Args:
        celltype (str): The celltype to load.
        timepoint (str): The timepoint to load.

    Returns:
        df_counts (pd.DataFrame): The counts for the celltype at the timepoint.
    """

    path = 'grn/data'
    df_counts = pd.read_csv(f"{path}/{celltype}/{celltype}_{timepoint}.csv", index_col=0)
    df_counts = df_counts.transpose()  # genes on y-axis, TFs on x-axis

    return df_counts


def plot_grn(celltype: str, timepoint: str) -> 'go.Heatmap':
    """Returns a plotly Heatmap of the GRN.

    Args:
        celltype (str): The celltype to plot.
        timepoint (str): The timepoint to plot.

    Returns:
        fig (go.Heatmap): The plotly figure.
    """

    # Create heatmap directly using graph_objects
    df_counts = load_data(celltype, timepoint)
    fig = go.Heatmap(
        z=df_counts.values,
        x=df_counts.columns,
        y=df_counts.index,
        colorscale='balance',
        zmin=-0.1, zmax=0.1
    )

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
