"""Master script for plotting GRNs.

Ben Iovino  08/09/24    CZ-Biohub
"""

import pandas as pd
import plotly.graph_objects as go
import streamlit as st


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

    # Plotting multiple time points for each cell type
    if control == 'timepoint':
        path += '/ct'  # this is very confusing
        df_counts = pd.read_csv(f"{path}/{celltype}/{celltype}_{timepoint}.csv", index_col=0)

    # Plotting multiple cell types for each time point
    if control == 'celltype':
        path += '/tp'  # still confusing
        df_counts = pd.read_csv(f"{path}/{timepoint}/{timepoint}_{celltype}.csv", index_col=0)
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

    # Create heatmap directly using graph_objects
    df_counts = load_data(path, celltype, timepoint, control)
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
