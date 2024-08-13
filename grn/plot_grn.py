"""Plot GRN cluster maps.

Ben Iovino  08/09/24    CZ-Biohub
"""

import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
import streamlit as st


def load_data(celltype: str, timepoint: str) -> 'tuple[pd.DataFrame, np.array, np.array]':
    """Returns counts for celltype at a timepoint, as well as row and column linkages.

    Args:
        celltype: The celltype to load.
        timepoint: The timepoint to load.

    Returns:
        df_counts (pd.DataFrame): The counts for the celltype at the timepoint.
        row_linkage (np.array): The row dendrogram linkage.
        col_linkage (np.array): The column dendrogram linkage.
    """

    path = 'grn/data'
    df_counts = pd.read_csv(f"{path}/{celltype}/{celltype}_{timepoint}.csv", index_col=0)
    row_linkage = np.load(f"{path}/{celltype}/{celltype}_row_linkage.npz")['arr_0']
    col_linkage = np.load(f"{path}/{celltype}/{celltype}_col_linkage.npz")['arr_0']

    return df_counts, row_linkage, col_linkage


def st_setup():
    """Initializes streamlit session.

    Args:
        celltypes: The list of cell types.
        timepoints: The list of timepoints.
    """

    st.set_page_config(layout="wide")
    st.title('Time-Resolved Gene Regulatory Networks (GRNs)')
    st.write('For each timepoint, we plot the CellOracle GRNs for each cell type.')
    st.sidebar.markdown('# Settings')


def plot_grn(celltype, timepoint) -> 'plotly.express.imshow':
    """Returns a plotly figure of the GRN clustermap.

    Args:
        celltype: The celltype to plot.
        timepoint: The timepoint to plot.

    Returns:
        fig (plotly.express.imshow): The plotly figure.
    """

    vmax, vmin = 0.1, -0.1
    df_counts, row_linkage, col_linkage = load_data(celltype, timepoint)
    g = sns.clustermap(df_counts, method='ward', metric='euclidean', 
                       cmap='coolwarm', standard_scale=None, 
                       row_cluster=True, col_cluster=True, 
                       xticklabels=df_counts.columns.tolist(),
                       yticklabels=df_counts.index.tolist(), 
                       vmax=vmax, vmin=vmin, 
                       row_linkage=row_linkage, col_linkage=col_linkage)

    # Adjust the font size for the labels
    g.ax_heatmap.tick_params(axis='x', labelsize=1)
    g.ax_heatmap.tick_params(axis='y', labelsize=1)

    # hide the dendrograms
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)

    # Convert seaborn plot to plotly and plot on page
    fig = px.imshow(g.data2d,
                    color_continuous_scale='balance')
    fig.update_layout(width=1000, height=800)

    return fig


def main():
    """
    """

    timepoints = {'TDR126': '0somites',
                'TDR127': '5somites',
                'TDR128': '10somites',
                'TDR118': '15 somites',
                'TDR125': '20somites',
                'TDR124': '30 somites'}
    celltypes = ['fast_muscle', 'neural_posterior', 'NMPs', 'PSM', 'somites', 'spinal_cord', 'tail_bud']

    st_setup()
    celltype = st.sidebar.selectbox(
        'Select a cell type to plot',
        celltypes
    )

    if "timepoints" not in st.session_state:
        st.session_state.selectboxes = [0]

    if st.sidebar.button("Add timepoint"):
        st.session_state.selectboxes.append(len(st.session_state.selectboxes))

    selected_timepoints = []
    for i, key in enumerate(st.session_state.selectboxes):
        st.sidebar.selectbox(
            'Select a timepoint to plot',
            list(timepoints.keys()),
            key=key
        )
        selected_timepoints.append(st.session_state[key])

    # Plot the clustermap
    fig = plot_grn(celltype, 'TDR126')
    st.plotly_chart(fig)


if __name__ == '__main__':
    main()
