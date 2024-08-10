"""Plot GRN cluster maps.

Ben Iovino  08/09/24    CZ-Biohub
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st



def main():
    """
    """

    timepoints = {'TDR126': '0somites',
                'TDR127': '5somites',
                'TDR128': '10somites',
                'TDR118': '15 somites',
                'TDR125': '20somites',
                'TDR124': '30 somites'}
    vmax, vmin = 0.1, -0.1

    # For loop to plot the GRNs over timepoints
    # Choose one celltype
    celltype = "PSM"
    timepoint = "TDR126"

    # Load celltype + counts and links
    df_counts = pd.read_csv(f"data/{celltype}/{celltype}_{timepoint}.csv", index_col=0)
    row_linkage = np.load(f"data/{celltype}/{celltype}_row_linkage.npz")['arr_0']
    col_linkage = np.load(f"data/{celltype}/{celltype}_col_linkage.npz")['arr_0']
    
    # plot the clustermap
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

    if not os.path.exists('output'):
        os.makedirs('output')

    # Save the figure
    plt.savefig(f"output/{celltype}_{timepoints[timepoint]}2.png", dpi=300)


if __name__ == '__main__':
    main()
