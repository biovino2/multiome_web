"""Loads the data from the csv file and plots CCAN on streamlit.

Ben Iovino  08/01/24    CZ-Biohub
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
import pandas as pd
import streamlit as st


def load_dfs(path: str, file_list: 'list[str]') -> 'list[pd.DataFrame]':
    """Returns a list of dataframes from the csv files in the given path.

    Args:
        path (str): The path to the csv files.
        file_list (list[str]): File names to load

    Returns:
        list[pd.Dataframe]: A list of dataframes.
    """

    df_list = []
    for file in file_list:
        df = pd.read_csv(f'{path}/{file}.csv')
        df_list.append(df)

    return df_list


def compute_start_end(df: pd.DataFrame):
    """Computes the chromosome, start, and end coordinates for a dataframe.

    Args:
        df (pd.DataFrame): The dataframe to compute the coordinates for.

    Returns:
        pd.DataFrame: The dataframe with the chromosome, start, and end columns.
    
    """

    df['chr'] = df['peak_id'].str.split('_').str[0].astype(str)
    df['start'] = df['peak_id'].str.split('_').str[1].astype(int)
    df['end'] = df['peak_id'].str.split('_').str[2].astype(int)

    return df


def plot_peaks(df: pd.DataFrame, ax: plt.Axes, y_base: int, color: str, label=None):
    """Plots peaks at a specified y-axis basis.

    Args:
        df (pd.DataFrame): The dataframe containing the peaks.
        ax (plt.Axes): The axis to plot the peaks on.
        y_base (int): The y-axis basis to plot the peaks.
        color (str): The color of the peaks.
        label (str): The label of the peaks.
    """

    for _, row in df.iterrows():
        rect = patches.Rectangle((row['start'], y_base), row['end'] - row['start'], 0.1,
                                  linewidth=1, edgecolor=color, facecolor=color, label=label)
        ax.add_patch(rect)
        label = None  # Set to None to prevent repeating the label


def plot_CCANs_genomic_loci(CCANs: 'list[pd.DataFrame]', timepoints: 'list[str]', gene_name: str, colordict: 'dict[str: np.array]', 
                            save_fig=False, figpath = None):
    """
    This function takes a list of CCANs (dataframes for each timepoint), and plots the genomic
      region with CCANs for each timepoint for "gene_name".
    
    Parameters:
    1) CCANs (list[pd.Dataframe]): A list of dataframes (one dataframe for each timepoint), i.e. [df1, df2, df3]
    2) timepoints (list(str)): A list of timepoint labels corresponding to each dataframe in CCANs
    3) gene_name (str): Name of the gene, i.e. "myf5", "her1"
    4) colordict (dict[str: np.array]): A dictionary of {timepoints:colors (viridis)}
    """
    
    if len(CCANs) != len(timepoints):
        raise ValueError("The number of CCANs dataframes and timepoints labels must be equal")
    
    # generate a figure object    
    fig, ax = plt.subplots(figsize=(10, 2))
    
    genomic_start, genomic_end = float('inf'), 0

    for index, (df, stage) in enumerate(zip(CCANs, timepoints)):
        # subset for the "gene_name" (gene of interest)
        df = df[df.gene_short_name == gene_name]
        df = compute_start_end(df)
        
        # # extract the chromosome and genomic coordinates
        # if index==0:
        #     chromosome = df.loc['chr'][0]
        
        genomic_start = min(genomic_start, df['start'].min())
        genomic_end = max(genomic_end, df['end'].max())

        plot_peaks(df, ax, 1 + index*0.1, colordict[stage], stage)
    
    ax.plot([genomic_start, genomic_end], [1, 1], color='grey', linewidth=2)
    
    ax.set_ylim(0.7, 1 + len(CCANs)*0.1 + 0.5)
    ax.set_yticks([])
    ax.set_xlabel('Genomic Coordinate')
    ax.set_title(f'Genomic Region Plot for {gene_name}')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    
    if save_fig==True:
        plt.savefig(figpath + "coverage_plot_CCANs_" + gene_name + ".png")
        plt.savefig(figpath + "coverage_plot_CCANs_" + gene_name + ".pdf")
    
    # Return plot
    return fig


def main():
    """
    """

    # Define the timepoints and load data in order
    timepoints = {"0budstage": 'TDR126', "5somites": 'TDR127', "10somites": 'TDR128',
                   "15somites": 'TDR118', "20somites": 'TDR125', "30somites": 'TDR124'}
    df_list = load_dfs('data', list(timepoints.values()))

    # Load the "viridis" colormap
    viridis = plt.cm.get_cmap('viridis', 256)

    # Select a subset of the colormap to ensure that "30 somites" is yellow
    # You can adjust the start and stop indices to shift the colors
    start = 50
    stop = 256
    colors = viridis(np.linspace(start/256, stop/256, len(timepoints)))

    # Create a dictionary to map timepoints to colors
    color_dict = dict(zip(timepoints, colors))
    
    # Get all gene names from the data and have user choose
    gene_names = df_list[0]['gene_short_name'].unique()
    option = st.selectbox(
        'Select a gene to plot',
        gene_names
    )

    # Create and plot figure
    fig = plot_CCANs_genomic_loci(df_list,
                                gene_name=option,
                                timepoints=list(timepoints.keys()),
                                colordict=color_dict)
    st.pyplot(fig)


if __name__ == "__main__":
    main()
