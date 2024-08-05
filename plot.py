"""Loads the data from the csv file and plots CCAN on streamlit.

Ben Iovino  08/01/24    CZ-Biohub
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import ScalarFormatter
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

    # Subset peak ID
    peak_id = df['peak_id']
    chrom = peak_id.str.split('_').str[0].astype(str)
    start = peak_id.str.split('_').str[1].astype(int)/1000
    end = peak_id.str.split('_').str[2].astype(int)/1000

    # Assign to DF
    df['chr'] = chrom
    df['start'] = start
    df['end'] = end

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


def get_xticks(beg: int, end: int) -> 'list[int]':
    """Returns xticks for the genomic region, depending on it's size.

    Args:
        beg (int): The beginning of the genomic region.
        end (int): The end of the genomic region.

    Returns:
        list[int]: The xticks for the genomic region plot.
    """

    st.write(beg, end)
    if end - beg <= 5.0:
        mid = (beg + end) / 2
        beg = mid-2.5
        end = mid+2.5
        st.write(beg, end, mid)
        xticks = np.arange(beg, end, 1, dtype=int)
    else:
        diff = end - beg
        xticks = np.arange(beg, end, diff/5)
        xticks = np.append(xticks, end)


    return xticks


def get_color_dict(timepoints: 'dict[str, str]') -> 'dict[str, np.array]':
    """Returns a dictionary mapping timepoints to colors.

    Args:
        timepoints (dict[str, str]): A dictionary mapping timepoints to labels.

    Returns:
        dict[str, np.array]: A dictionary mapping timepoints to colors.
    """

    # Load the "viridis" colormap
    viridis = plt.cm.get_cmap('viridis', 256)

    # Select a subset of the colormap to ensure that "30 somites" is yellow
    # You can adjust the start and stop indices to shift the colors
    start = 50
    stop = 256
    colors = viridis(np.linspace(start/256, stop/256, len(timepoints)))

    # Create a dictionary to map timepoints to colors
    color_dict = dict(zip(timepoints.keys(), colors))

    return color_dict


def plot_ccans_genomic_loci(dfs: 'list[pd.DataFrame]', timepoints: 'list[str]', gene_name: str,
                            colordict: 'dict[str: np.array]', save_fig=False, figpath = None):
    """
    This function takes a list of CCANs (dataframes for each timepoint), and plots the genomic
      region with CCANs for each timepoint for "gene_name".
    
    Parameters:
    1) dfs (list[pd.Dataframe]): A list of dataframes (one dataframe for each timepoint), i.e. [df1, df2, df3]
    2) timepoints (list(str)): A list of timepoint labels corresponding to each dataframe in CCANs
    3) gene_name (str): Name of the gene, i.e. "myf5", "her1"
    4) colordict (dict[str: np.array]): A dictionary of {timepoints:colors (viridis)}
    """

    if len(dfs) != len(timepoints):
        raise ValueError("The number of CCANs dataframes and timepoints labels must be equal")

    # generate a figure object
    fig, ax = plt.subplots(figsize=(10, 2))
    genomic_start, genomic_end = float('inf'), 0

     # one gene for each time point, extract genomic coordinates and plot peaks
    for index, (df, stage) in enumerate(zip(dfs, timepoints)):
        if df[df.gene_short_name == gene_name].empty:  # ignore if gene not present in timepoint
            continue
        df_gene = df[df.gene_short_name == gene_name]  # subset for the gene of interest
        df_gene = compute_start_end(df_gene)

        # update genomic start and end and plot peaks
        genomic_start = min(genomic_start, df_gene['start'].min())
        genomic_end = max(genomic_end, df_gene['end'].max())
        plot_peaks(df_gene, ax, 1 + index*0.1, colordict[stage], stage)

    ax.plot([genomic_start, genomic_end], [1, 1], color='grey', linewidth=2)

    # Set and format xticks
    ax.set_xticks(get_xticks(genomic_start, genomic_end))
    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)

    # Format rest of plot
    ax.set_ylim(0.7, 1 + len(dfs)*0.1 + 0.5)
    ax.set_yticks([])
    ax.set_xlabel('Genomic Coordinate (kbp)')
    ax.set_title(f'Genomic Region Plot for {gene_name} ({df_gene["chr"].iloc[0]})')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()

    if save_fig:
        plt.savefig(figpath + "coverage_plot_CCANs_" + gene_name + ".png")
        plt.savefig(figpath + "coverage_plot_CCANs_" + gene_name + ".pdf")

    # Return plot
    return fig


def main():
    """
    """

    # Define the timepoints and load data in order
    timepoints = {"0 budstage": 'TDR126', "5 somites": 'TDR127', "10 somites": 'TDR128',
                   "15 somites": 'TDR118', "20 somites": 'TDR125', "30 somites": 'TDR124'}
    df_list = load_dfs('data', list(timepoints.values()))
    color_dict = get_color_dict(timepoints)

    # Get all gene names from the data and have user choose
    gene_names = df_list[0]['gene_short_name'].unique()
    option = st.selectbox(
        'Select a gene to plot',
        gene_names
    )

    # Create and plot figure
    fig = plot_ccans_genomic_loci(df_list,
                                gene_name=option,
                                timepoints=list(timepoints.keys()),
                                colordict=color_dict)
    st.pyplot(fig)


if __name__ == "__main__":
    main()
