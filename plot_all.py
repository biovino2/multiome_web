"""Master script for plotting all figures on one page.

Ben Iovino  08/06/24    CZ-Biohub
"""

from plot_ccan import load_dfs, plot_ccans_genomic_loci, get_color_dict, get_gene_names
from plot_track import get_gene_info, define_config
from figeno import figeno_make
import streamlit as st
import polars as pl

def main():
    """
    """

    # Define the timepoints and load data in order
    timepoints = {"0 budstage": 'TDR126', "5 somites": 'TDR127', "10 somites": 'TDR128',
                   "15 somites": 'TDR118', "20 somites": 'TDR125', "30 somites": 'TDR124'}
    df_list = load_dfs('data', list(timepoints.values()))
    color_dict = get_color_dict(timepoints)
    gene_names = get_gene_names(df_list)

    # Take user input for gene name
    option = st.selectbox(
        'Select a gene to plot',
        gene_names
    )

    # Create and plot ccan plot
    fig = plot_ccans_genomic_loci(df_list,
                                gene_name=option,
                                timepoints=list(timepoints.keys()),
                                colordict=color_dict)
    st.pyplot(fig)

    # Get gene info and plot gene track
    df = pl.read_csv('data/GRCz11.csv')
    chrom, min, max = get_gene_info(df, option)
    config = define_config(chrom, min, max)
    figeno_make(config)
    st.image("figure.png")


if __name__ == "__main__":
    main()