"""Master script for plotting all figures on one page.

Ben Iovino  08/06/24    CZ-Biohub
"""

from plot_ccan import load_dfs, plot_ccans_genomic_loci, get_color_dict, get_gene_names
from plot_track import get_gene_info, define_config
from figeno import figeno_make
import streamlit as st
import polars as pl
from plot_both import combine_plot


def load_zfin_info() -> 'tuple[dict, dict]':
    """Loads ZFIN gene name and information from file.
    
    Returns:
        dict: A dictionary mapping gene names to ZFIN names.

    """

    # Gene name mapping to ZFIN ID
    mapping: dict[str, str] = {}
    with open('ccan/data/zfin/mapping.txt', 'r') as file:
        for line in file:
            line = line.strip().split()
            try:
                mapping[line[0]] = line[1]
            except IndexError:
                continue

    # Gene information from ZFIN
    info: dict[str, str] = {}
    with open('ccan/data/zfin/info.txt', 'r') as file:
        for line in file:
            line = line.split('\t')
            try:
                info[line[0]] = line[1]
            except IndexError:
                continue

    return mapping, info


def st_setup(gene_names: list[str]):
    """Initializes streamlit session.

    Args:
        gene_names (list): The list of gene names
    """

    st.set_page_config(layout="wide")
    st.title('Cis Co-Accessibility and Gene Track Plots')
    st.write('For each gene, we plot the cis co-accessible peaks at different timepoints, as well as the gene track.\
             The cis co-accessible peaks are represented along the chromosome on the x-axis, with the timepoints \
             stacked from earliest to latest on the y-axis. The gene track (made with figeno) shows the gene body \
             where the blue regions represent the exons and the blue arrow represents the direction of transcription.')
    st.write("We also provide a link to the ZFIN page for each gene, as well as ZFIN's annotation of the gene, if available.")
    st.sidebar.markdown('# Settings')
    option = st.sidebar.selectbox(
        'Select a gene to plot',
        gene_names
    )
    
    return option


def main():
    """
    """

    # Define the timepoints and load data in order
    timepoints = {"0 hours post fertilization ": 'TDR126', "5 hours post fertilization": 'TDR127',
                "10 hours post fertilization": 'TDR128', "15 hours post fertilization": 'TDR118',
                "20 hours post fertilization": 'TDR125', "30 hours post fertlization": 'TDR124'}
    df_list = load_dfs('ccan/data', list(timepoints.values()))
    color_dict = get_color_dict(timepoints)
    gene_names = get_gene_names(df_list)
    mapping, info = load_zfin_info()

    # Set up streamlit, get input
    option = st_setup(gene_names) 

    # Get gene info and plot gene track
    df = pl.read_csv('ccan/data/GRCz11.csv')
    chrom, min, max, strand = get_gene_info(df, option)

    # Create and plot ccan plot
    fig, ccan_start = plot_ccans_genomic_loci(df_list,
                                gene_name=option,
                                direction=strand,
                                timepoints=list(timepoints.keys()),
                                colordict=color_dict)
    config = define_config(chrom, min, max, option)
    st.markdown(f'# {option}')
    figeno_make(config)

    # Display
    try:
        st.markdown(f"[ZFIN](https://zfin.org/{mapping[option]}): {info[option]}")
    except KeyError:
        st.write("No ZFIN information available.")
    st.plotly_chart(combine_plot(option))
    st.markdown('## Cis Co-Accessibility')
    st.pyplot(fig)
    st.markdown('## Gene Track')
    st.image("ccan/data/figure.png")


if __name__ == "__main__":
    main()
