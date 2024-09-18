"""Master script for setting up streamlit session and plotting plotly figure.

Ben Iovino  08/06/24    CZ-Biohub
"""

import streamlit as st
import polars as pl
from plot_both import combined_plot, plot_legend


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
    st.write("For each gene, we plot the cis co-accessible peaks at different timepoints, as well as the gene track.\
             The cis co-accessible peaks are represented along the chromosome on the x-axis, with the timepoints \
             stacked from earliest to latest on the y-axis. The gene track shows the gene body where the blue regions \
             represent the exons and the arrow beneath represents the direction of transcription. \
            We also provide a link to the ZFIN page for each gene, as well as ZFIN's annotation of the gene, if available.")
    st.sidebar.markdown('# Settings')

    # Initialize drop down box
    gene_names = list(set((gene_names)))  # .unique() from polars messes up select box, set instead
    default = gene_names.index('myf5')
    option = st.sidebar.selectbox(
        'Select a gene to plot',
        gene_names,
        index=default
    )
    
    return option


def main():
    """
    """

    df = pl.read_csv('ccan/data/access.csv')
    gene_names = [gene[0] for gene in df.select('gene_name').rows()]
    mapping, info = load_zfin_info()

    # Set up streamlit, get input
    option = st_setup(gene_names) 

    # Get gene info and plot gene track
    df = pl.read_csv('ccan/data/GRCz11.csv')
    chrom = df.filter(pl.col('gene_name') == 'myf5')['seqname'][0]
    st.markdown(f'## {option} (chr{chrom})')

    # Display
    try:
        st.markdown(f"[(ZFIN) {mapping[option]}](https://zfin.org/{mapping[option]}): {info[option]}")
    except KeyError:
        st.write("No ZFIN information available.")
    st.write('The plot below is an explorable figure. You can zoom in and out, pan, and hover over the data points to see more information.')
    fig  = combined_plot(option)
    
    # Create a checkbox to show/hide the plot legend
    if 'hide_legend' not in st.session_state:
        fig = plot_legend(fig)  # Default is to show legend
        st.session_state.hide_legend = True
    hide_function = st.sidebar.checkbox('Hide Legend', value=False)
    if hide_function:  # If activated, plot is regenerated without adding legend
        if st.session_state.hide_legend:
            fig = combined_plot(option)
            st.session_state.hide_legend = False
    else:  # If deactivated, legend is added to the plot
        if not st.session_state.hide_legend:
            fig = plot_legend(fig)
            st.session_state.hide_legend = True

    st.plotly_chart(fig)


if __name__ == "__main__":
    main()
