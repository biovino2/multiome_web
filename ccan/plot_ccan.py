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
    st.write("For each gene, we plot the peaks and protein coding region. \
              We also provide a link to the ZFIN for each gene, as well as ZFIN's annotation, if available.") \
            
    st.sidebar.markdown('# Settings')

    # Initialize drop down box
    gene_names = list(set((gene_names)))  # .unique() from polars messes up select box, set instead
    default = gene_names.index('myf5')
    option = st.sidebar.selectbox(
        'Select a gene to plot',
        gene_names,
        index=default
    )

    # Remove extra space at top of the page
    st.markdown(
    """
    <style>
        /* Remove padding from main block */
        .block-container {
            padding-top: 2rem;
        }
    </style>
    """,
    unsafe_allow_html=True
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
    st.markdown(f'### {option} (chr{chrom})')

    # Display
    try:
        st.markdown(f"[(ZFIN) {mapping[option]}](https://zfin.org/{mapping[option]}): {info[option]}")
    except KeyError:
        st.write("No ZFIN information available.")
    fig  = combined_plot(option)
    
    # Create a checkbox to show/hide the plot legend
    if 'hide_legend' not in st.session_state:
        fig = plot_legend(fig)  # Default is to show legend
        st.session_state.hide_legend = True
    hide_function = st.sidebar.checkbox('Hide Legend', value=False)
    
    # If activated, plot is regenerated without adding legend
    if hide_function:
        if st.session_state.hide_legend:
            fig= combined_plot(option)
            st.session_state.hide_legend = False
    
    # If deactivated, legend is added to the plot
    else: 
        if not st.session_state.hide_legend:
            fig = plot_legend(fig)
            st.session_state.hide_legend = True

    st.plotly_chart(fig)


if __name__ == "__main__":
    main()
