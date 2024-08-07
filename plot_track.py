"""Loads the data from the csv file and plots a gene track plot on streamlit.

Ben Iovino  08/06/24    CZ-Biohub
"""

from figeno import figeno_make
import polars as pl
import streamlit as st
from plot_ccan import load_dfs, get_gene_names


def define_config(chrom: str, min: int, max: int) -> 'dict[str, dict]':
    """Returns a dictionary with the configuration for the gene rack plot.

    Args:
        chrom (str): The chromosome.
        min (int): The start coordinate.
        max (int): The end coordinate.

    Returns:
        dict[str, dict]: The configuration dictionary.
    """

    config = {}
    config["general"] = {
            "reference": "custom",
            "layout": "horizontal",
            "genes_file": "data/genes.gtf.gz"
        }

    config["output"] = {
            "file": "figure.png",
            "dpi": 300,
            "width": 180
        }

    config["regions"] = [
            {
                "chr": f"{chrom}",
                "start": f"{min}",
                "end": f"{max}",
                "color": "#f4a460"
            }
    ]
    
    config["tracks"] = [
            {
                "type": "genes",
                "height": 10,
                "margin_above": 1.5,
                "bounding_box": False,
                "fontscale": 1,
                "label": "",
                "style": "default",
                "collapsed": True,
			    "only_protein_coding": False,
			    "exon_color": "#2980b9",
			    "genes": "auto"
            },
            {
                "type": "chr_axis",
			    "height": 10,
			    "margin_above": 1.5,
			    "bounding_box": False,
			    "fontscale": 1,
			    "label": "",
			    "label_rotate": False,
			    "style": "default",
			    "unit": "kb",
			    "ticklabels_pos": "below",
			    "ticks_interval": "auto"
            }

    ]

    return config


def get_gene_info(df: pl.DataFrame, gene_name: str) -> 'tuple[str, int, int]':
    """Returns information about the gene for the gene track plot.
    
    Args:
        df (pl.DataFrame): The dataframe to subset.
        gene_name (str): The gene name to subset.
        
    Returns:
        tuple(str, int, int): The chromosome, start, and end coordinates for the gene.
    """

    gene = df.filter(pl.col('gene_name') == gene_name)
    chrom = gene['seqname'].to_list()[0]
    min = gene['start'].min()
    max = gene['end'].max()
    strand = gene['strand'].to_list()[0]

    return chrom, min, max, strand


def main():
    """
    """

    df = pl.read_csv('data/GRCz11.csv')
    timepoints = {"0 budstage": 'TDR126', "5 somites": 'TDR127', "10 somites": 'TDR128',
                   "15 somites": 'TDR118', "20 somites": 'TDR125', "30 somites": 'TDR124'}
    df_list = load_dfs('data', list(timepoints.values()))
    gene_names = get_gene_names(df_list)

    # Take user input for gene name
    option = st.selectbox(
        'Select a gene to plot',
        gene_names
    )
    chrom, min, max, _ = get_gene_info(df, option)
    config = define_config(chrom, min, max)
    
    fp = figeno_make(config)
    st.write(fp)
    st.image("figure.png")

    
if __name__ == "__main__":
    main()
