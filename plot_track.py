from figeno import figeno_make
import polars as pl
import streamlit as st
from plot_ccan import load_dfs


def main():
    """
    """

    df = pl.read_csv('data/GRCz11.csv')
    timepoints = {"0 budstage": 'TDR126', "5 somites": 'TDR127', "10 somites": 'TDR128',
                   "15 somites": 'TDR118', "20 somites": 'TDR125', "30 somites": 'TDR124'}
    df_list = load_dfs('data', list(timepoints.values()))

    # Get all gene names from the data and have user choose
    gene_names = set()
    for dfl in df_list:
        gene_names.update(dfl['gene_short_name'])
    option = st.selectbox(
        'Select a gene to plot',
        gene_names
    )

    gene = df.filter(pl.col('gene_name') == option)
    chrom = gene['seqname'].to_list()[0]
    if gene['strand'][0] == '+':
        min = gene['start'].min()
        max = gene['end'].max()
    else:
        min = gene['end'].max()
        max = gene['start'].min()

    config = {}
    config["general"] = {
            "reference": "custom",
            "layout": "horizontal",
            "genes_file": "data/GRCz11.gtf.gz"
        }

    config["output"] = {
            "file": "figure.png",
            "dpi": 400,
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
    
    figeno_make(config)
    st.image("figure.png")

    
if __name__ == "__main__":
    main()
