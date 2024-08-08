"""Loads the data from the csv file and plots a gene track plot on streamlit.

Ben Iovino  08/06/24    CZ-Biohub
"""

import polars as pl


def define_config(chrom: str, min: int, max: int, gene: str) -> 'dict[str, dict]':
    """Returns a dictionary with the configuration for the gene rack plot.

    Args:
        chrom (str): The chromosome.
        min (int): The start coordinate.
        max (int): The end coordinate.
        gene (str): The gene name.

    Returns:
        dict[str, dict]: The configuration dictionary.
    """

    config = {}
    config["general"] = {
            "reference": "custom",
            "layout": "horizontal",
            "genes_file": "ccan/data/GRCz11.gtf.gz"
        }

    config["output"] = {
            "file": "ccan/data/figure.png",
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
			    "genes": f"{gene}"
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
