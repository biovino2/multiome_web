"""Commands used for gathering and preprocessing data in data directory

Ben Iovino  08/07/24    CZ-Biohub
"""

from gtfparse import read_gtf
import polars as pl


def main():
    """
    """

    # Parse GTF for relevant info
    df = read_gtf('data/GRCz11.gtf.gz', usecols=['seqname', 'start', 'end', 'strand', 'gene_name'])
    df = df.filter(pl.col('seqname') != 'MT')
    df.write_csv('data/GRCz11.csv')


if __name__ == "__main__":
    main()
