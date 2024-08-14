"""Commands used for gathering and preprocessing data in data directory

Ben Iovino  08/07/24    CZ-Biohub
"""

import os
from gtfparse import read_gtf
import polars as pl
from urllib.request import urlretrieve


def parse_gtf(path: str, filename: str) -> 'list[str]':
    """Parses GTF file for relevant information.

    Args:
        path (str): The path to the GTF file.
        filename (str): The GTF filename.

    Returns:
        list[str]: A list of all genes in GTF file.
    """

    df = read_gtf(f'{path}/{filename}', usecols=['seqname', 'start', 'end', 'strand', 'gene_name'])
    df = df.filter(pl.col('seqname') != 'MT')
    df.write_csv(f'{path}/{filename.split(".")[0]}.csv')

    # Return gene names
    return df['gene_name'].unique().to_list()


def get_zfin_genes(path: str, gene_names: 'list[str]'):
    """Downloads ZFIN names and creates a mapping between gene names and ZFIN names.

    Args:
        path (str): The path to save the ZFIN names.
        gene_names (list): The list of gene names.
    """

    if not os.path.exists(f'{path}/zfin'):
        os.mkdir(f'{path}/zfin')
    url = 'https://zfin.org/downloads/all_rna_accessions.txt'
    urlretrieve(url, f'{path}/zfin/names.txt')

    # Match gene names to ZFIN names
    mapping = {name: '' for name in gene_names}
    with open(f'{path}/zfin/names.txt', 'r') as file:
        for line in file:
            line = line.strip().split()
            if line[2] in mapping:
                mapping[line[2]] = line[0]

    # Write to file
    with open(f'{path}/zfin/mapping.txt', 'w') as file:
        for key, value in mapping.items():
            key = key.replace('/', '-')  # some gene names have '/'
            file.write(f'{key}\t{value}\n')


def get_zfin_annotations(path: str):
    """Writes ZFINN annotations to file.

    Args:
        path (str): The path to save the ZFIN annotations.
        gene_names (list): The list of gene names.
    """

    if not os.path.exists(f'{path}/zfin/annotations'):
        os.mkdir(f'{path}/zfin/annotations')

    # Download ZFIN annotations (thousands of requests)
    downloaded = os.listdir(f'{path}/zfin/annotations')
    with open(f'{path}/zfin/mapping.txt', 'r') as file:
        for line in file:
            line = line.split()
            try:
                if f'{line[0]}.html' in downloaded:
                    continue
                url = f'https://zfin.org/{line[1]}'
                urlretrieve(url, f'{path}/zfin/annotations/{line[0]}.html')
            except IndexError:
                continue


def main():
    """
    """

    path = 'ccan/data'
    if not os.path.exists(path):
        os.mkdir(path)

    gene_names = parse_gtf(path, 'GRCz11.gtf.gz')
    get_zfin_genes(path, gene_names)
    get_zfin_annotations(path)


if __name__ == "__main__":
    main()
