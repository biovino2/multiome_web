"""Commands used for gathering and preprocessing data in data directory

Ben Iovino  08/07/24    CZ-Biohub
"""

from bs4 import BeautifulSoup
import os
from gtfparse import read_gtf
import polars as pl
import re
from urllib.request import urlretrieve


def parse_gtf(path: str, filename: str) -> 'list[str]':
    """Parses GTF file for relevant information.

    Args:
        path (str): The path to the GTF file.
        filename (str): The GTF filename.

    Returns:
        list[str]: A list of all genes in GTF file.
    """

    df = read_gtf(f'{path}/{filename}', usecols=['gene_name', 'seqname', 'start', 'end', 'strand', 'feature'])
    df = df.filter(pl.col('seqname') != 'MT')
    df = df.filter(pl.col('feature') == 'exon')
    df = df.drop('feature')
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


def parse_html(path: str):
    """Parses ZFIN annotations for gene information and writes to file.

    Args:
        path (str): The path to the ZFIN annotations.
    """

    info: dict[str, str] = {}
    key_words = ['Predicted', 'predicted',
                'Involved', 'involved',
                'Expressed', 'Is expressed',
                'Orthologous', 'orthologous',
                'Exhibits', 'exhibits']
    
    # Parse each HTML File
    regex = r'^(.*)\.html$'  # capture everything before .html extension
    for html in os.listdir(f'{path}/zfin/annotations'):
        with open(f'{path}/zfin/annotations/{html}', 'r') as file:
            html_content = file.read()

        # Parse for all dd tags
        soup = BeautifulSoup(html_content, 'html.parser')
        dd_tags = soup.find_all('dd', class_='col-sm-10')
        for tag in dd_tags:
            
            # If any key words are found in tag, save to dict
            if any(word in tag.text for word in key_words):
                gene = re.match(regex, html).group(1)
                info[gene] = tag.text.strip()
                break

    # Write info to file
    with open(f'{path}/zfin/info.txt', 'w') as file:
        for key, value in info.items():
            file.write(f'{key}\t{value}\n')


def main():
    """Two parts: reads GTF file and downloads ZFIN annotations.
    """

    path = 'ccan/data'
    if not os.path.exists(path):
        os.mkdir(path)

    gene_names = parse_gtf(path, 'GRCz11.gtf.gz')
    get_zfin_genes(path, gene_names)
    get_zfin_annotations(path)
    parse_html(path)


if __name__ == "__main__":
    main()
