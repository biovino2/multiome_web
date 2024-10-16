"""Commands used for gathering and preprocessing data in data directory

Ben Iovino  08/07/24    CZ-Biohub
"""

from bs4 import BeautifulSoup
import os
from gtfparse import read_gtf
import polars as pl
import re
from urllib.request import urlretrieve


def merge_exons(df: 'pl.DataFrame') -> 'pl.DataFrame':
    """Returns a DataFrame with no overlapping exons for each gene

    Args:
        df (pl.DataFrame): The DataFrame with all exons.

    Returns:
        pl.DataFrame: A DataFrame with no overlapping exons.
    """

    new_df = []
    genes = df['gene_name'].unique().to_list()
    for gene in genes:  # Sort exons for each gene
        gene_df = df.filter(pl.col('gene_name') == gene)
        
        # Merge overlapping exons by sorting and iterating
        rows = []
        gene_df = gene_df.sort('start')
        for row in gene_df.iter_rows():
            if not rows:  # initialize list of exons
                rows.append(list(row))
                continue

            # Merge if current start < previous end
            # If current end > previous end, update previous end
            if row[2] < rows[-1][3]:
                if rows[-1][3] < row[3]:
                    rows[-1][3] = row[3]

            # Otherwise, we have a new exon
            else:
                rows.append(list(row))
        for row in rows:
            new_df.append(row)

    # Create new DataFrame containing all genes
    return pl.DataFrame({
        'gene_name': pl.Series(values=[row[0] for row in new_df], dtype=pl.String),
        'seqname': pl.Series(values=[row[1] for row in new_df], dtype=pl.String),
        'start': pl.Series(values=[row[2] for row in new_df], dtype=pl.Int32),
        'end': pl.Series(values=[row[3] for row in new_df], dtype=pl.Int32),
        'strand': pl.Series(values=[row[4] for row in new_df], dtype=pl.String),
        'gene_id': pl.Series(values=[row[5] for row in new_df], dtype=pl.String)
    })


def parse_gtf(path: str, filename: str) -> 'list[str]':
    """Parses GTF file for relevant information.

    Args:
        path (str): The path to the GTF file.
        filename (str): The GTF filename.

    Returns:
        list[str]: A list of all genes in GTF file.
    """

    df = read_gtf(f'{path}/{filename}', usecols=['gene_name', 'seqname', 'start', 'end',
                                                  'strand', 'feature', 'gene_id'])
    df = df.filter(pl.col('seqname') != 'MT')  # Remove mitochondrial genes
    df = df.filter(pl.col('feature') == 'exon')  # Only keep exons
    df = df.drop('feature')
    df = merge_exons(df)
    df.write_csv(f'{path}/{filename.split(".")[0]}.csv')

    # Return gene names
    return df['gene_name'].unique().to_list()


def parse_atac(path: str):
    """Combines ATAC data into one file (path/access.csv).

    Args:
        path (str): The path to the ATAC files.
    """

    timepoints = ['TDR126', 'TDR127', 'TDR128', 'TDR118', 'TDR125', 'TDR124']

    # Read ATAC data for each timepoint
    rows: 'list[tuple[str, str, int, int]]' = []
    for file in timepoints:
        df = pl.read_csv(f'{path}/{file}.csv')
    
        # Extract gene, start, and end
        for row in df.iter_rows():
            gene = row[1]
            start = int(row[0].split('_')[1])
            end = int(row[0].split('_')[2])
            rows.append((gene, file, start, end))

    # Add rows to dataframe
    df = pl.DataFrame({
    'gene_name': pl.Series(values=[row[0] for row in rows], dtype=pl.String),
    'sample': pl.Series(values=[row[1] for row in rows], dtype=pl.String),
    'start': pl.Series(values=[row[2] for row in rows], dtype=pl.Int32),
    'end': pl.Series(values=[row[3] for row in rows], dtype=pl.Int32)
    })

    df.write_csv(f'{path}/access.csv')


def match_names(path: str, gene_names: 'list[str]'):
    """Writes matching gene names to file.

    Args:
        path (str): The path to save the gene names.
        gene_names (list): The list of gene names from the gtf file.
    """

    # Find shared names
    peaks = pl.read_csv(f'{path}/access.csv')
    peak_names = peaks['gene_name'].unique().to_list()
    shared_names = [name for name in sorted(gene_names) if name in peak_names]

    # Write to file
    df = pl.DataFrame({"names": shared_names})
    df.write_csv(f'{path}/gene_names.csv')


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
    """Three parts:
        1) Combines ATAC data into one file
        2) Parses GTF file for exons
        3) Downloads and parses ZFIN gene information

    Comment out any parts that you don't need (particularly part 3, lots of downloads).
    """

    path = 'src/dyn/data'
    if not os.path.exists(path):
        os.mkdir(path)

    parse_atac(path)  # part 1
    gene_names = parse_gtf(path, 'GRCz11.gtf.gz')  # part 2
    match_names(path, gene_names)
    #get_zfin_genes(path, gene_names)  # part 3
    #get_zfin_annotations(path)
    #parse_html(path)


if __name__ == "__main__":
    main()
