"""Plots peaks and gene track using plotly, along with ZFIN information.

Ben Iovino  08/27/24    CZ-Biohub
"""

import os
import streamlit as st
import plotly.graph_objects as go
import polars as pl
from util import get_timepoints


def get_data(path: str, option: str) -> 'tuple[pl.DataFrame, pl.DataFrame]':
    """Returns subset of data for the gene of interest.

    Args:
        path (str): The path to the data.
        option (str): The gene of interest.

    Returns:
        tuple: A tuple containing the gene data and ATAC data
    """

    genes = pl.read_csv(f'{path}/GRCz11.csv')
    access = pl.read_csv(f'{path}/access.csv')
    gene_data = genes.filter(pl.col('gene_name') == option)
    atac_data = access.filter(pl.col('gene_name') == option)

    return gene_data, atac_data


def plot_genomic_region(start: int, end: int) -> go.Figure:
    """Returns a plotly figure with a line plot representing the genomic region.

    Args:
        start (int): The start of the region.
        end (int): The end of the region.

    Returns:
        go.Figure: A plotly figure.
    """

    fig = go.Figure()

    # Set x-axis ticks and labels
    fig.update_xaxes(
        title="Genomic Coordinate (bp)",
        autorange=True,
        ticklen=10,
        tickcolor='black',
        showline=True,
        tickformat=",.0f",
    )

    return fig


def plot_gene(fig: go.Figure, gene_data: pl.DataFrame, start: int, end: int) -> go.Figure:
    """Plots the gene track on the figure.
    
    Args:
        fig (go.Figure): The plotly figure.
        gene_data (pl.DataFrame): The gene data.
        start (int): The start of the genomic region.
        end (int): The end of the genomic region.

    Returns:
        go.Figure: The plotly figure with the gene track.
    """

    exon_positions = gene_data[['start', 'end']].to_numpy()
    genomic_start = gene_data['start'].min()
    genomic_end = gene_data['end'].max()

    # Draw line representing gene region
    fig.add_trace(go.Scatter(
        x=[genomic_start, genomic_end],
        y=[0.1, 0.1],
        mode="lines",
        hoverinfo='none',
        line=dict(color="gray", width=2),
        showlegend=False
    ))

    # Draw rectangles representing exons
    for exon_start, exon_end in exon_positions:
        fig.add_shape(  # Blue rectangle for exon
            type="rect",
            x0=exon_start,
            y0=0.11,
            x1=exon_end,
            y1=0.09,
            line=dict(color="dodgerblue"),
            fillcolor="dodgerblue"
        )
        fig.add_trace(go.Scatter(  # Exon text
            x=[exon_start, exon_end, exon_end, exon_start, exon_start],
            y=[0.11, 0.11, 0.09, 0.09, 0.11],
            fill='toself',
            mode="lines",
            hoverinfo='text',
            text=f"<b>Start:</b> {exon_start}<br><b>End:</b> {exon_end}",
            showlegend=False
        ))

    # Draw arrow showing direction of transcription
    if gene_data['strand'][0] == '+':
        arrow = genomic_end
        ax = genomic_start
    else:
        arrow = genomic_start
        ax = genomic_end
    fig.add_annotation(
            x=arrow,
            y=0.075,
            ax=ax,
            ay=0.075,
            xref='x',
            yref='y',
            axref='x',
            ayref='y',
            showarrow=True,
            arrowhead=3,
            arrowsize=1,
            arrowwidth=2,
            arrowcolor="lightgray"
        )

    return fig


def plot_atac(fig: go.Figure, atac_data: pl.DataFrame) -> go.Figure:
    """Plots the ATAC data on the figure.

    Args:
        fig (go.Figure): The plotly figure.
        atac_data (pl.DataFrame): The ATAC data.
        timepoints (dict): A dictionary mapping timepoints to their names.
        colors (dict): A dictionary mapping timepoints to their colors.

    Returns:
        go.Figure: The plotly figure with the ATAC data.
    """

    timepoints = get_timepoints()
    heights = {'TDR126': 0.20, 'TDR127': 0.225, 'TDR128': 0.25,
                'TDR118': 0.275, 'TDR125': 0.30, 'TDR124': 0.325}
    colors = {'TDR126': '#440154', 'TDR127': '#414487', 'TDR128': '#2A788E',
               'TDR118': '#22A884', 'TDR125': '#7AD151', 'TDR124': '#FDE725'}

    # Plot each individual peak, color and height based on timepoint
    for row in atac_data.iter_rows():
        atac_start = row[2]
        atac_end = row[3]
        height = heights[row[1]]
        fig.add_shape(
            type="rect",
            x0=atac_start,
            y0=height-0.01,
            x1=atac_end,
            y1=height+0.01,
            line=dict(color=colors[row[1]]),
            fillcolor=colors[row[1]],
        )
        fig.add_trace(go.Scatter(
            x=[atac_start, atac_end, atac_end, atac_start, atac_start],
            y=[height+0.01, height+0.01, height-0.01, height-0.01, height+0.01],
            fill='toself',
            mode="lines",
            hoverinfo='text',
            text=
            f'<b>Sample:</b> {timepoints[row[1]]}<br>'
            f'<b>Start:</b> {atac_start}<br>'
            f'<b>End:</b> {atac_end}<br>',
            showlegend=False
        ))

    return fig


def plot_legend(fig: go.Figure) -> go.Figure:
    """Adds a legend for the transcription direction, exons, and ATAC peaks, individually.

    Args:
        fig (go.Figure): The plotly figure.

    Returns:
        go.Figure: The plotly figure with the legend.
    """

    timepoints = get_timepoints()
    colors = {'TDR124': '#FDE725', 'TDR125': '#7AD151', 'TDR118': '#22A884',
                    'TDR128': '#2A788E', 'TDR127': '#414487', 'TDR126': '#440154'}

    # Transcription direction
    fig.add_trace(go.Scatter(
        x=[None], y=[None],  # Draw the arrow line
        mode='lines+markers',
        line=dict(color='lightgray', width=3),  # Line for the arrow's body
        marker=dict(symbol='triangle-right', size=12, color='lightgray'),  # Arrowhead marker
        name='Transcription direction'
    ))

    # Exons
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(size=12, color='dodgerblue', symbol='square'),
        name='Exon'
    ))

    # ATAC peaks
    for key, color in colors.items():
        fig.add_trace(go.Scatter(
            x=[None], y=[0.15, 0.15],
            mode='markers',
            marker=dict(size=12, color=color, symbol='square'),
            name=f"ATAC peak - {timepoints[key]}"
        ))

    return fig


def combined_plot(path: str, option: str) -> go.Figure:
    """Returns a plotly figure with the gene track and ATAC data for the gene of interest.

    Args:
        path (str): The path to the data.
        option (str): The gene of interest.

    Returns:
        go.Figure: A plotly figure.
    """

    gene_data, atac_data = get_data(path, option)

    # Find minimum and maximum positions between both dataframes
    start = min(gene_data['start'].min(), atac_data['start'].min())
    end = max(gene_data['end'].max(), atac_data['end'].max())

    # Plot each component
    fig = plot_genomic_region(start, end)
    fig = plot_gene(fig, gene_data, start, end)
    fig = plot_atac(fig, atac_data)
    fig = plot_legend(fig)

    # Auto scale both axes to fit the data
    fig['layout']['yaxis'].update(autorange=True)
    fig['layout']['xaxis'].update(autorange=True)

    # Set y-axis limits and hide y-axis
    fig.update_yaxes(range=[-1, 1], visible=False, showgrid=False)

    # Remove background and gridlines
    fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
        xaxis=dict(showgrid=False),     # Hide gridlines
        yaxis=dict(showgrid=False)      # Hide gridlines
    )

    return fig


def load_zfin_info(path: str) -> 'tuple[dict, dict]':
    """Loads ZFIN gene name and information from file.

    Args:
        path (str): Path to the ZFIN data.
    
    Returns:
        dict: A dictionary mapping gene names to ZFIN names.

    """

    # Gene name mapping to ZFIN ID
    mapping: dict[str, str] = {}
    with open(f'{path}/zfin/mapping.txt', 'r') as file:
        for line in file:
            line = line.strip().split()
            try:
                mapping[line[0]] = line[1]
            except IndexError:
                continue

    # Gene information from ZFIN
    info: dict[str, str] = {}
    with open(f'{path}/zfin/info.txt', 'r') as file:
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
    st.write("For each gene, we plot the peaks and protein coding region. The peaks we show are \
             highly correlated with the transcription start site of the selected gene. For more \
             information about the gene and it's regulatory elements, we also provide links to \
             ZFIN and Ensembl, if they are available.")
            
    st.sidebar.markdown('# Settings')

    # Initialize drop down box
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


### Main

# Read in data
path = os.path.dirname(os.path.abspath(__file__))+'/data'
df = pl.read_csv(f'{path}/access.csv')
gene_names = pl.read_csv(f'{path}/gene_names.csv')['names'].to_list()
mapping, info = load_zfin_info(path)

# Set up streamlit, get input
option = st_setup(gene_names) 

# Get gene info and plot gene track
df = pl.read_csv(f'{path}/GRCz11.csv')
chrom = df.filter(pl.col('gene_name') == option)['seqname'][0]
st.markdown(f'### {option} (chr{chrom})')

# Display
try:
    st.markdown(f"[(ZFIN) {mapping[option]}](https://zfin.org/{mapping[option]}): {info[option]}")
except KeyError:
    st.write("No ZFIN information available.")
fig  = combined_plot(path, option)
st.plotly_chart(fig)
