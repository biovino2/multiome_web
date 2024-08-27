"""Creates single figure containing both ccan and gene track plot using plotly.

Ben Iovino  08/27/24    CZ-Biohub
"""

import plotly.graph_objects as go
import polars as pl


def get_data(option: str) -> 'tuple[pl.DataFrame, pl.DataFrame]':
    """Returns subset of data for the gene of interest.

    Args:
        option (str): The gene of interest.

    Returns:
        tuple: A tuple containing the gene data and ATAC data
    """

    genes = pl.read_csv('ccan/data/GRCz11.csv')
    access = pl.read_csv('ccan/data/access.csv')
    gene_data = genes.filter(pl.col('gene_name') == option).sort('end')
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

    # Create a line plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[start, end],
        y=[0, 0],
        mode="lines",
        line=dict(color="black", width=2),
        showlegend=False
    ))  

    # Add x-ticks
    tick_length = 0.015  # Increase this value to make ticks longer
    tick_width = 1.25   # Thickness of the ticks
    chromosome_ticks = list(range(start, end + 1, 5000))
    for tick in chromosome_ticks:
        fig.add_trace(go.Scatter(
            x=[tick, tick],
            y=[-tick_length / 2, tick_length / 2],  # Adjusted for longer ticks
            mode="lines",
            line=dict(color="black", width=tick_width),  # Thicker ticks
            showlegend=False
        ))

    # Set x-axis ticks and labels
    fig.update_xaxes(
        range=[start - 1, end + 1],
        title="Genomic Coordinate (kbp)",
        showgrid=False,
        tickvals=chromosome_ticks,
        ticktext=[f"{tick // 1000}" for tick in chromosome_ticks]
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
        line=dict(color="gray", width=2),
        showlegend=False
    ))

    # Draw rectangles representing exons
    exon_color = "dodgerblue"
    for exon_start, exon_end in exon_positions:
        fig.add_trace(go.Scatter(
            x=[exon_start, exon_end, exon_end, exon_start, exon_start],
            y=[0.075, 0.075, 0.125, 0.125, 0.075],
            fill="toself",
            fillcolor=exon_color,
            mode="lines",
            line=dict(color=exon_color, width=0.5),
            showlegend=False
        ))

    # Draw arrow beneath plot showing direction of transcription
    arrow_length = 0.025 * (end-start)
    if gene_data['strand'][0] == '+':
        arrow = genomic_end + arrow_length
    else:
        arrow = genomic_start - arrow_length
    fig.add_annotation(
            x=arrow,
            y=0.1,
            ax=genomic_end-25,
            ay=0.1,
            xref='x',
            yref='y',
            axref='x',
            ayref='y',
            showarrow=True,
            arrowhead=3,
            arrowsize=1,
            arrowwidth=2,
            arrowcolor="gray"
        )
    
    return fig


def plot_atac(fig: go.Figure, atac_data: pl.DataFrame) -> go.Figure:
    """Plots the ATAC data on the figure.

    Args:
        fig (go.Figure): The plotly figure.
        atac_data (pl.DataFrame): The ATAC data.

    Returns:
        go.Figure: The plotly figure with the ATAC data.
    """

    heights = {'TDR126': 0.20, 'TDR127': 0.225, 'TDR128': 0.25,
                'TDR118': 0.275, 'TDR125': 0.30, 'TDR124': 0.325}
    colors = {'TDR126': '#440154', 'TDR127': '#414487', 'TDR128': '#2A788E',
               'TDR118': '#22A884', 'TDR125': '#7AD151', 'TDR124': '#FDE725'}

    # Plot each individual peak, color and height based on timepoint
    for row in atac_data.iter_rows():
        atac_start = row[2]
        atac_end = row[3]
        height = heights[row[1]]
        fig.add_trace(go.Scatter(
            x=[atac_start, atac_end],
            y=[height, height],
            mode="lines",
            line=dict(color=colors[row[1]], width=10),
            showlegend=False
        ))

    return fig


def combined_plot(option: str) -> go.Figure:
    """Returns a plotly figure with the gene track and ATAC data for the gene of interest.

    Args:
        option (str): The gene of interest.

    Returns:
        go.Figure: A plotly figure.
    """

    gene_data, atac_data = get_data(option)

    # Find minimum and maximum positions between both dataframes
    start = min(gene_data['start'].min(), atac_data['start'].min())
    end = max(gene_data['end'].max(), atac_data['end'].max())

    # Plot each component
    fig = plot_genomic_region(start, end)
    fig = plot_gene(fig, gene_data, start, end)
    fig = plot_atac(fig, atac_data)

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
