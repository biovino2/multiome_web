import polars as pl


def combine_plot(option: str):

    genes = pl.read_csv('ccan/data/GRCz11.csv')
    access = pl.read_csv('ccan/data/access.csv')

    gene = option

    # Select the gene of interest
    gene_data = genes.filter(pl.col('gene_name') == gene)
    atac_data = access.filter(pl.col('gene_name') == gene)

    # Find minimum and maximum positions between both dataframes
    start = min(gene_data['start'].min(), atac_data['start'].min())
    end = max(gene_data['end'].max(), atac_data['end'].max())

    import plotly.graph_objects as go

    # Create a line plot
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=[start, end],
        y=[0, 0],
        mode="lines",
        line=dict(color="gray", width=2),
        showlegend=False
    ))  

    # Add x-ticks
    tick_length = 0.1  # Increase this value to make ticks longer
    tick_width = 2   # Thickness of the ticks
    chromosome_ticks = list(range(start, end + 1, 5000))

    for tick in chromosome_ticks:
        fig.add_trace(go.Scatter(
            x=[tick, tick],
            y=[-tick_length / 2, tick_length / 2],  # Adjusted for longer ticks
            mode="lines",
            line=dict(color="black", width=tick_width),  # Thicker ticks
            showlegend=False
        ))

    # Set y-axis limits and hide y-axis
    fig.update_yaxes(range=[-1, 1], visible=False, showgrid=False)

    # Set x-axis ticks
    fig.update_xaxes(
        range=[start - 1, end + 1],
        title="Genomic Coordinate (kbp)",
        showgrid=False,
        tickvals=chromosome_ticks,
        ticktext=[f"{tick // 1000}" for tick in chromosome_ticks]
    )

    # Set plot size and remove background color
    fig.update_layout(
        margin=dict(l=50, r=50, t=50, b=50),
        plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
        xaxis=dict(showgrid=False),     # Hide gridlines
        yaxis=dict(showgrid=False)      # Hide gridlines
    )

    exon_positions = gene_data[['start', 'end']].to_numpy()
    genomic_start = gene_data['start'].min()
    genomic_end = gene_data['end'].max()

    # Draw genomic region line
    fig.add_trace(go.Scatter(
        x=[genomic_start, genomic_end],
        y=[0.5, 0.5],
        mode="lines",
        line=dict(color="gray", width=2),
        showlegend=False
    ))

    # Draw exons as blue rectangles
    for exon_start, exon_end in exon_positions:  # For all exons except the last
        fig.add_trace(go.Scatter(
            x=[exon_start, exon_end, exon_end, exon_start, exon_start],
            y=[0.4, 0.4, 0.6, 0.6, 0.4],
            fill="toself",
            fillcolor="blue",
            mode="lines",
            line=dict(color="blue", width=0.5),
            showlegend=False
        ))

    heights = {'TDR126': 0.8, 'TDR127': 0.9, 'TDR128': 1, 'TDR118': 1.1, 'TDR125': 1.2, 'TDR124': 1.3}
    colors = {'TDR126': '#440154', 'TDR127': '#414487', 'TDR128': '#2A788E', 'TDR118': '#22A884', 'TDR125': '#7AD151', 'TDR124': '#FDE725'}

    # Overlay the ATAC bars on top
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

    fig['layout']['yaxis'].update(autorange=True)

    return fig
