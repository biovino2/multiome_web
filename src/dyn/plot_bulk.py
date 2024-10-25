"""Plot bulk gene activity and gene expression over time.

Ben Iovino  10/16/24    CZ-Biohub
"""

import os
import polars as pl
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import pearsonr
import streamlit as st


def plot_timepoints(gene: str, fig: go.Figure) -> go.Figure:
    """Returns a plotly figure (assumed to be a 1x2 subplots object) with two traces in the first
    subplot, one for gene expression and one for gene activity, both over time points.

    Args:
        gene (str): The gene name.
        fig (go.Figure): The plotly figure.

    Returns:
        go.Figure: The updated plotly figure.
    """

    # Subset data for gene expression
    x_tp = gex_values['timepoints']
    y_rna = gex_values[gene]
    yerr_rna = gex_errors[gene]

    # Add scatter plot for gene expression
    fig.add_trace(
        go.Scatter(
            x=x_tp,
            y=y_rna,
            mode='markers',
            marker=dict(color='dodgerblue'),
            legendgroup='group1',
            name='Gene Expression (RNA)',
            showlegend=True,
            hoverinfo='text',
            text=[f'<b>{tp}</b><br><b>Gene Expression:</b> {round(y, 5)}' for tp, y in zip(x_tp, y_rna)],
            error_y=dict(
                type='data',
                array=yerr_rna,
                visible=True,
                thickness=1,
                width=2)
        )
    )
    
    # Subset data for gene activity
    y_atac = atac_values[gene]
    yerr_atac = atac_errors[gene]

    # Add scatter plot for gene activity
    fig.add_trace(
        go.Scatter(
            x=x_tp,
            y=y_atac,
            mode='markers',
            marker=dict(color='orange'),
            legendgroup='group1',
            name='Gene Activity (ATAC)',
            showlegend=True,
            hoverinfo='text',
            text=[f'<b>{tp}</b><br><b>Gene Activity:</b> {round(y, 5)}' for tp, y in zip(x_tp, y_atac)],
            error_y=dict(
                type='data',
                array=yerr_atac,
                visible=True,
                thickness=1,
                width=2)
        )
    )

    return fig


def plot_gene(gene: str, fig: go.Figure) -> go.Figure:
    """Returns a plotly figure (assumed to be a 1x2 subplots object) with six traces in the second
    subplot, one dot with error bars for each timepoint.

    Args:
        gene (str): The gene name.
        fig (go.Figure): The plotly figure.

    Returns:
        go.Figure: The updated plotly figure 
    """

    colors = {'10hpf': '#440154', '12hpf': '#414487', '14hpf': '#2A788E',
               '16hpf': '#22A884', '19hpf': '#7AD151', '24hpf': '#FDE725'}
    timepoints = {'10hpf': '10 hours post fertilization', '12hpf': '12 hours post fertilization',
                    '14hpf': '14 hours post fertilization', '16hpf': '16 hours post fertilization',
                    '19hpf': '19 hours post fertilization', '24hpf': '24 hours post fertilization'}

    x_atac = atac_values[gene]
    y_rna = gex_values[gene]
    x_err = atac_errors[gene]
    y_err = gex_errors[gene]

    # Calculate correlation coefficient
    x_vals = atac_values[gene]
    y_vals = gex_values[gene]
    correlation, _ = pearsonr(x_vals, y_vals)

    # Plot each point with error bars (error bars take one color per trace, need to loop)
    for i in range(len(x_atac)):
        fig.add_trace(
            go.Scatter(
            x=[x_atac[i]],
            y=[y_rna[i]],
            mode='markers',
            marker=dict(color=list(colors.values())[i]),
            legendgroup='group2',
            name=timepoints[list(colors.keys())[i]],
            showlegend=True,
            hoverinfo='text',
            text=[f'<b>{list(timepoints.values())[i]}</b><br>'
                     f'<b>Gene Activity:</b> {round(x_atac[i], 5)}<br>'
                     f'<b>Gene Expression:</b> {round(y_rna[i], 5)}'],
            error_x=dict(
                type='data',
                array=[x_err[i]],
                visible=True,
                thickness=1,
                width=2
            ),
            error_y=dict(
                type='data',
                array=[y_err[i]],
                visible=True,
                thickness=1,
                width=2
            ),
        ), row=1, col=2
    )

    return fig


def st_setup(gene_names: list[str]):
    """Initializes streamlit session.

    Args:
        gene_names (list): The list of gene names
    """

    st.set_page_config(layout="wide")
    st.title('Pseudobulk Gene Activity and Gene Expression')
    st.write("For each gene, we present the change in gene activity (ATAC) and gene expression (RNA) \
              over time. Gene activity represents the sum of ATAC signals within the gene body, while \
              gene expression represents transcript counts (both are log-normalized). Each point on \
              the plot represents the average value across all cells at a given time point.")
    st.sidebar.markdown('# Settings')

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

    # Initialize drop down boxes
    if "selectboxes" not in st.session_state:
        st.session_state.selectboxes = [0]

    # Set up buttons for adding/removing genes
    with st.sidebar:
        add, reset = st.columns([1, 1])
        with add:
            if st.button('Add Gene'):
                st.session_state.selectboxes.append(len(st.session_state.selectboxes))
        with reset:
            if st.button('Reset'):
                st.session_state.selectboxes = [0]

    # Create selectboxes for adding genes
    selected_genes = []
    default = gene_names.index('myf5')
    for i, key in enumerate(st.session_state.selectboxes):
        st.sidebar.selectbox(
            f'Select gene {i+1}',
            gene_names,
            key=key,
            index=default
        )
        selected_genes.append(st.session_state[key])

    return selected_genes


# Main

# Load the csv files as polars dataframes
path = os.path.dirname(os.path.abspath(__file__))+'/data/bulk'
gex_values = pl.read_csv(f'{path}/gex_values.csv')
gex_errors = pl.read_csv(f'{path}/gex_errors.csv')
atac_values = pl.read_csv(f'{path}/atac_values.csv')
atac_errors = pl.read_csv(f'{path}/atac_errors.csv')

# Plot each gene
selected_genes = st_setup(gene_names=gex_values.columns[1:])
for gene in selected_genes:

    # Plot the figures (subplots)
    fig = make_subplots(rows=1, cols=2,
                    subplot_titles=('RNA and ATAC Over Time', 'RNA vs ATAC'),
                    vertical_spacing=0.15,
                    )
    fig = plot_timepoints(gene, fig)
    fig = plot_gene(gene, fig)

    #Update x and y axes for each subplot
    fig.update_xaxes(title_text='Timepoints', row=1, col=1)
    fig.update_yaxes(title_text='Gene Expression/Activity', row=1, col=1)
    fig.update_xaxes(title_text='Gene Activity', row=1, col=2)
    fig.update_yaxes(title_text='Gene Expression', row=1, col=2)

    # Update layout
    fig.update_layout(
            margin=dict(l=10, r=10, t=70, b=0),
            xaxis=dict(title_standoff=5)
        )

    st.markdown(f'### {gene}')
    st.plotly_chart(fig)
