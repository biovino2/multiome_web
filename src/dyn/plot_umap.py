"""Plot gene UMAP and show cluster enrichments.

Ben Iovino  12/04/24    CZ-Biohub
"""

import os
import pandas as pd
import plotly.graph_objects as go
import streamlit as st


def umap_color_dict() -> dict:
    """Returns dictionary mapping cluster to color.

    Returns:
        dict: Key is cluster ID, value is color.
    """

    color_dict = {
        0: '#1f77b4',
        1: '#ff7f0e',
        2: '#279e68',
        3: '#d62728',
        4: '#aa40fc',
        5: '#8c564b',
        6: '#e377c2',
        7: '#b5bd61',
        8: '#17becf',
        9: '#aec7e8',
        10: '#ffbb78',
        11: '#98df8a',
        12: '#ff9896',
        13: '#c5b0d5'
    }

    return color_dict


def st_setup():
    """Initializes streamlit session.
    """

    st.set_page_config(layout="wide")
    st.markdown("# Gene UMAP")
    st.write('')
    
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

    # Initialize drop down box
    cluster_ids = umap_color_dict().keys()
    option = st.sidebar.selectbox(
        'Select a cluster',
        cluster_ids,
    )

    return option


def plot_umap(df_umap: pd.DataFrame):
    """Plots gene UMAP.
    
    Args:
        df_umap (pd.DataFrame): The dataframe of UMAP coordinates and cluster assignments.
    """

    color_dict = umap_color_dict()
    fig = go.Figure()

    # Add all cells of same type to plot
    fig.add_trace(go.Scatter(
        x=df_umap['umap1'],
        y=df_umap['umap2'],
        mode='markers',
        marker=dict(
            size=5,
            color=df_umap['color'],
        ),
        hoverinfo='text',
        text=[f"Gene: {gene_name}<br>Cluster: {cluster_id}"
                for gene_name, cluster_id in zip(df_umap.index, df_umap['Leiden'])],
        showlegend=False,
    ))

    for i, cluster in enumerate(df_umap['Leiden'].unique()):
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(color=color_dict[i], size=10),
            showlegend=True,
            legendgroup='group1',
            name=i,
    ))

    fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False)
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False)
    fig.update_layout(margin=dict(l=10, r=10, t=20, b=0), plot_bgcolor='rgba(0,0,0,0)')

    return fig


def plot_enrichments(df_subset: pd.DataFrame, cluster: int):
    """
    """

    df_subset = df_subset.sort_values('Combined Score').tail(10)

    # Define color map matching the UMAP
    cluster_colors = {
        0: '#1f77b4',  # blue
        1: '#ff7f0e',  # orange
        2: '#2ca02c',  # green
        3: '#d62728',  # red
        4: '#9467bd',  # purple
        5: '#8c564b',  # brown
        6: '#e377c2',  # pink
        7: '#bcbd22',  # olive
        8: '#17becf',  # cyan
        9: '#aec7e8',  # light blue
        10: '#ffbb78', # light orange
        11: '#98df8a', # light green
        12: '#ff9896', # light red
        13: '#c5b0d5'  # light purple
    }
    nonsig_color = '#F5F5DC'  # light beige

    colors = [nonsig_color if p > 0.05 else cluster_colors[cluster] for p in df_subset['P-value']]

    # Create the horizontal bar plot
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=df_subset['Combined Score'],
            y=df_subset['Term'],
            orientation='h',
            marker=dict(color=colors, opacity=0.6),
            text=df_subset['Genes'],
            insidetextanchor='start',
        )
    )

    # Customize layout
    fig.update_layout(
        title=f"Cluster {cluster} Top Terms (FishEnrichr)",
        xaxis_title="Combined Score",
        yaxis_title="Term",
        plot_bgcolor='rgba(0,0,0,0)'
    )

    return fig


### Main

# Load data and set up streamlit
path = os.path.dirname(os.path.abspath(__file__))+'/data'
with open(f'{path}/umap.csv', 'r') as f:
    df_umap = pd.read_csv(f, index_col=0)
option = st_setup()

# Plot UMAP
fig = plot_umap(df_umap)
st.plotly_chart(fig, use_container_width=True)

# Display cluster enrichments
enrich_dir = "src/dyn/data/enrichr"
df = pd.read_csv(f"{enrich_dir}/cluster_{option}.txt", sep="\t")
df_subset = df[["Term", "Genes", "P-value", "Combined Score"]]
fig = plot_enrichments(df_subset, option)
st.plotly_chart(fig, use_container_width=True)
