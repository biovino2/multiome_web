"""Master script for plotting in-silico transcription factor knockouts.

Ben Iovino  09/16/24    CZ-Biohub
"""

import streamlit as st


def main():
    """Runs selected page.
    """
    
    pg = st.navigation([
        st.Page('plot_umap.py', title='Plots'),
        st.Page('plot_ps.py', title='Scores'),
        ])
    pg.run()


if __name__ == '__main__':
    main()