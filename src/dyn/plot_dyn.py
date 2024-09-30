"""Master script for plotting in-silico transcription factor knockouts.

Ben Iovino  09/24/24    CZ-Biohub
"""

import streamlit as st


def main():
    """Runs selected page.
    """
    
    pg = st.navigation([
        st.Page('plot_ccan.py', title='Peaks'),
        st.Page('plot_corr.py', title='Correlation'),
        ])
    pg.run()


if __name__ == '__main__':
    main()
