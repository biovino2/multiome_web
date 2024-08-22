# streamlit
streamlit stomping grounds

## Installation
You can install this repository with the following commands:

```
git clone https://github.com/biovino2/streamlit
cd streamlit/
```

Then install the following requirements (conda recommended):
- python==3.12.4
- streamlit==1.37.0
- figeno==1.4.5
- scanpy==1.10.2
- gtfparse==2.5.0
- plotly==5.23.0
- polars==0.20.31

### CCAN plots
Run the following command to run a streamlit instance that shows CCAN and gene track plots with your gene of choice:

```
streamlit run ccan/plot_all.py
```

### Correlation plots
Run the following command to run a streamlit instance that shows correlation scatterplots between ATAC and RNA expression data:

```
streamlit run corr/plot_corr.py
```

### GRN plots
Run the following command to run a streamlit instance that shows CellOracle GRNs as heatmaps for cell types and time points:

```
streamlit run grn/plot_grn.py
```

## Viewing a streamlit instance
You can connect to the instance via the URL streamlit gives you. If you are running from HPC, you must create a SSH tunnel via:

```
ssh -L <port>:localhost:<port> username@hpc_address
```

And then connect via browser. Streamlit defaults to port 8501, but you can set the port with the `-p` argument.
