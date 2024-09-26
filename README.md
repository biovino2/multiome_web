# multiome_web
This repository contains several interactive figures that can be run on any machine and accessed through a browser.

## Installation
You can install this repository with the following commands:

```
git clone https://github.com/biovino2/multiome_web
cd multiome_web/
conda env create -f env.yml
conda activate multiome_web
```

### Alternate installation

If creating a conda environment does not work for you, try installing each package one by one:
- python==3.12.4
- streamlit==1.37.0
- scanpy==1.10.2
- plotly==5.24.0
- polars==1.6.0

### Gene dynamics
Run the following command to run a streamlit instance that shows a track plot and associated peaks for every gene, as well as correlation scatter plots between gene activity and expression:

```
streamlit run dyn/plot_dyn.py
```

### GRN plots
Run the following command to run a streamlit instance that shows CellOracle GRNs as heatmaps for cell types and time points:

```
streamlit run grn/plot_grn.py
```

### TF Knockouts
Run the following command to run a streamlit instance that depicts the transition probabilities of metacells given a transcription factor knockout:

```
streamlit run tfko/plot_tfko.py
```

## Viewing a streamlit instance
You can connect to the instance via the URL streamlit gives you. If you are running from HPC, you must create a SSH tunnel via:

```
ssh -L <port>:localhost:<port> username@hpc_address
```

And then connect via browser. Streamlit defaults to port 8501, but you can set the port with the `-p` argument.
