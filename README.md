# streamlit
streamlit stomping grounds

## Installation
Not yet sure if you can run a streamlit instance on the HPC due to authorization problems. Best to run on local machine for now. All necessary data and scripts are contained in the repository.

You can install this repository with the following commands:

```
git clone https://github.com/biovino2/streamlit
cd streamlit/
```

Then install the following requirements (conda recommended):
- python==3.12.4
- streamlit==1.37.0
- figeno==1.4.5
- gtfparse==2.5.0

## CCAN plots with streamlit
Run the following command to run a streamlit instance that shows CCAN and gene track plots with your gene of choice:

```
streamlit run ccan/plot_all.py
```
