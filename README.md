# streamlit
streamlit stomping grounds

## Installation
Not yet sure if you can run a streamlit instance on the HPC due to authorization problems. Best to run on local machine for now. All necessary data and scripts are contained in the repository. It also appears that the environment that works on mac does not work on HPC, hence mac_env.yml.

You can install this repository with the following commands:

```
git clone https://github.com/biovino2/streamlit
cd streamlit/
```

Then install the following requirements:
- python==3.12.4
- streamlit==1.37.0
- figeno==1.4.4
- gtfparse==2.5.0 

## CCAN Plots
Run the following command to run a streamlit instance that shows CCAN plots for genes of your choice:

```
streamlit run ccan.py
```
