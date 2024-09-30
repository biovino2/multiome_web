## data explanation
These files contain the gene regulatory network (GRN) heatmaps as csv files. The heatmaps are hierachically clustered with respect to either cell type or time point. This affects how the GRN is displayed because in order to compare GRN's like we do in the plots, we must only show genes and TF's shared between all GRNs being compared.

By controlling for cell type, there are more genes and TF's shared between the same cell type at different time points, displaying a larger GRN. By controlling for time point, there are less genes and TF's shared between cell types at the same time point, displaying a smaller GRN.

### example for cell type clustering
To see GRN's clustered for spinal cord cells, go to `ct/spinal_cord`, where each csv gives the strength of interaction between genes and TF's for spinal cord cells that are shared among all time points (10, 12, 14, 16, 19, and 24 hours post fertilization).

### example for time point clustering
To see GRN's clustered for 10 hours post fertilization, go to `tp/10hpf`, where each csv gives the strength of interaction between genes and TF's for 10 hours post fertilization that are shared among all cell types (neural posterior, NMPs, PSM, somites, spinal cord, and tail bud).