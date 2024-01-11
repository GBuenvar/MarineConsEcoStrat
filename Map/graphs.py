#%%
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx
#%%

# Load data 
trajectories_df = pd.read_csv("../data/full_data_inds.csv.gz", sep=";", compression="gzip")
int_to_eez = {} # eez dictionary
with open("../data/eez_to_int.csv", "r") as f:  
    f.readline() # read header
    for line in f:
        elem = line.strip().split(";")
        int_to_eez[elem[1]] = elem[0]

# %%
"""
Para crear las red bipartita individuos-EEZs, basta con agrupar por individuo y EEZ, 
conservando la especie también para futuros análisis.
"""
unique_pairs = trajectories_df.groupby(["newid", "Species", "EEZ"]).sum().reset_index().loc[:, ["newid", "Species", "EEZ"]]
unique_pairs["EEZ"] = unique_pairs["EEZ"].apply(lambda x: int_to_eez[x])
# %%

# make a bipartite graph with newid and EEZ as nodes of each partition

B = nx.Graph()
B.add_nodes_from(unique_pairs["newid"].unique(), bipartite=0)
B.add_nodes_from(unique_pairs["EEZ"].unique(), bipartite=1)
B.add_edges_from(zip(unique_pairs["newid"], unique_pairs["EEZ"]))

pos = nx.bipartite_layout(B, nodes=unique_pairs["newid"].unique())
nx.draw_networkx(B, pos=pos)
# %%
