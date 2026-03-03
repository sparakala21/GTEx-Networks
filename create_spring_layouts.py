import pandas as pd
from sklearn.decomposition import PCA
from umap import UMAP
import os
import sys
import networkx as nx

def create_spring_layout(G):
    pos = nx.spring_layout(G, seed=42)
    return pos

def export_spring_layout(tissue_name, pos):
    df = pd.DataFrame(pos).T.reset_index()
    df.columns = ['node', 'x', 'y']
    df.to_csv(f"spring_layout/{tissue_name}.csv", index=False)
    
if __name__ == "__main__":
    folder = "tissue_networks"
    method = sys.argv[1] if len(sys.argv) > 1 else "umap"
    file = "Liver_network.gexf"
    print(f"Processing {file}...")
    G = nx.read_gexf(f"{folder}/{file}")
    print(f"Graph loaded with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    tissue_name = file.replace("_network.gexf", "").replace("_", " ")
    layout = create_spring_layout(G)
    print(f"Spring layout created for {tissue_name}.")
    export_spring_layout(tissue_name, layout)

            