import plotly.graph_objects as go
import networkx as nx
import sys
import pandas as pd

def load_graph(tissue_name, threshold=0.5):
    graph_path = f"tissue_networks/{tissue_name.replace(' ', '_')}_network.gexf"
    G = nx.read_gexf(graph_path)
    threshold = 0.5
    G_filtered = G.copy()
    edges_to_remove = [(u, v) for u, v, data in G_filtered.edges(data=True) 
                    if data.get('weight', 0) < threshold]
    G_filtered.remove_edges_from(edges_to_remove)
    return G_filtered

if __name__ == "__main__":
    tissue_name = sys.argv[1] if len(sys.argv) > 1 else "Subcutaneous_Adipose"
    threshold = sys.argv[2] if len(sys.argv) > 2 else "0.7"
    export = sys.argv[3]  if len(sys.argv) > 3 else "false"

    G = load_graph(tissue_name=tissue_name, threshold=0.7)

    #create a list of disjoint vertex sets
    communities = nx.community.louvain_communities(G, weight='weight', seed=42)


