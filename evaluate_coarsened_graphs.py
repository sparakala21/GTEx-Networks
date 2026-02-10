import os
import networkx as nx
import pandas as pd


def load_graph(tissue_name, threshold=0.5):
    graph_path = f"tissue_networks/{tissue_name.replace(' ', '_')}_network.gexf"
    G = nx.read_gexf(graph_path)
    G_filtered = G.copy()
    edges_to_remove = [(u, v) for u, v, data in G_filtered.edges(data=True) 
                    if data.get('weight', 0) < threshold]
    G_filtered.remove_edges_from(edges_to_remove)
    return G_filtered

def load_embeddings(tissue_name):
    embeddings = pd.read_csv(f"tissue_embeddings/2d-umap/{tissue_name}_embeddings_2d.csv")
    names = embeddings['node'].tolist()
    emb_matrix = embeddings.drop(columns=['node']).to_numpy()
    print(emb_matrix.shape)
    return names, emb_matrix

def evaluate_coarsened_graph_sizes(tissue_name, threshold=0.5):
    G = load_graph(tissue_name=tissue_name, threshold=threshold)
    print(f"Loading graph for {tissue_name}...")
    metrics = {
        "original": (G.number_of_nodes(), G.number_of_edges()),
        'coarsened': []
    }
    for coarsening_factor in [1, 2, 5, 10, 20]:
        print(f"Evaluating coarsened graph with factor {coarsening_factor}...")
        communities = nx.community.louvain_communities(G, weight='weight', seed=42, resolution=coarsening_factor)
        community_graph = nx.Graph()
        
        for i, community in enumerate(communities):
            community_size = len(community)
            community_graph.add_node(i, size=community_size)
            for j in range(i + 1, len(communities)):
                weight = sum(1 for u in community for v in communities[j] if G.has_edge(u, v))
                if weight > 0:
                    community_graph.add_edge(i, j, weight=weight)
        
        num_nodes = community_graph.number_of_nodes()
        num_edges = community_graph.number_of_edges()
        print(f"Coarsened graph with factor {coarsening_factor}: {num_nodes} nodes, {num_edges} edges")
        metrics['coarsened'].append((coarsening_factor, num_nodes, num_edges))
    
    return metrics

results = []
for file in os.listdir('tissue_networks'):
    if file.endswith("_network.gexf"):
        tissue_name = file.replace("_network.gexf", "")
        metrics = evaluate_coarsened_graph_sizes(tissue_name, threshold=0.0)
        
        row = {"Tissue": tissue_name}
        for factor, num_nodes, num_edges in metrics['coarsened']:
            row[f"num_nodes_{factor}"] = num_nodes
            row[f"num_edges_{factor}"] = num_edges
        results.append(row)

df = pd.DataFrame(results)
df.to_csv("coarsened_graph_metrics.csv", index=False)