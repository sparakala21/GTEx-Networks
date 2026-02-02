import plotly.graph_objects as go
import networkx as nx
import sys
import pandas as pd
from sklearn.decomposition import PCA
from umap import UMAP

def create_embeddings(tissue_name):
    embeddings = pd.read_csv(f"tissue_embeddings/{tissue_name}_embeddings.csv")
    names = embeddings['node'].tolist()
    emb_matrix = embeddings.drop(columns=['node']).to_numpy()
    print(emb_matrix.shape)

    pca = PCA(n_components=50, random_state=42)

    embeddings_50d = pca.fit_transform(emb_matrix)

    u = UMAP(n_components=2, random_state=42)

    embeddings_2d = u.fit_transform(embeddings_50d)
    return names, embeddings_2d

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
    
    tissue_name = sys.argv[1]
    print(f"creating 2D embeddings for {tissue_name}...")
    names, embeddings_2d = create_embeddings(tissue_name)
    print(f"created embeddings for {tissue_name}.")
    print(f"loading graph for {tissue_name}...")
    G = load_graph(tissue_name, threshold=0.5)
    print(f"loaded graph for {tissue_name} with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    pos = dict(zip(names, embeddings_2d))

    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = [pos[node][0] for node in G.nodes()]
    node_y = [pos[node][1] for node in G.nodes()]

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            size=5,
            color='blue',
            line_width=0.5))

    node_trace.text = list(G.nodes())

    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        title=f'{tissue_name.replace('_', " ")} Tissue Network(Node2Vec + UMAP layout)',
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=0,l=0,r=0,t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )

    fig.write_html(f"tissue_network_visualizations/{tissue_name}_graph.html")