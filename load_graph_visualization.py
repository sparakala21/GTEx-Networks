import plotly.graph_objects as go
import networkx as nx
import sys
import pandas as pd
from sklearn.decomposition import PCA
from umap import UMAP

tissue_name = sys.argv[1]

embeddings = pd.read_csv(f"{tissue_name}_embeddings.csv")
names = embeddings['node'].tolist()
emb_matrix = embeddings.drop(columns=['node']).to_numpy()
print(emb_matrix.shape)

pca = PCA(n_components=50, random_state=42)

embeddings_50d = pca.fit_transform(emb_matrix)

u = UMAP(n_components=2, random_state=42)

embeddings_2d = u.fit_transform(embeddings_50d)

tissue_name = "Adipose_Subcutaneous"
graph_path = f"tissue_networks/{tissue_name.replace(' ', '_')}_network.gexf"
G = nx.read_gexf(graph_path)
threshold = 0.5
G_filtered = G.copy()
edges_to_remove = [(u, v) for u, v, data in G_filtered.edges(data=True) 
                   if data.get('weight', 0) < threshold]
G_filtered.remove_edges_from(edges_to_remove)
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

# Create node trace
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

# Add node labels for hover
node_trace.text = list(G.nodes())

# Create figure
fig = go.Figure(data=[edge_trace, node_trace],
                layout=go.Layout(
                    title='Adipose Subcutaneous Tissue Network(Node2Vec + UMAP layout)',
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=0,l=0,r=0,t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )

fig.write_html("network_graph.html")