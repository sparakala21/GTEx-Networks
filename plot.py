#graphing an example tissue network
import matplotlib.pyplot as plt
import networkx as nx
from math import sqrt
import plotly.graph_objects as go
import plotly.express as px
import sys

k = 4
folder = "tissue_networks"
example_tissue = sys.argv[1] if len(sys.argv) > 1 else "Liver"
G = nx.read_gexf(f"{folder}/{example_tissue.replace(' ', '_')}_network.gexf")
core = nx.k_core(G, k=k)
print(f"{example_tissue} k-core (k={k}) has {core.number_of_nodes():,} nodes and {core.number_of_edges():,} edges.")
# graph using plotly
pos = nx.spring_layout(core, seed=42)

edge_x = []
edge_y = []
for edge in core.edges():
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

node_x = []
node_y = []
for node in core.nodes():
    x, y = pos[node]
    node_x.append(x)
    node_y.append(y)

node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode='markers',
    hoverinfo='text',
    text=list(core.nodes()),
    marker=dict(
        showscale=True,
        colorscale='YlGnBu',
        size=10,
        line_width=2))

fig = go.Figure(data=[edge_trace, node_trace],
                layout=go.Layout(
                    title=f"{example_tissue} Network (k-core={k})",
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                ))
fig.show()