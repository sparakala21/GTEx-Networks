import networkx as nx
import numpy as np

tissue = "Liver"  # or sys.argv[1]
G = nx.read_gexf(f"tissue_networks/{tissue.replace(' ', '_')}_network.gexf")

print(f"Total edges: {G.number_of_edges()}")

weights = [data.get('weight', 1.0) for u, v, data in G.edges(data=True)]
weights = [float(w) if isinstance(w, str) else w for w in weights]

print(f"Weight stats:")
print(f"  Min: {np.min(weights)}")
print(f"  Max: {np.max(weights)}")
print(f"  Mean: {np.mean(weights)}")
print(f"  NaN count: {np.sum(np.isnan(weights))}")
print(f"  Inf count: {np.sum(np.isinf(weights))}")
print(f"  Negative count: {np.sum(np.array(weights) < 0)}")
print(f"  Zero count: {np.sum(np.array(weights) == 0)}")