import networkx as nx
from node2vec import Node2Vec
import os
import sys

def embed_graph(tissue_name, folder="tissue_networks", dimensions=64, walk_length=10, num_walks=50, p=1, q=1):
    graph_path = f"{folder}/{tissue_name.replace(' ', '_')}_network.gexf"
    G = nx.read_gexf(graph_path)
    node2vec = Node2Vec(G, dimensions=dimensions, walk_length=walk_length, num_walks=num_walks, p=p, q=q, weight_key="weight", workers=os.cpu_count()//2)
    model = node2vec.fit(window=10, min_count=1, batch_words=4)
    return model
if __name__ == "__main__":
    tissue = sys.argv[1] if len(sys.argv) > 1 else "Liver"
    model = embed_graph(tissue)
    # model.wv.save_word2vec_format(f"{tissue.replace(' ', '_')}_embeddings.txt")
    # print(f"Embeddings for {tissue} saved to {tissue.replace(' ', '_')}_embeddings.txt")
    # save embeddings for each node to a file
    with open(f"{tissue.replace(' ', '_')}_embeddings.csv", "w") as f:
        f.write("node,embedding\n")
        for node in model.wv.index_to_key:
            embedding = ",".join([str(x) for x in model.wv[node]])
            f.write(f"{node},{embedding}\n")
    print(f"Embeddings for {tissue} saved to {tissue.replace(' ', '_')}_embeddings.csv")
    