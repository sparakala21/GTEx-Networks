import networkx as nx
from node2vec import Node2Vec
import os
import sys
import logging

def embed_graph(tissue_name, folder="tissue_networks", dimensions=64, walk_length=20, num_walks=100, p=1, q=1):
    logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

    graph_path = f"{folder}/{tissue_name.replace(' ', '_')}_network.gexf"
    G = nx.read_gexf(graph_path)
    node2vec = Node2Vec(G, dimensions=dimensions, walk_length=walk_length, num_walks=num_walks, p=p, q=q, weight_key="weight", workers=os.cpu_count()-2)
    print("walks complete. ")
    model = node2vec.fit(window=10, min_count=1, batch_words=4, workers=os.cpu_count()-2)
    return model
if __name__ == "__main__":
    tissue = sys.argv[1] if len(sys.argv) > 1 else "Liver"
    folder = "tissue_embeddings"
    model = embed_graph(tissue)
    # model.wv.save_word2vec_format(f"{tissue.replace(' ', '_')}_embeddings.txt")
    # print(f"Embeddings for {tissue} saved to {tissue.replace(' ', '_')}_embeddings.txt")
    # save embeddings for each node to a file
    with open(f"{folder}/{tissue.replace(' ', '_')}_embeddings.csv", "w") as f:
        f.write("node,d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23, d24, d25, d26, d27, d28, d29, d30, d31, d32, d33, d34, d35, d36, d37, d38, d39, d40, d41, d42, d43, d44, d45, d46, d47, d48, d49, d50, d51, d52, d53, d54, d55, d56, d57, d58, d59, d60, d61, d62, d63\n")
        for node in model.wv.index_to_key:
            embedding = ",".join([str(x) for x in model.wv[node]])
            f.write(f"{node},{embedding}\n")
    print(f"Embeddings for {tissue} saved to {tissue.replace(' ', '_')}_embeddings.csv")
    