import pandas as pd
import numpy as np
import sqlite3
from node2vec import Node2Vec
import networkx as nx
import pickle
import json
#this code burned through several ec2s until i realized ive been underprovisioning ram
#calculation for the ram usage:
# num_walks * walk_length * num_nodes * 40-80 bytes
# 40 * 15 * 25000 * 40-80 = 600,000,000 to 1.2 billion bytes which is roughly 4.4 to 9 GiB
def create_embeddings(G, dimensions=128, walk_length=15, num_walks=40, workers=1):
    node2vec = Node2Vec(
        G,
        dimensions=dimensions,
        walk_length=walk_length,
        num_walks=num_walks,
        p=0.5, q=2,
        workers=8,
        weight_key='weight'
    )
    model = node2vec.fit(window=10, min_count=1, batch_words=4)
    return model
if __name__ == "__main__":
    conn = sqlite3.connect('GTEx.db')
    df = pd.read_sql("SELECT * FROM GTEx_network_age_adjusted", conn)
    conn.close()

    df['weight'] = pd.to_numeric(df['weight'], errors='coerce')
    G = nx.from_pandas_edgelist(df, 'row', 'col', ['weight'])
    largest_cc = max(nx.connected_components(G), key=len)
    G = G.subgraph(largest_cc).copy()

    model = create_embeddings(G)

    embeddings = {node: model.wv[node].tolist() for node in G.nodes()}

    embedding_df = pd.DataFrame([
        {'gene': gene, 'embedding': json.dumps(emb)}
        for gene, emb in embeddings.items()
    ])

    conn = sqlite3.connect('GTEx.db')
    embedding_df.to_sql('gene_embeddings', conn, if_exists='replace', index=False)
    conn.close()

    conn = sqlite3.connect('GTEx.db')
    df = pd.read_sql("SELECT * FROM gene_embeddings", conn)
    df['embedding'] = df['embedding'].apply(json.loads)
    conn.close()
