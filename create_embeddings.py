import pandas as pd
import numpy as np
import sqlite3
from node2vec import Node2Vec
import networkx as nx
import pickle
import json

def create_embeddings(G, dimensions=128, walk_length=30, num_walks=200, workers=1):
    """Generate node2vec embeddings from the graph"""
    node2vec = Node2Vec(G, dimensions=dimensions, walk_length=walk_length, 
                       num_walks=num_walks, workers=workers)
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

    # Extract embeddings as dictionary
    embeddings = {node: model.wv[node].tolist() for node in G.nodes()}

    # Create DataFrame
    embedding_df = pd.DataFrame([
        {'gene': gene, 'embedding': json.dumps(emb)}
        for gene, emb in embeddings.items()
    ])

    # Save to database
    conn = sqlite3.connect('GTEx.db')
    embedding_df.to_sql('gene_embeddings', conn, if_exists='replace', index=False)
    conn.close()

    # Later, to retrieve:
    conn = sqlite3.connect('GTEx.db')
    df = pd.read_sql("SELECT * FROM gene_embeddings", conn)
    df['embedding'] = df['embedding'].apply(json.loads)
    conn.close()
