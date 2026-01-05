import pandas as pd
import numpy as np
import sqlite3
import networkx as nx
import json
from sklearn.metrics.pairwise import cosine_similarity
from collections import defaultdict
import matplotlib.pyplot as plt

def load_embeddings(db_path='GTEx.db'):
    """Load embeddings from database"""
    conn = sqlite3.connect(db_path)
    df = pd.read_sql("SELECT * FROM gene_embeddings", conn)
    conn.close()
    
    # Parse JSON strings to lists
    df['embedding'] = df['embedding'].apply(json.loads)
    
    # Convert to numpy array
    genes = df['gene'].values
    embeddings = np.array(df['embedding'].tolist())
    
    # Create gene to index mapping
    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
    
    return genes, embeddings, gene_to_idx

def load_graph(db_path='GTEx.db', weight_threshold=0.7):
    """Load the original co-expression graph"""
    conn = sqlite3.connect(db_path)
    df = pd.read_sql(f"""
        SELECT * FROM GTEx_network_age_adjusted 
        WHERE weight > {weight_threshold}
    """, conn)
    conn.close()
    
    df['weight'] = pd.to_numeric(df['weight'], errors='coerce')
    G = nx.from_pandas_edgelist(df, 'row', 'col', edge_attr='weight')
    
    # Get largest connected component
    largest_cc = max(nx.connected_components(G), key=len)
    G = G.subgraph(largest_cc).copy()
    
    return G

def get_knn_in_embedding_space(embeddings, gene_to_idx, target_gene, k=10):
    """Find k nearest neighbors in embedding space using cosine similarity"""
    if target_gene not in gene_to_idx:
        return []
    
    target_idx = gene_to_idx[target_gene]
    target_embedding = embeddings[target_idx].reshape(1, -1)
    
    # Calculate cosine similarity with all other genes
    similarities = cosine_similarity(target_embedding, embeddings)[0]
    
    # Get top k+1 (including the gene itself)
    top_indices = np.argsort(similarities)[::-1][1:k+1]  # Skip first (itself)
    
    genes = list(gene_to_idx.keys())
    return [genes[idx] for idx in top_indices]

def get_nth_order_neighbors(G, target_gene, n=2):
    """Find all neighbors within n hops in the graph"""
    if target_gene not in G:
        return set()
    
    neighbors = set()
    
    # BFS to find all nodes within n hops
    visited = {target_gene: 0}
    queue = [(target_gene, 0)]
    
    while queue:
        node, dist = queue.pop(0)
        
        if dist < n:
            for neighbor in G.neighbors(node):
                if neighbor not in visited:
                    visited[neighbor] = dist + 1
                    neighbors.add(neighbor)
                    queue.append((neighbor, dist + 1))
    
    return neighbors

def evaluate_embeddings(genes, embeddings, gene_to_idx, G, k_values=[5, 10, 20, 50], n_values=[1, 2, 3]):
    """
    Evaluate embeddings by comparing KNN in embedding space to graph neighbors
    
    Returns:
    - Precision: % of embedding KNN that are actual graph neighbors
    - Recall: % of actual graph neighbors that are in embedding KNN
    """
    results = []
    
    # Sample genes that exist in both embeddings and graph
    common_genes = [g for g in genes if g in G]
    
    # Sample for faster evaluation (adjust sample_size as needed)
    sample_size = min(500, len(common_genes))
    sampled_genes = np.random.choice(common_genes, sample_size, replace=False)
    
    print(f"Evaluating on {sample_size} genes...")
    
    for k in k_values:
        for n in n_values:
            precisions = []
            recalls = []
            overlaps = []
            
            for gene in sampled_genes:
                # Get KNN in embedding space
                knn_embedding = set(get_knn_in_embedding_space(embeddings, gene_to_idx, gene, k))
                
                # Get nth-order neighbors in graph
                graph_neighbors = get_nth_order_neighbors(G, gene, n)
                
                if len(graph_neighbors) == 0:
                    continue
                
                # Calculate overlap
                overlap = knn_embedding & graph_neighbors
                
                # Precision: what % of embedding KNN are graph neighbors?
                precision = len(overlap) / len(knn_embedding) if len(knn_embedding) > 0 else 0
                
                # Recall: what % of graph neighbors are in embedding KNN?
                recall = len(overlap) / len(graph_neighbors) if len(graph_neighbors) > 0 else 0
                
                precisions.append(precision)
                recalls.append(recall)
                overlaps.append(len(overlap))
            
            avg_precision = np.mean(precisions)
            avg_recall = np.mean(recalls)
            avg_overlap = np.mean(overlaps)
            f1_score = 2 * (avg_precision * avg_recall) / (avg_precision + avg_recall) if (avg_precision + avg_recall) > 0 else 0
            
            results.append({
                'k': k,
                'n': n,
                'precision': avg_precision,
                'recall': avg_recall,
                'f1_score': f1_score,
                'avg_overlap': avg_overlap
            })
            
            print(f"k={k}, n={n}: Precision={avg_precision:.3f}, Recall={avg_recall:.3f}, F1={f1_score:.3f}, Avg Overlap={avg_overlap:.1f}")
    
    return pd.DataFrame(results)

def visualize_results(results_df):
    """Create visualizations of evaluation results"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Precision heatmap
    pivot_precision = results_df.pivot(index='n', columns='k', values='precision')
    im1 = axes[0, 0].imshow(pivot_precision, cmap='YlGnBu', aspect='auto')
    axes[0, 0].set_title('Precision: % of Embedding KNN that are Graph Neighbors')
    axes[0, 0].set_xlabel('k (embedding neighbors)')
    axes[0, 0].set_ylabel('n (graph hops)')
    axes[0, 0].set_xticks(range(len(pivot_precision.columns)))
    axes[0, 0].set_xticklabels(pivot_precision.columns)
    axes[0, 0].set_yticks(range(len(pivot_precision.index)))
    axes[0, 0].set_yticklabels(pivot_precision.index)
    plt.colorbar(im1, ax=axes[0, 0])
    
    # Add values to heatmap
    for i in range(len(pivot_precision.index)):
        for j in range(len(pivot_precision.columns)):
            axes[0, 0].text(j, i, f'{pivot_precision.iloc[i, j]:.2f}',
                          ha="center", va="center", color="black")
    
    # Plot 2: Recall heatmap
    pivot_recall = results_df.pivot(index='n', columns='k', values='recall')
    im2 = axes[0, 1].imshow(pivot_recall, cmap='YlOrRd', aspect='auto')
    axes[0, 1].set_title('Recall: % of Graph Neighbors in Embedding KNN')
    axes[0, 1].set_xlabel('k (embedding neighbors)')
    axes[0, 1].set_ylabel('n (graph hops)')
    axes[0, 1].set_xticks(range(len(pivot_recall.columns)))
    axes[0, 1].set_xticklabels(pivot_recall.columns)
    axes[0, 1].set_yticks(range(len(pivot_recall.index)))
    axes[0, 1].set_yticklabels(pivot_recall.index)
    plt.colorbar(im2, ax=axes[0, 1])
    
    for i in range(len(pivot_recall.index)):
        for j in range(len(pivot_recall.columns)):
            axes[0, 1].text(j, i, f'{pivot_recall.iloc[i, j]:.2f}',
                          ha="center", va="center", color="black")
    
    # Plot 3: F1 Score heatmap
    pivot_f1 = results_df.pivot(index='n', columns='k', values='f1_score')
    im3 = axes[1, 0].imshow(pivot_f1, cmap='RdYlGn', aspect='auto')
    axes[1, 0].set_title('F1 Score (Harmonic Mean of Precision & Recall)')
    axes[1, 0].set_xlabel('k (embedding neighbors)')
    axes[1, 0].set_ylabel('n (graph hops)')
    axes[1, 0].set_xticks(range(len(pivot_f1.columns)))
    axes[1, 0].set_xticklabels(pivot_f1.columns)
    axes[1, 0].set_yticks(range(len(pivot_f1.index)))
    axes[1, 0].set_yticklabels(pivot_f1.index)
    plt.colorbar(im3, ax=axes[1, 0])
    
    for i in range(len(pivot_f1.index)):
        for j in range(len(pivot_f1.columns)):
            axes[1, 0].text(j, i, f'{pivot_f1.iloc[i, j]:.2f}',
                          ha="center", va="center", color="black")
    
    # Plot 4: Line plots
    for n in results_df['n'].unique():
        subset = results_df[results_df['n'] == n]
        axes[1, 1].plot(subset['k'], subset['f1_score'], marker='o', label=f'n={n} hops')
    
    axes[1, 1].set_title('F1 Score vs k for Different Graph Distances')
    axes[1, 1].set_xlabel('k (embedding neighbors)')
    axes[1, 1].set_ylabel('F1 Score')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('embedding_evaluation.png', dpi=300, bbox_inches='tight')
    print("\nVisualization saved as 'embedding_evaluation.png'")
    plt.savefig("weighted_embedding_evaluation.png")

def detailed_example(genes, embeddings, gene_to_idx, G, example_gene):
    """Show detailed results for a specific gene"""
    if example_gene not in gene_to_idx or example_gene not in G:
        print(f"Gene {example_gene} not found in embeddings or graph")
        return
    
    print(f"\n{'='*70}")
    print(f"Detailed Analysis for Gene: {example_gene}")
    print(f"{'='*70}\n")
    
    # Get embedding neighbors
    knn_10 = get_knn_in_embedding_space(embeddings, gene_to_idx, example_gene, k=10)
    print(f"Top 10 Nearest Neighbors in Embedding Space:")
    for i, gene in enumerate(knn_10, 1):
        print(f"  {i}. {gene}")
    
    # Get graph neighbors at different distances
    for n in [1, 2, 3]:
        graph_neighbors = get_nth_order_neighbors(G, example_gene, n)
        overlap = set(knn_10) & graph_neighbors
        
        print(f"\n{n}-hop Graph Neighbors: {len(graph_neighbors)} genes")
        print(f"Overlap with Top 10 Embedding KNN: {len(overlap)} genes")
        if overlap:
            print(f"Overlapping genes: {', '.join(list(overlap)[:5])}" + 
                  (f" ... and {len(overlap)-5} more" if len(overlap) > 5 else ""))

if __name__ == "__main__":
    print("Loading embeddings...")
    genes, embeddings, gene_to_idx = load_embeddings('GTEx.db')
    print(f"Loaded {len(genes)} gene embeddings with dimension {embeddings.shape[1]}")
    
    print("\nLoading graph...")
    G = load_graph('GTEx.db', weight_threshold=0)
    print(f"Loaded graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    
    print("\n" + "="*70)
    print("EMBEDDING EVALUATION")
    print("="*70)
    
    # Run evaluation
    results_df = evaluate_embeddings(
        genes, embeddings, gene_to_idx, G,
        k_values=[5, 10, 20, 50, 100],
        n_values=[1, 2, 3, 4]
    )
    
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print("\nBest configurations:")
    print("\nBy Precision (embedding KNN are graph neighbors):")
    print(results_df.nlargest(3, 'precision')[['k', 'n', 'precision', 'recall', 'f1_score']])
    
    print("\nBy Recall (graph neighbors are in embedding KNN):")
    print(results_df.nlargest(3, 'recall')[['k', 'n', 'precision', 'recall', 'f1_score']])
    
    print("\nBy F1 Score (balanced):")
    print(results_df.nlargest(3, 'f1_score')[['k', 'n', 'precision', 'recall', 'f1_score']])
    
    # Save results
    results_df.to_csv('embedding_evaluation_results.csv', index=False)
    print("\nResults saved to 'embedding_evaluation_results.csv'")
    
    # Visualize
    print("\nGenerating visualizations...")
    visualize_results(results_df)
    
    # Show detailed example for a specific gene
    common_genes = [g for g in genes if g in G]
    if common_genes:
        example_gene = common_genes[0]  # Or choose a specific gene like 'RPL31'
        detailed_example(genes, embeddings, gene_to_idx, G, example_gene)