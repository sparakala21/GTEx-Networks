import pandas as pd
from sklearn.decomposition import PCA
from umap import UMAP
import os
import sys
def create_embeddings(tissue_name, method='pca'):
    embeddings = pd.read_csv(f"tissue_embeddings/{tissue_name}_embeddings.csv")
    names = embeddings['node'].tolist()
    emb_matrix = embeddings.drop(columns=['node']).to_numpy()
    print(emb_matrix.shape)

    pca = PCA(n_components=50, random_state=42)

    embeddings_50d = pca.fit_transform(emb_matrix)
    if method == 'pca':
        pca = PCA(n_components=2, random_state=42)

        embeddings_2d = u.fit_transform(embeddings_50d)
        return names, embeddings_2d
    elif method == "umap":


        u = UMAP(n_components=2, random_state=42)

        embeddings_2d = u.fit_transform(embeddings_50d)
        return names, embeddings_2d
    else:
        print("Method not found")
if __name__ == "__main__":
    #create 2d embeddings for all 64d embedding files in tissue_embeddings folder
    folder = "tissue_embeddings"
    method = sys.argv[1] if len(sys.argv) > 1 else "umap"
    for file in os.listdir(folder):
        if file.endswith("_embeddings.csv"):
            tissue_name = file.replace("_embeddings.csv", "")
            print(f"creating 2D embeddings for {tissue_name}...")
            names, embeddings_2d = create_embeddings(tissue_name)
            print(f"created embeddings for {tissue_name}.")
            #save to csv
            df = pd.DataFrame(embeddings_2d, columns=['x', 'y'])
            df.insert(0, 'node', names)
            df.to_csv(f"{folder}/2d-umap/{tissue_name}_embeddings_2d.csv", index=False)
            print(f"saved 2D embeddings for {tissue_name} to {folder}/{tissue_name}_embeddings_2d.csv")