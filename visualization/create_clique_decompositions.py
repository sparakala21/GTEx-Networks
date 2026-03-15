import networkx as nx
import numpy as np
import pandas as pd
import psycopg2
from psycopg2.extras import execute_values
import json


def load_graph(tissue_name):
    graph_path = f"../tissue_networks/{tissue_name.replace(' ', '_')}_network.gexf"
    G = nx.read_gexf(graph_path)
    for u, v, data in G.edges(data=True):
        if 'weight' in data:
            data['weight'] = float(data['weight'])
    return G


def load_layout(tissue_name):
    embeddings = pd.read_csv(f"../spring_layout/{tissue_name}_embeddings_2d.csv")
    names = embeddings['node'].tolist()
    emb_matrix = embeddings.drop(columns=['node']).to_numpy(dtype=float)
    print(emb_matrix.shape)
    return names, emb_matrix


def build_layout_dict(G, names, emb_matrix):
    """Map node IDs that exist in G to (x, y) from the embedding matrix."""
    name_to_index = {name: idx for idx, name in enumerate(names)}
    layout = {}
    for node_id in G.nodes():
        if node_id in name_to_index:
            idx = name_to_index[node_id]
            layout[node_id] = (emb_matrix[idx, 0], emb_matrix[idx, 1])
        else:
            print(f"Warning: node {node_id} not found in layout, defaulting to (0, 0)")
            layout[node_id] = (0.0, 0.0)
    return layout


def export_nodes(G, layout, cur, conn):
    rows = []
    for node_id, (x, y) in layout.items():
        label = G.nodes[node_id].get("label", str(node_id))
        rows.append((str(node_id), float(x), float(y), label, None))

    execute_values(cur, """
        INSERT INTO nodes (id, x, y, label, parent_id)
        VALUES %s
        ON CONFLICT (id) DO UPDATE SET x=EXCLUDED.x, y=EXCLUDED.y
    """, rows)
    conn.commit()
    print(f"Inserted {len(rows)} nodes")


def export_edges(G, level, cur, conn):
    rows = []
    for source, target, data in G.edges(data=True):
        rows.append((
            str(source),
            str(target),
            level,
            False,
            float(data.get("weight", 1.0))
        ))

    execute_values(cur, """
        INSERT INTO edges (source_id, target_id, level, is_boundary, weight)
        VALUES %s
    """, rows)
    conn.commit()
    print(f"Inserted {len(rows)} edges at level {level}")


def find_and_collapse_cliques(G, level, cur, conn, min_size=3, max_size=4):
    clique_rows = []
    edge_rows = []
    node_to_clique = {}

    all_cliques = [
        c for c in nx.enumerate_all_cliques(G)
        if min_size <= len(c) <= max_size
    ]

    used_nodes = set()
    selected_cliques = []
    for clique in sorted(all_cliques, key=len, reverse=True):
        if not any(n in used_nodes for n in clique):
            selected_cliques.append(clique)
            used_nodes.update(clique)

    G_new = G.copy()

    for i, members in enumerate(selected_cliques):
        clique_id = f"clique_{level}_{i}"
        member_ids = [str(m) for m in members]

        xs = [float(G.nodes[m].get("x", 0.0)) for m in members]
        ys = [float(G.nodes[m].get("y", 0.0)) for m in members]
        cx, cy = sum(xs) / len(xs), sum(ys) / len(ys)

        bbox = {"minX": min(xs), "maxX": max(xs), "minY": min(ys), "maxY": max(ys)}
        clique_type = "K4" if len(members) == 4 else "K3"

        clique_rows.append((
            clique_id, clique_type, level, cx, cy,
            None, member_ids, json.dumps(bbox)
        ))

        internal = set(members)
        for m in members:
            for neighbor in G.neighbors(m):
                if neighbor not in internal:
                    edge_rows.append((clique_id, str(neighbor), level, True, 1.0))

        G_new.add_node(clique_id, x=cx, y=cy)
        for m in members:
            for neighbor in list(G_new.neighbors(m)):
                if neighbor not in internal and neighbor != clique_id:
                    G_new.add_edge(clique_id, neighbor)
        G_new.remove_nodes_from(members)

        for m in member_ids:
            node_to_clique[m] = clique_id

    execute_values(cur, """
        INSERT INTO cliques (id, clique_type, level, centroid_x, centroid_y,
                             parent_id, member_ids, bbox)
        VALUES %s
        ON CONFLICT (id) DO NOTHING
    """, clique_rows)

    if edge_rows:
        execute_values(cur, """
            INSERT INTO edges (source_id, target_id, level, is_boundary, weight)
            VALUES %s
        """, edge_rows)

    # Two separate execute calls — psycopg2 doesn't support multiple statements in one
    for member_id, clique_id in node_to_clique.items():
        cur.execute("UPDATE nodes SET parent_id = %s WHERE id = %s", (clique_id, member_id))
        cur.execute("UPDATE cliques SET parent_id = %s WHERE id = %s", (clique_id, member_id))

    conn.commit()
    print(f"Level {level}: {len(selected_cliques)} cliques, {G_new.number_of_nodes()} nodes remaining")
    return G_new, G_new.number_of_nodes()


def run_pipeline(tissue_name, cur, conn, depth):
    clear_database(cur, conn)
    G = load_graph(tissue_name)
    names, emb_matrix = load_layout(tissue_name)
    layout = build_layout_dict(G, names, emb_matrix)

    # Stamp x, y onto graph nodes so clique centroid math works
    for node_id, (x, y) in layout.items():
        G.nodes[node_id]["x"] = x
        G.nodes[node_id]["y"] = y

    export_nodes(G, layout, cur, conn)
    export_edges(G, level=0, cur=cur, conn=conn)

    level = 0
    prev_count = G.number_of_nodes()

    for i in range(depth):
        G, new_count = find_and_collapse_cliques(G, level=level, cur=cur, conn=conn)

        export_edges(G, level=level + 1, cur=cur, conn=conn)
        level += 1

        if new_count == prev_count:
            print(f"Converged at level {level} with {new_count} nodes")
            break
        prev_count = new_count
def clear_database(cur, conn):
    cur.execute("""
        TRUNCATE TABLE edges;
        TRUNCATE TABLE nodes CASCADE;
        TRUNCATE TABLE cliques CASCADE;
    """)
    conn.commit()
    print("Database cleared")

# --- Entry point ---
conn = psycopg2.connect(
    host="localhost",
    dbname="viz_db",
    user="viz_user",
    password="viz_password"
)
cur = conn.cursor()

run_pipeline("Liver", cur, conn, 20)

cur.close()
conn.close()