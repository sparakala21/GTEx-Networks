"""
Gene Co-expression Network — Hierarchy Exporter
================================================
Builds the full clique-coarsening hierarchy and exports a single JSON file
that a frontend (Dash, React, etc.) can use to render the interactive graph.

Usage
-----
    python export_hierarchy.py [tissue_name] [threshold] [output_path]

    tissue_name  : e.g. "Liver"  (default: "Liver")
    threshold    : float 0-1     (default: 0.5)
    output_path  : path to write (default: "graph_export.json")

Output JSON schema
------------------
{
  "meta": {
    "tissue":    str,
    "threshold": float,
    "n_levels":  int           // total hierarchy levels (0 = coarsest)
  },

  "levels": [                  // ordered 0 = coarsest … N = individual genes
    {
      "level_idx":   int,
      "label":       str,      // e.g. "Coarsening Level 3" / "Individual Nodes"
      "is_full_res": bool,
      "nodes": {
        "<node_id>": {
          "x":             float,
          "y":             float,
          "size":          int,   // number of genes collapsed into this node
          "clique_weight": float,
          "members":       [str, ...],   // gene names inside this node
          "label":         str,
          "degree":        int           // only present on full-res nodes
        },
        ...
      },
      "edges": [
        {"source": str, "target": str, "weight": float},
        ...
      ]
    },
    ...
  ],

  "initial_visible": [
    {"node": str, "level": int},
    ...
  ],

  // Pre-computed children map so the frontend never has to recompute it.
  // children_map[level_idx][node_id] = [{"node": str, "level": int}, ...]
  "children_map": {
    "<level_idx>": {
      "<node_id>": [{"node": str, "level": int}, ...],
      ...
    },
    ...
  }
}

Requires: networkx, numpy, pandas
"""

import sys
import json
import numpy as np
import networkx as nx
import pandas as pd
from itertools import combinations
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional


# ══════════════════════════════════════════════════════════════════════════
# Clique-coarsening core
# ══════════════════════════════════════════════════════════════════════════

def _clique_weight_sum(G, clique):
    return sum(
        G[u][v].get("weight", 1)
        for u, v in combinations(clique, 2)
        if G.has_edge(u, v)
    )


def _find_target_cliques(G):
    cliques = set()
    for clique in nx.find_cliques(G):
        size = len(clique)
        if size >= 3:
            if size >= 4:
                for sub in combinations(clique, 4):
                    cliques.add(frozenset(sub))
            for sub in combinations(clique, 3):
                cliques.add(frozenset(sub))
    return list(cliques)


def _assign_nodes_to_cliques(G, cliques):
    clique_weights = {i: _clique_weight_sum(G, c) for i, c in enumerate(cliques)}
    node_to_clique = {}
    for node in G.nodes():
        best_idx, best_w = None, -1
        for i, clique in enumerate(cliques):
            if node in clique:
                w = clique_weights[i]
                if w > best_w:
                    best_w, best_idx = w, i
        if best_idx is not None:
            node_to_clique[node] = best_idx
    return node_to_clique, clique_weights


def _average_coords(members_flat, coord_cache):
    coords = [coord_cache[m] for m in members_flat if m in coord_cache]
    if not coords:
        return (0.0, 0.0)
    return (
        sum(c[0] for c in coords) / len(coords),
        sum(c[1] for c in coords) / len(coords),
    )


def _resolve_members(raw_members, previous_members):
    if previous_members is None:
        return set(raw_members)
    resolved = set()
    for m in raw_members:
        resolved |= previous_members.get(m, {m})
    return resolved


def _coarsen_one_pass(G, previous_members=None, coord_cache=None):
    cliques = _find_target_cliques(G)

    if not cliques:
        flat = {n: _resolve_members({n}, previous_members) for n in G.nodes()}
        return G.copy(), {n: n for n in G.nodes()}, {}, flat

    node_to_clique, clique_weights = _assign_nodes_to_cliques(G, cliques)

    supernode_members_raw = defaultdict(set)
    for node, cidx in node_to_clique.items():
        supernode_members_raw[f"C{cidx}"].add(node)
    for node in G.nodes():
        if node not in node_to_clique:
            supernode_members_raw[f"S_{node}"].add(node)

    node_map = {
        node: sn
        for sn, members in supernode_members_raw.items()
        for node in members
    }
    supernode_members_flat = {
        sn: _resolve_members(members, previous_members)
        for sn, members in supernode_members_raw.items()
    }

    CG = nx.Graph()
    for sn, members in supernode_members_raw.items():
        cidx = int(sn[1:]) if sn.startswith("C") else None
        flat = supernode_members_flat[sn]
        x, y = _average_coords(flat, coord_cache) if coord_cache else (0.0, 0.0)
        CG.add_node(sn,
                    members=supernode_members_flat[sn],
                    clique_weight=clique_weights.get(cidx, 0.0),
                    size=len(supernode_members_flat[sn]),
                    x=x, y=y)

    inter_weights = defaultdict(float)
    for u, v, data in G.edges(data=True):
        su, sv = node_map[u], node_map[v]
        if su != sv:
            key = tuple(sorted([su, sv]))
            inter_weights[key] += data.get("weight", 1.0)
    for (su, sv), w in inter_weights.items():
        CG.add_edge(su, sv, weight=w)

    supernode_info = {
        sn: {
            "members_flat":  supernode_members_flat[sn],
            "clique_weight": CG.nodes[sn]["clique_weight"],
            "size":          CG.nodes[sn]["size"],
            "x":             CG.nodes[sn]["x"],
            "y":             CG.nodes[sn]["y"],
        }
        for sn in supernode_members_raw
    }
    return CG, node_map, supernode_info, supernode_members_flat


def coarsen_iterative(G, coord_cache=None, max_iterations=20, min_nodes=5, verbose=True):
    current_G        = G
    previous_members = None
    levels           = []
    iteration = 0
    while nx.find_cliques(G):
        n_before = current_G.number_of_nodes()
        if verbose:
            print(f"[coarsen] iteration {iteration}: "
                  f"{n_before} nodes, {current_G.number_of_edges()} edges")

        CG, node_map, info, members_flat = _coarsen_one_pass(
            current_G, previous_members=previous_members, coord_cache=coord_cache)
        n_after = CG.number_of_nodes()

        levels.append({"iteration": iteration, "graph": CG,
                        "node_map": node_map, "info": info, "members": members_flat})

        if verbose:
            print(f"          -> {n_after} nodes  (delta {n_before - n_after})")

        if n_after >= n_before:
            if verbose: print("          No further reduction -- stopping.")
            break
        if n_after <= min_nodes:
            if verbose: print(f"          Reached min_nodes={min_nodes} -- stopping.")
            break

        current_G        = CG
        previous_members = members_flat
        iteration += 1

    return levels


# ══════════════════════════════════════════════════════════════════════════
# Data structure
# ══════════════════════════════════════════════════════════════════════════

@dataclass
class CoarseLevel:
    level_idx:   int
    label:       str
    G:           nx.Graph
    pos:         Dict
    node_meta:   Dict
    is_full_res: bool = False


# ══════════════════════════════════════════════════════════════════════════
# GraphHierarchy — builds levels, no UI dependencies
# ══════════════════════════════════════════════════════════════════════════

class GraphHierarchy:
    """
    Precomputes the full clique-coarsening hierarchy.

    UI levels stored in self.levels:
        0         = most coarse (last coarsening iteration, fewest nodes)
        1 ... N-1 = progressively finer coarse levels
        N         = full resolution (individual genes)
    """

    def __init__(self, tissue_name: str, threshold: float = 0.5):
        self.tissue_name = tissue_name
        self.threshold   = threshold

        print(f"[hierarchy] Loading '{tissue_name}' @ threshold={threshold} ...")
        self.G_full = load_graph(tissue_name, threshold)
        names, emb  = load_embeddings(tissue_name)
        self.coord_cache = {
            str(n): (float(emb[i, 0]), float(emb[i, 1]))
            for i, n in enumerate(names)
        }
        print(f"[hierarchy] {self.G_full.number_of_nodes()} nodes, "
              f"{self.G_full.number_of_edges()} edges")

        self.levels: List[CoarseLevel] = []
        self._build_all_levels()

    # ── builders ──────────────────────────────────────────────────────────

    def _build_all_levels(self):
        raw_levels = coarsen_iterative(
            self.G_full, coord_cache=self.coord_cache,
            max_iterations=20, min_nodes=5, verbose=True)

        # Reverse: UI level 0 = most coarse (last iteration)
        for ui_idx, raw in enumerate(reversed(raw_levels)):
            lvl = self._raw_to_coarse_level(ui_idx, raw)
            self.levels.append(lvl)
            print(f"[hierarchy] UI level {ui_idx} '{lvl.label}': "
                  f"{lvl.G.number_of_nodes()} nodes")

        full_idx = len(self.levels)
        self.levels.append(self._build_full_resolution(full_idx))
        print(f"[hierarchy] UI level {full_idx} 'Individual Nodes': "
              f"{self.levels[-1].G.number_of_nodes()} nodes")

    def _raw_to_coarse_level(self, ui_idx: int, raw: dict) -> CoarseLevel:
        CG, info = raw["graph"], raw["info"]
        label    = f"Coarsening Level {raw['iteration']}"
        pos, node_meta = {}, {}

        for node in CG.nodes():
            nd  = CG.nodes[node]
            pos[node] = (float(nd.get("x", 0.0)), float(nd.get("y", 0.0)))

            if node in info:
                d = info[node]
                members       = list(d.get("members_flat", {node}))
                clique_weight = float(d.get("clique_weight", 0.0))
                size          = int(d.get("size", 1))
            else:
                raw_members = raw["members"].get(node, {node})   # fully-resolved gene names
                members     = list(raw_members)
                size        = int(nd.get("size", len(members)))
                clique_weight = float(nd.get("clique_weight", 0.0))
                size          = int(nd.get("size", 1))

            node_meta[node] = {
                "size":          size,
                "clique_weight": clique_weight,
                "members":       members,
                "label":         str(node),
            }

        return CoarseLevel(ui_idx, label, CG, pos, node_meta)

    def _build_full_resolution(self, idx: int) -> CoarseLevel:
        G   = nx.Graph(self.G_full)
        pos = {n: (float(self.coord_cache[n][0]), float(self.coord_cache[n][1]))
               for n in G.nodes() if n in self.coord_cache}
        meta = {n: {"size": 1, "clique_weight": 0.0,
                    "members": [n], "label": str(n),
                    "degree": G.degree(n)}
                for n in G.nodes()}
        return CoarseLevel(idx, "Individual Nodes", G, pos, meta, is_full_res=True)

    # ── hierarchy queries ─────────────────────────────────────────────────

    @property
    def n_levels(self) -> int:
        return len(self.levels)

    def can_expand(self, ui_level: int) -> bool:
        return ui_level < self.n_levels - 1

    def initial_visible(self) -> List[dict]:
        """All nodes from the coarsest level — the starting expansion state."""
        return [{"node": str(n), "level": 0}
                for n in self.levels[0].G.nodes()]

    def get_children(self, ui_level: int, node_id: str) -> List[dict]:
        """
        Returns the one-level-finer constituents of node_id.

        Looks at ui_level+1 and returns all nodes whose gene-membership
        is a strict subset of node_id's gene-membership.
        Falls back to returning the node unchanged if no children found.
        """
        next_level = ui_level + 1
        if next_level >= self.n_levels:
            return []

        parent_members = set(self.levels[ui_level].node_meta[node_id]["members"])
        next_lvl       = self.levels[next_level]

        children = [
            {"node": str(node), "level": next_level}
            for node in next_lvl.G.nodes()
            if set(next_lvl.node_meta[node]["members"]) <= parent_members
        ]
        return children if children else [{"node": node_id, "level": ui_level}]


# ══════════════════════════════════════════════════════════════════════════
# I/O helpers
# ══════════════════════════════════════════════════════════════════════════

def load_graph(tissue_name, threshold=0.5):
    path = f"../tissue_networks/{tissue_name.replace(' ', '_')}_network.gexf"
    G    = nx.read_gexf(path)
    G.remove_edges_from([
        (u, v) for u, v, d in G.edges(data=True)
        if d.get("weight", 0) < threshold
    ])
    return G


def load_embeddings(tissue_name):
    df         = pd.read_csv(f"../node2vec_layout/{tissue_name}_embeddings_2d.csv")
    names      = df["node"].tolist()
    emb_matrix = df.drop(columns=["node"]).to_numpy()
    print(emb_matrix.shape)
    return names, emb_matrix


# ══════════════════════════════════════════════════════════════════════════
# Export
# ══════════════════════════════════════════════════════════════════════════

def _level_edges(lvl: CoarseLevel, G_full: nx.Graph) -> List[dict]:
    """
    Reconstruct edges for a mixed/coarse level from G_full.
    Each node at this level 'owns' a set of genes; we draw an edge between
    two nodes if any of their genes are connected in G_full.
    """
    # Build gene → supernode map for this level
    gene_to_node: Dict[str, str] = {}
    for nid, meta in lvl.node_meta.items():
        for gene in meta["members"]:
            gene_to_node[gene] = nid

    edge_weights: Dict[tuple, float] = defaultdict(float)
    for u, v, data in G_full.edges(data=True):
        su = gene_to_node.get(u)
        sv = gene_to_node.get(v)
        if su and sv and su != sv:
            key = (min(su, sv), max(su, sv))
            edge_weights[key] += data.get("weight", 1.0)

    return [
        {"source": s, "target": t, "weight": round(w, 6)}
        for (s, t), w in edge_weights.items()
    ]


def _serialise_level(lvl: CoarseLevel, G_full: nx.Graph) -> dict:
    nodes_out = {}
    for nid, meta in lvl.node_meta.items():
        x, y = lvl.pos.get(nid, (0.0, 0.0))
        entry = {
            "x":             round(x, 6),
            "y":             round(y, 6),
            "size":          meta["size"],
            "clique_weight": round(meta["clique_weight"], 6),
            "members":       list(meta["members"]),
            "label":         meta["label"],
        }
        if lvl.is_full_res:
            entry["degree"] = meta.get("degree", 0)
        nodes_out[nid] = entry

    return {
        "level_idx":   lvl.level_idx,
        "label":       lvl.label,
        "is_full_res": lvl.is_full_res,
        "nodes":       nodes_out,
        "edges":       _level_edges(lvl, G_full),
    }


def build_children_map(h: GraphHierarchy) -> dict:
    """
    Pre-compute the full children map so the frontend does zero graph math.
    Structure: { level_idx (str) -> { node_id -> [child dicts] } }
    Only includes levels that can be expanded (all except the last).
    """
    children_map = {}
    for lvl in h.levels[:-1]:           # skip the last (full-res) level
        level_key = str(lvl.level_idx)
        children_map[level_key] = {}
        for nid in lvl.G.nodes():
            children_map[level_key][nid] = h.get_children(lvl.level_idx, nid)
    return children_map


def export_hierarchy(tissue_name: str, threshold: float, output_path: str):
    h = GraphHierarchy(tissue_name, threshold)

    print("[export] Serialising levels ...")
    levels_out = [_serialise_level(lvl, h.G_full) for lvl in h.levels]

    print("[export] Building children map ...")
    children_map = build_children_map(h)

    payload = {
        "meta": {
            "tissue":    tissue_name,
            "threshold": threshold,
            "n_levels":  h.n_levels,
        },
        "levels":          levels_out,
        "initial_visible": h.initial_visible(),
        "children_map":    children_map,
    }

    print(f"[export] Writing → {output_path}")
    with open(output_path, "w") as f:
        json.dump(payload, f, separators=(",", ":"), indent=2)   # compact — no extra whitespace

    # Human-readable summary
    print("\n── Export summary ─────────────────────────────────────────")
    print(f"  Tissue    : {tissue_name}")
    print(f"  Threshold : {threshold}")
    print(f"  Levels    : {h.n_levels}  (0 = coarsest, {h.n_levels-1} = individual genes)")
    for lvl in h.levels:
        n_nodes = len(lvl.node_meta)
        n_edges = len(_level_edges(lvl, h.G_full))
        flag = " ← initial view" if lvl.level_idx == 0 else ""
        print(f"    [{lvl.level_idx}] {lvl.label:<30} {n_nodes:>5} nodes  {n_edges:>6} edges{flag}")
    print(f"\n  Output    : {output_path}")
    print("───────────────────────────────────────────────────────────\n")


# ── Entry point ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    tissue_name = sys.argv[1] if len(sys.argv) > 1 else "Liver"
    threshold   = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    output_path = sys.argv[3] if len(sys.argv) > 3 else "graph_export.json"

    export_hierarchy(tissue_name, threshold, output_path)