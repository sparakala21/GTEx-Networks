"""
Gene Co-expression Network — K-Core Hierarchy Exporter
=======================================================
Builds a k-core decomposition hierarchy and exports a single JSON file
that a frontend (Dash, React, etc.) can use to render the interactive graph.

K-core levels are ordered from most restrictive (fewest nodes, densest core)
to least restrictive (all nodes), matching the original "coarsest-first" convention:

    Level 0  →  k_max-core  (densest subgraph, initial view)
    Level 1  →  (k_max-1)-core
    ...
    Level N  →  1-core  (all nodes with at least one edge)

Usage
-----
    python export_hierarchy_kcore.py [tissue_name] [threshold] [output_path]

    tissue_name  : e.g. "Liver"  (default: "Liver")
    threshold    : float 0-1     (default: 0.5)
    output_path  : path to write (default: "graph_export_kcore.json")

Output JSON schema
------------------
Identical to the clique-coarsening exporter, with two additions in "meta":

{
  "meta": {
    "tissue":        str,
    "threshold":     float,
    "n_levels":      int,
    "decomposition": "kcore",       // NEW — tells the frontend which mode
    "k_values":      [int, ...]     // NEW — k value for each level index
  },
  ...
}

children_map notes
------------------
Because k-cores are *nested subgraphs* (not supernode hierarchies), the
expansion semantics differ slightly from clique-coarsening:

  children_map[level_idx][node_id]
      → [{"node": node_id, "level": level_idx + 1}]
        The same node persists into the next (less restrictive) k-core.

  children_map[level_idx]["__new__"]
      → [{"node": str, "level": level_idx + 1}, ...]
        Nodes that first *appear* at level_idx + 1 (their coreness equals
        exactly k at that level). These have no single parent; the frontend
        should surface them as newly visible when the user steps outward.

Requires: networkx, numpy, pandas
"""

import sys
import json
import numpy as np
import networkx as nx
import pandas as pd
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ══════════════════════════════════════════════════════════════════════════
# K-core decomposition
# ══════════════════════════════════════════════════════════════════════════

def kcore_decomposition(G: nx.Graph, verbose: bool = True) -> List[dict]:
    """
    Compute every distinct k-core of G.

    Returns a list of dicts ordered from the *highest* k (densest) down to
    k = 1, skipping any k whose subgraph is identical in size to the next:

        [
          {"k": k_max, "graph": <subgraph>},
          {"k": k_max - 1, "graph": <subgraph>},
          ...
          {"k": 1, "graph": <subgraph>},
        ]

    Empty k-cores (no nodes) are omitted.
    """
    if G.number_of_nodes() == 0:
        return []

    coreness: Dict[str, int] = nx.core_number(G)
    k_max = max(coreness.values())
    k_min = 1  # 0-core is the whole graph incl. isolates — not useful here

    if verbose:
        print(f"[kcore] max coreness = {k_max}, min = {k_min}")

    levels = []
    seen_sizes = set()

    for k in range(k_max, k_min - 1, -1):
        subgraph = nx.k_core(G, k, core_number=coreness)
        n = subgraph.number_of_nodes()
        if n == 0:
            continue
        # Skip k values that produce a subgraph identical in node-set to the
        # previous level (can happen when no node has coreness exactly k).
        node_key = frozenset(subgraph.nodes())
        if node_key in seen_sizes:
            if verbose:
                print(f"[kcore]   k={k}: {n} nodes (duplicate — skipped)")
            continue
        seen_sizes.add(node_key)
        if verbose:
            print(f"[kcore]   k={k}: {n} nodes, {subgraph.number_of_edges()} edges")
        levels.append({"k": k, "graph": subgraph})

    return levels


# ══════════════════════════════════════════════════════════════════════════
# Data structure
# ══════════════════════════════════════════════════════════════════════════

@dataclass
class CoreLevel:
    level_idx:   int
    k:           int           # the k value for this core
    label:       str
    G:           nx.Graph
    pos:         Dict
    node_meta:   Dict
    is_full_res: bool = False  # True only for the outermost (1-core) level


# ══════════════════════════════════════════════════════════════════════════
# GraphHierarchy — builds k-core levels, no UI dependencies
# ══════════════════════════════════════════════════════════════════════════

class GraphHierarchy:
    """
    Precomputes the full k-core hierarchy.

    self.levels ordering:
        0         = k_max-core  (densest, fewest nodes — initial view)
        1 ... N-1 = progressively less restrictive cores
        N         = 1-core  (most nodes; flagged is_full_res=True)
    """

    def __init__(self, tissue_name: str, threshold: float = 0.5):
        self.tissue_name = tissue_name
        self.threshold   = threshold

        print(f"[hierarchy] Loading '{tissue_name}' @ threshold={threshold} ...")
        self.G_full = load_graph(tissue_name, threshold)
        names, emb  = load_embeddings(tissue_name)
        self.coord_cache: Dict[str, Tuple[float, float]] = {
            str(n): (float(emb[i, 0]), float(emb[i, 1]))
            for i, n in enumerate(names)
        }
        print(f"[hierarchy] {self.G_full.number_of_nodes()} nodes, "
              f"{self.G_full.number_of_edges()} edges")

        self.levels:   List[CoreLevel] = []
        self.k_values: List[int]       = []   # k value per level index
        self._build_all_levels()

    # ── builders ──────────────────────────────────────────────────────────

    def _build_all_levels(self):
        raw_levels = kcore_decomposition(self.G_full, verbose=True)

        if not raw_levels:
            raise ValueError("Graph has no k-cores (possibly empty or fully disconnected).")

        for ui_idx, raw in enumerate(raw_levels):
            is_last = (ui_idx == len(raw_levels) - 1)
            lvl = self._raw_to_core_level(ui_idx, raw, is_full_res=is_last)
            self.levels.append(lvl)
            self.k_values.append(raw["k"])
            print(f"[hierarchy] UI level {ui_idx} '{lvl.label}': "
                  f"{lvl.G.number_of_nodes()} nodes  (k={raw['k']})")

    def _raw_to_core_level(
        self,
        ui_idx: int,
        raw: dict,
        is_full_res: bool = False,
    ) -> CoreLevel:
        k: int         = raw["k"]
        SG: nx.Graph   = raw["graph"]
        label          = f"{k}-core" if not is_full_res else f"1-core (all nodes)"
        coreness       = nx.core_number(SG)

        pos, node_meta = {}, {}
        for node in SG.nodes():
            node_s = str(node)
            x, y   = self.coord_cache.get(node_s, (0.0, 0.0))
            pos[node_s] = (x, y)
            node_meta[node_s] = {
                "size":          1,
                "clique_weight": float(coreness.get(node, 0)),  # reuse field for coreness
                "members":       [node_s],
                "label":         node_s,
                "degree":        SG.degree(node),
            }

        return CoreLevel(
            level_idx=ui_idx,
            k=k,
            label=label,
            G=SG,
            pos=pos,
            node_meta=node_meta,
            is_full_res=is_full_res,
        )

    # ── hierarchy queries ─────────────────────────────────────────────────

    @property
    def n_levels(self) -> int:
        return len(self.levels)

    def can_expand(self, ui_level: int) -> bool:
        return ui_level < self.n_levels - 1

    def initial_visible(self) -> List[dict]:
        """All nodes from the densest core — the starting view."""
        return [{"node": str(n), "level": 0}
                for n in self.levels[0].G.nodes()]

    def get_children(self, ui_level: int, node_id: str) -> List[dict]:
        """
        For k-core hierarchies the node itself persists into the next
        (less restrictive) level, so the 'child' is the same node one
        level out.  Returns [] if already at the outermost level.
        """
        next_level = ui_level + 1
        if next_level >= self.n_levels:
            return []
        # The node always exists in the next level (k-cores are nested).
        return [{"node": node_id, "level": next_level}]

    def get_new_nodes_at(self, ui_level: int) -> List[dict]:
        """
        Nodes that *first appear* at ui_level — i.e. their coreness equals
        exactly the k of the previous level minus one.
        For ui_level == 0 this is the full node set of that level.
        """
        if ui_level == 0:
            return [{"node": str(n), "level": 0}
                    for n in self.levels[0].G.nodes()]

        prev_nodes = set(self.levels[ui_level - 1].G.nodes())
        curr_nodes = set(self.levels[ui_level].G.nodes())
        new_nodes  = curr_nodes - prev_nodes
        return [{"node": str(n), "level": ui_level} for n in sorted(new_nodes)]


# ══════════════════════════════════════════════════════════════════════════
# I/O helpers
# ══════════════════════════════════════════════════════════════════════════

def load_graph(tissue_name: str, threshold: float = 0.5) -> nx.Graph:
    path = f"../tissue_networks/{tissue_name.replace(' ', '_')}_network.gexf"
    G    = nx.read_gexf(path)
    G.remove_edges_from([
        (u, v)
        for u, v, d in G.edges(data=True)
        if d.get("weight", 0) < threshold
    ])
    return G


def load_embeddings(tissue_name: str):
    df         = pd.read_csv(f"../node2vec_layout/{tissue_name}_embeddings_2d.csv")
    names      = df["node"].tolist()
    emb_matrix = df.drop(columns=["node"]).to_numpy()
    print(emb_matrix.shape)
    return names, emb_matrix


# ══════════════════════════════════════════════════════════════════════════
# Serialisation
# ══════════════════════════════════════════════════════════════════════════

def _level_edges(lvl: CoreLevel, G_full: nx.Graph) -> List[dict]:
    """
    Edges for a k-core level — only edges whose *both* endpoints are present
    in this level's subgraph, using weights from G_full.
    """
    node_set = set(lvl.G.nodes())
    edges_out = []
    for u, v, data in G_full.edges(data=True):
        if u in node_set and v in node_set:
            edges_out.append({
                "source": str(u),
                "target": str(v),
                "weight": round(float(data.get("weight", 1.0)), 6),
            })
    return edges_out


def _serialise_level(lvl: CoreLevel, G_full: nx.Graph) -> dict:
    nodes_out = {}
    for nid, meta in lvl.node_meta.items():
        x, y = lvl.pos.get(nid, (0.0, 0.0))
        nodes_out[nid] = {
            "x":             round(x, 6),
            "y":             round(y, 6),
            "size":          meta["size"],
            "clique_weight": round(meta["clique_weight"], 6),  # stores coreness here
            "members":       list(meta["members"]),
            "label":         meta["label"],
            "degree":        meta.get("degree", 0),
            "k":             lvl.k,
        }
    return {
        "level_idx":   lvl.level_idx,
        "label":       lvl.label,
        "k":           lvl.k,
        "is_full_res": lvl.is_full_res,
        "nodes":       nodes_out,
        "edges":       _level_edges(lvl, G_full),
    }


def build_children_map(h: GraphHierarchy) -> dict:
    """
    Pre-compute children map for every non-outermost level.

    children_map[level_idx][node_id]
        → list of child dicts for that node (same node, next level)

    children_map[level_idx]["__new__"]
        → list of node dicts that *first appear* at level_idx + 1
          (i.e. their coreness means they enter the picture one step out)
    """
    children_map: Dict[str, dict] = {}

    for lvl in h.levels[:-1]:   # skip outermost; nothing to expand into
        level_key = str(lvl.level_idx)
        children_map[level_key] = {}

        # Per-node children (same node persists one level out)
        for nid in lvl.G.nodes():
            children_map[level_key][str(nid)] = h.get_children(lvl.level_idx, str(nid))

        # New nodes that appear at the next level
        next_idx = lvl.level_idx + 1
        new_nodes = h.get_new_nodes_at(next_idx)
        if new_nodes:
            children_map[level_key]["__new__"] = new_nodes

    return children_map


# ══════════════════════════════════════════════════════════════════════════
# Top-level export
# ══════════════════════════════════════════════════════════════════════════

def export_hierarchy(tissue_name: str, threshold: float, output_path: str):
    h = GraphHierarchy(tissue_name, threshold)

    print("[export] Serialising levels ...")
    levels_out = [_serialise_level(lvl, h.G_full) for lvl in h.levels]

    print("[export] Building children map ...")
    children_map = build_children_map(h)

    payload = {
        "meta": {
            "tissue":        tissue_name,
            "threshold":     threshold,
            "n_levels":      h.n_levels,
            "decomposition": "kcore",      # signals k-core mode to the frontend
            "k_values":      h.k_values,   # [k_max, k_max-1, ..., 1] (one per level)
        },
        "levels":          levels_out,
        "initial_visible": h.initial_visible(),
        "children_map":    children_map,
    }

    print(f"[export] Writing → {output_path}")
    with open(output_path, "w") as f:
        json.dump(payload, f, separators=(",", ":"), indent=2)

    # ── Human-readable summary ────────────────────────────────────────────
    print("\n── Export summary ─────────────────────────────────────────")
    print(f"  Tissue    : {tissue_name}")
    print(f"  Threshold : {threshold}")
    print(f"  Levels    : {h.n_levels}  (0 = densest k-core, {h.n_levels - 1} = 1-core)")
    for lvl in h.levels:
        n_nodes  = len(lvl.node_meta)
        n_edges  = len(_level_edges(lvl, h.G_full))
        new_here = len(h.get_new_nodes_at(lvl.level_idx))
        flag     = " ← initial view" if lvl.level_idx == 0 else ""
        print(
            f"    [{lvl.level_idx}] k={lvl.k:<4} {lvl.label:<25} "
            f"{n_nodes:>5} nodes  {n_edges:>6} edges  "
            f"+{new_here} new{flag}"
        )
    print(f"\n  Output    : {output_path}")
    print("───────────────────────────────────────────────────────────\n")


# ── Entry point ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    tissue_name = sys.argv[1] if len(sys.argv) > 1 else "Liver"
    threshold   = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    output_path = sys.argv[3] if len(sys.argv) > 3 else "graph_export_kcore.json"

    export_hierarchy(tissue_name, threshold, output_path)