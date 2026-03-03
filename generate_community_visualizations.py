"""
Usage
-----
    python dynamic_graph_visualizer.py

Requires: dash, plotly, networkx, numpy (and your load_graph / load_embeddings helpers)
"""


import numpy as np
import networkx as nx
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output, State, callback_context
import dash_bootstrap_components as dbc   # pip install dash-bootstrap-components
from functools import lru_cache
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
import colorsys
import pandas as pd
import sys


# ── Palette & constants ────────────────────────────────────────────────────
LEVELS = [
    # (louvain_resolution, label, min_zoom_span)
    # zoom_span = visible x-range width; smaller → more zoomed in
    (0.3,  "Macro Communities",   float("inf")),   # level 0 – most coarse
    (0.7,  "Communities",         0.60),
    (1.5,  "Sub-communities",     0.30),
    (3.0,  "Fine clusters",       0.12),
    (None, "Individual Nodes",    0.05),            # level 4 – full resolution
]

COLORSCALES = [
    "#e63946", "#f4a261", "#2a9d8f", "#457b9d", "#6a0572",
    "#e9c46a", "#264653", "#a8dadc", "#f1faee", "#ff6b6b",
    "#48cae4", "#90e0ef", "#ade8f4", "#caf0f8", "#023e8a",
]

BACKGROUND   = "#f0fff0"
EDGE_COLOR   = "rgba(100,180,255,0.18)"
GRID_COLOR   = "#1c2333"
TEXT_COLOR   = "#000080"
ACCENT_COLOR = "#58a6ff"

def load_graph(tissue_name, threshold=0.5):
    graph_path = f"tissue_networks/{tissue_name.replace(' ', '_')}_network.gexf"
    G = nx.read_gexf(graph_path)
    G_filtered = G.copy()
    edges_to_remove = [(u, v) for u, v, data in G_filtered.edges(data=True) 
                    if data.get('weight', 0) < threshold]
    G_filtered.remove_edges_from(edges_to_remove)
    return G_filtered

def load_embeddings(tissue_name):
    embeddings = pd.read_csv(f"node2vec_layout/{tissue_name}_embeddings_2d.csv")
    names = embeddings['node'].tolist()
    emb_matrix = embeddings.drop(columns=['node']).to_numpy()
    print(emb_matrix.shape)
    return names, emb_matrix


# ── Data structures ────────────────────────────────────────────────────────
@dataclass
class CoarseLevel:
    level_idx:  int
    label:      str
    G:          nx.Graph                        # coarsened graph
    pos:        Dict[int, Tuple[float, float]]  # node → (x, y)
    node_meta:  Dict[int, dict]                 # size, community_id, members, ...


# ── Hierarchy builder ──────────────────────────────────────────────────────
class GraphHierarchy:
    """Pre-computes the full coarsening hierarchy once at startup."""

    def __init__(self, tissue_name: str, threshold: float = 0.5):
        self.tissue_name = tissue_name
        self.threshold   = threshold

        print(f"[hierarchy] Loading graph '{tissue_name}' @ threshold={threshold} …")
        self.G_full = load_graph(tissue_name, threshold)
        self.names, self.emb = load_embeddings(tissue_name)
        self.name_to_emb = {n: self.emb[i] for i, n in enumerate(self.names)}

        print(f"[hierarchy] {self.G_full.number_of_nodes()} nodes, "
              f"{self.G_full.number_of_edges()} edges")

        self.levels: List[CoarseLevel] = []
        self._build_all_levels()

    # ── level builders ────────────────────────────────────────────────────

    def _build_all_levels(self):
        for idx, (resolution, label, _) in enumerate(LEVELS):
            if resolution is None:
                lvl = self._build_full_resolution(idx, label)
            else:
                lvl = self._build_coarse_level(idx, label, resolution)
            self.levels.append(lvl)
            print(f"[hierarchy] Level {idx} '{label}': {lvl.G.number_of_nodes()} nodes")

    def _build_coarse_level(self, idx: int, label: str, resolution: float) -> CoarseLevel:
        communities = nx.community.louvain_communities(
            self.G_full, weight="weight", seed=42, resolution=resolution
        )

        CG   = nx.Graph()
        pos  = {}
        meta = {}

        for i, community in enumerate(communities):
            members = list(community)
            embs    = [self.name_to_emb[n] for n in members if n in self.name_to_emb]
            center  = np.mean(embs, axis=0) if embs else np.array([0.0, 0.0])
            CG.add_node(i)
            pos[i]  = (float(center[0]), float(center[1]))
            meta[i] = {
                "size":         len(members),
                "community_id": i,
                "members":      members,
                "label":        f"Community {i}",
            }

        for i, ci in enumerate(communities):
            for j in range(i + 1, len(communities)):
                w = sum(1 for u in ci for v in communities[j] if self.G_full.has_edge(u, v))
                if w > 0:
                    CG.add_edge(i, j, weight=w)

        return CoarseLevel(idx, label, CG, pos, meta)

    def _build_full_resolution(self, idx: int, label: str) -> CoarseLevel:
        G   = nx.Graph(self.G_full)
        pos = {n: (float(self.name_to_emb[n][0]), float(self.name_to_emb[n][1]))
               for n in G.nodes() if n in self.name_to_emb}
        meta = {n: {"size": 1, "community_id": -1,
                    "members": [n], "label": n,
                    "degree": G.degree(n)}
                for n in G.nodes()}
        return CoarseLevel(idx, label, G, pos, meta)

    # ── level selection ───────────────────────────────────────────────────

    def level_for_zoom(self, x_span: float) -> int:
        """Return the coarsest level whose min_zoom_span > x_span."""
        for i, (_, _, min_span) in reversed(list(enumerate(LEVELS))):
            if x_span <= min_span:
                return i
        return 0

    # ── figure builder ────────────────────────────────────────────────────

    def build_figure(self,
                     level_idx:    int,
                     x_range:      Optional[Tuple[float, float]] = None,
                     y_range:      Optional[Tuple[float, float]] = None) -> go.Figure:
        lvl  = self.levels[level_idx]
        G    = lvl.G
        pos  = lvl.pos
        meta = lvl.node_meta

        # ── edges ───────────────────────────────────────────────────────
        ex, ey = [], []
        for u, v in G.edges():
            if u in pos and v in pos:
                x0, y0 = pos[u]
                x1, y1 = pos[v]
                ex.extend([x0, x1, None])
                ey.extend([y0, y1, None])

        edge_trace = go.Scatter(
            x=ex, y=ey,
            mode="lines",
            line=dict(width=0.6, color=EDGE_COLOR),
            hoverinfo="none",
        )

        # ── nodes ───────────────────────────────────────────────────────
        nx_list   = []
        ny_list   = []
        sizes     = []
        colors    = []
        texts     = []
        hover     = []

        is_full = (level_idx == len(LEVELS) - 1)

        for node in G.nodes():
            if node not in pos:
                continue
            x, y   = pos[node]
            m      = meta[node]
            nx_list.append(x)
            ny_list.append(y)

            if is_full:
                deg  = G.degree(node)
                sz   = max(4, min(12, 4 + deg * 0.4))
                col  = deg
                txt  = str(node)
                htxt = f"<b>{node}</b><br>Degree: {deg}"
            else:
                sz   = max(8, min(60, 8 + m["size"] ** 0.55 * 3))
                col  = m["community_id"]
                txt  = m["label"] if m["size"] > 5 else ""
                htxt = (f"<b>Community {m['community_id']}</b><br>"
                        f"Genes: {m['size']}<br>"
                        f"Edges out: {G.degree(node)}<br>"
                        f"Top members: {', '.join(m['members'][:5])}"
                        + (" …" if len(m["members"]) > 5 else ""))

            sizes.append(sz)
            colors.append(col)
            texts.append(txt)
            hover.append(htxt)

        node_trace = go.Scatter(
            x=nx_list, y=ny_list,
            mode="markers+text" if not is_full else "markers",
            text=texts,
            textposition="top center",
            textfont=dict(size=9, color=TEXT_COLOR, family="'JetBrains Mono', monospace"),
            hovertext=hover,
            hoverinfo="text",
            marker=dict(
                size=sizes,
                sizemode="diameter",
                color=colors,
                colorscale="Turbo",
                showscale=True,
                colorbar=dict(
                    title=dict(text="Degree" if is_full else "Community",
                               font=dict(color=TEXT_COLOR, size=11)),
                    tickfont=dict(color=TEXT_COLOR),
                    bgcolor="rgba(0,0,0,0)",
                    bordercolor=GRID_COLOR,
                ),
                line=dict(width=0.8, color="rgba(255,255,255,0.15)"),
                opacity=0.90,
            ),
        )

        # ── layout ──────────────────────────────────────────────────────
        layout = go.Layout(
            paper_bgcolor=BACKGROUND,
            plot_bgcolor=BACKGROUND,
            margin=dict(l=0, r=0, t=52, b=0),
            showlegend=False,
            hovermode="closest",
            title=dict(
                text=(f"<span style='font-size:15px;letter-spacing:2px;'>"
                      f"{self.tissue_name.upper()} NETWORK</span>"
                      f"<span style='font-size:11px;color:{ACCENT_COLOR};'>"
                      f"  ·  {lvl.label}  ·  "
                      f"{G.number_of_nodes()} nodes</span>"),
                font=dict(color=TEXT_COLOR, family="'JetBrains Mono', monospace"),
                x=0.02, xanchor="left",
            ),
            xaxis=dict(
                showgrid=False, zeroline=False, showticklabels=False,
                range=x_range,
            ),
            yaxis=dict(
                showgrid=False, zeroline=False, showticklabels=False,
                range=y_range, scaleanchor="x",
            ),
            transition=dict(duration=350, easing="cubic-in-out"),
        )

        return go.Figure(data=[edge_trace, node_trace], layout=layout)


# ── Dash app ───────────────────────────────────────────────────────────────

def build_app(tissue_name: str = "Liver", threshold: float = 0.5) -> Dash:
    hierarchy = GraphHierarchy(tissue_name, threshold)
    initial_fig = hierarchy.build_figure(level_idx=0)

    app = Dash(
        __name__,
        external_stylesheets=[dbc.themes.DARKLY,
                               "https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;600&display=swap"],
        title=f"{tissue_name} Network Explorer",
    )

    # ── UI ────────────────────────────────────────────────────────────────
    level_labels = [l[1] for l in LEVELS]

    app.layout = dbc.Container(fluid=True, style={"background": BACKGROUND,
                                                   "minHeight": "100vh",
                                                   "padding": "0"}, children=[

        # top bar
        dbc.Row(style={"background": "#010409",
                       "borderBottom": f"1px solid {GRID_COLOR}",
                       "padding": "10px 20px", "margin": "0"}, children=[
            dbc.Col(html.H5("Gene Co-expression Explorer",
                            style={"color": ACCENT_COLOR,
                                   "fontFamily": "'JetBrains Mono', monospace",
                                   "margin": "0", "letterSpacing": "1px"}), width=4),
            dbc.Col([
                html.Span("TISSUE ", style={"color": "#555", "fontSize": "11px",
                                            "fontFamily": "'JetBrains Mono'"}),
                dcc.Input(id="tissue-input", value=tissue_name, debounce=True,
                          style={"background": "#161b22", "color": TEXT_COLOR,
                                 "border": f"1px solid {GRID_COLOR}",
                                 "borderRadius": "4px", "padding": "4px 10px",
                                 "fontFamily": "'JetBrains Mono'", "marginRight": "16px"}),
                html.Span("THRESHOLD ", style={"color": "#555", "fontSize": "11px",
                                               "fontFamily": "'JetBrains Mono'"}),
                dcc.Input(id="threshold-input", value=str(threshold),
                          type="number", min=0, max=1, step=0.05, debounce=True,
                          style={"background": "#161b22", "color": TEXT_COLOR,
                                 "border": f"1px solid {GRID_COLOR}",
                                 "borderRadius": "4px", "padding": "4px 10px",
                                 "fontFamily": "'JetBrains Mono'", "width": "90px",
                                 "marginRight": "16px"}),
                dbc.Button("RELOAD", id="reload-btn", color="primary", size="sm",
                           style={"fontFamily": "'JetBrains Mono'", "letterSpacing": "1px"}),
            ], width=5),
            dbc.Col([
                html.Span("DETAIL LEVEL ", style={"color": "#555", "fontSize": "11px",
                                                   "fontFamily": "'JetBrains Mono'"}),
                html.Span(id="level-label", children=level_labels[0],
                          style={"color": ACCENT_COLOR,
                                 "fontFamily": "'JetBrains Mono'",
                                 "fontSize": "12px", "marginLeft": "6px"}),
                html.Br(),
                dcc.Slider(id="level-slider", min=0, max=len(LEVELS)-1,
                           step=1, value=0,
                           marks={i: {"label": str(i),
                                      "style": {"color": "#555", "fontSize": "10px"}}
                                  for i in range(len(LEVELS))},
                           tooltip={"placement": "bottom"},
                           updatemode="drag"),
            ], width=3),
        ]),

        # main graph
        dbc.Row(style={"margin": "0"}, children=[
            dbc.Col(
                dcc.Graph(
                    id="main-graph",
                    figure=initial_fig,
                    style={"height": "calc(100vh - 70px)"},
                    config={"scrollZoom": True,
                            "displayModeBar": True,
                            "modeBarButtonsToRemove": ["select2d", "lasso2d"],
                            "toImageButtonOptions": {"format": "svg", "filename": tissue_name}},
                ), width=12
            ),
        ]),

        # hidden store for hierarchy state
        dcc.Store(id="zoom-store", data={"x_range": None, "y_range": None,
                                          "level": 0, "auto": True}),
        dcc.Store(id="hierarchy-store", data={"tissue": tissue_name,
                                               "threshold": threshold}),
    ])

    # ── callbacks ─────────────────────────────────────────────────────────

    @app.callback(
        Output("zoom-store",   "data"),
        Output("level-label",  "children"),
        Output("level-slider", "value"),
        Input("main-graph",    "relayoutData"),
        Input("level-slider",  "value"),
        State("zoom-store",    "data"),
        prevent_initial_call=True,
    )
    def on_zoom_or_slider(relayout_data, slider_value, store):
        ctx         = callback_context
        trigger_id  = ctx.triggered[0]["prop_id"].split(".")[0]
        auto        = store.get("auto", True)

        x_range = store.get("x_range")
        y_range = store.get("y_range")

        # ── user dragged the manual slider ──────────────────────────────
        if trigger_id == "level-slider":
            store["level"] = slider_value
            store["auto"]  = False      # manual override — don't auto-switch
            return store, level_labels[slider_value], slider_value

        # ── user zoomed / panned the graph ──────────────────────────────
        if relayout_data:
            if "xaxis.range[0]" in relayout_data:
                x0 = relayout_data["xaxis.range[0]"]
                x1 = relayout_data["xaxis.range[1]"]
                y0 = relayout_data.get("yaxis.range[0]")
                y1 = relayout_data.get("yaxis.range[1]")
                x_range = [x0, x1]
                y_range = [y0, y1] if y0 is not None else y_range
                store["x_range"] = x_range
                store["y_range"] = y_range
                store["auto"]    = True  # zoom resets to auto mode

            elif "xaxis.autorange" in relayout_data or "autosize" in relayout_data:
                # full zoom-out / reset
                x_range = None
                y_range = None
                store["x_range"] = None
                store["y_range"] = None
                store["auto"]    = True

        # ── determine level from zoom span ──────────────────────────────
        if store.get("auto", True) and x_range is not None:
            x_span    = abs(x_range[1] - x_range[0])
            new_level = hierarchy.level_for_zoom(x_span)
        elif store.get("auto", True):
            new_level = 0
        else:
            new_level = store["level"]

        store["level"] = new_level
        return store, level_labels[new_level], new_level

    @app.callback(
        Output("main-graph",   "figure"),
        Input("zoom-store",    "data"),
        Input("reload-btn",    "n_clicks"),
        State("tissue-input",  "value"),
        State("threshold-input","value"),
    )
    def update_figure(store, _reload_clicks, tissue_input, threshold_input):
        ctx = callback_context
        if ctx.triggered and ctx.triggered[0]["prop_id"] == "reload-btn.n_clicks":
            nonlocal hierarchy
            t  = tissue_input or tissue_name
            th = float(threshold_input or threshold)
            hierarchy = GraphHierarchy(t, th)
            store     = {"x_range": None, "y_range": None, "level": 0, "auto": True}

        level_idx = store.get("level", 0)
        x_range   = store.get("x_range")
        y_range   = store.get("y_range")

        return hierarchy.build_figure(
            level_idx=level_idx,
            x_range=x_range,
            y_range=y_range,
        )

    return app


# ── entry point ───────────────────────────────────────────────────────────
if __name__ == "__main__":
    tissue_name = sys.argv[1] if len(sys.argv) > 1 else "Liver"
    threshold   = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
    app = build_app(tissue_name=tissue_name, threshold=threshold)
    app.run(debug=True, port=8050)