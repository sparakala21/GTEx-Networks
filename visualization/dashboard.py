"""
Gene Co-expression Network — Multi-Mode Hierarchy Viewer
=========================================================
Visualises JSON produced by either:
    • export_hierarchy.py          (k3-clique / supernode coarsening)
    • export_hierarchy_kcore.py    (k-core decomposition)

The viewer auto-detects the decomposition mode from meta.decomposition
and adapts labels, colourbars, and sidebar content accordingly.

Features
--------
  • Level selector — browse all hierarchy levels
  • Mode-aware colouring (clique weight vs coreness)
  • k-core: highlights newly appearing nodes at each level
  • k-core: "__new__" node count badge per level transition
  • Click a node to inspect members / gene list
  • Top-N label slider
  • PNG export

Usage
-----
    python visualize_hierarchy.py [graph_export.json] [--port 8050] [--debug]

Dependencies
------------
    pip install dash plotly dash-bootstrap-components
"""

import sys
import json
import math
import argparse
import numpy as np
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output, State, callback_context, no_update
import dash_bootstrap_components as dbc


# ══════════════════════════════════════════════════════════════════════════
# Constants
# ══════════════════════════════════════════════════════════════════════════

MODE_KCORE  = "kcore"
MODE_CLIQUE = "clique"   # default / original export (no decomposition key)

# Colour palettes per mode
PALETTE = {
    MODE_KCORE:  dict(
        node_scale="Plasma",
        edge_color="rgba(130,100,220,0.18)",
        accent="#a07cf5",
        new_node_border="#f5c842",
        new_node_border_width=2.5,
        colorbar_title="Coreness (k)",
    ),
    MODE_CLIQUE: dict(
        node_scale="Viridis",
        edge_color="rgba(120,190,230,0.18)",
        accent="#7ad4cc",
        new_node_border="#ffffff",
        new_node_border_width=0.5,
        colorbar_title="Clique Weight",
    ),
}


# ══════════════════════════════════════════════════════════════════════════
# Load & detect mode
# ══════════════════════════════════════════════════════════════════════════

def load_hierarchy(path: str) -> dict:
    with open(path) as f:
        data = json.load(f)

    meta = data["meta"]
    mode = meta.get("decomposition", MODE_CLIQUE)
    if mode not in (MODE_KCORE, MODE_CLIQUE):
        mode = MODE_CLIQUE   # graceful fallback

    # Index levels by level_idx for O(1) lookup
    level_index = {lvl["level_idx"]: lvl for lvl in data["levels"]}

    # Pre-extract "new nodes at level N" for k-core mode
    new_at_level: dict[int, set] = {}
    if mode == MODE_KCORE and "children_map" in data:
        for li_str, mapping in data["children_map"].items():
            li = int(li_str)
            next_li = li + 1
            new_nodes = {d["node"] for d in mapping.get("__new__", [])}
            if new_nodes:
                new_at_level[next_li] = new_nodes

    return {
        "meta":         meta,
        "mode":         mode,
        "levels":       data["levels"],
        "level_index":  level_index,
        "new_at_level": new_at_level,
        "k_values":     meta.get("k_values", []),  # empty for clique mode
        "n_levels":     meta["n_levels"],
    }


# ══════════════════════════════════════════════════════════════════════════
# Position normalisation
# ══════════════════════════════════════════════════════════════════════════

def normalise_positions(nodes: dict) -> dict:
    xs = [n["x"] for n in nodes.values()]
    ys = [n["y"] for n in nodes.values()]
    if not xs:
        return nodes
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    span = max(max(xs) - min(xs), max(ys) - min(ys)) or 1.0
    return {
        nid: {**node,
              "nx": (node["x"] - cx) / span * 2,
              "ny": (node["y"] - cy) / span * 2}
        for nid, node in nodes.items()
    }


# ══════════════════════════════════════════════════════════════════════════
# Figure builder
# ══════════════════════════════════════════════════════════════════════════

def _degree_map(level: dict) -> dict:
    deg = {nid: 0 for nid in level["nodes"]}
    for e in level.get("edges", []):
        if e["source"] in deg: deg[e["source"]] += 1
        if e["target"] in deg: deg[e["target"]] += 1
    return deg


def make_figure(
    hierarchy:      dict,
    level_idx:      int,
    highlight_id:   str  = None,
    label_top_n:    int  = 10,
) -> go.Figure:

    level  = hierarchy["level_index"][level_idx]
    meta   = hierarchy["meta"]
    mode   = hierarchy["mode"]
    pal    = PALETTE[mode]
    nodes  = normalise_positions(level["nodes"])
    edges  = level.get("edges", [])
    new_at = hierarchy["new_at_level"].get(level_idx, set())
    k_vals = hierarchy["k_values"]

    # ── top-N label selection ─────────────────────────────────────────
    degree      = _degree_map(level)
    label_top_n = min(label_top_n, len(nodes))
    ranked      = sorted(degree.items(), key=lambda x: x[1], reverse=True)
    central_ids = {nid for nid, _ in ranked[:label_top_n]}

    # ── edge trace ────────────────────────────────────────────────────
    edge_x, edge_y = [], []
    for e in edges:
        src, tgt = e["source"], e["target"]
        if src not in nodes or tgt not in nodes:
            continue
        edge_x += [nodes[src]["nx"], nodes[tgt]["nx"], None]
        edge_y += [nodes[src]["ny"], nodes[tgt]["ny"], None]

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        mode="lines",
        line=dict(color=pal["edge_color"], width=1),
        hoverinfo="none",
        showlegend=False,
    )

    # ── node arrays ───────────────────────────────────────────────────
    nids       = list(nodes.keys())
    nx_vals    = [nodes[n]["nx"] for n in nids]
    ny_vals    = [nodes[n]["ny"] for n in nids]
    cq_vals    = [nodes[n]["clique_weight"] for n in nids]   # coreness in kcore mode
    sizes_raw  = [nodes[n]["size"] for n in nids]
    max_size   = max(sizes_raw) or 1

    def node_px(s):
        return 8 + math.sqrt(s / max_size) * 38

    node_sizes = [node_px(s) for s in sizes_raw]

    # ── hover text ────────────────────────────────────────────────────
    def _hover(i):
        n    = nodes[nids[i]]
        lbl  = n.get("label", nids[i])
        mems = n.get("members", [])
        is_new = nids[i] in new_at

        if mode == MODE_KCORE:
            k_here = n.get("k", "?")
            core_line = f"Coreness: {cq_vals[i]:.0f}"
            new_badge = "  🆕 NEW at this level<br>" if is_new else ""
            return (
                f"<b>{lbl}</b><br>"
                f"{new_badge}"
                f"k (this level): {k_here}<br>"
                f"{core_line}<br>"
                f"Degree (within core): {n.get('degree', '?')}"
            )
        else:
            mems_str = ", ".join(mems[:10]) + ("…" if len(mems) > 10 else "")
            return (
                f"<b>{lbl}</b><br>"
                f"Genes collapsed: {n['size']}<br>"
                f"Clique weight: {cq_vals[i]:.4f}<br>"
                f"Members: {mems_str}"
            )

    hover_texts = [_hover(i) for i in range(len(nids))]

    # ── per-node border (highlight selected; gold ring for new nodes) ──
    border_colors = []
    border_widths = []
    for nid in nids:
        if highlight_id and nid == highlight_id:
            border_colors.append("#ffffff")
            border_widths.append(3)
        elif nid in new_at:
            border_colors.append(pal["new_node_border"])
            border_widths.append(pal["new_node_border_width"])
        else:
            border_colors.append("rgba(255,255,255,0.15)")
            border_widths.append(0.5)

    # ── text labels ───────────────────────────────────────────────────
    text_labels = [
        nodes[n].get("label", n) if n in central_ids else ""
        for n in nids
    ]

    node_trace = go.Scatter(
        x=nx_vals, y=ny_vals,
        mode="markers+text",
        marker=dict(
            size=node_sizes,
            color=cq_vals,
            colorscale=pal["node_scale"],
            showscale=True,
            colorbar=dict(
                title=dict(text=pal["colorbar_title"], side="right"),
                thickness=14,
                len=0.55,
                bgcolor="rgba(18,18,28,0.75)",
                bordercolor="rgba(255,255,255,0.1)",
                tickfont=dict(color="#aacccc", size=10),
            ),
            line=dict(color=border_colors, width=border_widths),
            sizemode="diameter",
            opacity=0.87,
        ),
        text=text_labels,
        textposition="top center",
        textfont=dict(size=8, color="rgba(200,230,230,0.75)", family="'IBM Plex Mono', monospace"),
        hovertext=hover_texts,
        hoverinfo="text",
        customdata=nids,
    )

    # ── title ─────────────────────────────────────────────────────────
    tissue    = meta.get("tissue", "?")
    threshold = meta.get("threshold", "?")
    lvl_label = level["label"]

    if mode == MODE_KCORE and k_vals and level_idx < len(k_vals):
        k_tag = f"k = {k_vals[level_idx]}"
        new_count = len(new_at)
        new_tag = (
            f"  <span style='color:#f5c842'>+{new_count} new gene{'s' if new_count!=1 else ''}</span>"
            if new_count else ""
        )
    else:
        k_tag, new_tag = "", ""

    subtitle = (
        f"{lvl_label}  {k_tag}  ·  "
        f"{len(nodes)} nodes, {len(edges)} edges"
        + new_tag
    )

    mode_badge = (
        "k-core decomposition"
        if mode == MODE_KCORE
        else "k₃-clique coarsening"
    )

    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        title=dict(
            text=(
                f"<b>Gene Co-expression</b>  ·  "
                f"<span style='color:{pal['accent']}'>{tissue}</span>  "
                f"· thr {threshold}  · {mode_badge}<br>"
                f"<span style='font-size:12px;color:#7aadad'>{subtitle}</span>"
            ),
            font=dict(color="#c8e8e8", size=15, family="'IBM Plex Mono', monospace"),
            x=0.5, xanchor="center",
        ),
        paper_bgcolor="#0b0f1a",
        plot_bgcolor="#0b0f1a",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, scaleanchor="x"),
        hovermode="closest",
        dragmode="pan",
        margin=dict(l=20, r=20, t=90, b=20),
        font=dict(color="#c8e8e8"),
        height=760,
        transition=dict(duration=350, easing="cubic-in-out"),
    )
    return fig


# ══════════════════════════════════════════════════════════════════════════
# Sidebar helpers
# ══════════════════════════════════════════════════════════════════════════

SIDEBAR_STYLE = {
    "padding": "20px 16px",
    "height": "100vh",
    "overflowY": "auto",
    "borderLeft": "1px solid rgba(120,180,180,0.15)",
    "backgroundColor": "#0b0f1a",
    "fontFamily": "'IBM Plex Mono', monospace",
}

HEADING_STYLE = {
    "color": "#7ad4cc",
    "borderBottom": "1px solid rgba(120,180,180,0.25)",
    "paddingBottom": "6px",
    "marginBottom": "10px",
    "fontSize": "14px",
    "letterSpacing": "0.06em",
    "textTransform": "uppercase",
}


def _row(label, value, accent=False):
    colour = "#f5c842" if accent else "#9abebe"
    return html.P(
        [html.B(f"{label}: "), str(value)],
        style={"color": colour, "fontSize": "12px", "margin": "3px 0"},
    )


def _hr():
    return html.Hr(style={"borderColor": "rgba(120,180,180,0.15)", "margin": "14px 0"})


def _placeholder(mode):
    tip = (
        "Click a node for coreness & degree stats."
        if mode == MODE_KCORE
        else "Click a node to inspect its collapsed genes."
    )
    return html.P(tip, style={"color": "#3a6060", "fontStyle": "italic", "fontSize": "12px"})


def _node_detail(nid: str, node: dict, mode: str):
    members = node.get("members", [])
    if mode == MODE_KCORE:
        return [
            _row("Gene",      nid),
            _row("Coreness",  f"{node['clique_weight']:.0f}", accent=True),
            _row("k (level)", node.get("k", "?")),
            _row("Degree",    node.get("degree", "?")),
        ]
    else:
        rows = [
            _row("Node ID",       nid),
            _row("Label",         node.get("label", nid)),
            _row("Genes inside",  node["size"]),
            _row("Clique weight", f"{node['clique_weight']:.6f}", accent=True),
            html.P(html.B(f"Members ({len(members)}):"),
                   style={"color": "#9abebe", "fontSize": "12px", "margin": "6px 0 2px"}),
            html.Div(
                ", ".join(members),
                style={
                    "maxHeight": "260px", "overflowY": "auto",
                    "backgroundColor": "rgba(255,255,255,0.03)",
                    "padding": "8px", "borderRadius": "4px",
                    "fontSize": "10px", "lineHeight": "1.7",
                    "color": "#88cccc", "wordBreak": "break-word",
                }
            ),
        ]
        return rows


def _level_options(hierarchy: dict) -> list:
    mode   = hierarchy["mode"]
    k_vals = hierarchy["k_values"]
    opts   = []
    for lvl in hierarchy["levels"]:
        li  = lvl["level_idx"]
        lbl = lvl["label"]
        n   = len(lvl["nodes"])
        if mode == MODE_KCORE and k_vals and li < len(k_vals):
            k = k_vals[li]
            opts.append({"label": f"[{li}] {lbl}  ({n} genes)", "value": li})
        else:
            opts.append({"label": f"[{li}] {lbl}  ({n} nodes)", "value": li})
    return opts


# ══════════════════════════════════════════════════════════════════════════
# App builder
# ══════════════════════════════════════════════════════════════════════════

def build_app(path: str) -> Dash:
    hierarchy = load_hierarchy(path)
    meta      = hierarchy["meta"]
    mode      = hierarchy["mode"]
    pal       = PALETTE[mode]
    n_levels  = hierarchy["n_levels"]

    app = Dash(
        __name__,
        external_stylesheets=[
            dbc.themes.CYBORG,
            "https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&display=swap",
        ],
        title="Co-expression Hierarchy",
    )

    level_opts  = _level_options(hierarchy)
    max_nodes_l = max(len(lvl["nodes"]) for lvl in hierarchy["levels"])

    # ── Mode badge ────────────────────────────────────────────────────
    mode_label = (
        "k-core Decomposition"
        if mode == MODE_KCORE
        else "k₃-Clique Coarsening"
    )
    mode_colour = pal["accent"]

    app.layout = dbc.Container(
        fluid=True,
        style={"backgroundColor": "#0b0f1a", "minHeight": "100vh", "padding": 0},
        children=[
            dcc.Store(id="selected-node"),
            dcc.Store(id="current-level", data=0),

            # ── Top bar ───────────────────────────────────────────────
            dbc.Row([
                dbc.Col(
                    html.Div(
                        style={
                            "display": "flex", "alignItems": "center",
                            "gap": "24px", "padding": "10px 20px",
                            "borderBottom": f"1px solid {mode_colour}22",
                            "backgroundColor": "#090d17",
                        },
                        children=[
                            html.Span(
                                mode_label,
                                style={
                                    "color": mode_colour,
                                    "fontFamily": "'IBM Plex Mono', monospace",
                                    "fontWeight": "600",
                                    "fontSize": "13px",
                                    "letterSpacing": "0.08em",
                                    "textTransform": "uppercase",
                                    "borderRight": f"1px solid {mode_colour}44",
                                    "paddingRight": "24px",
                                }
                            ),
                            html.Div([
                                html.Span("Level: ",
                                          style={"color": "#4a7070", "fontSize": "12px"}),
                                dcc.Dropdown(
                                    id="level-select",
                                    options=level_opts,
                                    value=0,
                                    clearable=False,
                                    style={
                                        "width": "380px",
                                        "backgroundColor": "#12192e",
                                        "color": "#c8e8e8",
                                        "border": f"1px solid {mode_colour}44",
                                        "fontFamily": "'IBM Plex Mono', monospace",
                                        "fontSize": "12px",
                                    },
                                ),
                            ], style={"display": "flex", "alignItems": "center", "gap": "8px"}),

                            # Prev / Next buttons
                            dbc.ButtonGroup([
                                dbc.Button("◀", id="btn-prev", n_clicks=0, size="sm",
                                           style={"backgroundColor": "#12192e",
                                                  "borderColor": f"{mode_colour}44",
                                                  "color": mode_colour, "fontFamily": "monospace"}),
                                dbc.Button("▶", id="btn-next", n_clicks=0, size="sm",
                                           style={"backgroundColor": "#12192e",
                                                  "borderColor": f"{mode_colour}44",
                                                  "color": mode_colour, "fontFamily": "monospace"}),
                            ]),

                            # New-node badge (k-core only)
                            html.Div(id="new-node-badge",
                                     style={"marginLeft": "auto", "color": "#f5c842",
                                            "fontSize": "12px",
                                            "fontFamily": "'IBM Plex Mono', monospace"}),
                        ]
                    ),
                    width=12,
                )
            ], style={"margin": 0}),

            # ── Main content ──────────────────────────────────────────
            dbc.Row([
                dbc.Col([
                    dcc.Graph(
                        id="network-graph",
                        figure=make_figure(hierarchy, 0),
                        config={
                            "scrollZoom": True,
                            "displayModeBar": True,
                            "modeBarButtonsToRemove": ["lasso2d", "select2d"],
                            "toImageButtonOptions": {
                                "format": "png",
                                "filename": f"coexpression_{meta.get('tissue', 'network')}",
                                "scale": 2,
                            },
                        },
                        style={"height": "calc(100vh - 52px)"},
                    )
                ], width=8, style={"padding": 0}),

                dbc.Col([
                    html.Div(
                        style=SIDEBAR_STYLE,
                        children=[

                            # Network info
                            html.H6("Network Info", style=HEADING_STYLE),
                            _row("Tissue",    meta.get("tissue", "—")),
                            _row("Threshold", meta.get("threshold", "—")),
                            _row("Levels",    n_levels),
                            *(
                                [_row("k range",
                                      f"{hierarchy['k_values'][-1]} – {hierarchy['k_values'][0]}"
                                      if hierarchy["k_values"] else "—")]
                                if mode == MODE_KCORE
                                else [_row("Mode", "k₃-clique coarsening")]
                            ),
                            _hr(),

                            # Current level info
                            html.H6("Current Level", style=HEADING_STYLE),
                            html.Div(id="level-info",
                                     style={"fontSize": "12px", "color": "#9abebe"}),
                            _hr(),

                            # Selected node
                            html.H6("Selected Node", style=HEADING_STYLE),
                            html.Div(id="node-detail",
                                     style={"fontSize": "12px", "color": "#9abebe"},
                                     children=_placeholder(mode)),
                            _hr(),

                            # Label settings
                            html.H6("Label Settings", style=HEADING_STYLE),
                            html.P(
                                "Label top-N nodes by degree:",
                                style={"color": "#6a8e8e", "fontSize": "11px",
                                       "margin": "0 0 4px"}
                            ),
                            dcc.Slider(
                                id="label-top-n",
                                min=0,
                                max=min(60, max_nodes_l),
                                step=1,
                                value=10,
                                marks={
                                    0:  {"label": "0",  "style": {"color": "#4a7070"}},
                                    10: {"label": "10", "style": {"color": "#4a7070"}},
                                    30: {"label": "30", "style": {"color": "#4a7070"}},
                                    60: {"label": "60", "style": {"color": "#4a7070"}},
                                },
                                tooltip={"placement": "bottom", "always_visible": True},
                            ),
                            _hr(),

                            # Legend for k-core new-node highlighting
                            *(
                                [
                                    html.H6("Legend", style=HEADING_STYLE),
                                    html.Div([
                                        html.Span("●", style={"color": "#f5c842",
                                                              "fontSize": "16px",
                                                              "marginRight": "6px"}),
                                        html.Span("Newly added at this level",
                                                  style={"fontSize": "11px",
                                                         "color": "#9abebe"}),
                                    ], style={"display": "flex", "alignItems": "center",
                                              "margin": "4px 0"}),
                                    _hr(),
                                ]
                                if mode == MODE_KCORE
                                else []
                            ),

                            # Tips
                            html.H6("Tips", style=HEADING_STYLE),
                            html.Ul([
                                html.Li("Scroll to zoom"),
                                html.Li("Drag to pan"),
                                html.Li("Click node for details"),
                                html.Li("▶ / ◀ to step levels"),
                            ], style={"color": "#3a6060", "fontSize": "11px",
                                      "lineHeight": "1.9", "paddingLeft": "16px"}),
                        ]
                    )
                ], width=4, style={"padding": 0}),
            ], style={"margin": 0}),
        ]
    )

    # ── Callbacks ─────────────────────────────────────────────────────────

    # Sync level via prev/next buttons → level-select dropdown
    @app.callback(
        Output("level-select", "value"),
        Input("btn-prev", "n_clicks"),
        Input("btn-next", "n_clicks"),
        State("level-select", "value"),
        prevent_initial_call=True,
    )
    def step_level(prev_clicks, next_clicks, current_level):
        triggered = callback_context.triggered[0]["prop_id"]
        if "btn-prev" in triggered:
            return max(0, current_level - 1)
        if "btn-next" in triggered:
            return min(n_levels - 1, current_level + 1)
        return current_level

    # Main callback: level change + node click + label slider
    @app.callback(
        Output("network-graph", "figure"),
        Output("node-detail",   "children"),
        Output("level-info",    "children"),
        Output("new-node-badge","children"),
        Output("selected-node", "data"),
        Input("level-select",   "value"),
        Input("network-graph",  "clickData"),
        Input("label-top-n",    "value"),
        State("selected-node",  "data"),
    )
    def update(level_idx, click_data, label_top_n, current_sel):
        ctx       = callback_context
        triggered = ctx.triggered[0]["prop_id"] if ctx.triggered else ""
        level_idx = level_idx or 0
        label_top_n = label_top_n or 0

        level = hierarchy["level_index"][level_idx]
        new_at = hierarchy["new_at_level"].get(level_idx, set())

        # ── current level info sidebar ────────────────────────────────
        n_nodes = len(level["nodes"])
        n_edges = len(level.get("edges", []))
        if mode == MODE_KCORE and hierarchy["k_values"]:
            k = hierarchy["k_values"][level_idx] if level_idx < len(hierarchy["k_values"]) else "?"
            lvl_info_rows = [
                _row("Label",    level["label"]),
                _row("k value",  k, accent=True),
                _row("Nodes",    n_nodes),
                _row("Edges",    n_edges),
                _row("New here", len(new_at), accent=bool(new_at)),
            ]
        else:
            lvl_info_rows = [
                _row("Label",  level["label"]),
                _row("Nodes",  n_nodes),
                _row("Edges",  n_edges),
            ]

        # ── new-node badge ────────────────────────────────────────────
        badge = (
            f"+ {len(new_at)} new gene{'s' if len(new_at)!=1 else ''} at this level"
            if new_at else ""
        )

        # ── handle node click ─────────────────────────────────────────
        selected = current_sel
        detail   = _placeholder(mode)

        if "clickData" in triggered and click_data:
            nid = click_data["points"][0].get("customdata")
            if nid == current_sel:
                selected = None
                detail   = _placeholder(mode)
            elif nid and nid in level["nodes"]:
                selected = nid
                detail   = _node_detail(nid, level["nodes"][nid], mode)
            # If clicked node not in new level after level switch, clear
        elif "level-select" in triggered:
            # Reset selection on level change
            selected = None
            detail   = _placeholder(mode)

        # Validate selection still in current level
        if selected and selected not in level["nodes"]:
            selected = None
            detail   = _placeholder(mode)

        fig = make_figure(
            hierarchy, level_idx,
            highlight_id=selected,
            label_top_n=label_top_n,
        )

        return fig, detail, lvl_info_rows, badge, selected

    return app


# ══════════════════════════════════════════════════════════════════════════
# Entry point
# ══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Visualise k-core or k3-clique hierarchy JSON in Dash."
    )
    parser.add_argument(
        "json_path", nargs="?", default="graph_export_kcore.json",
        help="Path to hierarchy JSON (default: graph_export.json)"
    )
    parser.add_argument("--port",  type=int,          default=8050)
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    print(f"[viz] Loading: {args.json_path}")
    app = build_app(args.json_path)
    detected = app.server.config.get("mode", "")
    print(f"[viz] Detected mode: {detected or '(auto from JSON)'}")
    print(f"[viz] Serving on http://localhost:{args.port}")
    app.run(debug=args.debug, port=args.port)