-- Leaf nodes (original graph)
CREATE TABLE nodes (
  id          TEXT PRIMARY KEY,
  x           FLOAT,
  y           FLOAT,
  label       TEXT,
  parent_id   TEXT REFERENCES cliques(id)  -- which clique absorbed this node
);

-- Cliques at every decomposition level
CREATE TABLE cliques (
  id           TEXT PRIMARY KEY,
  clique_type  TEXT CHECK (clique_type IN ('K3','K4')),
  level        INT,           -- 0 = leaves, higher = more decomposed
  centroid_x   FLOAT,
  centroid_y   FLOAT,
  parent_id    TEXT REFERENCES cliques(id),  -- null if top-level
  member_ids   TEXT[],        -- node or clique IDs contained within
  bbox         JSONB          -- {minX, minY, maxX, maxY} for viewport culling
);

-- All edges, tagged by what level they belong to
CREATE TABLE edges (
  id           BIGSERIAL PRIMARY KEY,
  source_id    TEXT,          -- node or clique id
  target_id    TEXT,
  level        INT,           -- which decomposition level this edge is visible at
  is_boundary  BOOL,          -- crosses a clique boundary (out_edge)?
  weight       FLOAT
);

-- Index heavily
CREATE INDEX idx_cliques_parent   ON cliques(parent_id);
CREATE INDEX idx_cliques_level    ON cliques(level);
CREATE INDEX idx_edges_source     ON edges(source_id);
CREATE INDEX idx_edges_level      ON edges(level);
```

---

## Decomposition Pipeline

Run this **offline** as a preprocessing step before the app is used:
```
Raw graph (level 0)
  → Find all K4 subgraphs → collapse → store as cliques
  → Find all K3 subgraphs → collapse → store as cliques
  → Compute out_edges (boundary edges) for each clique
  → Record new node_count at this level
  → Repeat until node_count stabilizes
  → Store max_level in a config table
```

Each iteration produces a new set of clique rows at `level = N+1`. The termination condition (node count unchanged) is your convergence check. Tools: **NetworkX** in Python handles K3/K4 subgraph detection natively.

---

## Backend API

**FastAPI** (Python) or **Express** (Node) — two key endpoints:
```
GET /graph/top
  → Returns cliques at max_level with their centroid x,y and out_edges
  → This is what the UI loads first

GET /graph/expand/:clique_id
  → Returns the direct children of this clique (one level down)
  → Returns their internal edges + boundary edges
  → Called on user click
```

Lazy loading is critical here — never send the full 20K node graph to the browser.

---

## Frontend

**Cytoscape.js** is better than D3 for this use case — it's built for graph interaction, handles large node counts better, and has a built-in compound node model that maps directly to your hierarchy.
```
State: Map<cliqueId, 'collapsed' | 'expanded'>

On click(cliqueId):
  fetch /graph/expand/:cliqueId
  replace that node with its children in the Cytoscape instance
  animate layout transition
```

Use a **viewport-aware renderer** — only render nodes within the current viewport bounds (use the `bbox` column for culling).

---

## Caching Layer

Add **Redis** between the API and database:

- Cache `expand/:clique_id` responses (they never change after preprocessing)
- TTL: indefinite or very long — this data is static
- The top-level graph view should be prewarmed on startup

---

## Full Stack Summary
```
Python (NetworkX)          PostgreSQL
preprocessing pipeline  →  nodes + cliques + edges tables
                                    ↓
                           FastAPI backend
                           + Redis cache
                                    ↓
                           Cytoscape.js frontend
                           (lazy expand on click)