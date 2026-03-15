from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional
import psycopg2
import psycopg2.extras
import os

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- DB Connection ---
def get_conn():
    return psycopg2.connect(
        host=os.getenv("DB_HOST", "localhost"),
        dbname=os.getenv("DB_NAME", "viz_db"),
        user=os.getenv("DB_USER", "viz_user"),
        password=os.getenv("DB_PASSWORD", "viz_password"),
        cursor_factory=psycopg2.extras.RealDictCursor
    )

# --- Response Models ---
class NodeOut(BaseModel):
    id: str
    x: float
    y: float
    label: Optional[str]
    is_clique: bool
    clique_type: Optional[str]
    member_count: Optional[int]

class EdgeOut(BaseModel):
    source_id: str
    target_id: str
    weight: float
    is_boundary: bool

class GraphResponse(BaseModel):
    nodes: list[NodeOut]
    edges: list[EdgeOut]


# --- Endpoints ---

@app.get("/graph/top", response_model=GraphResponse)
def get_top_graph():
    """
    Returns the most collapsed view of the graph — cliques at max level
    whose parent_id is NULL (i.e. have not been absorbed into anything).
    """
    conn = get_conn()
    cur = conn.cursor()

    try:
        # Get the top-level clique nodes (never absorbed)
        cur.execute("""
            SELECT id, centroid_x AS x, centroid_y AS y, clique_type,
                   array_length(member_ids, 1) AS member_count
            FROM cliques
            WHERE parent_id IS NULL
        """)
        clique_rows = cur.fetchall()

        if not clique_rows:
            raise HTTPException(status_code=404, detail="No top-level cliques found")

        clique_ids = [r["id"] for r in clique_rows]

        nodes = [
            NodeOut(
                id=r["id"],
                x=r["x"],
                y=r["y"],
                label=r["id"],
                is_clique=True,
                clique_type=r["clique_type"],
                member_count=r["member_count"]
            )
            for r in clique_rows
        ]

        # Get edges between top-level cliques only
        cur.execute("""
            SELECT source_id, target_id, weight, is_boundary
            FROM edges
            WHERE source_id = ANY(%s)
              AND target_id = ANY(%s)
        """, (clique_ids, clique_ids))
        edge_rows = cur.fetchall()

        edges = [
            EdgeOut(
                source_id=r["source_id"],
                target_id=r["target_id"],
                weight=r["weight"],
                is_boundary=r["is_boundary"]
            )
            for r in edge_rows
        ]

        return GraphResponse(nodes=nodes, edges=edges)

    finally:
        cur.close()
        conn.close()


@app.get("/graph/expand/{clique_id}", response_model=GraphResponse)
def expand_clique(clique_id: str):
    """
    Returns the direct children of a clique (one level down).
    Children can be either cliques or leaf nodes.
    Also returns:
      - internal edges between children
      - boundary edges from children to already-visible external nodes
    """
    conn = get_conn()
    cur = conn.cursor()

    try:
        # Fetch the clique and its member_ids
        cur.execute("""
            SELECT member_ids, level FROM cliques WHERE id = %s
        """, (clique_id,))
        row = cur.fetchone()

        if not row:
            raise HTTPException(status_code=404, detail=f"Clique {clique_id} not found")

        member_ids = row["member_ids"]
        level = row["level"]

        # Separate member_ids into cliques vs leaf nodes
        cur.execute("""
            SELECT id, centroid_x AS x, centroid_y AS y, clique_type,
                   array_length(member_ids, 1) AS member_count
            FROM cliques
            WHERE id = ANY(%s)
        """, (member_ids,))
        child_cliques = cur.fetchall()
        child_clique_ids = {r["id"] for r in child_cliques}

        cur.execute("""
            SELECT id, x, y, label
            FROM nodes
            WHERE id = ANY(%s)
        """, (member_ids,))
        child_nodes = cur.fetchall()

        nodes = [
            NodeOut(
                id=r["id"],
                x=r["x"],
                y=r["y"],
                label=r["id"],
                is_clique=True,
                clique_type=r["clique_type"],
                member_count=r["member_count"]
            )
            for r in child_cliques
        ] + [
            NodeOut(
                id=r["id"],
                x=r["x"],
                y=r["y"],
                label=r["label"] or r["id"],
                is_clique=False,
                clique_type=None,
                member_count=None
            )
            for r in child_nodes
        ]

        # Internal edges — between members at this level
        cur.execute("""
            SELECT source_id, target_id, weight, is_boundary
            FROM edges
            WHERE level = %s
              AND source_id = ANY(%s)
              AND target_id = ANY(%s)
        """, (level, member_ids, member_ids))
        internal_edges = cur.fetchall()

        # Boundary edges — from members outward to already-visible nodes
        cur.execute("""
            SELECT source_id, target_id, weight, is_boundary
            FROM edges
            WHERE level = %s
              AND is_boundary = TRUE
              AND source_id = ANY(%s)
        """, (level, member_ids))
        boundary_edges = cur.fetchall()

        edges = [
            EdgeOut(
                source_id=r["source_id"],
                target_id=r["target_id"],
                weight=r["weight"],
                is_boundary=r["is_boundary"]
            )
            for r in internal_edges + boundary_edges
        ]

        return GraphResponse(nodes=nodes, edges=edges)

    finally:
        cur.close()
        conn.close()