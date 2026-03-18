from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional
import psycopg2
import psycopg2.extras
import os
import traceback
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

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
        os.environ["DATABASE_URL"],
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
    expression: Optional[float] = None

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
    conn = None
    cur = None
    try:
        logger.debug("Connecting to database...")
        conn = get_conn()
        cur = conn.cursor()
        logger.debug("Connected. Fetching top-level cliques...")

        cur.execute("""
            SELECT id, centroid_x AS x, centroid_y AS y, clique_type,
                   array_length(member_ids, 1) AS member_count,
                   (SELECT COALESCE(AVG(COALESCE(expression, 0)), 0) 
                    FROM nodes WHERE id = ANY(cliques.member_ids)) AS avg_expression
            FROM cliques
            WHERE parent_id IS NULL
        """)
        clique_rows = cur.fetchall()
        logger.debug(f"environment variables {os.environ}")
        logger.debug(f"Fetched {len(clique_rows)} clique rows")

        if not clique_rows:
            raise HTTPException(status_code=404, detail="No top-level cliques found")

        clique_ids = [r["id"] for r in clique_rows]
        logger.debug(f"Clique IDs: {clique_ids[:5]}...")  # first 5 only

        nodes = [
            NodeOut(
                id=r["id"], x=r["x"], y=r["y"], label=r["id"],
                is_clique=True, clique_type=r["clique_type"],
                member_count=r["member_count"], expression=r["avg_expression"]
            )
            for r in clique_rows
        ]

        logger.debug("Fetching edges...")
        cur.execute("""
            SELECT source_id, target_id, weight, is_boundary
            FROM edges
            WHERE source_id = ANY(%s)
              AND target_id = ANY(%s)
        """, (clique_ids, clique_ids))
        edge_rows = cur.fetchall()
        logger.debug(f"Fetched {len(edge_rows)} edges")

        edges = [
            EdgeOut(
                source_id=r["source_id"], target_id=r["target_id"],
                weight=r["weight"], is_boundary=r["is_boundary"]
            )
            for r in edge_rows
        ]

        return GraphResponse(nodes=nodes, edges=edges)

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error in get_top_graph: {e}")
        logger.error(traceback.format_exc())
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        if cur: cur.close()
        if conn: conn.close()

@app.get("/graph/expand/{clique_id}", response_model=GraphResponse)
def expand_clique(clique_id: str):
    conn = get_conn()
    cur = conn.cursor()

    try:
        # 1. Fetch the clique info
        cur.execute("SELECT member_ids, level FROM cliques WHERE id = %s", (clique_id,))
        row = cur.fetchone()

        if not row:
            raise HTTPException(status_code=404, detail=f"Clique {clique_id} not found")

        member_ids = row["member_ids"]
        level = row["level"]

        # 2. Fetch child cliques (Calculates average of THEIR members)
        cur.execute("""
            SELECT id, centroid_x AS x, centroid_y AS y, clique_type,
                   array_length(member_ids, 1) AS member_count,
                   (SELECT COALESCE(AVG(COALESCE(expression, 0)), 0) 
                    FROM nodes WHERE id = ANY(cliques.member_ids)) AS avg_expression
            FROM cliques
            WHERE id = ANY(%s)
        """, (member_ids,))
        child_cliques = cur.fetchall()

        # 3. Fetch leaf nodes (Directly gets their expression value)
        cur.execute("""
            SELECT id, x, y, label, 
                   COALESCE(expression, 0) AS expression  -- Added this column
            FROM nodes
            WHERE id = ANY(%s)
        """, (member_ids,))
        child_nodes = cur.fetchall()

        # 4. Map to Pydantic objects
        nodes = [
            NodeOut(
                id=r["id"], x=r["x"], y=r["y"], label=r["id"],
                is_clique=True, clique_type=r["clique_type"],
                member_count=r["member_count"],
                expression=r["avg_expression"]
            )
            for r in child_cliques
        ] + [
            NodeOut(
                id=r["id"], x=r["x"], y=r["y"], label=r["label"] or r["id"],
                is_clique=False, clique_type=None, member_count=None,
                expression=r["expression"] # Use "expression", NOT "avg_expression"
            )
            for r in child_nodes
        ]

        # ... (Edge logic remains the same)
        # Internal edges
        cur.execute("""
            SELECT source_id, target_id, weight, is_boundary
            FROM edges
            WHERE level = %s AND source_id = ANY(%s) AND target_id = ANY(%s)
        """, (level, member_ids, member_ids))
        internal_edges = cur.fetchall()

        # Boundary edges
        cur.execute("""
            SELECT source_id, target_id, weight, is_boundary
            FROM edges
            WHERE level = %s AND is_boundary = TRUE AND source_id = ANY(%s)
        """, (level, member_ids))
        boundary_edges = cur.fetchall()

        edges = [
            EdgeOut(
                source_id=r["source_id"], target_id=r["target_id"],
                weight=r["weight"], is_boundary=r["is_boundary"]
            )
            for r in internal_edges + boundary_edges
        ]

        return GraphResponse(nodes=nodes, edges=edges)

    finally:
        cur.close()
        conn.close()

@app.get("/graph/parent/{node_id}")
def get_parent_clique(node_id: str):
    conn = get_conn()
    cur = conn.cursor()
    try:
        # Find which clique contains this node as a member
        cur.execute("""
            SELECT id, centroid_x, centroid_y, clique_type, member_ids,
                   (SELECT AVG(expression) FROM nodes WHERE id = ANY(member_ids)) as avg_expr
            FROM cliques 
            WHERE %s = ANY(member_ids) 
            LIMIT 1
        """, (node_id,))
        r = cur.fetchone()
        
        if not r:
            raise HTTPException(status_code=404, detail="Parent not found")

        return {
            "clique": {
                "id": r["id"], "x": r["centroid_x"], "y": r["centroid_y"],
                "clique_type": r["clique_type"], "member_count": len(r["member_ids"]),
                "expression": r["avg_expr"]
            },
            "member_ids": r["member_ids"]
        }
    finally:
        cur.close()
        conn.close()
#note to self: python3 -m uvicorn api:app --reload to run