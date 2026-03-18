const API = 'http://localhost:8000';

let cy;
let expanding = false;

function showError(msg) {
  const bar = document.getElementById('error-bar');
  bar.textContent = msg;
  bar.style.display = 'block';
}

function showToast(msg) {
  const t = document.getElementById('expand-toast');
  t.textContent = msg;
  t.classList.add('visible');
  setTimeout(() => t.classList.remove('visible'), 1200);
}

function updateStats() {
  document.getElementById('stat-nodes').textContent = cy.nodes().length;
  document.getElementById('stat-edges').textContent = cy.edges().length;
}

function buildElements(data) {
  const nodes = data.nodes.map(n => ({
    data: {
      id: n.id,
      label: n.label || n.id,
      x: n.x,
      y: n.y,
      is_clique: n.is_clique,
      clique_type: n.clique_type,
      member_count: n.member_count,
      expression: n.expression,
    },
    position: { x: n.x * 1000, y: n.y * 1000 }
  }));

  const edges = data.edges.map(e => ({
    data: {
      id: `${e.source_id}__${e.target_id}`,
      source: e.source_id,
      target: e.target_id,
      weight: e.weight,
      is_boundary: e.is_boundary,
    }
  }));

  return [...nodes, ...edges];
}

function clearSearchHighlights() {
  cy.elements().removeClass('search-hit search-clique-hit search-dim');
}

async function searchGene(query) {
  query = query.trim();
  if (!query) {
    clearSearchHighlights();
    document.getElementById('search-status').textContent = '';
    return;
  }

  const q = query.toLowerCase();

  const directHit = cy.nodes().filter(n => {
    return n.id().toLowerCase() === q || (n.data('label') || '').toLowerCase() === q;
  });

  if (directHit.length) {
    clearSearchHighlights();
    cy.elements().addClass('search-dim');
    directHit.removeClass('search-dim').addClass('search-hit');
    directHit.connectedEdges().removeClass('search-dim');
    cy.animate({ fit: { eles: directHit, padding: 120 }, duration: 400 });
    renderDetail(directHit.first());
    document.getElementById('search-status').textContent =
      `found: ${directHit.first().id()}`;
    return;
  }

  document.getElementById('search-status').textContent = 'searching…';
  try {
    const res = await fetch(`${API}/graph/parent/${encodeURIComponent(query)}`);
    if (!res.ok) {
      document.getElementById('search-status').textContent = 'not found';
      clearSearchHighlights();
      return;
    }
    const data = await res.json();
    const cliqueId = data.clique?.id;
    if (!cliqueId) {
      document.getElementById('search-status').textContent = 'not found';
      clearSearchHighlights();
      return;
    }

    const cliqueNode = cy.getElementById(cliqueId);
    if (!cliqueNode.length) {
      document.getElementById('search-status').textContent =
        `inside clique: ${cliqueId.slice(0, 24)}… (not visible)`;
      clearSearchHighlights();
      return;
    }

    clearSearchHighlights();
    cy.elements().addClass('search-dim');
    cliqueNode.removeClass('search-dim').addClass('search-clique-hit');
    cliqueNode.connectedEdges().removeClass('search-dim');
    cy.animate({ center: { eles: cliqueNode }, zoom: Math.max(cy.zoom(), 2), duration: 400 });
    renderDetail(cliqueNode);
    document.getElementById('search-status').textContent =
      `"${query}" is inside this cluster`;
  } catch (err) {
    document.getElementById('search-status').textContent = `error: ${err.message}`;
  }
}

async function expandClique(cliqueNode) {
  if (expanding) return;
  expanding = true;

  const cliqueId = cliqueNode.id();
  const cliquePos = cliqueNode.position();

  cliqueNode.addClass('expanding');
  showToast('expanding…');

  try {
    const res = await fetch(`${API}/graph/expand/${encodeURIComponent(cliqueId)}`);
    if (!res.ok) throw new Error(`expand returned ${res.status}`);
    const data = await res.json();

    const memberIds = new Set(data.nodes.map(n => n.id));
    const externalEdges = [];

    cliqueNode.connectedEdges().forEach(e => {
      const srcId = e.data('source');
      const tgtId = e.data('target');
      const otherId = srcId === cliqueId ? tgtId : srcId;
      if (!memberIds.has(otherId)) {
        externalEdges.push({
          otherId,
          weight: e.data('weight'),
          is_boundary: e.data('is_boundary'),
        });
      }
    });

    cy.remove(cliqueNode);

    const newNodeEls = data.nodes
      .filter(n => !cy.getElementById(n.id).length)
      .map(n => ({
        data: {
          id: n.id,
          label: n.label || n.id,
          x: n.x,
          y: n.y,
          is_clique: n.is_clique,
          clique_type: n.clique_type,
          member_count: n.member_count,
          expression: n.expression,
        },
        position: { x: n.x * 1000, y: n.y * 1000 }
      }));

    cy.add(newNodeEls);

    const internalEdgeEls = data.edges
      .filter(e => !e.is_boundary)
      .filter(e => !cy.getElementById(`${e.source_id}__${e.target_id}`).length)
      .map(e => ({
        data: {
          id: `${e.source_id}__${e.target_id}`,
          source: e.source_id,
          target: e.target_id,
          weight: e.weight,
          is_boundary: false,
        }
      }));

    cy.add(internalEdgeEls);

    const boundaryEdgeEls = data.edges
      .filter(e => e.is_boundary)
      .filter(e => cy.getElementById(e.target_id).length)
      .filter(e => !cy.getElementById(`${e.source_id}__${e.target_id}`).length)
      .map(e => ({
        data: {
          id: `${e.source_id}__${e.target_id}`,
          source: e.source_id,
          target: e.target_id,
          weight: e.weight,
          is_boundary: true,
        }
      }));

    cy.add(boundaryEdgeEls);

    externalEdges.forEach(({ otherId, weight, is_boundary }) => {
      const otherNode = cy.getElementById(otherId);
      if (!otherNode.length) return;

      const otherPos = otherNode.position();
      let closest = null;
      let closestDist = Infinity;
      data.nodes.forEach(n => {
        const nEl = cy.getElementById(n.id);
        if (!nEl.length) return;
        const p = nEl.position();
        const d = Math.hypot(p.x - otherPos.x, p.y - otherPos.y);
        if (d < closestDist) { closestDist = d; closest = n.id; }
      });

      if (!closest) return;

      const edgeId = `${closest}__${otherId}`;
      const edgeIdRev = `${otherId}__${closest}`;
      if (cy.getElementById(edgeId).length || cy.getElementById(edgeIdRev).length) return;

      cy.add({ data: { id: edgeId, source: closest, target: otherId, weight, is_boundary } });
    });

    updateStats();
    showToast(`+${data.nodes.length} nodes`);

  } catch (err) {
    showError(`expand failed: ${err.message}`);
    cliqueNode.removeClass('expanding');
  } finally {
    expanding = false;
  }
}

async function repackClique(childNode) {
  if (expanding) return;

  const nodeId = childNode.id();
  showToast('repacking...');

  try {
    const res = await fetch(`${API}/graph/parent/${encodeURIComponent(nodeId)}`);
    if (!res.ok) throw new Error("Could not find parent clique");
    const data = await res.json();

    const parentClique = data.clique;
    const memberIds = data.member_ids;

    const nodesToRemove = cy.nodes().filter(n => memberIds.includes(n.id()));
    const memberIdsArray = nodesToRemove.map(n => n.id());
    const externalConnections = [];

    nodesToRemove.connectedEdges().forEach(e => {
      const src = e.source().id();
      const tgt = e.target().id();
      const srcIsMember = memberIdsArray.includes(src);
      const tgtIsMember = memberIdsArray.includes(tgt);

      if (srcIsMember !== tgtIsMember) {
        const other = srcIsMember ? tgt : src;
        externalConnections.push({
          otherId: other,
          weight: e.data('weight'),
          is_boundary: e.data('is_boundary')
        });
      }
    });

    cy.remove(nodesToRemove);

    cy.add({
      group: 'nodes',
      data: {
        id: parentClique.id,
        label: parentClique.label || parentClique.id,
        x: parentClique.x,
        y: parentClique.y,
        is_clique: true,
        clique_type: parentClique.clique_type,
        member_count: parentClique.member_count,
        expression: parentClique.expression
      },
      position: { x: parentClique.x * 1000, y: parentClique.y * 1000 }
    });

    externalConnections.forEach(conn => {
      const edgeId = `${parentClique.id}__${conn.otherId}`;
      if (!cy.getElementById(edgeId).length) {
        cy.add({
          group: 'edges',
          data: {
            id: edgeId,
            source: parentClique.id,
            target: conn.otherId,
            weight: conn.weight,
            is_boundary: conn.is_boundary
          }
        });
      }
    });

    updateStats();
    showToast('repacked');

  } catch (err) {
    showError(`repack failed: ${err.message}`);
  }
}

async function loadTop() {
  const res = await fetch(`${API}/graph/top`);
  if (!res.ok) throw new Error(`/graph/top returned ${res.status}`);
  return res.json();
}

function renderDetail(node) {
  const n = node.data();
  const detail = document.getElementById('detail-content');
  const formatExpr = (val) => (val !== null && val !== undefined) ? Number(val).toFixed(2) : '—';

  if (n.is_clique) {
    detail.innerHTML = `
      <div class="clique-badge">${n.clique_type || 'CLIQUE'}</div>
      <div class="row"><span>id</span><span class="val">${n.id.slice(0, 18)}</span></div>
      <div class="row"><span>members</span><span class="val">${n.member_count ?? '—'}</span></div>
      <div class="row"><span>avg expr</span><span class="val">${formatExpr(n.expression)}</span></div>
      <div class="row"><span>x</span><span class="val">${Number(n.x).toFixed(4)}</span></div>
      <div class="row"><span>y</span><span class="val">${Number(n.y).toFixed(4)}</span></div>
      <div class="row"><span>degree</span><span class="val">${node.degree()}</span></div>
      <div style="margin-top:8px; font-size:9px; color: var(--accent); letter-spacing:0.1em;">
        ▸ double-click to expand
      </div>
    `;
  } else {
    detail.innerHTML = `
      <div class="row"><span>id</span><span class="val">${n.id.slice(0, 18)}</span></div>
      <div class="row"><span>label</span><span class="val">${(n.label || '—').slice(0, 16)}</span></div>
      <div class="row"><span>expression</span><span class="val">${formatExpr(n.expression)}</span></div>
      <div class="row"><span>x</span><span class="val">${Number(n.x).toFixed(4)}</span></div>
      <div class="row"><span>y</span><span class="val">${Number(n.y).toFixed(4)}</span></div>
      <div class="row"><span>degree</span><span class="val">${node.degree()}</span></div>
    `;
  }
}

function initCytoscape(elements) {
  cy = cytoscape({
    container: document.getElementById('cy'),
    elements,
    layout: { name: 'preset' },

    style: [
      {
        selector: 'node',
        style: {
          'width': 10,
          'height': 10,
          'background-color': '#4a9eff',
          'border-width': 0,
          'label': '',
          'transition-property': 'background-color, border-color, width, height, opacity',
          'transition-duration': '150ms',
          'cursor': 'default',
        }
      },
      {
        selector: 'node[?is_clique]',
        style: {
          'width': 'mapData(member_count, 3, 10, 14, 26)',
          'height': 'mapData(member_count, 3, 10, 14, 26)',
          'background-color': '#00e5ff',
          'border-width': 1.5,
          'border-color': '#00e5ff66',
          'shadow-blur': 8,
          'shadow-color': '#00e5ff',
          'shadow-opacity': 0.6,
          'shadow-offset-x': 0,
          'shadow-offset-y': 0,
          'cursor': 'pointer',
        }
      },
      {
        selector: 'node.expanding',
        style: {
          'opacity': 0.3,
        }
      },
      {
        selector: 'node:selected',
        style: {
          'border-width': 2,
          'border-color': '#ffffff',
          'background-color': '#ffffff',
        }
      },
      {
        selector: 'node:hover',
        style: {
          'border-width': 1.5,
          'border-color': '#ffffff88',
        }
      },
      {
        selector: 'edge',
        style: {
          'width': 'mapData(weight, 0, 5, 0.4, 2)',
          'line-color': '#1e3a4a',
          'opacity': 0.7,
          'curve-style': 'straight',
        }
      },
      {
        selector: 'edge[?is_boundary]',
        style: {
          'line-color': '#ff6b35',
          'opacity': 0.35,
          'line-style': 'dashed',
          'line-dash-pattern': [4, 4],
        }
      },
      // --- Search highlight styles ---
      {
        selector: 'node.search-hit',
        style: {
          'background-color': '#ffe040',
          'border-width': 2.5,
          'border-color': '#ffcc00',
          'width': 16,
          'height': 16,
          'z-index': 999,
        }
      },
      {
        selector: 'node.search-clique-hit',
        style: {
          'background-color': '#ffe040',
          'border-width': 3,
          'border-color': '#ff9900',
          'shadow-blur': 14,
          'shadow-color': '#ffcc00',
          'shadow-opacity': 0.9,
          'shadow-offset-x': 0,
          'shadow-offset-y': 0,
          'z-index': 999,
        }
      },
      {
        selector: 'node.search-dim, edge.search-dim',
        style: {
          'opacity': 0.12,
        }
      },
    ],

    wheelSensitivity: 0.3,
    minZoom: 0.01,
    maxZoom: 20,
  });

  updateStats();

  cy.on('tap', 'node', function(evt) {
    renderDetail(evt.target);
  });

  cy.on('cxttap', 'node', function(evt) {
    const n = evt.target;
    if (!n.data('is_clique')) {
      repackClique(n);
    }
  });

  cy.on('dblclick tap', 'node[?is_clique]', function(evt) {
    evt.stopPropagation();
    expandClique(evt.target);
  });

  cy.on('tap', function(evt) {
    if (evt.target === cy) {
      document.getElementById('detail-content').textContent = 'click a node';
      clearSearchHighlights();
      document.getElementById('search-status').textContent = '';
    }
  });
}

(async () => {
  try {
    const data = await loadTop();
    const elements = buildElements(data);
    initCytoscape(elements);
    document.getElementById('loading').style.display = 'none';

    const searchInput = document.getElementById('search-input');
    const searchBtn = document.getElementById('search-btn');
    const searchClear = document.getElementById('search-clear');

    const doSearch = () => searchGene(searchInput.value);

    searchBtn.addEventListener('click', doSearch);
    searchInput.addEventListener('keydown', e => {
      if (e.key === 'Enter') doSearch();
    });
    searchClear.addEventListener('click', () => {
      searchInput.value = '';
      clearSearchHighlights();
      document.getElementById('search-status').textContent = '';
      searchInput.focus();
    });

  } catch (err) {
    document.getElementById('loading').style.display = 'none';
    showError(`failed to load: ${err.message}`);
  }
})();