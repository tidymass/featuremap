# ============================================================
# 47_synthesis_html.R
# Interactive hypothesis-level synthesis network HTML
# ============================================================

.syn_he <- function(s) {
  s <- as.character(s %||% "")
  s <- gsub("&",  "&amp;",  s, fixed = TRUE)
  s <- gsub("<",  "&lt;",   s, fixed = TRUE)
  s <- gsub(">",  "&gt;",   s, fixed = TRUE)
  s <- gsub('"',  "&quot;", s, fixed = TRUE)
  s
}

.syn_scalar <- function(x, default = "") {
  x <- x %||% default
  x <- as.character(x)
  if (length(x) == 0L) return(as.character(default)[1L])
  x <- trimws(x[1L])
  if (is.na(x) || !nzchar(x) || toupper(x) %in% c("NA", "NULL")) {
    return(as.character(default)[1L])
  }
  x
}

.syn_get_synthesis <- function(x) {
  if (!is.null(x$module_synthesis)) return(x$module_synthesis)
  x
}

.syn_json_records <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(list())
  unname(lapply(seq_len(nrow(df)), function(i) {
    row <- as.list(df[i, , drop = FALSE])
    lapply(row, function(x) .syn_scalar(x, ""))
  }))
}

# ── CSS ───────────────────────────────────────────────────────
.syn_css <- function() {
'*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Helvetica Neue",Arial,sans-serif;
     background:#f4f5f7;color:#1c1f26;font-size:14px;line-height:1.55}

.page-header{background:#519d9e;color:#fff;padding:20px 30px 18px}
.page-header h1{font-size:1.15rem;font-weight:600;letter-spacing:-.01em}
.page-header .sub{font-size:.8rem;color:rgba(255,255,255,.82);margin-top:4px}

.layout{display:grid;grid-template-columns:1fr 300px;gap:16px;
        max-width:1360px;margin:0 auto;padding:18px 18px 0}
.panel{background:#fff;border:1px solid #dbe5f2;border-radius:10px;overflow:hidden}

/* toolbar + legend */
.toolbar{padding:10px 14px;display:flex;align-items:center;gap:10px;flex-wrap:wrap;
         border-bottom:1px solid #eef0f4}
.toolbar label{font-size:.78rem;color:#4b5563;display:flex;align-items:center;gap:5px;
               cursor:pointer}
.toolbar input[type=text]{border:1px solid #d8dde6;border-radius:6px;padding:5px 9px;
                          min-width:210px;font-size:.82rem;outline:none}
.toolbar input[type=text]:focus{border-color:#3966a2}
.legend{padding:8px 14px;display:flex;gap:13px;flex-wrap:wrap;color:#6b7280;font-size:.73rem;
        border-bottom:1px solid #eef0f4;align-items:center}
.legend span{display:inline-flex;align-items:center;gap:4px}
.leg-circle{width:11px;height:11px;border-radius:50%;background:#667085;flex-shrink:0}
.leg-edge-biochem{width:18px;height:0;border-top:2px dashed #2f9e6a;flex-shrink:0}
.leg-edge-biology{width:18px;height:0;border-top:2px dashed #c1602e;flex-shrink:0}
.leg-edge-both{width:18px;height:2px;background:#3966a2;flex-shrink:0}
.leg-sep{width:1px;height:11px;background:#d1d5db}
#cluster-legend{display:flex;gap:8px;flex-wrap:wrap}
#cluster-legend span{display:inline-flex;align-items:center;gap:4px;font-size:.72rem}
.cluster-dot{width:9px;height:9px;border-radius:50%;flex-shrink:0}

/* network SVG */
.network-wrap{position:relative;overflow:hidden;height:600px}
svg{width:100%;height:100%;display:block;background:#fbfcfe}
.edge{stroke-opacity:.65}
.edge.biochemical{stroke:#2f9e6a;stroke-width:1.8;stroke-dasharray:5 4}
.edge.biology{stroke:#c1602e;stroke-width:1.8;stroke-dasharray:5 4}
.edge.both{stroke:#3966a2;stroke-width:2.6}
.node{cursor:default}
.node text{font-size:10.5px;fill:#1f2937;paint-order:stroke;stroke:#fbfcfe;
           stroke-width:2.5px;user-select:none;pointer-events:none}
.isolate-divider{stroke:#e4e7ed;stroke-width:1;stroke-dasharray:3 5}
.isolate-heading{font-size:10px;font-weight:600;fill:#8b909e;text-anchor:middle;
                 text-transform:uppercase;letter-spacing:0}

/* hover tooltip — lightweight, title only, no pointer-events */
#tooltip{position:absolute;pointer-events:none;display:none;z-index:30;
         background:rgba(28,31,38,.82);color:#fff;
         padding:5px 10px;border-radius:6px;
         font-size:.76rem;font-weight:500;white-space:nowrap;
         max-width:260px;overflow:hidden;text-overflow:ellipsis;
         line-height:1.4}

/* pinned card — click to open, scrollable, pointer-events:auto */
#pinned{position:absolute;display:none;z-index:40;
        background:#fff;border:1px solid #d1d5db;border-radius:10px;
        box-shadow:0 6px 24px rgba(0,0,0,.16);
        width:270px;max-height:280px;flex-direction:column;overflow:hidden}

/* study-context pills in page header */
.ctx-pill{display:inline-flex;align-items:center;
          font-size:.75rem;padding:2px 10px;border-radius:20px;
          background:rgba(255,255,255,.16);border:1px solid rgba(255,255,255,.3);
          color:#fff;margin-right:5px}
#pinned .pin-hdr{padding:11px 13px 8px;border-bottom:1px solid #eef0f4;
                 display:flex;align-items:flex-start;gap:6px;flex-shrink:0}
#pinned .pin-id {font-size:.68rem;font-weight:700;color:#3966a2;
                 background:#eef4fb;padding:1px 6px;border-radius:3px;flex-shrink:0}
#pinned .pin-pills{display:flex;flex-wrap:wrap;gap:3px;margin-top:5px}
#pinned .pin-pill{font-size:.67rem;padding:1px 6px;border-radius:3px;
                  background:#f1f3f5;color:#4b5563;border:1px solid #e4e7ed}
#pinned .pin-pill.accepted {background:#eaf7ef;color:#21643a;border-color:#bddfc9}
#pinned .pin-close{margin-left:auto;flex-shrink:0;cursor:pointer;
                   font-size:1rem;color:#9ca3af;line-height:1;
                   padding:0 2px;user-select:none}
#pinned .pin-close:hover{color:#374151}
#pinned .pin-body{padding:10px 13px 13px;overflow-y:auto;flex:1}
#pinned .pin-title{font-size:.85rem;font-weight:600;color:#1c1f26;
                   line-height:1.45;margin-bottom:7px}
#pinned .pin-section{margin-top:9px;padding-top:9px;border-top:1px solid #f0f2f5}
#pinned .pin-section:first-child{margin-top:0;padding-top:0;border-top:none}
#pinned .pin-lbl{font-size:.63rem;font-weight:700;letter-spacing:.08em;
                 text-transform:uppercase;color:#9ca3af;margin-bottom:3px}
#pinned .pin-txt{font-size:.79rem;color:#4b5563;line-height:1.6;white-space:pre-wrap}

/* right summary panel */
.sum-panel{padding:16px;display:flex;flex-direction:column;gap:0}
.sum-panel h2{font-size:.9rem;font-weight:600;margin-bottom:8px;color:#1c1f26}
.sum-label{font-size:.63rem;font-weight:700;letter-spacing:.09em;text-transform:uppercase;
           color:#8b909e;margin-bottom:6px}
.sum-text{font-size:.81rem;color:#374151;line-height:1.65;white-space:pre-wrap}
.sum-empty{font-size:.82rem;color:#9ca3af;font-style:italic}

/* bottom table tabs */
.table-wrap{max-width:1360px;margin:14px auto 24px;padding:0 18px}
.tab-bar{display:flex;border-bottom:2px solid #e4e7ed;margin-bottom:0}
.tab-btn{border:none;background:none;padding:8px 18px;font-size:.8rem;font-weight:600;
         color:#6b7280;cursor:pointer;border-bottom:2px solid transparent;
         margin-bottom:-2px;transition:color .15s}
.tab-btn.active{color:#3966a2;border-bottom-color:#3966a2}
.tab-btn:hover:not(.active){color:#374151}
.tab-panel{display:none}
.tab-panel.active{display:block}

.node-tbl{width:100%;border-collapse:collapse;background:#fff;
          border:1px solid #dbe5f2;border-radius:0 0 10px 10px;overflow:hidden;
          font-size:.78rem;box-shadow:0 1px 3px rgba(0,0,0,.04)}
.node-tbl thead tr{background:#f8f9fb;border-bottom:2px solid #dbe5f2}
.node-tbl th{padding:9px 12px;text-align:left;font-size:.66rem;font-weight:700;
             letter-spacing:.07em;text-transform:uppercase;color:#6b7280}
.node-tbl td{padding:8px 12px;border-bottom:1px solid #f0f2f5;color:#374151;
             vertical-align:top}
.node-tbl tr:last-child td{border-bottom:none}
.node-tbl tr:hover td{background:#fafbff}
.tbl-id{font-weight:600;color:#3966a2;white-space:nowrap}
.tbl-mod{white-space:nowrap}
.tbl-pill{display:inline-block;font-size:.66rem;padding:1px 6px;border-radius:3px;
          font-weight:600;white-space:nowrap}
.tbl-biochemical{background:#eaf7ef;color:#21643a;border:1px solid #bddfc9}
.tbl-biology{background:#faece7;color:#993c1d;border:1px solid #f0997b}
.tbl-both{background:#eef4fb;color:#1c3957;border:1px solid #cfe0f2}
.tbl-conf{color:#6b7280;font-size:.74rem}
.tbl-rat{font-size:.76rem;color:#6b7280;line-height:1.55}

/* hypothesis type pills — pinned card + table */
.type-endo {background:#eaf7ef;color:#21643a;border:1px solid #bddfc9}
.type-micro {background:#fff6e6;color:#80551b;border:1px solid #f2d29a}
.type-diet  {background:#eef4fb;color:#3966a2;border:1px solid #cfe0f2}
.type-xeno  {background:#faece7;color:#993c1d;border:1px solid #f0997b}
.type-mixed {background:#eef1f5;color:#44546b;border:1px solid #d7dee7}
.tbl-type{display:inline-block;font-size:.66rem;padding:1px 7px;border-radius:3px;
          font-weight:500;white-space:nowrap;border:1px solid #dbe5f2;
          background:#f1f3f5;color:#4b5563}

@media(max-width:960px){.layout{grid-template-columns:1fr}
  .network-wrap{height:480px}}'
}

# ── JavaScript ────────────────────────────────────────────────
.syn_js <- function() {
'const data = GRAPH_DATA;
const svg    = document.getElementById("network");
const gGuides = document.getElementById("guides");
const gEdges = document.getElementById("edges");
const gNodes = document.getElementById("nodes");
const tooltip = document.getElementById("tooltip");
const pinned  = document.getElementById("pinned");
const search  = document.getElementById("search");
const controls = {
  biochemical: document.getElementById("fBiochemical"),
  biology:     document.getElementById("fBiology"),
  both:        document.getElementById("fBoth")
};
let pinnedNode = null;

/* ── colour / size maps ── */
const PALETTE = [
  "#3966a2","#519d9e","#c1602e","#2f9e6a",
  "#b5892f","#7a5ea0","#4a7c8c","#a24a6b",
  "#5c8a3f","#8a5a3f","#3f6b8a","#6b6f4a"
];
const UNCLUSTERED = "Unclustered";
const clusterId = n => (n.cluster_id && n.cluster_id.trim()) ? n.cluster_id : UNCLUSTERED;
const naturalCollator = new Intl.Collator(undefined, {numeric: true, sensitivity: "base"});
function clusterCompare(a, b){
  if (a === UNCLUSTERED) return b === UNCLUSTERED ? 0 : 1;
  if (b === UNCLUSTERED) return -1;
  return naturalCollator.compare(a, b);
}
const clusterRows = data.clusters || [];
const clusterTheme = Object.fromEntries(
  clusterRows.map(c => [c.cluster_id, c.cluster_theme || c.cluster_id])
);
const clusterIds = [...new Set(data.nodes.map(clusterId))]
  .sort(clusterCompare);
const clusterColor = Object.fromEntries(
  clusterIds.map((c, i) => [c, c === UNCLUSTERED ? "#667085" : PALETTE[i % PALETTE.length]])
);
const SIZE = {high: 17, medium: 13, low: 10};
const DEF_SIZE = 12;

function nodeSize(n) { return SIZE[n.confidence_label] || DEF_SIZE; }
function nodeColor(n){ return clusterColor[clusterId(n)] || "#667085"; }
function edgeType(e){ return e.relation_type || e.relation_class || ""; }

/* ── cluster colour legend ── */
const clusterLeg = document.getElementById("cluster-legend");
clusterIds.forEach(c => {
  const sp = document.createElement("span");
  const label = c === UNCLUSTERED ? c : `${c}: ${clusterTheme[c] || c}`;
  sp.innerHTML = `<i class="cluster-dot" style="background:${clusterColor[c]}"></i>${esc(label)}`;
  clusterLeg.appendChild(sp);
});

/* ── node / edge data ── */
let nodes = data.nodes.map(d => ({...d, x: 0, y: 0}));
const nodeById = Object.fromEntries(nodes.map(d => [d.node_id, d]));
let edges = data.edges
  .map(e => ({...e, source: nodeById[e.node_1], target: nodeById[e.node_2]}))
  .filter(e => e.source && e.target);
const structurallyConnectedIds = new Set();
edges.forEach(e => {
  structurallyConnectedIds.add(e.node_1);
  structurallyConnectedIds.add(e.node_2);
});
nodes.forEach(n => { n.isIsolated = !structurallyConnectedIds.has(n.node_id); });
let isolateRegion = null;

/* ── helpers ── */
function esc(s){
  return String(s ?? "")
    .replace(/&/g,"&amp;")
    .replace(/</g,"&lt;")
    .replace(/>/g,"&gt;")
    .replace(/"/g,"&quot;");
}

function typeCls(t){
  if (!t) return "";
  const tl = t.toLowerCase();
  if (tl.includes("host") && tl.includes("microbiome")) return "type-mixed";
  if (tl.includes("microbiome") || tl.includes("microbial"))  return "type-micro";
  if (tl.includes("diet") || tl.includes("food"))             return "type-diet";
  if (tl.includes("xenobiotic") || tl.includes("exogenous"))  return "type-xeno";
  if (tl.includes("endogenous") || tl.includes("host"))       return "type-endo";
  return "";
}

function filters(){
  return {
    biochemical: controls.biochemical.checked,
    biology:     controls.biology.checked,
    both:        controls.both.checked,
    q: search.value.trim().toLowerCase()
  };
}

function nodeVisible(n, f){
  if (f.q){
    const hay = [n.node_id,n.module_id,n.hypothesis_id,n.hypothesis_type,
                 n.cluster_id,n.title,n.functional_summary,n.phenotype_relationship]
                .join(" ").toLowerCase();
    if (!hay.includes(f.q)) return false;
  }
  return true;
}

function edgeVisible(e, visIds, f){
  if (!visIds.has(e.node_1) || !visIds.has(e.node_2)) return false;
  const t = edgeType(e);
  if (t === "biochemical") return f.biochemical;
  if (t === "biology")     return f.biology;
  if (t === "both")        return f.both;
  return false;
}

/* ── hit test: all nodes are circles ── */
function isOver(n, mx, my){
  const sz = nodeSize(n);
  const dx = n.x - mx, dy = n.y - my;
  return dx*dx + dy*dy <= sz*sz;
}

/* ── hit test: point-to-segment distance for edge tooltips ── */
function ptSegDist(px, py, x1, y1, x2, y2){
  const dx = x2-x1, dy = y2-y1;
  const lenSq = dx*dx + dy*dy;
  if (lenSq === 0) return Math.hypot(px-x1, py-y1);
  const t = Math.max(0, Math.min(1, ((px-x1)*dx + (py-y1)*dy) / lenSq));
  return Math.hypot(px-(x1+t*dx), py-(y1+t*dy));
}

/* ── hover tooltip — title only, follows mouse, no pointer-events ── */
function showTooltip(ev, n){
  tooltip.textContent = n.title || n.node_id;
  tooltip.style.display = "block";
  positionTooltip(ev);
}

function positionTooltip(ev){
  const wrap = document.querySelector(".network-wrap");
  const r    = wrap.getBoundingClientRect();
  let tx = ev.clientX - r.left + 14;
  let ty = ev.clientY - r.top  - 30;
  if (tx + 270 > wrap.offsetWidth) tx = ev.clientX - r.left - 270;
  if (ty < 4) ty = ev.clientY - r.top + 14;
  tooltip.style.left = Math.max(4, tx) + "px";
  tooltip.style.top  = Math.max(4, ty) + "px";
}

/* ── pinned card — click to open, scrollable ── */
function pinSection(label, text){
  if (!text || !text.trim()) return "";
  return `<div class="pin-section">
    <div class="pin-lbl">${esc(label)}</div>
    <div class="pin-txt">${esc(text)}</div>
  </div>`;
}

function showPin(ev, n){
  pinnedNode = n.node_id;
  const confLabel = n.confidence_label ? n.confidence_label + " confidence" : "";
  const cl = clusterId(n);
  const clLabel = cl === UNCLUSTERED ? cl : `${cl}: ${clusterTheme[cl] || cl}`;
  pinned.innerHTML =
    `<div class="pin-hdr">
       <div>
         <span class="pin-id">${esc(n.node_id)}</span>
         <div class="pin-pills">
           ${confLabel ? `<span class="pin-pill">${esc(confLabel)}</span>` : ""}
           <span class="pin-pill">${esc(clLabel)}</span>
           ${n.hypothesis_type ? `<span class="pin-pill ${typeCls(n.hypothesis_type)}">${esc(n.hypothesis_type)}</span>` : ""}
         </div>
       </div>
       <span class="pin-close" id="pin-close-btn">&#x2715;</span>
     </div>
     <div class="pin-body">
       ${n.title ? `<div class="pin-title">${esc(n.title)}</div>` : ""}
       ${pinSection("Functional summary", n.functional_summary)}
     </div>`;

  /* position near node */
  const wrap = document.querySelector(".network-wrap");
  const wW = wrap.offsetWidth, wH = wrap.offsetHeight;
  const cardW = 270, cardH = Math.min(280, wH - 20);
  let px = n.x + nodeSize(n) + 10;
  let py = n.y - 60;
  if (px + cardW > wW - 8) px = n.x - nodeSize(n) - cardW - 10;
  if (py + cardH > wH - 8) py = wH - cardH - 8;
  if (py < 8) py = 8;
  if (px < 8) px = 8;

  pinned.style.left    = px + "px";
  pinned.style.top     = py + "px";
  pinned.style.display = "flex";

  /* close button */
  document.getElementById("pin-close-btn")
    .addEventListener("click", ev2 => { ev2.stopPropagation(); hidePin(); });
}

function hidePin(){
  pinnedNode = null;
  pinned.style.display = "none";
  render();
}

/* ── SVG-level mouse tracking (survives DOM re-renders every tick) ── */
svg.addEventListener("mousemove", ev => {
  const r  = svg.getBoundingClientRect();
  const mx = ev.clientX - r.left;
  const my = ev.clientY - r.top;
  const f  = filters();
  const vis    = nodes.filter(n => nodeVisible(n, f));
  const visIds = new Set(vis.map(n => n.node_id));
  let hit = null;
  for (let i = vis.length - 1; i >= 0; i--) {
    if (isOver(vis[i], mx, my)) { hit = vis[i]; break; }
  }
  if (hit && hit.node_id !== pinnedNode) {
    /* node hover: show title */
    showTooltip(ev, hit);
  } else if (!hit) {
    /* edge hover: find nearest visible edge within 6px */
    const visE = edges.filter(e => edgeVisible(e, visIds, f));
    let hitEdge = null;
    for (const e of visE) {
      if (ptSegDist(mx, my,
                    e.source.x, e.source.y,
                    e.target.x, e.target.y) < 6) {
        hitEdge = e; break;
      }
    }
    if (hitEdge) {
      tooltip.textContent = edgeType(hitEdge);
      tooltip.style.display = "block";
      positionTooltip(ev);
    } else {
      tooltip.style.display = "none";
    }
  } else {
    tooltip.style.display = "none";
  }
});

svg.addEventListener("mouseleave", () => {
  tooltip.style.display = "none";
});

/* click on SVG background dismisses pin */
svg.addEventListener("click", ev => {
  const r  = svg.getBoundingClientRect();
  const mx = ev.clientX - r.left;
  const my = ev.clientY - r.top;
  const f  = filters();
  const vis = nodes.filter(n => nodeVisible(n, f));
  const hit = vis.find(n => isOver(n, mx, my));
  if (hit) {
    tooltip.style.display = "none";
    showPin(ev, hit);
    render();  /* re-render to highlight selected node */
  } else {
    hidePin();
  }
});

/* ── static layout ── */
function computeLayout(){
  let rect = svg.getBoundingClientRect();
  const W = rect.width || 860;
  const connectedNodes = nodes.filter(n => !n.isIsolated);
  const isolatedNodes = nodes.filter(n => n.isIsolated);

  const railWidth = isolatedNodes.length > 0
    ? Math.min(420, Math.max(250, W * 0.34))
    : 0;
  const graphW = isolatedNodes.length > 0 ? Math.max(260, W - railWidth - 16) : W;
  const railX = graphW + 16;
  const railInnerWidth = Math.max(1, W - railX - 20);
  const isoCols = isolatedNodes.length > 0
    ? Math.max(1, Math.floor(railInnerWidth / 72))
    : 1;
  const isolatedGroups = clusterIds
    .map(id => ({
      id,
      nodes: isolatedNodes
        .filter(n => clusterId(n) === id)
        .sort((a, b) => naturalCollator.compare(a.node_id, b.node_id))
    }))
    .filter(group => group.nodes.length > 0);
  const isoRows = isolatedGroups.reduce(
    (sum, group) => sum + Math.ceil(group.nodes.length / isoCols), 0
  );

  const wrap = svg.closest(".network-wrap");
  const baseHeight = window.matchMedia("(max-width:960px)").matches ? 480 : 600;
  const requiredHeight = isolatedNodes.length > 0 ? 58 + isoRows * 52 : baseHeight;
  wrap.style.height = Math.max(baseHeight, requiredHeight) + "px";
  rect = svg.getBoundingClientRect();
  const H = rect.height || baseHeight;
  const cx = graphW / 2, cy = H / 2;

  /* group nodes by module */
  const groups = {};
  connectedNodes.forEach(n => { (groups[n.module_id] = groups[n.module_id] || []).push(n); });
  const mids = Object.keys(groups).sort((a, b) => naturalCollator.compare(a, b));
  const nMod = mids.length;

  /* place module centres in a circle */
  const Rmod = Math.min(cx, cy) * (nMod <= 2 ? 0.30 : nMod <= 4 ? 0.44 : 0.52);
  mids.forEach((mid, mi) => {
    const angle = 2 * Math.PI * mi / nMod - Math.PI / 2;
    const mcx   = cx + Rmod * Math.cos(angle);
    const mcy   = cy + Rmod * Math.sin(angle);
    const grp   = groups[mid];
    const nN    = grp.length;
    /* spread nodes in a small circle around module centre */
    const Rnode = nN === 1 ? 0 : Math.max(28, nN * 13);
    grp.forEach((n, ni) => {
      const a = 2 * Math.PI * ni / nN - Math.PI / 2;
      n.x = mcx + Rnode * Math.cos(a);
      n.y = mcy + Rnode * Math.sin(a);
    });
  });

  /* iterative overlap resolution — no animation, just converge */
  for (let iter = 0; iter < 200; iter++){
    /* repulsion between all pairs */
    for (let i = 0; i < connectedNodes.length; i++){
      for (let j = i+1; j < connectedNodes.length; j++){
        const a = connectedNodes[i], b = connectedNodes[j];
        let dx = a.x - b.x, dy = a.y - b.y;
        const d2 = dx*dx + dy*dy + 0.01;
        const d  = Math.sqrt(d2);
        const minD = nodeSize(a) + nodeSize(b) + 44;
        if (d < minD){
          const push = (minD - d) * 0.45;
          dx /= d; dy /= d;
          a.x += dx * push; a.y += dy * push;
          b.x -= dx * push; b.y -= dy * push;
        }
      }
    }
    /* gentle spring along edges (pull connected nodes closer) */
    edges.forEach(e => {
      const a = e.source, b = e.target;
      let dx = b.x - a.x, dy = b.y - a.y;
      const d = Math.sqrt(dx*dx + dy*dy) + 0.01;
      const target = edgeType(e) === "both" ? 140 : 180;
      if (d > target){
        const pull = (d - target) * 0.04;
        dx /= d; dy /= d;
        a.x += dx * pull; a.y += dy * pull;
        b.x -= dx * pull; b.y -= dy * pull;
      }
    });
    /* soft pull back to module centre */
    mids.forEach(mid => {
      const grp = groups[mid];
      const mcx = grp.reduce((s, n) => s + n.x, 0) / grp.length;  /* current centroid */
      const mcy = grp.reduce((s, n) => s + n.y, 0) / grp.length;
    });
    /* boundary clamp */
    connectedNodes.forEach(n => {
      const pad = nodeSize(n) + 16;
      n.x = Math.max(pad, Math.min(graphW - pad, n.x));
      n.y = Math.max(pad, Math.min(H - pad, n.y));
    });
  }

  if (isolatedNodes.length > 0){
    const top = 42;
    const bottom = H - 18;
    const rowStep = isoRows > 0 ? (bottom - top) / isoRows : 0;
    const slotWidth = railInnerWidth / isoCols;
    let rowCursor = 0;

    isolatedGroups.forEach(group => {
      const rowsInGroup = Math.ceil(group.nodes.length / isoCols);
      for (let r = 0; r < rowsInGroup; r++){
        const rowNodes = group.nodes.slice(r * isoCols, (r + 1) * isoCols);
        const rowWidth = rowNodes.length * slotWidth;
        const startX = railX + 10 + (railInnerWidth - rowWidth) / 2;
        rowNodes.forEach((n, col) => {
          n.x = startX + (col + 0.5) * slotWidth;
          n.y = top + (rowCursor + r + 0.5) * rowStep;
        });
      }
      rowCursor += rowsInGroup;
    });
    isolateRegion = {dividerX: railX - 8, centreX: railX + railInnerWidth / 2 + 10, height: H};
  } else {
    isolateRegion = null;
  }
}

/* ── render ── */
function render(){
  const f      = filters();
  const vis    = nodes.filter(n => nodeVisible(n, f));
  const visIds = new Set(vis.map(n => n.node_id));
  const visE   = edges.filter(e => edgeVisible(e, visIds, f));

  /* stable side area for nodes with no relationships */
  gGuides.innerHTML = "";
  if (isolateRegion){
    const divider = document.createElementNS("http://www.w3.org/2000/svg", "line");
    divider.setAttribute("x1", isolateRegion.dividerX);
    divider.setAttribute("x2", isolateRegion.dividerX);
    divider.setAttribute("y1", 16);
    divider.setAttribute("y2", isolateRegion.height - 16);
    divider.setAttribute("class", "isolate-divider");
    gGuides.appendChild(divider);

    const heading = document.createElementNS("http://www.w3.org/2000/svg", "text");
    heading.setAttribute("x", isolateRegion.centreX);
    heading.setAttribute("y", 20);
    heading.setAttribute("class", "isolate-heading");
    heading.textContent = "Isolated nodes";
    gGuides.appendChild(heading);
  }

  /* edges */
  gEdges.innerHTML = "";
  for (const e of visE){
    const line = document.createElementNS("http://www.w3.org/2000/svg","line");
    line.setAttribute("x1", e.source.x); line.setAttribute("y1", e.source.y);
    line.setAttribute("x2", e.target.x); line.setAttribute("y2", e.target.y);
    line.setAttribute("class", `edge ${edgeType(e)}`);
    gEdges.appendChild(line);
  }

  /* nodes */
  gNodes.innerHTML = "";
  for (const n of vis){
    const grp = document.createElementNS("http://www.w3.org/2000/svg","g");
    grp.setAttribute("class", "node");
    grp.setAttribute("transform", `translate(${n.x},${n.y})`);

    const sz   = nodeSize(n);
    const fill = nodeColor(n);
    let shape;

    shape = document.createElementNS("http://www.w3.org/2000/svg","circle");
    shape.setAttribute("r", sz);
    const isPinned = n.node_id === pinnedNode;
    shape.setAttribute("fill", fill);
    shape.setAttribute("stroke", isPinned ? "#1c1f26" : "#fff");
    shape.setAttribute("stroke-width", isPinned ? "3.5" : "2.5");

    const t = document.createElementNS("http://www.w3.org/2000/svg","text");
    if (n.isIsolated){
      t.setAttribute("x", 0);
      t.setAttribute("y", sz + 12);
      t.setAttribute("text-anchor", "middle");
    } else {
      t.setAttribute("x", sz + 4);
      t.setAttribute("y", 4);
    }
    t.textContent = n.node_id;

    grp.appendChild(shape);
    grp.appendChild(t);
    gNodes.appendChild(grp);
  }
}

/* ── node table ── */
const tbody = document.getElementById("tbl-body");
[...nodes]
  .sort((a, b) => clusterCompare(clusterId(a), clusterId(b)) ||
    naturalCollator.compare(a.node_id, b.node_id))
  .forEach(n => {
    const tr = document.createElement("tr");
    tr.innerHTML =
      `<td class="tbl-id">${esc(n.node_id)}</td>` +
      `<td class="tbl-mod">${esc(n.module_id)}</td>` +
      `<td class="tbl-mod">${esc(clusterId(n))}</td>` +
      `<td><span class="tbl-type ${typeCls(n.hypothesis_type)}">${esc(n.hypothesis_type || "—")}</span></td>` +
      `<td class="tbl-conf">${esc(n.confidence_label || "—")}</td>` +
      `<td>${esc(n.title || "—")}</td>`;
    tbody.appendChild(tr);
  });

/* ── edge table ── */
const etbody = document.getElementById("etbl-body");
data.edges.forEach(e => {
  const tr = document.createElement("tr");
  const rel = edgeType(e);
  const pillCls = rel === "both" ? "tbl-both" :
    rel === "biochemical" ? "tbl-biochemical" : "tbl-biology";
  tr.innerHTML =
    `<td class="tbl-id">${esc(e.node_1)}</td>` +
    `<td class="tbl-id">${esc(e.node_2)}</td>` +
    `<td><span class="tbl-pill ${pillCls}">${esc(rel)}</span></td>` +
    `<td class="tbl-rat">${esc(e.biochemical_evidence_detail || "—")}</td>` +
    `<td class="tbl-rat">${esc(e.biology_rationale || "—")}</td>`;
  etbody.appendChild(tr);
});

/* ── tab switching ── */
document.querySelectorAll(".tab-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    document.querySelectorAll(".tab-btn").forEach(b => b.classList.remove("active"));
    document.querySelectorAll(".tab-panel").forEach(p => p.classList.remove("active"));
    btn.classList.add("active");
    document.getElementById("tab-" + btn.dataset.tab).classList.add("active");
  });
});

/* ── summary ── */
const sumEl = document.getElementById("sum-text");
if (data.summary && data.summary.trim()) {
  sumEl.textContent = data.summary;
} else {
  sumEl.innerHTML = `<span class="sum-empty">No integrated summary available.</span>`;
}

/* ── wire controls ── */
for (const el of Object.values(controls)) el.addEventListener("change", render);
search.addEventListener("input", render);

window.addEventListener("load", function(){
  computeLayout();
  render();
});'
}

# ── R public function ─────────────────────────────────────────

#' Generate an interactive hypothesis-level synthesis network HTML
#'
#' @param synthesis A \code{res$module_synthesis} list, or a featuremap result
#'   object containing \code{module_synthesis}.
#' @param study_context Optional list with \code{species}, \code{biospecimen},
#'   \code{phenotype} to display in the page header.
#' @param output_dir Directory for the output HTML.
#' @param output_file File name for the output HTML.
#' @param show_same_module Ignored (same-module links are never rendered).
#' @return Invisibly, the written HTML path.
#' @export
generate_synthesis_network_html <- function(synthesis,
                                            study_context    = NULL,
                                            output_dir       = getwd(),
                                            output_file      = "synthesis_network.html",
                                            show_same_module = FALSE) {
  synthesis <- .syn_get_synthesis(synthesis)
  if (is.null(synthesis$nodes) || !is.data.frame(synthesis$nodes))
    stop("generate_synthesis_network_html() requires synthesis$nodes.")

  nodes <- synthesis$nodes
  relationships <- synthesis$relationships %||% synthesis$relationships_sparse
  if (is.null(relationships) || !is.data.frame(relationships)) {
    relationships <- data.frame(
      node_1 = character(), node_2 = character(),
      relation_type = character(), relation_class = character(), subtype = character(),
      rationale = character(), stringsAsFactors = FALSE
    )
  }
  if (!"relation_type" %in% colnames(relationships)) {
    relationships$relation_type <- relationships$relation_class %||% NA_character_
  }

  nodes <- nodes[nodes$verdict_status == "accepted", , drop = FALSE]
  vn    <- nodes$node_id
  relationships <- relationships[
    relationships$node_1 %in% vn &
    relationships$node_2 %in% vn &
    relationships$relation_type %in% c("biochemical", "biology", "both"),
    , drop = FALSE
  ]

  graph_json <- jsonlite::toJSON(
    list(nodes   = .syn_json_records(nodes),
         edges   = .syn_json_records(relationships),
         clusters = .syn_json_records(synthesis$clusters),
         summary = .syn_scalar(synthesis$summary, "")),
    auto_unbox = TRUE, dataframe = "rows", na = "null"
  )
  js <- sub("GRAPH_DATA", graph_json, .syn_js(), fixed = TRUE)

  html <- sprintf(
'<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Synthesis Network</title>
<style>%s</style>
</head>
<body>

<div class="layout">

  <!-- network -->
  <div class="panel">
    <div class="toolbar">
      <label><input type="checkbox" id="fBiochemical" checked> Biochemical links</label>
      <label><input type="checkbox" id="fBiology" checked> Biology links</label>
      <label><input type="checkbox" id="fBoth" checked> Both</label>
      <input type="text" id="search" placeholder="Search module, cluster, hypothesis, title&hellip;">
    </div>
    <div class="legend">
      <span><i class="leg-circle" style="background:#667085"></i>Node color&nbsp;&rarr;&nbsp;cluster</span>
      <i class="leg-sep"></i>
      <span>Size&nbsp;&rarr;&nbsp;confidence</span>
      <i class="leg-sep"></i>
      <span><i class="leg-edge-biochem"></i>Biochemical</span>
      <span><i class="leg-edge-biology"></i>Biology</span>
      <span><i class="leg-edge-both"></i>Both</span>
      <i class="leg-sep"></i>
      <div id="cluster-legend"></div>
    </div>
    <div class="network-wrap">
      <svg id="network"><g id="guides"></g><g id="edges"></g><g id="nodes"></g></svg>
      <div id="tooltip"></div>
      <div id="pinned"></div>
    </div>
  </div>

  <!-- summary -->
  <div class="panel sum-panel">
    <div class="sum-label">Integrated Summary</div>
    <div class="sum-text" id="sum-text"></div>
  </div>

</div>

<!-- node / edge tables with tabs -->
<div class="table-wrap">
  <div class="tab-bar">
    <button class="tab-btn active" data-tab="nodes">Nodes</button>
    <button class="tab-btn" data-tab="edges">Edges</button>
  </div>

  <div class="tab-panel active" id="tab-nodes">
    <table class="node-tbl">
      <thead>
        <tr>
          <th>Hypothesis</th>
          <th>Module</th>
          <th>Cluster</th>
          <th>Type</th>
          <th>Confidence</th>
          <th>Title</th>
        </tr>
      </thead>
      <tbody id="tbl-body"></tbody>
    </table>
  </div>

  <div class="tab-panel" id="tab-edges">
    <table class="node-tbl">
      <thead>
        <tr>
          <th>Node 1</th>
          <th>Node 2</th>
          <th>Relation</th>
          <th>Biochemical Evidence</th>
          <th>Biology Rationale</th>
        </tr>
      </thead>
      <tbody id="etbl-body"></tbody>
    </table>
  </div>
</div>

<script>%s</script>
</body>
</html>',
    .syn_css(), js
  )

  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(output_dir, output_file)
  writeLines(html, out_path, useBytes = FALSE)
  message(sprintf("[HTML] Synthesis network -> %s", out_path))
  invisible(out_path)
}

#' Export an interactive synthesis network from a featuremap result
#'
#' @param res A featuremap result containing \code{module_synthesis}.
#' @param html_dir Directory for the output HTML.
#' @param output_file File name for the output HTML.
#' @param show_same_module Ignored (same-module links are never rendered).
#' @param output_dir Deprecated alias for \code{html_dir}.
#' @return Invisibly, the written HTML path.
#' @export
export_synthesis_network_html <- function(res,
                                          html_dir    = getwd(),
                                          output_file = "synthesis_network.html",
                                          show_same_module = FALSE,
                                          output_dir  = html_dir) {
  if (is.null(res$module_synthesis))
    stop("res$module_synthesis is missing. Run synthesize_modules() first.")
  if (!missing(output_dir)) html_dir <- output_dir
  generate_synthesis_network_html(
    synthesis        = res$module_synthesis,
    study_context    = res$study_context,
    output_dir       = html_dir,
    output_file      = output_file,
    show_same_module = show_same_module
  )
}
