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

.page-header{background:#1c1f26;color:#fff;padding:20px 30px 18px}
.page-header h1{font-size:1.15rem;font-weight:600;letter-spacing:-.01em}
.page-header .sub{font-size:.8rem;color:#a8adba;margin-top:4px}

.layout{display:grid;grid-template-columns:1fr 300px;gap:16px;
        max-width:1360px;margin:0 auto;padding:18px 18px 0}
.panel{background:#fff;border:1px solid #e4e7ed;border-radius:10px;overflow:hidden}

/* toolbar + legend */
.toolbar{padding:10px 14px;display:flex;align-items:center;gap:10px;flex-wrap:wrap;
         border-bottom:1px solid #eef0f4}
.toolbar label{font-size:.78rem;color:#4b5563;display:flex;align-items:center;gap:5px;
               cursor:pointer}
.toolbar input[type=text]{border:1px solid #d8dde6;border-radius:6px;padding:5px 9px;
                          min-width:210px;font-size:.82rem;outline:none}
.toolbar input[type=text]:focus{border-color:#6741d9}
.legend{padding:8px 14px;display:flex;gap:13px;flex-wrap:wrap;color:#6b7280;font-size:.73rem;
        border-bottom:1px solid #eef0f4;align-items:center}
.legend span{display:inline-flex;align-items:center;gap:4px}
.leg-circle{width:11px;height:11px;border-radius:50%;background:#667085;flex-shrink:0}
.leg-square{width:11px;height:11px;border-radius:2px;background:#667085;flex-shrink:0}
.leg-edge-d{width:18px;height:2px;background:#516f9f;flex-shrink:0}
.leg-edge-i{width:18px;height:0;border-top:2px dashed #9a8064;flex-shrink:0}
.leg-sep{width:1px;height:11px;background:#d1d5db}
#module-legend{display:flex;gap:8px;flex-wrap:wrap}
#module-legend span{display:inline-flex;align-items:center;gap:4px;font-size:.72rem}
.mod-dot{width:9px;height:9px;border-radius:50%;flex-shrink:0}

/* network SVG */
.network-wrap{position:relative;overflow:hidden;height:600px}
svg{width:100%;height:100%;display:block;background:#fbfcfe}
.edge{stroke-opacity:.65}
.edge.direct  {stroke:#516f9f;stroke-width:2}
.edge.indirect{stroke:#9a8064;stroke-width:1.6;stroke-dasharray:5 4}
.node{cursor:default}
.node text{font-size:10.5px;fill:#1f2937;paint-order:stroke;stroke:#fbfcfe;
           stroke-width:2.5px;user-select:none;pointer-events:none}

/* hover tooltip card */
#tooltip{position:absolute;pointer-events:none;display:none;z-index:30;
         background:#fff;border:1px solid #d1d5db;border-radius:9px;
         box-shadow:0 4px 18px rgba(0,0,0,.13);padding:12px 14px;
         max-width:280px;min-width:180px;line-height:1.5}
#tooltip .tt-id   {font-size:.68rem;font-weight:700;color:#6741d9;
                   background:#ede9fb;padding:1px 6px;border-radius:3px;
                   display:inline-block;margin-bottom:6px}
#tooltip .tt-pills{display:flex;flex-wrap:wrap;gap:4px;margin-bottom:7px}
#tooltip .tt-pill {font-size:.68rem;padding:1px 6px;border-radius:4px;
                   background:#f1f3f5;color:#4b5563;border:1px solid #e4e7ed}
#tooltip .tt-pill.accepted {background:#ecfdf5;color:#065f46;border-color:#a7f3d0}
#tooltip .tt-pill.uncertain{background:#fffbeb;color:#92400e;border-color:#fcd34d}
#tooltip .tt-title  {font-size:.83rem;font-weight:600;color:#1c1f26;margin-bottom:5px}
#tooltip .tt-summary{font-size:.78rem;color:#4b5563;line-height:1.55;
                     max-height:110px;overflow-y:auto}

/* right summary panel */
.sum-panel{padding:16px;display:flex;flex-direction:column;gap:0}
.sum-panel h2{font-size:.9rem;font-weight:600;margin-bottom:8px;color:#1c1f26}
.sum-label{font-size:.63rem;font-weight:700;letter-spacing:.09em;text-transform:uppercase;
           color:#8b909e;margin-bottom:6px}
.sum-text{font-size:.81rem;color:#374151;line-height:1.65;white-space:pre-wrap}
.sum-empty{font-size:.82rem;color:#9ca3af;font-style:italic}

/* bottom table */
.table-wrap{max-width:1360px;margin:14px auto 24px;padding:0 18px}
.node-tbl{width:100%;border-collapse:collapse;background:#fff;
          border:1px solid #e4e7ed;border-radius:10px;overflow:hidden;
          font-size:.78rem;box-shadow:0 1px 3px rgba(0,0,0,.04)}
.node-tbl thead tr{background:#f8f9fb;border-bottom:2px solid #e4e7ed}
.node-tbl th{padding:9px 12px;text-align:left;font-size:.66rem;font-weight:700;
             letter-spacing:.07em;text-transform:uppercase;color:#6b7280}
.node-tbl td{padding:8px 12px;border-bottom:1px solid #f0f2f5;color:#374151;
             vertical-align:top}
.node-tbl tr:last-child td{border-bottom:none}
.node-tbl tr:hover td{background:#fafbff}
.tbl-id{font-weight:600;color:#6741d9;white-space:nowrap}
.tbl-mod{white-space:nowrap}
.tbl-pill{display:inline-block;font-size:.66rem;padding:1px 6px;border-radius:3px;
          font-weight:600;white-space:nowrap}
.tbl-acc{background:#ecfdf5;color:#065f46}
.tbl-unc{background:#fffbeb;color:#92400e}
.tbl-conf{color:#6b7280;font-size:.74rem}

@media(max-width:960px){.layout{grid-template-columns:1fr}
  .network-wrap{height:480px}}'
}

# ── JavaScript ────────────────────────────────────────────────
.syn_js <- function() {
'const data = GRAPH_DATA;
const svg    = document.getElementById("network");
const gEdges = document.getElementById("edges");
const gNodes = document.getElementById("nodes");
const tooltip = document.getElementById("tooltip");
const search  = document.getElementById("search");
const controls = {
  direct:   document.getElementById("fDirect"),
  indirect: document.getElementById("fIndirect")
};

/* ── colour / size maps ── */
const PALETTE = [
  "#395f9d","#8a5a9d","#2d7a55","#b55a2a",
  "#7a6a2d","#2f7280","#8a3f58","#5b6f91",
  "#4a7c59","#7c4a6b","#6b4a7c","#4a6b7c"
];
const moduleIds   = [...new Set(data.nodes.map(d => d.module_id))].sort();
const moduleColor = Object.fromEntries(
  moduleIds.map((m, i) => [m, PALETTE[i % PALETTE.length]])
);
const SIZE = {high: 17, medium: 13, low: 10};
const DEF_SIZE = 12;

function nodeSize(n) { return SIZE[n.confidence_label] || DEF_SIZE; }
function nodeColor(n){ return moduleColor[n.module_id] || "#667085"; }

/* ── module colour legend ── */
const modLeg = document.getElementById("module-legend");
moduleIds.forEach(m => {
  const sp = document.createElement("span");
  sp.innerHTML = `<i class="mod-dot" style="background:${nodeColor({module_id:m})}"></i>${esc(m)}`;
  modLeg.appendChild(sp);
});

/* ── physics state ── */
const W0 = 900, H0 = 600;
let nodes = data.nodes.map((d, i) => {
  const cols = Math.max(1, Math.ceil(Math.sqrt(data.nodes.length * 1.4)));
  const col  = i % cols, row = Math.floor(i / cols);
  return {
    ...d,
    x: 80 + col * (W0 - 160) / Math.max(cols - 1, 1) + (Math.random() - .5) * 60,
    y: 80 + row * (H0 - 160) / Math.max(Math.ceil(data.nodes.length / cols) - 1, 1) + (Math.random() - .5) * 60,
    vx: 0, vy: 0
  };
});
const nodeById = Object.fromEntries(nodes.map(d => [d.node_id, d]));
let edges = data.edges
  .map(e => ({...e, source: nodeById[e.node_1], target: nodeById[e.node_2]}))
  .filter(e => e.source && e.target);

/* ── helpers ── */
function esc(s){
  return String(s ?? "")
    .replace(/&/g,"&amp;")
    .replace(/</g,"&lt;")
    .replace(/>/g,"&gt;")
    .replace(/"/g,"&quot;");
}

function filters(){
  return {
    direct:   controls.direct.checked,
    indirect: controls.indirect.checked,
    q: search.value.trim().toLowerCase()
  };
}

function nodeVisible(n, f){
  if (f.q){
    const hay = [n.node_id,n.module_id,n.hypothesis_id,
                 n.title,n.functional_summary,n.phenotype_relationship]
                .join(" ").toLowerCase();
    if (!hay.includes(f.q)) return false;
  }
  return true;
}

function edgeVisible(e, visIds, f){
  if (!visIds.has(e.node_1) || !visIds.has(e.node_2)) return false;
  if (e.relation_class === "direct")   return f.direct;
  if (e.relation_class === "indirect") return f.indirect;
  return false;
}

/* ── hit test (works for circle AND rect nodes) ── */
function isOver(n, mx, my){
  const sz = nodeSize(n);
  const dx = n.x - mx, dy = n.y - my;
  if (n.verdict_status === "accepted") return dx*dx + dy*dy <= sz*sz;
  return Math.abs(dx) <= sz && Math.abs(dy) <= sz;
}

/* ── tooltip via SVG-level mouse tracking (survives DOM re-renders) ── */
function showTooltip(ev, n){
  const statusCls = n.verdict_status === "accepted" ? "accepted" : "uncertain";
  tooltip.innerHTML =
    `<div class="tt-id">${esc(n.node_id)}</div>` +
    `<div class="tt-pills">` +
      `<span class="tt-pill ${statusCls}">${esc(n.verdict_status)}</span>` +
      (n.confidence_label ? `<span class="tt-pill">${esc(n.confidence_label)} confidence</span>` : "") +
      `<span class="tt-pill">${esc(n.module_id)}</span>` +
    `</div>` +
    (n.title    ? `<div class="tt-title">${esc(n.title)}</div>` : "") +
    (n.functional_summary ? `<div class="tt-summary">${esc(n.functional_summary)}</div>` : "");
  tooltip.style.display = "block";
  positionTooltip(ev);
}

function positionTooltip(ev){
  const wrap = document.querySelector(".network-wrap");
  const r    = wrap.getBoundingClientRect();
  let tx = ev.clientX - r.left + 16;
  let ty = ev.clientY - r.top  + 16;
  // flip left if overflowing right
  if (tx + 290 > wrap.offsetWidth)  tx = ev.clientX - r.left - 294;
  // flip up if overflowing bottom
  if (ty + tooltip.offsetHeight + 10 > wrap.offsetHeight)
    ty = ev.clientY - r.top - tooltip.offsetHeight - 10;
  tooltip.style.left = Math.max(4, tx) + "px";
  tooltip.style.top  = Math.max(4, ty) + "px";
}

svg.addEventListener("mousemove", ev => {
  const r  = svg.getBoundingClientRect();
  const mx = ev.clientX - r.left;
  const my = ev.clientY - r.top;
  const f  = filters();
  const vis = nodes.filter(n => nodeVisible(n, f));

  // check nodes in reverse paint order so top-most wins
  let hit = null;
  for (let i = vis.length - 1; i >= 0; i--) {
    if (isOver(vis[i], mx, my)) { hit = vis[i]; break; }
  }
  if (hit) showTooltip(ev, hit);
  else      tooltip.style.display = "none";
});

svg.addEventListener("mouseleave", () => {
  tooltip.style.display = "none";
});

/* ── physics ── */
function tick(){
  const rect = svg.getBoundingClientRect();
  const W = rect.width  || W0;
  const H = rect.height || H0;
  const f  = filters();
  const vis    = nodes.filter(n => nodeVisible(n, f));
  const visIds = new Set(vis.map(n => n.node_id));
  const visE   = edges.filter(e => edgeVisible(e, visIds, f));

  /* centering */
  for (const n of vis){
    n.vx *= .84; n.vy *= .84;
    n.vx += (W/2 - n.x) * .0007;
    n.vy += (H/2 - n.y) * .0007;
  }
  /* repulsion – stronger to keep nodes spread */
  for (let i = 0; i < vis.length; i++){
    for (let j = i+1; j < vis.length; j++){
      const a = vis[i], b = vis[j];
      let dx = a.x - b.x, dy = a.y - b.y;
      const d2 = dx*dx + dy*dy + 1;
      const d  = Math.sqrt(d2);
      const minD = nodeSize(a) + nodeSize(b) + 55;
      /* hard push when too close */
      const push = minD*minD*3 / (d2 + 200);
      dx /= d; dy /= d;
      a.vx += dx*push; a.vy += dy*push;
      b.vx -= dx*push; b.vy -= dy*push;
    }
  }
  /* spring attraction along edges */
  for (const e of visE){
    const a = e.source, b = e.target;
    let dx = b.x - a.x, dy = b.y - a.y;
    const d = Math.sqrt(dx*dx + dy*dy) + .01;
    const target = e.relation_class === "direct" ? 170 : 220;
    const k = .004;
    const force = (d - target) * k;
    dx /= d; dy /= d;
    a.vx += dx*force; a.vy += dy*force;
    b.vx -= dx*force; b.vy -= dy*force;
  }
  /* boundary */
  for (const n of vis){
    const pad = nodeSize(n) + 12;
    n.x = Math.max(pad, Math.min(W-pad, n.x + n.vx));
    n.y = Math.max(pad, Math.min(H-pad, n.y + n.vy));
  }
  render();
  requestAnimationFrame(tick);
}

/* ── render ── */
function render(){
  const f      = filters();
  const vis    = nodes.filter(n => nodeVisible(n, f));
  const visIds = new Set(vis.map(n => n.node_id));
  const visE   = edges.filter(e => edgeVisible(e, visIds, f));

  /* edges */
  gEdges.innerHTML = "";
  for (const e of visE){
    const line = document.createElementNS("http://www.w3.org/2000/svg","line");
    line.setAttribute("x1", e.source.x); line.setAttribute("y1", e.source.y);
    line.setAttribute("x2", e.target.x); line.setAttribute("y2", e.target.y);
    line.setAttribute("class", `edge ${e.relation_class}`);
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
    const isAcc = n.verdict_status === "accepted";
    let shape;

    if (isAcc){
      shape = document.createElementNS("http://www.w3.org/2000/svg","circle");
      shape.setAttribute("r", sz);
    } else {
      shape = document.createElementNS("http://www.w3.org/2000/svg","rect");
      shape.setAttribute("x", -sz); shape.setAttribute("y", -sz);
      shape.setAttribute("width", sz*2); shape.setAttribute("height", sz*2);
      shape.setAttribute("rx", "3");
    }
    shape.setAttribute("fill", fill);
    shape.setAttribute("fill-opacity", isAcc ? "0.9" : "0.6");
    shape.setAttribute("stroke", "#fff");
    shape.setAttribute("stroke-width", "2");

    const t = document.createElementNS("http://www.w3.org/2000/svg","text");
    t.setAttribute("x", sz + 4);
    t.setAttribute("y", 4);
    t.textContent = n.node_id;

    grp.appendChild(shape);
    grp.appendChild(t);
    gNodes.appendChild(grp);
  }
}

/* ── bottom table ── */
const tbody = document.getElementById("tbl-body");
nodes.forEach(n => {
  const isAcc = n.verdict_status === "accepted";
  const tr = document.createElement("tr");
  tr.innerHTML =
    `<td class="tbl-id">${esc(n.node_id)}</td>` +
    `<td class="tbl-mod">${esc(n.module_id)}</td>` +
    `<td><span class="tbl-pill ${isAcc ? "tbl-acc" : "tbl-unc"}">${esc(n.verdict_status)}</span></td>` +
    `<td class="tbl-conf">${esc(n.confidence_label || "—")}</td>` +
    `<td>${esc(n.title || "—")}</td>`;
  tbody.appendChild(tr);
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

tick();'
}

# ── R public function ─────────────────────────────────────────

#' Generate an interactive hypothesis-level synthesis network HTML
#'
#' @param synthesis A \code{res$module_synthesis} list, or a featuremap result
#'   object containing \code{module_synthesis}.
#' @param output_dir Directory for the output HTML.
#' @param output_file File name for the output HTML.
#' @param show_same_module Ignored (same-module links are never rendered).
#' @return Invisibly, the written HTML path.
#' @export
generate_synthesis_network_html <- function(synthesis,
                                            output_dir  = getwd(),
                                            output_file = "synthesis_network.html",
                                            show_same_module = FALSE) {
  synthesis <- .syn_get_synthesis(synthesis)
  if (is.null(synthesis$nodes) || !is.data.frame(synthesis$nodes))
    stop("generate_synthesis_network_html() requires synthesis$nodes.")

  nodes <- synthesis$nodes
  relationships <- synthesis$relationships %||% synthesis$relationships_sparse
  if (is.null(relationships) || !is.data.frame(relationships)) {
    relationships <- data.frame(
      node_1 = character(), node_2 = character(),
      relation_class = character(), subtype = character(),
      rationale = character(), stringsAsFactors = FALSE
    )
  }

  nodes <- nodes[nodes$verdict_status %in% c("accepted","uncertain"), , drop = FALSE]
  vn    <- nodes$node_id
  relationships <- relationships[
    relationships$node_1 %in% vn &
    relationships$node_2 %in% vn &
    relationships$relation_class %in% c("direct","indirect"),
    , drop = FALSE
  ]

  graph_json <- jsonlite::toJSON(
    list(nodes   = .syn_json_records(nodes),
         edges   = .syn_json_records(relationships),
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

<div class="page-header">
  <h1>&#x1F9EC; Hypothesis-level Synthesis Network</h1>
  <div class="sub">Accepted &#x25CF; and uncertain &#x25A0; hypotheses &middot;
    direct and indirect biological relationships &middot;
    node size = confidence &middot; node colour = source module</div>
</div>

<div class="layout">

  <!-- network -->
  <div class="panel">
    <div class="toolbar">
      <label><input type="checkbox" id="fDirect"   checked> Direct edges</label>
      <label><input type="checkbox" id="fIndirect" checked> Indirect edges</label>
      <input type="text" id="search" placeholder="Search module, hypothesis, title&hellip;">
    </div>
    <div class="legend">
      <span><i class="leg-circle" style="background:#667085"></i>Accepted</span>
      <span><i class="leg-square" style="background:#667085;opacity:.6"></i>Uncertain</span>
      <i class="leg-sep"></i>
      <span>Size&nbsp;&rarr;&nbsp;confidence</span>
      <i class="leg-sep"></i>
      <span><i class="leg-edge-d"></i>Direct</span>
      <span><i class="leg-edge-i"></i>Indirect</span>
      <i class="leg-sep"></i>
      <div id="module-legend"></div>
    </div>
    <div class="network-wrap">
      <svg id="network"><g id="edges"></g><g id="nodes"></g></svg>
      <div id="tooltip"></div>
    </div>
  </div>

  <!-- summary -->
  <div class="panel sum-panel">
    <div class="sum-label">Integrated Summary</div>
    <div class="sum-text" id="sum-text"></div>
  </div>

</div>

<!-- node table -->
<div class="table-wrap">
  <table class="node-tbl">
    <thead>
      <tr>
        <th>Hypothesis</th>
        <th>Module</th>
        <th>Status</th>
        <th>Confidence</th>
        <th>Title</th>
      </tr>
    </thead>
    <tbody id="tbl-body"></tbody>
  </table>
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
#' @inheritParams generate_synthesis_network_html
#' @return Invisibly, the written HTML path.
#' @export
export_synthesis_network_html <- function(res,
                                          output_dir  = getwd(),
                                          output_file = "synthesis_network.html",
                                          show_same_module = FALSE) {
  if (is.null(res$module_synthesis))
    stop("res$module_synthesis is missing. Run synthesize_modules() first.")
  generate_synthesis_network_html(
    synthesis        = res$module_synthesis,
    output_dir       = output_dir,
    output_file      = output_file,
    show_same_module = show_same_module
  )
}
