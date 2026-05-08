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
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return(list())
  }
  unname(lapply(seq_len(nrow(df)), function(i) {
    row <- as.list(df[i, , drop = FALSE])
    lapply(row, function(x) .syn_scalar(x, ""))
  }))
}

.syn_css <- function() {
'*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Helvetica Neue",Arial,sans-serif;
     background:#f4f5f7;color:#1c1f26;font-size:14px;line-height:1.55}

/* header */
.page-header{background:#1c1f26;color:#fff;padding:20px 30px 18px}
.page-header h1{font-size:1.15rem;font-weight:600;letter-spacing:-.01em}
.page-header .sub{font-size:.8rem;color:#a8adba;margin-top:4px}

/* layout */
.layout{display:grid;grid-template-columns:minmax(520px,1fr) 340px;gap:16px;
        max-width:1380px;margin:0 auto;padding:18px}
.panel{background:#fff;border:1px solid #e4e7ed;border-radius:10px;overflow:hidden}

/* toolbar */
.toolbar{padding:11px 14px;display:flex;align-items:center;gap:10px;flex-wrap:wrap;
         border-bottom:1px solid #eef0f4}
.toolbar label{font-size:.78rem;color:#4b5563;display:flex;align-items:center;gap:5px}
.toolbar input[type="text"]{border:1px solid #d8dde6;border-radius:6px;padding:5px 9px;
                            min-width:200px;font-size:.82rem;outline:none}
.toolbar input[type="text"]:focus{border-color:#6741d9}

/* legend */
.legend{padding:9px 14px;display:flex;gap:14px;flex-wrap:wrap;color:#6b7280;font-size:.74rem;
        border-bottom:1px solid #eef0f4;align-items:center}
.legend span{display:inline-flex;align-items:center;gap:5px}
.leg-circle{width:11px;height:11px;border-radius:50%;background:#667085;display:inline-block}
.leg-square{width:11px;height:11px;border-radius:2px;background:#667085;display:inline-block}
.leg-edge-d{width:20px;height:2px;background:#516f9f;display:inline-block}
.leg-edge-i{width:20px;height:0;border-top:2px dashed #9a8064;display:inline-block}
.leg-sep{width:1px;height:12px;background:#d1d5db;display:inline-block}
#module-legend{display:flex;gap:10px;flex-wrap:wrap}
#module-legend span{display:inline-flex;align-items:center;gap:4px}
.mod-dot{width:9px;height:9px;border-radius:50%;display:inline-block;flex-shrink:0}

/* network */
.network-wrap{height:640px;position:relative;overflow:hidden}
svg{width:100%;height:100%;display:block;background:#fbfcfe}
.edge{stroke-opacity:.7}
.edge.direct{stroke:#516f9f;stroke-width:2}
.edge.indirect{stroke:#9a8064;stroke-width:1.6;stroke-dasharray:5 4}
.node{cursor:pointer}
.node text{font-size:10.5px;fill:#1f2937;paint-order:stroke;stroke:#fbfcfe;
           stroke-width:2.5px;user-select:none}

/* tooltip */
#tooltip{position:absolute;background:rgba(28,31,38,.88);color:#fff;
         padding:8px 11px;border-radius:7px;font-size:.76rem;max-width:230px;
         pointer-events:none;line-height:1.45;display:none;z-index:20;
         backdrop-filter:blur(2px)}
#tooltip strong{font-weight:600;display:block;margin-bottom:3px;font-size:.79rem}
#tooltip .tt-sum{opacity:.88;font-size:.74rem}

/* side panel */
.side{padding:0;display:flex;flex-direction:column;max-height:720px}
.side-top{padding:14px 16px 0;flex-shrink:0}
.side h2{font-size:.93rem;font-weight:600;margin-bottom:5px}
.side .hint{color:#8b909e;font-size:.8rem;margin-bottom:0}
.detail-scroll{flex:1;overflow-y:auto;padding:0 16px 16px}
.detail-section{border-top:1px solid #e4e7ed;padding-top:10px;margin-top:10px}
.detail-section:first-child{border-top:none;padding-top:0;margin-top:8px}
.detail-label{font-size:.65rem;font-weight:700;letter-spacing:.08em;text-transform:uppercase;
              color:#8b909e;margin-bottom:4px}
.detail-text{font-size:.82rem;color:#374151;white-space:pre-wrap;line-height:1.6}
.pill-row{display:flex;flex-wrap:wrap;gap:5px;margin-top:7px;margin-bottom:2px}
.pill{font-size:.72rem;padding:2px 8px;border-radius:4px;background:#f1f3f5;color:#4b5563;
      border:1px solid #e4e7ed}
.pill.accepted{background:#ecfdf5;color:#065f46;border-color:#a7f3d0}
.pill.uncertain{background:#fffbeb;color:#92400e;border-color:#fcd34d}

/* integrated summary */
.summary-panel{flex-shrink:0;margin:0 16px 16px;background:#f8f9fb;
               border:1px solid #e4e7ed;border-radius:8px;padding:11px 13px}
.summary-panel .detail-label{margin-bottom:6px}
.summary-panel .detail-text{font-size:.81rem;color:#374151}

@media(max-width:980px){.layout{grid-template-columns:1fr}
  .side{max-height:none}.network-wrap{height:520px}}'
}

.syn_js <- function() {
'const data = GRAPH_DATA;
const svg    = document.getElementById("network");
const gEdges = document.getElementById("edges");
const gNodes = document.getElementById("nodes");
const detail  = document.getElementById("detail");
const tooltip = document.getElementById("tooltip");
const search  = document.getElementById("search");
const controls = {
  direct:   document.getElementById("fDirect"),
  indirect: document.getElementById("fIndirect")
};

/* ── palette ── */
const modulePalette = [
  "#395f9d","#8a5a9d","#2d7a55","#b55a2a",
  "#7a6a2d","#2f7280","#8a3f58","#5b6f91",
  "#4a7c59","#7c4a6b","#6b4a7c","#4a6b7c"
];
const moduleIds    = [...new Set(data.nodes.map(d => d.module_id))].sort();
const moduleColor  = Object.fromEntries(
  moduleIds.map((m, i) => [m, modulePalette[i % modulePalette.length]])
);
const sizeByConf   = {high: 16, medium: 13, low: 10};
const defaultSize  = 12;

function nodeSize(n) { return sizeByConf[n.confidence_label] || defaultSize; }
function nodeColor(n) { return moduleColor[n.module_id] || "#667085"; }

/* ── build module legend ── */
const modLeg = document.getElementById("module-legend");
moduleIds.forEach(m => {
  const sp = document.createElement("span");
  sp.innerHTML = `<i class="mod-dot" style="background:${moduleColor[m]}"></i>${esc(m)}`;
  modLeg.appendChild(sp);
});

/* ── physics state ── */
let nodes = data.nodes.map((d, i) => ({
  ...d,
  x: 160 + (i % 7) * 85 + Math.random() * 40,
  y: 110 + Math.floor(i / 7) * 85 + Math.random() * 40,
  vx: 0, vy: 0
}));
let nodeById = Object.fromEntries(nodes.map(d => [d.node_id, d]));
let edges = data.edges
  .map(e => ({...e, source: nodeById[e.node_1], target: nodeById[e.node_2]}))
  .filter(e => e.source && e.target && e.relation_class !== "same_module");
let selectedId = null;

/* ── filter helpers ── */
function filters() {
  return {
    direct:   controls.direct.checked,
    indirect: controls.indirect.checked,
    q: search.value.trim().toLowerCase()
  };
}

function nodeVisible(n, f) {
  if (f.q) {
    const hay = [n.node_id, n.module_id, n.hypothesis_id,
                 n.title, n.functional_summary, n.phenotype_relationship]
                .join(" ").toLowerCase();
    if (!hay.includes(f.q)) return false;
  }
  return true;
}

function edgeVisible(e, visibleIds, f) {
  if (!visibleIds.has(e.node_1) || !visibleIds.has(e.node_2)) return false;
  if (e.relation_class === "direct")   return f.direct;
  if (e.relation_class === "indirect") return f.indirect;
  return false;
}

/* ── detail panel ── */
function showNode(n) {
  selectedId = n.node_id;
  const statusCls = n.verdict_status === "accepted" ? "accepted" : "uncertain";
  detail.innerHTML = `
    <div class="pill-row">
      <span class="pill ${statusCls}">${esc(n.verdict_status)}</span>
      <span class="pill">${esc(n.confidence_label || "confidence N/A")}</span>
      <span class="pill">${esc(n.module_id)}</span>
    </div>
    ${section("Title", n.title)}
    ${section("Functional summary", n.functional_summary)}
    ${section("Phenotype relationship", n.phenotype_relationship)}
    ${section("Confidence rationale", n.rationale)}
    ${section("Core drivers", n.core_drivers)}
    ${section("Top KEGG pathways", n.top_kegg_pathways)}
    ${section("Top HMDB pathways", n.top_hmdb_pathways)}
    ${section("Phenotype literature support",
      [n.phenotype_literature_support, n.phenotype_evidence_note].filter(Boolean).join("\\n"))}
  `;
  render();
}

function showEdge(e) {
  selectedId = null;
  detail.innerHTML = `
    <div class="pill-row">
      <span class="pill">${esc(e.relation_class)}</span>
      <span class="pill">${esc(e.subtype)}</span>
    </div>
    ${section("From → To", e.node_1 + " → " + e.node_2)}
    ${section("Rationale", e.rationale)}
  `;
  render();
}

function section(label, text) {
  if (!text) return "";
  return `<div class="detail-section">
    <div class="detail-label">${esc(label)}</div>
    <div class="detail-text">${esc(text)}</div>
  </div>`;
}

function esc(s) {
  return String(s ?? "").replace(/[&<>"]/g,
    c => ({"&":"&amp;","<":"&lt;",">":"&gt;",\'"\':"&quot;"}[c]));
}

/* ── tooltip ── */
const wrap = document.querySelector(".network-wrap");

function showTooltip(ev, n) {
  const title   = n.title || n.node_id;
  const summary = n.functional_summary || "";
  tooltip.innerHTML =
    `<strong>${esc(title)}</strong>` +
    (summary ? `<div class="tt-sum">${esc(summary.length > 110 ? summary.slice(0,110)+"…" : summary)}</div>` : "");
  tooltip.style.display = "block";
  moveTooltip(ev);
}

function moveTooltip(ev) {
  const r  = wrap.getBoundingClientRect();
  let tx = ev.clientX - r.left + 14;
  let ty = ev.clientY - r.top  + 14;
  if (tx + 240 > wrap.offsetWidth)  tx = ev.clientX - r.left - 244;
  if (ty + 90  > wrap.offsetHeight) ty = ev.clientY - r.top  - 94;
  tooltip.style.left = tx + "px";
  tooltip.style.top  = ty + "px";
}

/* ── physics ── */
function tick() {
  const rect = svg.getBoundingClientRect();
  const w = rect.width || 900, h = rect.height || 640;
  const f = filters();
  const visible = nodes.filter(n => nodeVisible(n, f));
  const visibleIds = new Set(visible.map(n => n.node_id));
  const visibleEdges = edges.filter(e => edgeVisible(e, visibleIds, f));

  for (const n of visible) {
    n.vx *= .85; n.vy *= .85;
    n.vx += (w / 2 - n.x) * .0008;
    n.vy += (h / 2 - n.y) * .0008;
  }
  for (let i = 0; i < visible.length; i++) {
    for (let j = i + 1; j < visible.length; j++) {
      const a = visible[i], b = visible[j];
      let dx = a.x - b.x, dy = a.y - b.y, d2 = dx*dx + dy*dy + 0.01;
      const rsum = nodeSize(a) + nodeSize(b) + 28;
      const force = Math.max(0, rsum*rsum*1.1 - d2) / (d2 + 0.01) * .8;
      const d = Math.sqrt(d2); dx /= d; dy /= d;
      a.vx += dx * force; a.vy += dy * force;
      b.vx -= dx * force; b.vy -= dy * force;
    }
  }
  for (const e of visibleEdges) {
    const a = e.source, b = e.target;
    let dx = b.x - a.x, dy = b.y - a.y, d = Math.sqrt(dx*dx + dy*dy) + 0.01;
    const target = e.relation_class === "direct" ? 160 : 210;
    const k = .005;
    const force = (d - target) * k; dx /= d; dy /= d;
    a.vx += dx * force; a.vy += dy * force;
    b.vx -= dx * force; b.vy -= dy * force;
  }
  for (const n of visible) {
    const pad = nodeSize(n) + 8;
    n.x = Math.max(pad, Math.min(w - pad, n.x + n.vx));
    n.y = Math.max(pad, Math.min(h - pad, n.y + n.vy));
  }
  render();
  requestAnimationFrame(tick);
}

/* ── render ── */
function render() {
  const f = filters();
  const visible = nodes.filter(n => nodeVisible(n, f));
  const visibleIds = new Set(visible.map(n => n.node_id));
  const visibleEdges = edges.filter(e => edgeVisible(e, visibleIds, f));

  /* edges */
  gEdges.innerHTML = "";
  for (const e of visibleEdges) {
    const line = document.createElementNS("http://www.w3.org/2000/svg","line");
    line.setAttribute("x1", e.source.x); line.setAttribute("y1", e.source.y);
    line.setAttribute("x2", e.target.x); line.setAttribute("y2", e.target.y);
    line.setAttribute("class", `edge ${e.relation_class}`);
    line.style.cursor = "pointer";
    line.addEventListener("click", () => showEdge(e));
    gEdges.appendChild(line);
  }

  /* nodes */
  gNodes.innerHTML = "";
  for (const n of visible) {
    const grp = document.createElementNS("http://www.w3.org/2000/svg","g");
    grp.setAttribute("class", `node ${selectedId === n.node_id ? "selected" : ""}`);
    grp.setAttribute("transform", `translate(${n.x},${n.y})`);

    const sz    = nodeSize(n);
    const fill  = nodeColor(n);
    const isAcc = n.verdict_status === "accepted";

    /* shape: circle = accepted, rect = uncertain */
    let shape;
    if (isAcc) {
      shape = document.createElementNS("http://www.w3.org/2000/svg","circle");
      shape.setAttribute("r", sz);
    } else {
      shape = document.createElementNS("http://www.w3.org/2000/svg","rect");
      shape.setAttribute("x", -sz); shape.setAttribute("y", -sz);
      shape.setAttribute("width", sz * 2); shape.setAttribute("height", sz * 2);
      shape.setAttribute("rx", "3"); shape.setAttribute("ry", "3");
    }
    shape.setAttribute("fill", fill);
    shape.setAttribute("fill-opacity", isAcc ? "0.88" : "0.6");
    shape.setAttribute("stroke", selectedId === n.node_id ? "#1c1f26" : "#fff");
    shape.setAttribute("stroke-width", selectedId === n.node_id ? "3" : "2");

    /* label */
    const t = document.createElementNS("http://www.w3.org/2000/svg","text");
    t.setAttribute("x", sz + 4); t.setAttribute("y", 4);
    t.textContent = n.node_id;

    grp.appendChild(shape); grp.appendChild(t);

    /* events */
    grp.addEventListener("click",      ev => { ev.stopPropagation(); showNode(n); });
    grp.addEventListener("mouseenter", ev => showTooltip(ev, n));
    grp.addEventListener("mousemove",  ev => moveTooltip(ev));
    grp.addEventListener("mouseleave", ()  => { tooltip.style.display = "none"; });

    gNodes.appendChild(grp);
  }
}

/* ── init ── */
for (const el of Object.values(controls)) el.addEventListener("change", render);
search.addEventListener("input", render);
svg.addEventListener("click", () => { selectedId = null; render(); });

if (data.summary) {
  document.getElementById("summary-text").textContent = data.summary;
} else {
  document.getElementById("summary-panel").style.display = "none";
}
tick();'
}

#' Generate an interactive hypothesis-level synthesis network HTML
#'
#' @param synthesis A \code{res$module_synthesis} list, or a featuremap result
#'   object containing \code{module_synthesis}.
#' @param output_dir Directory for the output HTML.
#' @param output_file File name for the output HTML.
#' @param show_same_module Logical; ignored (same-module links not shown).
#' @return Invisibly, the written HTML path.
#' @export
generate_synthesis_network_html <- function(synthesis,
                                            output_dir = getwd(),
                                            output_file = "synthesis_network.html",
                                            show_same_module = FALSE) {
  synthesis <- .syn_get_synthesis(synthesis)
  if (is.null(synthesis$nodes) || !is.data.frame(synthesis$nodes)) {
    stop("generate_synthesis_network_html() requires synthesis$nodes.")
  }

  nodes <- synthesis$nodes
  relationships <- synthesis$relationships
  if (is.null(relationships) || !is.data.frame(relationships)) {
    relationships <- synthesis$relationships_sparse
  }
  if (is.null(relationships) || !is.data.frame(relationships)) {
    relationships <- data.frame(
      node_1 = character(), node_2 = character(), relation_class = character(),
      subtype = character(), rationale = character(),
      stringsAsFactors = FALSE
    )
  }

  nodes <- nodes[nodes$verdict_status %in% c("accepted", "uncertain"), , drop = FALSE]
  valid_nodes <- nodes$node_id
  # exclude same_module edges entirely
  relationships <- relationships[
    relationships$node_1 %in% valid_nodes &
      relationships$node_2 %in% valid_nodes &
      relationships$relation_class %in% c("direct", "indirect"),
    ,
    drop = FALSE
  ]

  graph <- list(
    nodes   = .syn_json_records(nodes),
    edges   = .syn_json_records(relationships),
    summary = .syn_scalar(synthesis$summary, "")
  )
  graph_json <- jsonlite::toJSON(graph, auto_unbox = TRUE, dataframe = "rows", na = "null")
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
  <div class="sub">Accepted and uncertain hypotheses &middot; direct and indirect biological relationships</div>
</div>
<div class="layout">

  <!-- network panel -->
  <div class="panel">
    <div class="toolbar">
      <label><input type="checkbox" id="fDirect" checked> Direct edges</label>
      <label><input type="checkbox" id="fIndirect" checked> Indirect edges</label>
      <input type="text" id="search" placeholder="Search module, hypothesis, title&hellip;">
    </div>
    <div class="legend">
      <span><i class="leg-circle"></i> Accepted</span>
      <span><i class="leg-square"></i> Uncertain</span>
      <span style="opacity:.4"><i class="leg-sep"></i></span>
      <span>Size &rarr; confidence (large = high)</span>
      <span style="opacity:.4"><i class="leg-sep"></i></span>
      <span><i class="leg-edge-d"></i> Direct</span>
      <span><i class="leg-edge-i"></i> Indirect</span>
      <span style="opacity:.4"><i class="leg-sep"></i></span>
      <span>Color &rarr; module:</span>
      <div id="module-legend"></div>
    </div>
    <div class="network-wrap">
      <svg id="network"><g id="edges"></g><g id="nodes"></g></svg>
      <div id="tooltip"></div>
    </div>
  </div>

  <!-- detail panel -->
  <div class="panel side">
    <div class="side-top">
      <h2>Details</h2>
      <p class="hint">Click a node or edge to inspect.</p>
    </div>
    <div class="detail-scroll" id="detail"></div>
    <div class="summary-panel" id="summary-panel">
      <div class="detail-label">Integrated summary</div>
      <div class="detail-text" id="summary-text"></div>
    </div>
  </div>

</div>
<script>%s</script>
</body>
</html>',
    .syn_css(), js
  )

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
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
                                          output_dir = getwd(),
                                          output_file = "synthesis_network.html",
                                          show_same_module = FALSE) {
  if (is.null(res$module_synthesis)) {
    stop("res$module_synthesis is missing. Run synthesize_modules() first.")
  }
  generate_synthesis_network_html(
    synthesis = res$module_synthesis,
    output_dir = output_dir,
    output_file = output_file,
    show_same_module = show_same_module
  )
}
