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
.page-header{background:#1c1f26;color:#fff;padding:20px 30px 18px}
.page-header h1{font-size:1.15rem;font-weight:600;letter-spacing:-.01em}
.page-header .sub{font-size:.8rem;color:#a8adba;margin-top:4px}
.layout{display:grid;grid-template-columns:minmax(520px,1fr) 360px;gap:16px;
        max-width:1380px;margin:0 auto;padding:18px}
.panel{background:#fff;border:1px solid #e4e7ed;border-radius:10px}
.toolbar{padding:12px 14px;display:flex;align-items:center;gap:10px;flex-wrap:wrap}
.toolbar label{font-size:.78rem;color:#4b5563;display:flex;align-items:center;gap:5px}
.toolbar input[type="text"]{border:1px solid #d8dde6;border-radius:6px;padding:6px 8px;
                            min-width:190px;font-size:.82rem}
.legend{padding:0 14px 12px;display:flex;gap:12px;flex-wrap:wrap;color:#6b7280;font-size:.75rem}
.legend span{display:inline-flex;align-items:center;gap:5px}
.dot{width:10px;height:10px;border-radius:50%;display:inline-block}
.network-wrap{height:650px;border-top:1px solid #e4e7ed;position:relative;overflow:hidden}
svg{width:100%;height:100%;display:block;background:#fbfcfe}
.edge{stroke-opacity:.75}
.edge.direct{stroke:#516f9f;stroke-width:2.2}
.edge.indirect{stroke:#9a8064;stroke-width:1.8;stroke-dasharray:5 5}
.edge.same_module{stroke:#c4c9d2;stroke-width:1.2;stroke-dasharray:3 5}
.node circle{cursor:pointer}
.node text{font-size:11px;fill:#1f2937;paint-order:stroke;stroke:#fff;stroke-width:3px}
.node.selected circle{stroke-width:4px}
.side{padding:15px 16px;min-height:650px}
.side h2{font-size:.95rem;margin-bottom:8px}
.side .hint{color:#8b909e;font-size:.82rem;margin-bottom:12px}
.detail-section{border-top:1px solid #e4e7ed;padding-top:10px;margin-top:10px}
.detail-label{font-size:.66rem;font-weight:700;letter-spacing:.08em;text-transform:uppercase;
              color:#8b909e;margin-bottom:4px}
.detail-text{font-size:.83rem;color:#374151;white-space:pre-wrap}
.pill-row{display:flex;flex-wrap:wrap;gap:5px;margin-top:6px}
.pill{font-size:.72rem;padding:2px 7px;border-radius:4px;background:#f1f3f5;color:#4b5563}
.summary{font-size:.83rem;color:#374151;background:#f8f9fb;border:1px solid #e4e7ed;
         padding:10px 12px;border-radius:8px;margin-top:10px}
.node-table{font-size:.76rem;margin-top:10px;color:#4b5563}
.node-table div{padding:5px 0;border-bottom:1px solid #eef0f4}
@media(max-width:980px){.layout{grid-template-columns:1fr}.side{min-height:auto}.network-wrap{height:560px}}'
}

.syn_js <- function() {
'const data = GRAPH_DATA;
const svg = document.getElementById("network");
const gEdges = document.getElementById("edges");
const gNodes = document.getElementById("nodes");
const detail = document.getElementById("detail");
const search = document.getElementById("search");
const controls = {
  accepted: document.getElementById("fAccepted"),
  uncertain: document.getElementById("fUncertain"),
  direct: document.getElementById("fDirect"),
  indirect: document.getElementById("fIndirect"),
  same: document.getElementById("fSame")
};

const modulePalette = ["#395f9d","#8a5a9d","#2d7a55","#a7652a","#7a6a2d","#2f7280","#8a3f58","#5b6f91"];
const moduleIds = [...new Set(data.nodes.map(d => d.module_id))];
const moduleColor = Object.fromEntries(moduleIds.map((m,i) => [m, modulePalette[i % modulePalette.length]]));
const fillByConfidence = {high:"#d1fae5", medium:"#fef9c3", low:"#fee2e2", "":"#eef2f7"};

let nodes = data.nodes.map((d,i) => ({
  ...d,
  x: 160 + (i % 6) * 95 + Math.random() * 30,
  y: 130 + Math.floor(i / 6) * 90 + Math.random() * 30,
  vx: 0, vy: 0
}));
let nodeById = Object.fromEntries(nodes.map(d => [d.node_id, d]));
let edges = data.edges.map(e => ({...e, source: nodeById[e.node_1], target: nodeById[e.node_2]}))
                      .filter(e => e.source && e.target);
let selectedId = null;

function filters(){
  return {
    accepted: controls.accepted.checked,
    uncertain: controls.uncertain.checked,
    direct: controls.direct.checked,
    indirect: controls.indirect.checked,
    same: controls.same.checked,
    q: search.value.trim().toLowerCase()
  };
}

function nodeVisible(n, f){
  if (n.verdict_status === "accepted" && !f.accepted) return false;
  if (n.verdict_status === "uncertain" && !f.uncertain) return false;
  if (f.q) {
    const hay = [n.node_id,n.module_id,n.hypothesis_id,n.title,n.functional_summary,n.phenotype_relationship].join(" ").toLowerCase();
    if (!hay.includes(f.q)) return false;
  }
  return true;
}

function edgeVisible(e, visibleIds, f){
  if (!visibleIds.has(e.node_1) || !visibleIds.has(e.node_2)) return false;
  if (e.same_module === "TRUE" || e.same_module === true || e.relation_class === "same_module") return f.same;
  if (e.relation_class === "direct") return f.direct;
  if (e.relation_class === "indirect") return f.indirect;
  return false;
}

function wrapText(txt, max=68){
  txt = txt || "";
  if (txt.length <= max) return txt;
  return txt.slice(0, max - 1) + "…";
}

function showNode(n){
  selectedId = n.node_id;
  detail.innerHTML = `
    <h2>${esc(n.node_id)}</h2>
    <div class="pill-row">
      <span class="pill">${esc(n.verdict_status)}</span>
      <span class="pill">${esc(n.confidence_label || "confidence NA")}</span>
      <span class="pill">${esc(n.module_id)}</span>
    </div>
    ${section("Title", n.title)}
    ${section("Functional summary", n.functional_summary)}
    ${section("Phenotype relationship", n.phenotype_relationship)}
    ${section("Confidence rationale", n.rationale)}
    ${section("Core drivers", n.core_drivers)}
    ${section("Top KEGG pathways", n.top_kegg_pathways)}
    ${section("Top HMDB pathways", n.top_hmdb_pathways)}
    ${section("Phenotype literature support", [n.phenotype_literature_support, n.phenotype_evidence_note].filter(Boolean).join("\\n"))}
  `;
  render();
}

function showEdge(e){
  selectedId = null;
  detail.innerHTML = `
    <h2>${esc(e.node_1)} → ${esc(e.node_2)}</h2>
    <div class="pill-row">
      <span class="pill">${esc(e.relation_class)}</span>
      <span class="pill">${esc(e.subtype)}</span>
    </div>
    ${section("Rationale", e.rationale || "Same source module, not inferred biological relationship.")}
  `;
  render();
}

function section(label, text){
  if (!text) return "";
  return `<div class="detail-section"><div class="detail-label">${esc(label)}</div><div class="detail-text">${esc(text)}</div></div>`;
}
function esc(s){return String(s ?? "").replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;","\\"":"&quot;"}[c]));}

function tick(){
  const rect = svg.getBoundingClientRect();
  const w = rect.width || 900, h = rect.height || 650;
  const f = filters();
  const visible = nodes.filter(n => nodeVisible(n, f));
  const visibleIds = new Set(visible.map(n => n.node_id));
  const visibleEdges = edges.filter(e => edgeVisible(e, visibleIds, f));

  for (const n of visible) {
    n.vx *= .86; n.vy *= .86;
    n.vx += (w/2 - n.x) * .0009;
    n.vy += (h/2 - n.y) * .0009;
  }
  for (let i=0;i<visible.length;i++){
    for (let j=i+1;j<visible.length;j++){
      const a=visible[i], b=visible[j];
      let dx=a.x-b.x, dy=a.y-b.y, d2=dx*dx+dy*dy+0.01;
      const force = 950 / d2;
      const d=Math.sqrt(d2); dx/=d; dy/=d;
      a.vx += dx*force; a.vy += dy*force;
      b.vx -= dx*force; b.vy -= dy*force;
    }
  }
  for (const e of visibleEdges) {
    const a=e.source, b=e.target;
    let dx=b.x-a.x, dy=b.y-a.y, d=Math.sqrt(dx*dx+dy*dy)+0.01;
    const target = (e.relation_class === "same_module") ? 110 : (e.relation_class === "direct" ? 150 : 190);
    const k = (e.relation_class === "same_module") ? .002 : .006;
    const force=(d-target)*k; dx/=d; dy/=d;
    a.vx += dx*force; a.vy += dy*force;
    b.vx -= dx*force; b.vy -= dy*force;
  }
  for (const n of visible) {
    n.x = Math.max(30, Math.min(w-30, n.x+n.vx));
    n.y = Math.max(30, Math.min(h-30, n.y+n.vy));
  }
  render();
  requestAnimationFrame(tick);
}

function render(){
  const f = filters();
  const visible = nodes.filter(n => nodeVisible(n, f));
  const visibleIds = new Set(visible.map(n => n.node_id));
  const visibleEdges = edges.filter(e => edgeVisible(e, visibleIds, f));

  gEdges.innerHTML = "";
  for (const e of visibleEdges) {
    const line = document.createElementNS("http://www.w3.org/2000/svg","line");
    line.setAttribute("x1", e.source.x); line.setAttribute("y1", e.source.y);
    line.setAttribute("x2", e.target.x); line.setAttribute("y2", e.target.y);
    line.setAttribute("class", `edge ${e.relation_class}`);
    line.addEventListener("click", () => showEdge(e));
    gEdges.appendChild(line);
  }

  gNodes.innerHTML = "";
  for (const n of visible) {
    const grp = document.createElementNS("http://www.w3.org/2000/svg","g");
    grp.setAttribute("class", `node ${selectedId === n.node_id ? "selected" : ""}`);
    grp.setAttribute("transform", `translate(${n.x},${n.y})`);
    const c = document.createElementNS("http://www.w3.org/2000/svg","circle");
    c.setAttribute("r", n.verdict_status === "accepted" ? 13 : 11);
    c.setAttribute("fill", fillByConfidence[n.confidence_label] || fillByConfidence[""]);
    c.setAttribute("stroke", moduleColor[n.module_id] || "#667085");
    c.setAttribute("stroke-width", n.verdict_status === "accepted" ? 3 : 2);
    if (n.verdict_status === "uncertain") c.setAttribute("stroke-dasharray","4 3");
    c.addEventListener("click", () => showNode(n));
    const t = document.createElementNS("http://www.w3.org/2000/svg","text");
    t.setAttribute("x", 17); t.setAttribute("y", 4);
    t.textContent = n.node_id;
    grp.appendChild(c); grp.appendChild(t);
    gNodes.appendChild(grp);
  }
}

for (const el of Object.values(controls)) el.addEventListener("change", render);
search.addEventListener("input", render);
document.getElementById("nodeList").innerHTML = nodes.map(n => `<div><strong>${esc(n.node_id)}</strong><br>${esc(wrapText(n.title, 82))}</div>`).join("");
if (data.summary) document.getElementById("summary").innerHTML = esc(data.summary);
tick();'
}

#' Generate an interactive hypothesis-level synthesis network HTML
#'
#' @param synthesis A \code{res$module_synthesis} list, or a featuremap result
#'   object containing \code{module_synthesis}.
#' @param output_dir Directory for the output HTML.
#' @param output_file File name for the output HTML.
#' @param show_same_module Logical; whether same-module links are shown by
#'   default. They can also be toggled in the HTML.
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
      subtype = character(), rationale = character(), same_module = logical(),
      stringsAsFactors = FALSE
    )
  }

  nodes <- nodes[nodes$verdict_status %in% c("accepted", "uncertain"), , drop = FALSE]
  valid_nodes <- nodes$node_id
  relationships <- relationships[
    relationships$node_1 %in% valid_nodes &
      relationships$node_2 %in% valid_nodes &
      relationships$relation_class %in% c("direct", "indirect", "same_module"),
    ,
    drop = FALSE
  ]

  graph <- list(
    nodes = .syn_json_records(nodes),
    edges = .syn_json_records(relationships),
    summary = .syn_scalar(synthesis$summary, ""),
    show_same_module = isTRUE(show_same_module)
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
  <h1>Hypothesis-level synthesis network</h1>
  <div class="sub">Accepted and uncertain hypotheses; direct and indirect biological relationships.</div>
</div>
<div class="layout">
  <div class="panel">
    <div class="toolbar">
      <label><input type="checkbox" id="fAccepted" checked> accepted</label>
      <label><input type="checkbox" id="fUncertain" checked> uncertain</label>
      <label><input type="checkbox" id="fDirect" checked> direct</label>
      <label><input type="checkbox" id="fIndirect" checked> indirect</label>
      <label><input type="checkbox" id="fSame" %s> same-module links</label>
      <input type="text" id="search" placeholder="Search module, hypothesis, title...">
    </div>
    <div class="legend">
      <span><i class="dot" style="background:#d1fae5"></i>high confidence</span>
      <span><i class="dot" style="background:#fef9c3"></i>medium</span>
      <span><i class="dot" style="background:#fee2e2"></i>low</span>
      <span>solid edge: direct</span>
      <span>dashed edge: indirect</span>
      <span>node outline: source module</span>
    </div>
    <div class="network-wrap">
      <svg id="network"><g id="edges"></g><g id="nodes"></g></svg>
    </div>
  </div>
  <div class="panel side">
    <h2>Details</h2>
    <div class="hint">Click a node or edge to inspect the synthesis evidence.</div>
    <div id="detail">
      <div class="detail-section">
        <div class="detail-label">Integrated summary</div>
        <div id="summary" class="summary"></div>
      </div>
      <div class="detail-section">
        <div class="detail-label">Nodes</div>
        <div id="nodeList" class="node-table"></div>
      </div>
    </div>
  </div>
</div>
<script>%s</script>
</body>
</html>',
    .syn_css(),
    if (isTRUE(show_same_module)) "checked" else "",
    js
  )

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(output_dir, output_file)
  writeLines(html, out_path, useBytes = FALSE)
  message(sprintf("[HTML] Synthesis network → %s", out_path))
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
