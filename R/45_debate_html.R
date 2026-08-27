# ============================================================
# 45_debate_html.R
# Debate HTML visualizer — academic chat-style layout
#
# Sourced automatically by 45_advanced_mode.R.
# Public export: generate_debate_html()
# ============================================================


# ── Internal helpers ──────────────────────────────────────────

.he <- function(s) {
  s <- as.character(s %||% "")
  s <- gsub("&",  "&amp;",  s, fixed = TRUE)
  s <- gsub("<",  "&lt;",   s, fixed = TRUE)
  s <- gsub(">",  "&gt;",   s, fixed = TRUE)
  s <- gsub('"',  "&quot;", s, fixed = TRUE)
  s
}

#' Coerce any list/vector field to a single character scalar (prevents
#' "condition has length > 1" when JSON parsing returns unexpected vectors)
.s <- function(x, default = "") {
  x <- x %||% default
  x <- as.character(x)
  if (length(x) == 0L) return(as.character(default)[1L])
  x <- trimws(x[1L])
  if (is.na(x) || !nzchar(x) || toupper(x) %in% c("NA", "NULL")) {
    return(as.character(default)[1L])
  }
  x
}

.is_nonempty <- function(x) {
  x <- .s(x, "")
  nzchar(x)
}

#' Truncate long string(s), appending ellipsis.  Vectorised — safe for
#' character vectors (e.g. when called from .met_chips() with multiple names).
.trunc <- function(s, n) {
  s <- as.character(s %||% "")
  ifelse(nchar(s) > n, paste0(substr(s, 1L, n - 1L), "…"), s)
}

#' Render all metabolite chips
.met_chips <- function(names) {
  names <- unique(trimws(as.character(names %||% character(0))))
  names <- names[nzchar(names) & !is.na(names) & toupper(names) != "NA"]
  if (length(names) == 0L) return('<span class="muted">—</span>')
  chips <- paste(
    sprintf('<span class="met-chip" title="%s">%s</span>',
            .he(names), .he(names)),
    collapse = "")
  chips
}

#' Render top class profile from precomputed Class-level structure
.class_bars <- function(module_obj, n = 3L) {
  df <- module_obj$other_metabolites
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L ||
      !"class" %in% colnames(df) || !"n_metabolites" %in% colnames(df)) {
    return('<span class="muted">—</span>')
  }

  cls <- trimws(as.character(df$class))
  keep <- !is.na(cls) & nzchar(cls) & !toupper(cls) %in% c("NA", "NULL", "UNKNOWN")
  df <- df[keep, , drop = FALSE]
  cls <- cls[keep]
  if (nrow(df) == 0L) return('<span class="muted">—</span>')
  df$class <- cls
  df$n_metabolites <- suppressWarnings(as.numeric(df$n_metabolites))
  df <- df[is.finite(df$n_metabolites) & df$n_metabolites > 0, , drop = FALSE]
  if (nrow(df) == 0L) return('<span class="muted">—</span>')

  df <- df[order(-df$n_metabolites), , drop = FALSE]
  rows <- head(df, n)
  bars <- paste(vapply(seq_len(nrow(rows)), function(i) {
    pct <- round(rows$n_metabolites[[i]] * 100)
    lbl <- rows$class[[i]]
    sprintf(
      '<div class="cls-row">
        <div class="cls-name">%s</div>
        <div class="cls-track"><div class="cls-fill" style="width:%d%%"></div></div>
        <div class="cls-pct">%d%%</div>
      </div>',
      .he(lbl), pct, pct)
  }, character(1L)), collapse = "\n")

  sprintf('<div class="cls-block">%s</div>', bars)
}

#' Render anchor / tension metabolite tags (subtle, dot-prefixed)
.met_tags <- function(csv, variant) {
  # variant: "anc" | "ten"
  vals <- trimws(strsplit(as.character(csv %||% ""), ",")[[1L]])
  vals <- vals[nzchar(vals)]
  if (length(vals) == 0L) return("")
  label <- if (identical(variant, "anc")) "Core anchors" else "Tension metabolites"
  tags  <- paste(
    sprintf('<span class="mtag mtag-%s">%s</span>', variant, .he(vals)),
    collapse = "")
  sprintf('<div class="mtag-row"><span class="mtag-label">%s</span>%s</div>',
          label, tags)
}

#' Render one hypothesis card
.hyp_card <- function(h) {
  anc <- paste(as.character(h$anchor_metabolites  %||% character(0)), collapse = ", ")
  ten <- paste(as.character(h$tension_metabolites %||% character(0)), collapse = ", ")
  typ <- .s(h$type, "")
  type_html <- if (nzchar(typ))
    sprintf('<span class="type-pill">%s</span>', .he(typ)) else ""
  sprintf(
    '<div class="hyp-card">
      <div class="hyp-header">
        <span class="hyp-id">%s</span>
        <span class="hyp-title">%s</span>
        %s
      </div>
      <div class="hyp-summary">%s</div>
      %s%s
    </div>',
    .he(.s(h$hypothesis_id)),
    .he(.s(h$title)),
    type_html,
    .he(.s(h$summary)),
    .met_tags(anc, "anc"),
    .met_tags(ten, "ten")
  )
}

#' Render one Skeptic challenge
.ch_item <- function(ch) {
  sev   <- .s(ch$severity,       "minor")
  ctype <- .s(ch$challenge_type, "")
  if (identical(sev, "no_challenge")) return("")
  sev_cls <- if (identical(sev, "major")) "sev-major" else "sev-minor"
  sprintf(
    '<div class="ch-item ch-%s">
      <div class="ch-meta">
        <span class="sev-pill %s">%s</span>
        <span class="ch-type">%s</span>
        <span class="ch-target">%s</span>
      </div>
      <div class="ch-text">%s</div>
    </div>',
    sev, sev_cls, .he(toupper(sev)),
    .he(ctype),
    .he(.s(ch$target_hypothesis_id)),
    .he(.s(ch$challenge))
  )
}

#' Render one Proposer revision
.rev_item <- function(rev, hypotheses_after) {
  hid    <- .s(rev$hypothesis_id)
  action <- toupper(.s(rev$action, "keep"))
  just   <- .s(rev$justification)
  h_new  <- Filter(function(x) x$hypothesis_id == hid, hypotheses_after)

  act_cls <- switch(action, MODIFY = "act-modify", DROP = "act-drop", "act-keep")

  detail <- ""
  if (identical(action, "MODIFY") && length(h_new) > 0L) {
    h   <- h_new[[1L]]
    anc <- paste(as.character(h$anchor_metabolites  %||% character(0)), collapse = ", ")
    ten <- paste(as.character(h$tension_metabolites %||% character(0)), collapse = ", ")
    typ <- .s(h$type, "")
    type_html <- if (nzchar(typ))
      sprintf('<div class="rev-new-type">%s</div>', .he(typ)) else ""
    detail <- sprintf(
      '<div class="rev-detail">
        %s
        <div class="rev-new-title">%s</div>
        <div class="rev-new-summary">%s</div>
        %s%s
      </div>',
      type_html,
      .he(.s(h$title)),
      .he(.s(h$summary)),
      .met_tags(anc, "anc"),
      .met_tags(ten, "ten")
    )
  }
  sprintf(
    '<div class="rev-item">
      <div class="rev-header">
        <span class="rev-hid">%s</span>
        <span class="act-pill %s">%s</span>
      </div>
      <div class="rev-just">%s</div>
      %s
    </div>',
    .he(hid), act_cls, action, .he(just), detail
  )
}

#' Render one Arbiter verdict card
.verdict_card <- function(v, hypotheses_final) {
  hid    <- .s(v$hypothesis_id)
  status <- .s(v$status,               "rejected")
  conf   <- .s(v$confidence_label,     "")
  vrat   <- .s(v$verdict_rationale,    "")
  crat   <- .s(v$confidence_rationale, "")
  st_cls <- switch(status,
                   accepted  = "st-accepted",
                   "st-rejected")
  h_fin       <- Filter(function(x) x$hypothesis_id == hid, hypotheses_final %||% list())
  fin_title   <- if (length(h_fin) > 0L) .s(h_fin[[1L]]$title)   else ""
  fin_summary <- if (length(h_fin) > 0L) .s(h_fin[[1L]]$summary) else ""
  fin_type    <- if (length(h_fin) > 0L) .s(h_fin[[1L]]$type)    else ""

  status_pill <- sprintf(
    '<span class="status-pill %s">%s</span>', st_cls, .he(toupper(status)))
  conf_pill <- if (nzchar(conf))
    sprintf('<span class="conf-pill conf-%s">%s confidence</span>', .he(conf), .he(conf))
  else ""

  title_html <- if (nzchar(fin_title))
    sprintf('<div class="verdict-title">%s</div>', .he(fin_title)) else ""
  type_html <- if (nzchar(fin_type))
    sprintf('<div class="verdict-type">%s</div>', .he(fin_type)) else ""
  summary_html <- if (nzchar(fin_summary))
    sprintf('<div class="verdict-fin-summary">%s</div>', .he(fin_summary)) else ""

  # Verdict rationale block
  vrat_html <- if (nzchar(vrat)) sprintf(
    '<div class="rat-block">
      <div class="rat-label">Verdict rationale</div>
      <div class="rat-text">%s</div>
    </div>', .he(vrat)) else ""

  # Confidence rationale block — with divider when verdict block also present
  crat_html <- if (nzchar(crat)) sprintf(
    '%s<div class="rat-block">
      <div class="rat-label">Confidence rationale
        %s
      </div>
      <div class="rat-text">%s</div>
    </div>',
    if (nzchar(vrat)) '<div class="rat-divider"></div>' else "",
    conf_pill,
    .he(crat)
  ) else ""

  sprintf(
    '<div class="verdict-card %s">
      <div class="verdict-hdr">
        <span class="verdict-hid">%s</span>
        %s
      </div>
      %s%s%s
      %s%s
    </div>',
    st_cls,
    .he(hid), status_pill,
    type_html, title_html, summary_html,
    vrat_html, crat_html
  )
}

#' Render one PubMed phenotype-evidence card
.pubmed_card <- function(hid, rel = NULL, evidence = NULL) {
  support <- .s(rel$literature_support, "insufficient")
  relationship <- .s(rel$phenotype_relationship, "")
  note <- .s(rel$evidence_note, "")

  pmids <- as.character(evidence$pmids %||% character(0))
  pmids <- unique(trimws(pmids))
  pmids <- pmids[nzchar(pmids) & !is.na(pmids) & toupper(pmids) != "NA"]

  pmid_html <- ""
  if (length(pmids) > 0L) {
    pmid_html <- paste(
      sprintf('<span class="pmid-chip">PMID:%s</span>', .he(pmids)),
      collapse = ""
    )
  } else {
    pmid_html <- '<span class="pm-empty">No PubMed records retrieved.</span>'
  }

  sprintf(
    '<div class="pubmed-card">
      <div class="pubmed-hdr">
        <span class="pubmed-hid">%s</span>
        <span class="lit-pill lit-%s">%s literature support</span>
      </div>
      %s
      %s
      <div class="pmid-chips">%s</div>
    </div>',
    .he(.s(hid)),
    .he(support),
    .he(support),
    if (.is_nonempty(relationship))
      sprintf('<div class="pm-relationship"><strong>Phenotype relationship:</strong> %s</div>', .he(relationship))
    else "",
    if (.is_nonempty(note))
      sprintf('<div class="pm-note"><strong>Evidence note:</strong> %s</div>', .he(note))
    else "",
    pmid_html
  )
}

#' Count revision actions (safe for NULL/empty)
.count_action <- function(revs, action) {
  if (is.null(revs) || length(revs) == 0L) return(0L)
  sum(vapply(revs, function(r) identical(r$action, action), logical(1L)))
}

#' Wrap HTML in a chat bubble
.bubble <- function(role, avatar, label, content, extra_cls = "") {
  sprintf(
    '<div class="msg msg-%s %s">
      <div class="av av-%s" title="%s">%s</div>
      <div class="bub bub-%s">%s</div>
    </div>',
    role, extra_cls,
    role, .he(label), .he(avatar),
    role, content
  )
}

#' Full inline CSS
.debate_css <- function() {
'/* reset */
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}

/* page */
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Helvetica Neue",Arial,sans-serif;
     background:#f4f5f7;color:#1c1f26;font-size:14px;line-height:1.65}

/* header */
.page-header{background:#519d9e;color:#fff;padding:20px 32px 18px}
.page-header h1{font-size:1.15rem;font-weight:600;letter-spacing:-.01em}
.page-header .sub{font-size:.78rem;color:rgba(255,255,255,.82);margin-top:3px;display:flex;gap:14px;flex-wrap:wrap}
.page-header .sub span{display:flex;align-items:center;gap:4px}

/* layout */
.chat-body{max-width:900px;margin:0 auto;padding:28px 20px;
           display:flex;flex-direction:column;gap:18px}

/* bubbles */
.msg{display:flex;gap:12px;align-items:flex-start}
.msg-left {flex-direction:row}
.msg-right{flex-direction:row-reverse}
.msg-center{flex-direction:column;align-items:stretch}

.av{width:36px;height:36px;border-radius:50%;display:flex;align-items:center;
    justify-content:center;font-size:15px;flex-shrink:0;
    border:1px solid rgba(0,0,0,.08)}
.av-left  {background:#eef4fb}
.av-right {background:#e1f5ee}
.av-center{background:#f6f9fc;width:36px;height:36px;border-radius:50%;
           display:flex;align-items:center;justify-content:center;
           font-size:15px;flex-shrink:0;margin:0 auto 10px;
           border:1px solid rgba(0,0,0,.08)}

.bub{background:#fff;border-radius:12px;padding:14px 16px;
     max-width:84%;border:1px solid #dbe5f2}
.bub-left {border-top-left-radius:3px}
.bub-right{border-top-right-radius:3px}
.bub-center{max-width:100%;border-radius:12px}
.msg-center .bub{max-width:100%}

.skeptic .bub{background:#fdf7f5;border-color:#f0997b;border-left:3px solid #c1602e}
.arbiter .bub{background:#f6fafa;border-color:#bcdad9;border-left:3px solid #519d9e}

/* role label */
.role-label{font-size:.68rem;font-weight:600;letter-spacing:.07em;
            text-transform:uppercase;color:#8b909e;margin-bottom:8px}
.skeptic .role-label{color:#993c1d}
.arbiter .role-label{color:#2d6e6e}
.proposer-r .role-label{color:#1c3957}

/* summary line */
.summary-line{font-size:.88rem;color:#3a3f4b;margin-bottom:10px}
.summary-line strong{font-weight:600}
.no-challenge-note{font-size:.88rem;color:#21643a;font-style:italic}

/* ── query bubble ── */
.query-grid{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-top:10px}
.query-section{padding:10px 12px;background:#f8f9fb;
               border:1px solid #dbe5f2;border-radius:8px}
.section-label{font-size:.65rem;font-weight:600;letter-spacing:.08em;
               text-transform:uppercase;color:#8b909e;margin-bottom:6px}
.met-chips{display:flex;flex-wrap:wrap;gap:5px;margin-top:2px}
.met-chip{display:inline-block;font-size:.76rem;padding:2px 8px;
          background:#eef4fb;color:#2a4662;border-radius:4px;
          border:1px solid #cfe0f2;white-space:normal}
.muted{color:#9ca3af;font-size:.82rem}

/* class bars */
.cls-block{margin-top:2px;display:flex;flex-direction:column;gap:5px}
.cls-row{display:grid;grid-template-columns:1fr 60px 36px;
         align-items:center;gap:6px;font-size:.75rem}
.cls-name{color:#3a3f4b;white-space:normal;line-height:1.35}
.cls-track{height:6px;background:#e4e7ed;border-radius:3px;overflow:hidden}
.cls-fill {height:100%;background:#6191d3;border-radius:3px}
.cls-pct  {color:#8b909e;text-align:right;white-space:nowrap}

/* context pills (phenotype / species / biospecimen) */
.context-pills{display:flex;flex-wrap:wrap;gap:5px;margin-top:2px}
.ctx-pill{font-size:.75rem;padding:2px 9px;background:#f6f9fc;
          color:#52677e;border-radius:20px;border:1px solid #dbe5f2}

/* ── hypothesis cards ── */
.hyp-card{background:#f8f9fb;border:1px solid #dbe5f2;border-radius:8px;
          padding:11px 13px;margin-top:9px}
.hyp-header{display:flex;align-items:baseline;gap:8px;margin-bottom:5px}
.hyp-id   {font-size:.68rem;font-weight:700;color:#3966a2;
           background:#eef4fb;padding:1px 6px;border-radius:3px;flex-shrink:0}
.hyp-title{font-weight:600;font-size:.88rem;color:#1c1f26}
.type-pill{font-size:.68rem;color:#2d6e6e;background:#e1f5ee;border:1px solid #bcdad9;
           padding:2px 6px;border-radius:999px;white-space:nowrap}
.hyp-summary{font-size:.82rem;color:#4b5563;line-height:1.6;margin-top:2px}

/* metabolite tags inside hypothesis cards */
.mtag-row  {margin-top:7px;display:flex;flex-wrap:wrap;gap:4px;align-items:center}
.mtag-label{font-size:.65rem;font-weight:600;letter-spacing:.05em;
            text-transform:uppercase;color:#9ca3af;margin-right:3px;flex-shrink:0}
.mtag      {font-size:.72rem;padding:1px 7px;border-radius:3px}
.mtag-anc  {background:#eef4fb;color:#1c3957;border:1px solid #cfe0f2}
.mtag-ten  {background:#fff6e6;color:#80551b;border:1px solid #f2d29a}

/* ── challenge items ── */
.ch-item  {border-radius:6px;padding:9px 11px;margin-top:8px;
           border-left:3px solid transparent}
.ch-major {background:#fdeeee;border-left-color:#8c3030}
.ch-minor {background:#fff6e6;border-left-color:#c1602e}
.ch-meta  {display:flex;align-items:center;gap:7px;flex-wrap:wrap;margin-bottom:5px}
.sev-pill {font-size:.65rem;font-weight:700;letter-spacing:.04em;
           padding:2px 7px;border-radius:3px}
.sev-major{background:#8c3030;color:#fff}
.sev-minor{background:#c1602e;color:#fff}
.ch-type  {font-size:.72rem;color:#6b7280;font-style:italic}
.ch-target{font-size:.75rem;font-weight:600;color:#374151;
           background:#f1f3f5;padding:1px 6px;border-radius:3px}
.ch-text  {font-size:.82rem;color:#4b5563;line-height:1.6}

/* ── revision items ── */
.rev-item  {padding:9px 11px;margin-top:8px;border-radius:6px;
            background:#f8f9fb;border:1px solid #dbe5f2;
            border-left:3px solid #6191d3}
.act-modify{border-left-color:#3966a2}
.act-drop  {border-left-color:#8c3030}
.act-keep  {border-left-color:#21643a}
.rev-header{display:flex;align-items:center;gap:7px;margin-bottom:5px}
.act-pill  {font-size:.65rem;font-weight:700;letter-spacing:.04em;
            padding:2px 7px;border-radius:3px}
.act-modify .act-pill,.act-pill.act-modify{background:#3966a2;color:#fff}
.act-drop   .act-pill,.act-pill.act-drop  {background:#8c3030;color:#fff}
.act-keep   .act-pill,.act-pill.act-keep  {background:#21643a;color:#fff}
.rev-hid   {font-size:.8rem;font-weight:600;color:#374151}
.rev-just  {font-size:.82rem;color:#4b5563;line-height:1.6}
.rev-detail{margin-top:9px;padding:9px 11px;background:#fff;
            border:1px solid #dbe5f2;border-radius:6px}
.rev-new-type   {display:inline-block;font-size:.68rem;color:#2d6e6e;background:#e1f5ee;
                 border:1px solid #bcdad9;padding:2px 6px;border-radius:999px;margin-bottom:5px}
.rev-new-title  {font-size:.85rem;font-weight:600;color:#1c1f26;margin-bottom:4px}
.rev-new-summary{font-size:.8rem;color:#4b5563;line-height:1.55}

/* ── verdict cards ── */
.verdict-card{border-radius:8px;padding:14px 15px;margin-top:10px;
              border:1px solid #dbe5f2;background:#fff}
.st-accepted {background:#f3faf6;border-color:#bddfc9}
.st-rejected {background:#fdf3f3;border-color:#efc0c0}

.verdict-hdr{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin-bottom:9px}
.verdict-hid{font-size:.8rem;font-weight:700;color:#374151;
             background:#f1f3f5;padding:2px 8px;border-radius:3px}
.status-pill{font-size:.65rem;font-weight:700;letter-spacing:.05em;
             padding:2px 8px;border-radius:3px}
.st-accepted  .status-pill{background:#2f9e6a;color:#fff}
.st-rejected  .status-pill{background:#c1602e;color:#fff}

.verdict-type       {display:inline-block;font-size:.68rem;color:#2d6e6e;background:#e1f5ee;
                    border:1px solid #bcdad9;padding:2px 6px;border-radius:999px;margin-bottom:5px}
.verdict-title      {font-size:.9rem;font-weight:600;color:#1c1f26;margin-bottom:5px}
.verdict-fin-summary{font-size:.82rem;color:#4b5563;line-height:1.6;
                     margin-bottom:10px;
                     padding-bottom:10px;border-bottom:1px solid #dbe5f2}

/* rationale sections */
.rat-block {margin-top:8px}
.rat-label {font-size:.65rem;font-weight:600;letter-spacing:.07em;
            text-transform:uppercase;color:#8b909e;margin-bottom:4px;
            display:flex;align-items:center;gap:8px}
.rat-text  {font-size:.82rem;color:#4b5563;line-height:1.65}
.rat-divider{height:1px;background:#dbe5f2;margin:10px 0}

/* confidence pill inside rationale label */
.conf-pill{font-size:.62rem;font-weight:600;padding:1px 7px;
           border-radius:20px;letter-spacing:.03em}
.conf-high  {background:#eaf7ef;color:#21643a}
.conf-medium{background:#fff6e6;color:#80551b}
.conf-low   {background:#fdeeee;color:#8c3030}

/* ambiguity note */
.ambiguity-note{margin-top:12px;padding:9px 12px;background:#fff6e6;
                border-left:3px solid #c1602e;border-radius:0 6px 6px 0;
                font-size:.82rem;color:#4b5563}
.ambiguity-note strong{color:#1c1f26}

/* PubMed evidence */
.pubmed .bub{background:#f3faf6;border-color:#bddfc9;border-left:3px solid #21643a}
.pubmed .role-label{color:#21643a}
.pubmed-card{background:#fff;border:1px solid #bddfc9;border-radius:8px;
             padding:11px 13px;margin-top:9px}
.pubmed-hdr{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin-bottom:7px}
.pubmed-hid{font-size:.78rem;font-weight:700;color:#21643a;
            background:#eaf7ef;padding:2px 8px;border-radius:3px}
.lit-pill{font-size:.65rem;font-weight:700;letter-spacing:.04em;
          padding:2px 8px;border-radius:3px;text-transform:uppercase}
.lit-direct{background:#2f9e6a;color:#fff}
.lit-indirect{background:#bddfc9;color:#174d27}
.lit-insufficient{background:#e5e7eb;color:#4b5563}
.pm-relationship,.pm-note{font-size:.82rem;color:#4b5563;line-height:1.6;margin-top:5px}
.pmid-chips{display:flex;flex-wrap:wrap;gap:5px;margin-top:9px}
.pmid-chip{display:inline-block;font-size:.74rem;padding:2px 8px;
           background:#eaf7ef;color:#21643a;border-radius:4px;
           border:1px solid #bddfc9}
.pm-empty{font-size:.8rem;color:#8b909e;font-style:italic}

/* round divider */
.round-divider{display:flex;align-items:center;gap:10px;
               font-size:.65rem;font-weight:600;letter-spacing:.1em;
               text-transform:uppercase;color:#9ca3af;margin:4px 0}
.round-divider::before,.round-divider::after{
  content:"";flex:1;height:1px;background:#dbe5f2}
'
}


# ── Public function ───────────────────────────────────────────

#' Generate a chat-style debate HTML for one module
#'
#' Called automatically at the end of annotate_one_module_advanced().
#' Writes debate_<module_id>.html to output_dir (default: getwd()).
#'
#' @param module_id       Module identifier string
#' @param module_obj      Original module object (topN_metabolites used)
#' @param study_context   List with species, biospecimen, phenotype
#' @param hypotheses_v0   P0 hypotheses (list of hypothesis objects)
#' @param s1_challenges   S1 parsed challenges (list)
#' @param p1_revisions    P1 revision decisions (list)
#' @param hypotheses_v1   Hypotheses after P1 revisions (list)
#' @param s2_challenges   S2 challenges, or NULL if Round 2 did not run
#' @param p2_revisions    P2 revision decisions, or NULL
#' @param hypotheses_final Final hypotheses after all revisions
#' @param arbiter_result  List with $verdicts and $residual_ambiguity
#' @param phenotype_relationships Named list of PubMed phenotype summaries
#' @param pubmed_per_hyp Named list of retrieved PubMed records per hypothesis
#' @param output_dir      Directory for the output file (default: getwd())
#' @return Invisibly, the path of the written HTML file
generate_debate_html <- function(module_id,
                                  module_obj,
                                  study_context,
                                  hypotheses_v0,
                                  s1_challenges,
                                  p1_revisions,
                                  hypotheses_v1,
                                  s2_challenges          = NULL,
                                  p2_revisions           = NULL,
                                  hypotheses_final       = list(),
                                  arbiter_result         = list(verdicts = list(),
                                                                residual_ambiguity = NA),
                                  phenotype_relationships = NULL,
                                  pubmed_per_hyp         = NULL,
                                  output_dir             = getwd()) {

  # ── Extract module metadata ─────────────────────────────────
  top_names <- character(0)
  if (!is.null(module_obj$topN_metabolites) &&
      "Compound_name" %in% colnames(module_obj$topN_metabolites))
    top_names <- unique(as.character(module_obj$topN_metabolites$Compound_name))

  phenotype   <- .s(study_context$phenotype,   "—")
  species     <- .s(study_context$species,     "—")
  biospecimen <- .s(study_context$biospecimen, "—")
  # ── 1. Query bubble ─────────────────────────────────────────
  ctx_pills <- paste(
    sprintf('<span class="ctx-pill">%s</span>',
            .he(c(species, biospecimen, phenotype))),
    collapse = "")

  user_content <- sprintf(
    '<div class="role-label">\U0001f50d Query</div>
    <div class="summary-line">
      <strong>%s</strong> is a detected dysregulated functional module.
      Please interpret its biological significance.
    </div>
    <div class="query-section">
      <div class="section-label">Core metabolites</div>
      <div class="met-chips">%s</div>
    </div>
    <div class="query-section" style="margin-top:8px">
      <div class="section-label">Dominant chemical classes</div>
      %s
    </div>
    <div class="query-section" style="margin-top:8px">
      <div class="section-label">Study context</div>
      <div class="context-pills">%s</div>
    </div>',
    .he(module_id),
    .met_chips(top_names),
    .class_bars(module_obj, 3L),
    ctx_pills
  )
  bubbles <- .bubble("left", "\U0001f50d", "Query", user_content)

  # ── 2. P0 Proposer ──────────────────────────────────────────
  p0_cards <- paste(vapply(hypotheses_v0, .hyp_card, character(1L)),
                    collapse = "\n")
  p0_content <- sprintf(
    '<div class="role-label proposer-r">\U0001f52c Proposer &mdash; Initial</div>
    <div class="summary-line">
      I propose <strong>%d</strong> candidate hypothesis%s based on the module data:
    </div>
    %s',
    length(hypotheses_v0),
    if (length(hypotheses_v0) != 1L) "es" else "",
    p0_cards
  )
  bubbles <- paste0(bubbles, "\n",
                    .bubble("right", "\U0001f52c", "Proposer", p0_content, "proposer-r"))

  # ── 3. S1 Skeptic ───────────────────────────────────────────
  real_s1 <- Filter(function(c) !identical(c$severity, "no_challenge"), s1_challenges)
  if (length(real_s1) > 0L) {
    s1_items  <- paste(vapply(real_s1, .ch_item, character(1L)), collapse = "\n")
    s1_content <- sprintf(
      '<div class="role-label">\U0001f9d0 Skeptic &mdash; Round 1</div>
      <div class="summary-line">I raise <strong>%d</strong> challenge%s:</div>
      %s',
      length(real_s1),
      if (length(real_s1) != 1L) "s" else "",
      s1_items
    )
  } else {
    s1_content <- sprintf(
      '<div class="role-label">\U0001f9d0 Skeptic &mdash; Round 1</div>
      <div class="no-challenge-note">%s</div>',
      "No substantive challenges. The proposed hypotheses are well-supported by the data."
    )
  }
  bubbles <- paste0(bubbles, "\n",
                    .bubble("left", "\U0001f9d0", "Skeptic", s1_content, "skeptic"))

  # ── 4. P1 Proposer revision ─────────────────────────────────
  p1_items    <- paste(vapply(p1_revisions, .rev_item,
                              character(1L), hypotheses_after = hypotheses_v1),
                       collapse = "\n")
  modified_r1 <- .count_action(p1_revisions, "modify")
  dropped_r1  <- .count_action(p1_revisions, "drop")
  kept_r1     <- length(p1_revisions) - modified_r1 - dropped_r1
  p1_content  <- sprintf(
    '<div class="role-label proposer-r">✏️ Proposer &mdash; Revision 1</div>
    <div class="summary-line">
      Addressing Round 1:
      <strong>%d modified</strong>, <strong>%d dropped</strong>, %d kept.
    </div>
    %s',
    modified_r1, dropped_r1, kept_r1, p1_items
  )
  bubbles <- paste0(bubbles, "\n",
                    .bubble("right", "✏️", "Proposer", p1_content, "proposer-r"))

  # ── 5. S2 Skeptic (conditional) ─────────────────────────────
  if (!is.null(s2_challenges) && length(s2_challenges) > 0L) {
    real_s2 <- Filter(function(c) !identical(c$severity, "no_challenge"), s2_challenges)
    bubbles <- paste0(bubbles,
                      '\n<div class="round-divider">Round 2</div>')
    if (length(real_s2) > 0L) {
      s2_items   <- paste(vapply(real_s2, .ch_item, character(1L)), collapse = "\n")
      s2_content <- sprintf(
        '<div class="role-label">\U0001f9d0 Skeptic &mdash; Round 2</div>
        <div class="summary-line">After revisions, I raise <strong>%d</strong> further challenge%s:</div>
        %s',
        length(real_s2),
        if (length(real_s2) != 1L) "s" else "",
        s2_items
      )
      bubbles <- paste0(bubbles, "\n",
                        .bubble("left", "\U0001f9d0", "Skeptic", s2_content, "skeptic"))
    } else {
      s2_content <- sprintf(
        '<div class="role-label">\U0001f9d0 Skeptic &mdash; Round 2</div>
        <div class="no-challenge-note">%s</div>',
        "No substantive challenges after reviewing the revision."
      )
      bubbles <- paste0(bubbles, "\n",
                        .bubble("left", "\U0001f9d0", "Skeptic", s2_content, "skeptic"))
    }
  }

  # ── 6. P2 Proposer revision (conditional) ───────────────────
  if (!is.null(p2_revisions) && length(p2_revisions) > 0L) {
    p2_items    <- paste(vapply(p2_revisions, .rev_item,
                                character(1L), hypotheses_after = hypotheses_final),
                         collapse = "\n")
    modified_r2 <- .count_action(p2_revisions, "modify")
    dropped_r2  <- .count_action(p2_revisions, "drop")
    kept_r2     <- length(p2_revisions) - modified_r2 - dropped_r2
    p2_content  <- sprintf(
      '<div class="role-label proposer-r">✏️ Proposer &mdash; Revision 2</div>
      <div class="summary-line">
        Addressing Round 2:
        <strong>%d modified</strong>, <strong>%d dropped</strong>, %d kept.
      </div>
      %s',
      modified_r2, dropped_r2, kept_r2, p2_items
    )
    bubbles <- paste0(bubbles, "\n",
                      .bubble("right", "✏️", "Proposer", p2_content, "proposer-r"))
  }

  # ── 7. PubMed phenotype evidence ────────────────────────────
  pubmed_ids <- unique(c(names(phenotype_relationships %||% list()),
                         names(pubmed_per_hyp %||% list())))
  pubmed_ids <- pubmed_ids[nzchar(pubmed_ids)]
  if (length(pubmed_ids) > 0L) {
    pm_cards <- paste(vapply(pubmed_ids, function(hid) {
      .pubmed_card(
        hid = hid,
        rel = (phenotype_relationships %||% list())[[hid]],
        evidence = (pubmed_per_hyp %||% list())[[hid]]
      )
    }, character(1L)), collapse = "\n")
    pm_content <- sprintf(
      '<div class="role-label">📚 PubMed phenotype evidence</div>
      <div class="summary-line">
        Retrieved literature and phenotype relationship summaries used by the Arbiter:
      </div>
      %s',
      pm_cards
    )
    bubbles <- paste0(bubbles, "\n",
                      .bubble("left", "📚", "PubMed", pm_content, "pubmed"))
  }

  # ── 8. A0 Arbiter ───────────────────────────────────────────
  verdicts <- arbiter_result$verdicts %||% list()
  n_acc    <- sum(vapply(verdicts, function(v) identical(v$status, "accepted"),  logical(1L)))
  n_rej    <- sum(vapply(verdicts, function(v) identical(v$status, "rejected"),  logical(1L)))

  vcards <- paste(vapply(verdicts, .verdict_card,
                         character(1L), hypotheses_final = hypotheses_final),
                  collapse = "\n")

  amb      <- .s(arbiter_result$residual_ambiguity, "")
  amb_html <- if (.is_nonempty(amb))
    sprintf('<div class="ambiguity-note"><strong>Residual ambiguity:</strong> %s</div>',
            .he(amb))
  else ""

  arbiter_content <- sprintf(
    '<div class="role-label">⚖ Arbiter &mdash; Final Verdict</div>
    <div class="summary-line">
      Reviewing the full debate:
      <strong>%d accepted</strong> &middot;
      <strong>%d rejected</strong>.
    </div>
    %s%s',
    n_acc, n_rej, vcards, amb_html
  )
  bubbles <- paste0(bubbles, "\n",
                    .bubble("center", "⚖", "Arbiter",
                            arbiter_content, "arbiter"))

  # ── Assemble page ────────────────────────────────────────────
  html <- sprintf(
'<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Module Debate &mdash; %s</title>
<style>%s</style>
</head>
<body>
<div class="chat-body">
%s
</div>
</body>
</html>',
    .he(module_id),
    .debate_css(),
    bubbles
  )

  fname    <- gsub("[^A-Za-z0-9._-]", "_", module_id)
  out_path <- file.path(output_dir, sprintf("debate_%s.html", fname))
  writeLines(html, out_path, useBytes = FALSE)
  message(sprintf("[HTML] Debate visualization → %s", out_path))
  invisible(out_path)
}
