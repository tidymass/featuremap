# ============================================================
# 45_debate_html.R
# Debate HTML visualizer for the advanced multi-agent pipeline
#
# Sourced automatically by 45_advanced_mode.R.
# Exports one public function:
#   generate_debate_html(module_id, module_obj, study_context,
#                        hypotheses_v0, s1_challenges, p1_revisions,
#                        hypotheses_v1, s2_challenges, p2_revisions,
#                        hypotheses_final, arbiter_result, output_dir)
# ============================================================


# ── Internal helpers (not exported) ──────────────────────────

#' Escape special HTML characters
.he <- function(s) {
  s <- as.character(s %||% "")
  s <- gsub("&",  "&amp;",  s, fixed = TRUE)
  s <- gsub("<",  "&lt;",   s, fixed = TRUE)
  s <- gsub(">",  "&gt;",   s, fixed = TRUE)
  s <- gsub('"',  "&quot;", s, fixed = TRUE)
  s
}

#' Render metabolite tag spans from a comma-separated string
.tag_row <- function(csv, cls, label) {
  if (!nzchar(trimws(csv))) return("")
  tags <- paste(
    sprintf('<span class="tag %s">%s</span>', cls,
            .he(trimws(strsplit(csv, ",")[[1]]))),
    collapse = "")
  sprintf('<div class="tag-row"><span class="tag-label">%s</span> %s</div>',
          label, tags)
}

#' Format the top chemical classes for the query bubble
.top_class_text <- function(module_obj, n = 3L) {
  class_df <- module_obj$other_metabolites
  if (is.null(class_df) || !is.data.frame(class_df) || nrow(class_df) == 0 ||
      !"class" %in% colnames(class_df)) {
    return("—")
  }

  if ("n_metabolites" %in% colnames(class_df)) {
    class_df <- class_df[order(-as.numeric(class_df$n_metabolites)), , drop = FALSE]
  }
  rows <- head(class_df, n)
  vals <- vapply(seq_len(nrow(rows)), function(i) {
    cls <- as.character(rows$class[[i]] %||% "unknown")
    prop <- if ("n_metabolites" %in% colnames(rows)) {
      x <- suppressWarnings(as.numeric(rows$n_metabolites[[i]]))
      if (!is.na(x)) sprintf(" %.0f%%", x * 100) else ""
    } else {
      ""
    }
    paste0(cls, prop)
  }, character(1))
  paste(vals, collapse = "; ")
}

#' Render one hypothesis card (P0 / P1 / P2 bubbles)
.hyp_card <- function(h, badge = NULL) {
  anc <- paste(h$anchor_metabolites  %||% character(0), collapse = ", ")
  ten <- paste(h$tension_metabolites %||% character(0), collapse = ", ")
  badge_html <- if (!is.null(badge))
    sprintf('<span class="badge badge-%s">%s</span>', .he(badge), .he(toupper(badge)))
  else ""
  sprintf(
    '<div class="hyp-card">
      <div class="hyp-id">%s %s</div>
      <div class="hyp-title">%s</div>
      <div class="hyp-summary">%s</div>
      %s%s
    </div>',
    .he(h$hypothesis_id), badge_html,
    .he(h$title   %||% ""),
    .he(h$summary %||% ""),
    .tag_row(anc, "anc", "Anchors:"),
    .tag_row(ten, "ten", "Unexplained/conflicting:")
  )
}

#' Render one Skeptic challenge item
.ch_item <- function(ch) {
  sev   <- ch$severity       %||% "minor"
  ctype <- ch$challenge_type %||% ""
  if (identical(sev, "no_challenge")) {
    return("")
  }
  sev_cls <- if (identical(sev, "major")) "sev-major" else "sev-minor"
  sprintf(
    '<div class="ch-item">
      <span class="sev-badge %s">%s</span>
      <span class="ch-type">%s</span> → <strong>%s</strong>
      <div class="ch-text">%s</div>
    </div>',
    sev_cls, .he(toupper(sev)),
    .he(ctype), .he(ch$target_hypothesis_id),
    .he(ch$challenge %||% "")
  )
}

#' Render one Proposer revision item
.rev_item <- function(rev, hypotheses_after) {
  hid     <- rev$hypothesis_id %||% ""
  action  <- toupper(rev$action %||% "keep")
  just    <- rev$justification %||% ""
  act_cls <- switch(action, MODIFY = "act-modify", DROP = "act-drop", "act-keep")
  h_new   <- Filter(function(x) x$hypothesis_id == hid, hypotheses_after)

  detail <- ""
  if (action == "MODIFY" && length(h_new) > 0) {
    h   <- h_new[[1]]
    anc <- paste(h$anchor_metabolites  %||% character(0), collapse = ", ")
    ten <- paste(h$tension_metabolites %||% character(0), collapse = ", ")
    detail <- sprintf(
      '<div class="rev-detail">
        <div><em>Revised title:</em> %s</div>
        <div class="rev-summary"><em>Revised summary:</em> %s</div>
        %s%s
      </div>',
      .he(h$title %||% ""),
      .he(h$summary %||% ""),
      .tag_row(anc, "anc", "Anchors:"),
      .tag_row(ten, "ten", "Unexplained/conflicting:")
    )
  }
  sprintf(
    '<div class="rev-item">
      <span class="act-badge %s">%s</span> <strong>%s</strong>
      <div class="rev-just">%s</div>
      %s
    </div>',
    act_cls, action, .he(hid), .he(just), detail
  )
}

#' Render one Arbiter verdict card
.verdict_card <- function(v, hypotheses_lookup) {
  hid    <- v$hypothesis_id    %||% ""
  status <- v$status           %||% "rejected"
  conf   <- v$confidence_label %||% ""
  rat    <- v$verdict_rationale %||% ""
  conf_rat <- v$confidence_rationale %||% ""
  st_cls <- switch(status,
                   accepted  = "st-accepted",
                   uncertain = "st-uncertain",
                   "st-rejected")
  h_fin  <- Filter(function(x) x$hypothesis_id == hid, hypotheses_lookup %||% list())
  final_title <- if (length(h_fin) > 0) h_fin[[1]]$title %||% "" else ""
  final_summary <- if (length(h_fin) > 0) h_fin[[1]]$summary %||% "" else ""
  sprintf(
    '<div class="verdict-card %s">
      <div class="verdict-header">
        <span class="verdict-id">%s</span>
        <span class="verdict-status %s">%s</span>
        %s
      </div>
      %s
      %s
      %s
      %s
    </div>',
    st_cls, .he(hid), st_cls, .he(toupper(status)),
    if (nzchar(conf))
      sprintf('<span class="conf-badge conf-%s">%s confidence</span>', .he(conf), .he(conf))
    else "",
    if (nzchar(final_title))
      sprintf('<div class="verdict-title">%s</div>', .he(final_title))
    else "",
    if (nzchar(final_summary))
      sprintf('<div class="verdict-summary">%s</div>', .he(final_summary))
    else "",
    if (nzchar(rat))
      sprintf('<div class="verdict-rat"><strong>Arbiter rationale:</strong> %s</div>', .he(rat))
    else "",
    if (nzchar(conf_rat))
      sprintf('<div class="verdict-rat"><strong>Confidence rationale:</strong> %s</div>', .he(conf_rat))
    else ""
  )
}

#' Count revision actions safely for empty revision lists
.count_action <- function(revisions, action) {
  if (is.null(revisions) || length(revisions) == 0) return(0L)
  sum(vapply(revisions, function(r) identical(r$action, action), logical(1)))
}

#' Render PubMed phenotype relationship summaries
.pubmed_item <- function(hid, rel) {
  support <- rel$literature_support %||% "insufficient"
  relationship <- rel$phenotype_relationship %||% ""
  note <- rel$evidence_note %||% ""
  sprintf(
    '<div class="pubmed-item">
      <div><strong>%s</strong> <span class="lit-badge lit-%s">Literature support: %s</span></div>
      <div class="pubmed-rel"><strong>Phenotype relationship:</strong> %s</div>
      %s
    </div>',
    .he(hid),
    .he(support),
    .he(toupper(support)),
    .he(relationship),
    if (nzchar(note)) sprintf('<div class="pubmed-note"><strong>Evidence note:</strong> %s</div>', .he(note)) else ""
  )
}

#' Wrap content in a chat message bubble
.bubble <- function(role, avatar, avatar_label, content_html, extra_class = "") {
  sprintf(
    '<div class="msg msg-%s %s">
      <div class="avatar avatar-%s" title="%s">%s</div>
      <div class="bubble bubble-%s">%s</div>
    </div>',
    role, extra_class,
    role, .he(avatar_label), .he(avatar),
    role, content_html
  )
}

#' Inline CSS for the debate HTML page
.debate_css <- function() {
'*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
       background: #f0f2f5; color: #1a1a2e; min-height: 100vh; }
.chat-header { background: linear-gradient(135deg,#1a1a2e,#16213e);
               color:#fff; padding:24px 32px; }
.chat-header h1 { font-size:1.4rem; font-weight:700; }
.chat-header p  { font-size:.85rem; color:#8892b0; margin-top:4px; }
.chat-body { max-width:860px; margin:0 auto; padding:28px 20px;
             display:flex; flex-direction:column; gap:20px; }
.msg { display:flex; gap:14px; align-items:flex-start; }
.msg-left   { flex-direction:row; }
.msg-right  { flex-direction:row-reverse; }
.msg-center { flex-direction:column; align-items:stretch; }
.avatar { width:44px; height:44px; border-radius:50%; display:flex;
          align-items:center; justify-content:center; font-size:1.25rem;
          flex-shrink:0; box-shadow:0 2px 6px rgba(0,0,0,.15); }
.avatar-left   { background:#e8eaf6; }
.avatar-right  { background:#e3f2fd; }
.avatar-center { width:44px; height:44px; border-radius:50%; display:flex;
                 align-items:center; justify-content:center; font-size:1.25rem;
                 flex-shrink:0; background:#f3e5f5; margin:0 auto 8px; }
.bubble { background:#fff; border-radius:16px; padding:16px 18px;
          box-shadow:0 2px 8px rgba(0,0,0,.08); max-width:82%; line-height:1.55; }
.bubble-left  { border-top-left-radius:4px; }
.bubble-right { border-top-right-radius:4px; }
.bubble-center { max-width:100%; border-radius:16px; }
.msg-center .bubble { max-width:100%; }
.skeptic .bubble { background:#fff8f5; border-left:3px solid #ff7043; }
.arbiter .bubble { background:#fdf6ff; border:2px solid #ab47bc; }
.msg-role { font-size:.78rem; font-weight:700; letter-spacing:.05em;
            text-transform:uppercase; color:#5c6bc0; margin-bottom:8px; }
.skeptic .msg-role { color:#e64a19; }
.arbiter .msg-role { color:#8e24aa; }
p { font-size:.9rem; margin-bottom:10px; }
.meta-row { font-size:.8rem; color:#666; display:flex; flex-wrap:wrap;
            gap:8px; margin-top:6px; }
.meta-row span { background:#eee; padding:2px 8px; border-radius:20px; }
.no-issue { color:#2e7d32; font-style:italic; }
.hyp-card { background:#f8faff; border:1px solid #dde3f0;
            border-radius:10px; padding:12px 14px; margin-top:10px; }
.hyp-id      { font-size:.72rem; font-weight:700; color:#5c6bc0; margin-bottom:4px; }
.hyp-title   { font-weight:600; font-size:.9rem; margin-bottom:6px; }
.hyp-summary { font-size:.83rem; color:#444; line-height:1.5; white-space:pre-wrap; }
.tag-row  { margin-top:7px; display:flex; flex-wrap:wrap; gap:5px; align-items:center; }
.tag-label{ font-size:.72rem; color:#888; }
.tag      { font-size:.72rem; padding:2px 8px; border-radius:20px; }
.tag.anc  { background:#e3f2fd; color:#1565c0; }
.tag.ten  { background:#fff3e0; color:#e65100; }
.ch-item { background:#fff5f2; border-radius:8px; padding:10px 12px;
           margin-top:8px; font-size:.85rem; }
.ch-item.no-challenge { background:#f1f8e9; }
.ch-text { margin-top:6px; color:#555; }
.ch-type { font-size:.78rem; color:#888; }
.sev-badge { display:inline-block; font-size:.68rem; font-weight:700;
             padding:2px 7px; border-radius:20px; margin-right:6px; }
.sev-major { background:#ff5252; color:#fff; }
.sev-minor { background:#ffd740; color:#333; }
.sev-none  { background:#aed581; color:#33691e; }
.rev-item { border-left:3px solid #90caf9; padding:8px 12px;
            margin-top:10px; font-size:.85rem; }
.act-badge  { display:inline-block; font-size:.7rem; font-weight:700;
              padding:2px 8px; border-radius:4px; margin-right:6px; }
.act-modify { background:#42a5f5; color:#fff; }
.act-drop   { background:#ef5350; color:#fff; }
.act-keep   { background:#66bb6a; color:#fff; }
.rev-just   { color:#555; margin-top:4px; font-size:.82rem; }
.rev-detail { margin-top:8px; font-size:.82rem; background:#f0f7ff;
              padding:8px 10px; border-radius:6px; }
.rev-summary { margin-top:6px; line-height:1.5; white-space:pre-wrap; }
.verdict-card { border-radius:10px; padding:14px 16px; margin-top:10px; font-size:.85rem; }
.st-accepted  { background:#e8f5e9; border:1px solid #81c784; }
.st-uncertain { background:#fff8e1; border:1px solid #ffd54f; }
.st-rejected  { background:#ffebee; border:1px solid #ef9a9a; }
.verdict-header { display:flex; align-items:center; gap:10px; flex-wrap:wrap;
                  margin-bottom:8px; }
.verdict-id     { font-weight:700; color:#333; }
.verdict-status { font-size:.72rem; font-weight:700; padding:2px 9px; border-radius:20px; }
.st-accepted  .verdict-status { background:#2e7d32; color:#fff; }
.st-uncertain .verdict-status { background:#f9a825; color:#fff; }
.st-rejected  .verdict-status { background:#c62828; color:#fff; }
.conf-badge  { font-size:.7rem; padding:2px 8px; border-radius:20px;
               background:#ede7f6; color:#6a1b9a; }
.conf-high   { background:#ce93d8; color:#4a148c; }
.conf-medium { background:#e1bee7; color:#6a1b9a; }
.conf-low    { background:#f3e5f5; color:#6a1b9a; }
.verdict-title   { font-weight:600; font-size:.9rem; margin-bottom:6px; }
.verdict-summary { font-size:.83rem; color:#444; margin-bottom:8px; line-height:1.5; white-space:pre-wrap; }
.verdict-rat     { font-size:.82rem; color:#555; }
.pubmed-item { background:#f1f8f4; border-left:3px solid #43a047;
               border-radius:8px; padding:10px 12px; margin-top:8px; font-size:.85rem; }
.pubmed-rel { color:#444; margin-top:6px; }
.pubmed-note { color:#666; margin-top:6px; font-size:.8rem; font-style:italic; }
.lit-badge { display:inline-block; font-size:.68rem; font-weight:700;
             padding:2px 7px; border-radius:20px; margin-left:6px; }
.lit-direct { background:#2e7d32; color:#fff; }
.lit-indirect { background:#81c784; color:#1b5e20; }
.lit-insufficient { background:#c8e6c9; color:#33691e; }
.ambiguity-note  { margin-top:14px; background:#fffde7; border-left:3px solid #f9a825;
                   padding:10px 14px; border-radius:6px; font-size:.85rem; }
.round-divider { text-align:center; font-size:.75rem; color:#aaa;
                 letter-spacing:.1em; text-transform:uppercase; margin:4px 0; }
.round-divider::before, .round-divider::after {
  content:""; display:inline-block; width:30%; height:1px;
  background:#ddd; vertical-align:middle; margin:0 10px; }
.badge          { font-size:.65rem; padding:1px 6px; border-radius:20px;
                  margin-left:5px; vertical-align:middle; }
.badge-major    { background:#ff5252; color:#fff; }
.badge-minor    { background:#ffd740; color:#333; }
.badge-accepted { background:#2e7d32; color:#fff; }
.badge-uncertain{ background:#f9a825; color:#fff; }
.badge-rejected { background:#c62828; color:#fff; }'
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
#' @param output_dir      Directory for the output file (default: getwd())
#' @return Invisibly, the path of the written HTML file
#' @export
generate_debate_html <- function(module_id,
                                  module_obj,
                                  study_context,
                                  hypotheses_v0,
                                  s1_challenges,
                                  p1_revisions,
                                  hypotheses_v1,
                                  s2_challenges    = NULL,
                                  p2_revisions     = NULL,
                                  hypotheses_final = list(),
                                  arbiter_result   = list(verdicts = list(),
                                                          residual_ambiguity = NA),
                                  phenotype_relationships = NULL,
                                  output_dir       = getwd()) {

  # ── Study context ────────────────────────────────────────────
  top_names <- character(0)
  if (!is.null(module_obj$topN_metabolites) &&
      "Compound_name" %in% colnames(module_obj$topN_metabolites))
    top_names <- unique(as.character(module_obj$topN_metabolites$Compound_name))
  top_str     <- if (length(top_names) > 0) paste(head(top_names, 10), collapse = ", ") else "—"
  class_str   <- .top_class_text(module_obj, n = 3L)
  phenotype   <- study_context$phenotype   %||% "—"
  species     <- study_context$species     %||% "—"
  biospecimen <- study_context$biospecimen %||% "—"

  # ── 1. User query bubble ─────────────────────────────────────
  user_content <- sprintf(
    '<div class="msg-role">\U0001f50d Research Query</div>
    <p><strong>%s</strong> is a detected dysregulated functional module.
       Please interpret its biological significance.</p>
    <div class="meta-row">
      <span>Top metabolites: <em>%s</em></span>
      <span>Top classes: <em>%s</em></span>
      <span>Phenotype: <em>%s</em></span>
      <span>Species: <em>%s</em></span>
      <span>Biospecimen: <em>%s</em></span>
    </div>',
    .he(module_id),
    .he(top_str), .he(class_str), .he(phenotype), .he(species), .he(biospecimen)
  )
  bubbles <- .bubble("left", "\U0001f50d", "User", user_content)

  # ── 2. P0 Proposer ───────────────────────────────────────────
  p0_cards <- paste(sapply(hypotheses_v0, .hyp_card), collapse = "\n")
  p0_content <- sprintf(
    '<div class="msg-role">\U0001f52c Proposer</div>
    <p>Based on the module data, I propose <strong>%d</strong> candidate hypothesis%s:</p>
    %s',
    length(hypotheses_v0),
    if (length(hypotheses_v0) != 1) "es" else "",
    p0_cards
  )
  bubbles <- paste0(bubbles, "\n", .bubble("right", "\U0001f52c", "Proposer", p0_content))

  # ── 3. S1 Skeptic ────────────────────────────────────────────
  real_s1 <- Filter(function(c) !identical(c$severity, "no_challenge"), s1_challenges)
  if (length(real_s1) > 0) {
    s1_items <- paste(vapply(real_s1, .ch_item, character(1)), collapse = "\n")
    s1_content <- sprintf(
      '<div class="msg-role">\U0001f914 Skeptic — Round 1</div>
      <p>I raise <strong>%d</strong> challenge%s against the proposed hypotheses:</p>
      %s',
      length(real_s1),
      if (length(real_s1) != 1) "s" else "",
      s1_items
    )
  } else {
    s1_content <- '<div class="msg-role">\U0001f914 Skeptic — Round 1</div>
    <p class="no-issue">No substantive challenges. The proposed hypotheses are well-supported.</p>'
  }
  bubbles <- paste0(bubbles, "\n",
                    .bubble("left", "\U0001f914", "Skeptic", s1_content, "skeptic"))

  # ── 4. P1 Proposer revision ──────────────────────────────────
  p1_items    <- paste(sapply(p1_revisions, .rev_item, hypotheses_after = hypotheses_v1),
                       collapse = "\n")
  dropped_r1  <- .count_action(p1_revisions, "drop")
  modified_r1 <- .count_action(p1_revisions, "modify")
  p1_content  <- sprintf(
    '<div class="msg-role">✏️ Proposer — Round 1 Revision</div>
    <p>Addressing Round 1: <strong>%d modified</strong>, <strong>%d dropped</strong>, %d kept.</p>
    %s',
    modified_r1, dropped_r1,
    length(p1_revisions) - modified_r1 - dropped_r1,
    p1_items
  )
  bubbles <- paste0(bubbles, "\n",
                    .bubble("right", "✏️", "Proposer", p1_content))

  # ── 5. S2 Skeptic (conditional) ──────────────────────────────
  if (!is.null(s2_challenges) && length(s2_challenges) > 0) {
    real_s2   <- Filter(function(c) !identical(c$severity, "no_challenge"), s2_challenges)
    if (length(real_s2) > 0) {
      bubbles <- paste0(bubbles,
                        '\n<div class="round-divider">Round 2</div>')
      s2_items  <- paste(vapply(real_s2, .ch_item, character(1)), collapse = "\n")
      s2_content <- sprintf(
        '<div class="msg-role">\U0001f914 Skeptic — Round 2</div>
        <p>After Round 1 revisions, I raise <strong>%d</strong> further challenge%s:</p>
        %s',
        length(real_s2),
        if (length(real_s2) != 1) "s" else "",
        s2_items
      )
      bubbles <- paste0(bubbles, "\n",
                        .bubble("left", "\U0001f914", "Skeptic", s2_content, "skeptic"))
    }
  }

  # ── 6. P2 Proposer revision (conditional) ────────────────────
  if (!is.null(p2_revisions) && length(p2_revisions) > 0) {
    p2_items    <- paste(sapply(p2_revisions, .rev_item, hypotheses_after = hypotheses_final),
                         collapse = "\n")
    dropped_r2  <- .count_action(p2_revisions, "drop")
    modified_r2 <- .count_action(p2_revisions, "modify")
    p2_content  <- sprintf(
      '<div class="msg-role">✏️ Proposer — Round 2 Revision</div>
      <p>Addressing Round 2: <strong>%d modified</strong>, <strong>%d dropped</strong>, %d kept.</p>
      %s',
      modified_r2, dropped_r2,
      length(p2_revisions) - modified_r2 - dropped_r2,
      p2_items
    )
    bubbles <- paste0(bubbles, "\n",
                      .bubble("right", "✏️", "Proposer", p2_content))
  }

  # ── 7. PubMed relationship summaries ─────────────────────────
  if (!is.null(phenotype_relationships) && length(phenotype_relationships) > 0) {
    pm_items <- paste(
      vapply(names(phenotype_relationships), function(hid) {
        .pubmed_item(hid, phenotype_relationships[[hid]])
      }, character(1)),
      collapse = "\n"
    )
    pm_content <- sprintf(
      '<div class="msg-role">📚 PubMed Relationship Summary</div>
      <p>Phenotype relationship evidence summarized before final arbitration:</p>
      %s',
      pm_items
    )
    bubbles <- paste0(bubbles, "\n",
                      .bubble("left", "📚", "PubMed Summary", pm_content))
  }

  # ── 8. A0 Arbiter ────────────────────────────────────────────
  verdicts     <- arbiter_result$verdicts %||% list()
  n_acc        <- sum(sapply(verdicts, function(v) identical(v$status, "accepted")))
  n_unc        <- sum(sapply(verdicts, function(v) identical(v$status, "uncertain")))
  n_rej        <- sum(sapply(verdicts, function(v) identical(v$status, "rejected")))
  hypotheses_lookup <- c(
    hypotheses_final %||% list(),
    hypotheses_v0 %||% list()
  )
  if (length(hypotheses_lookup) > 0) {
    hypothesis_ids <- vapply(hypotheses_lookup, function(h) h$hypothesis_id %||% "", character(1))
    hypotheses_lookup <- hypotheses_lookup[!duplicated(hypothesis_ids)]
  }
  verdict_cards <- paste(sapply(verdicts, .verdict_card,
                                hypotheses_lookup = hypotheses_lookup),
                         collapse = "\n")
  amb      <- arbiter_result$residual_ambiguity %||% NA
  amb_html <- if (!is.na(amb) && nzchar(trimws(amb)))
    sprintf('<div class="ambiguity-note"><strong>Residual ambiguity:</strong> %s</div>', .he(amb))
  else ""
  arbiter_content <- sprintf(
    '<div class="msg-role">⚖️ Arbiter — Final Verdict</div>
    <p>After reviewing the full debate:
       <strong>%d accepted</strong>,
       <strong>%d uncertain</strong>,
       <strong>%d rejected</strong>.</p>
    %s%s',
    n_acc, n_unc, n_rej, verdict_cards, amb_html
  )
  bubbles <- paste0(bubbles, "\n",
                    .bubble("center", "⚖️", "Arbiter",
                            arbiter_content, "arbiter"))

  # ── Assemble final HTML ──────────────────────────────────────
  html <- sprintf(
'<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Debate — Module %s</title>
  <style>%s</style>
</head>
<body>
  <div class="chat-header">
    <h1>\U0001f9ec Module Debate: %s</h1>
    <p>Multi-agent metabolic annotation pipeline &nbsp;·&nbsp; %s</p>
  </div>
  <div class="chat-body">
%s
  </div>
</body>
</html>',
    .he(module_id),
    .debate_css(),
    .he(module_id),
    .he(as.character(Sys.Date())),
    bubbles
  )

  # ── Write file ───────────────────────────────────────────────
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  fname    <- gsub("[^A-Za-z0-9._-]", "_", module_id)
  out_path <- file.path(output_dir, sprintf("debate_%s.html", fname))
  writeLines(html, out_path, useBytes = FALSE)
  message(sprintf("[HTML] Debate visualization → %s", out_path))
  invisible(out_path)
}
