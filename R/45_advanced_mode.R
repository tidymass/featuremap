# ============================================================
# Advanced Mode Pipeline
#
# Multi-agent debate for metabolic module annotation.
#
# Agents (per module):
#   P0  Proposer initial      – generates 1-3 candidate biochemical hypotheses
#   S1  Skeptic round 1       – challenges each hypothesis
#   P1  Proposer revision 1   – keep / modify / drop
#   S2  Skeptic round 2       – (conditional) re-challenges revised set
#   P2  Proposer revision 2   – (conditional) final keep / modify / drop
#   A0  Arbiter               – accepted / rejected + confidence
#
# Information visible per agent:
#   P0  : full module data (metabolites, class, pathway, Signal_concentration)
#          → outputs: hypothesis_id, type, title, summary, anchor_metabolites, tension_metabolites
#          → NO claim field (removed; full summary serves as the mechanistic statement)
#
#   S1/S2: module evidence (top metabolites, class, pathway coverage)
#          + biological context: origin_profile, species, biospecimen (NO phenotype)
#          + current hypotheses: type, title, FULL biochemical summary, anchor, tension metabolites
#
#   P1/P2: module evidence + biological context (same as S1/S2, no phenotype)
#          + current full hypothesis text + Skeptic challenges
#          → provides biochemical grounding for rebuttal decisions
#
#   A0  : raw module evidence (top metabolites, class, pathway coverage)
#         + surviving hypotheses (full) + debate history digest
#         + origin_profile + species + biospecimen + phenotype
#         + per-hypothesis PubMed phenotype summaries
#         → does NOT receive global_direction (too complex to use as mechanistic constraint)
#
# PubMed: two-layer search (deduplicated anchors + per-hypothesis title, each with phenotype),
#   then anchor evidence is assigned back to the hypotheses that use those anchors.
#   Search (fetch) and summarization (LLM) are separate steps.
#
# global_direction: computed from module data, stored in output row ONLY.
# It is NOT passed to any agent.
# ============================================================


# ============================================================
# SECTION 1 – JSON parsing utilities (top-level, advanced mode)
# ============================================================

adv_normalize_parsed <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x) && nrow(x) >= 1) return(as.list(x[1, , drop = FALSE]))
  if (is.list(x) && !is.null(names(x))) return(x)
  if (is.list(x) && length(x) >= 1 && is.list(x[[1]]) && !is.null(names(x[[1]])))
    return(x[[1]])
  NULL
}

adv_parse_json_response <- function(txt) {
  if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)

  direct <- tryCatch(
    jsonlite::fromJSON(txt, simplifyVector = FALSE),
    error = function(e) NULL
  )
  direct <- adv_normalize_parsed(direct)
  if (!is.null(direct)) return(direct)

  json_txt <- tryCatch(extract_json_object(txt), error = function(e) NULL)
  if (!is.null(json_txt)) {
    obj <- tryCatch(
      jsonlite::fromJSON(json_txt, simplifyVector = FALSE),
      error = function(e) NULL
    )
    obj <- adv_normalize_parsed(obj)
    if (!is.null(obj)) return(obj)
  }

  fence <- regmatches(txt, regexpr("```\\s*[\\s\\S]*?\\s*```", txt, perl = TRUE))
  if (length(fence) > 0 && nzchar(fence[1])) {
    inner <- sub("^```[a-zA-Z]*\\s*", "", fence[1], perl = TRUE)
    inner <- sub("\\s*```$", "", inner, perl = TRUE)
    fobj <- tryCatch(
      jsonlite::fromJSON(inner, simplifyVector = FALSE),
      error = function(e) NULL
    )
    fobj <- adv_normalize_parsed(fobj)
    if (!is.null(fobj)) return(fobj)
  }
  NULL
}

adv_validate_required_keys <- function(x, required_keys) {
  if (is.null(x)) return(NULL)
  ok <- vapply(required_keys, function(k) {
    v <- x[[k]]
    !is.null(v) && length(v) >= 1
  }, logical(1))
  if (all(ok)) x else NULL
}

adv_scalar <- function(x, default = "") {
  x <- x %||% default
  x <- as.character(x)
  if (length(x) == 0L) return(as.character(default)[1L])
  trimws(x[1L])
}

# Single LLM call with 2-attempt retry.
# Returns list(parsed = <named list or NULL>, raw_text = <character>)
adv_request_json <- function(messages,
                              required_keys,
                              schema        = NULL,
                              retry_instr   = "Return only valid JSON.",
                              provider,
                              llm_api_key,
                              model,
                              max_tokens    = 2000,
                              temperature   = 0) {

  parsed   <- NULL
  raw_text <- NA_character_
  cur_max  <- if (identical(provider, "gemini")) max(4000L, max_tokens) else max_tokens

  for (attempt in 1:2) {
    res_call <- safe_chat_call(
      messages        = messages,
      provider        = provider,
      api_key         = llm_api_key,
      model           = model,
      max_tokens      = cur_max,
      temperature     = temperature,
      response_schema = schema
    )

    if (!is.null(res_call) &&
        is.character(res_call$content) &&
        nzchar(res_call$content)) {

      raw_text <- res_call$content
      parsed   <- adv_validate_required_keys(
        adv_parse_json_response(raw_text), required_keys
      )
      if (!is.null(parsed)) break

      if (identical(provider, "gemini") &&
          !is.null(res_call$finish_reason) &&
          identical(toupper(as.character(res_call$finish_reason)), "MAX_TOKENS")) {
        cur_max <- min(max(3000L, cur_max * 2L), 8192L)
      }
    }

    if (attempt == 1)
      messages <- c(messages, list(list(role = "user", content = retry_instr)))
  }

  list(parsed = parsed, raw_text = raw_text)
}


# ============================================================
# SECTION 2 – Module summary builders
# ============================================================

# Build module evidence used by S1/S2.
# Does NOT include phenotype or global_direction.
build_skeptic_module_evidence <- function(module_obj) {

  scalar <- function(x, na = "NA") {
    if (is.null(x) || length(x) == 0) return(na)
    x <- trimws(as.character(x[1]))
    if (!nzchar(x) || is.na(x) || toupper(x) %in% c("NA", "NULL", "UNKNOWN")) na else x
  }

  # Top metabolites (all annotation-selected top metabolites).
  top_df <- module_obj$topN_metabolites
  if (!is.null(top_df) && nrow(top_df) > 0) {
    top_lines <- apply(top_df, 1, function(r) {
      sprintf("  %s. %s [%s; class: %s]",
              scalar(r[["top_rank"]]),
              scalar(r[["Compound_name"]]),
              scalar(r[["priority_level"]]),
              scalar(r[["class"]], na = "not annotated"))
    })
    top_text <- paste(top_lines, collapse = "\n")
  } else {
    top_text <- "  (none)"
  }

  # Class structure
  class_df <- module_obj$other_metabolites
  if (!is.null(class_df) && nrow(class_df) > 0) {
    cls_lines <- apply(class_df, 1, function(r) {
      sprintf("  %s: %.0f%% [%s]",
              r[["class"]],
              as.numeric(r[["n_metabolites"]]) * 100,
              r[["priority_level"]])
    })
    cls_text <- paste(cls_lines, collapse = "\n")
  } else {
    cls_text <- "  (none)"
  }

  # Pathway coverage
  pc <- module_obj$pathway_coverage
  fmt_pw <- function(df, label) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0)
      return(paste0("  ", label, ": none matched"))
    rows <- apply(head(df, 3), 1, function(r) {
      sprintf("  %s (%d/%d, %.0f%%)",
              r[["pathway_name"]],
              as.integer(r[["n_hit"]]),
              as.integer(r[["n_metabolites"]]),
              as.numeric(r[["coverage"]]) * 100)
    })
    paste0("  ", label, ":\n", paste(rows, collapse = "\n"))
  }
  pw_text <- paste(
    fmt_pw(pc$kegg, "KEGG"),
    fmt_pw(pc$hmdb, "HMDB"),
    sep = "\n"
  )

  glue::glue(
    "--- Module evidence for biochemical review ---
Top metabolites (rank, name, priority):
{top_text}

Chemical class structure:
{cls_text}

Pathway coverage:
{pw_text}"
  )
}


# Build biological context block (no phenotype).
# Used by S1/S2/P1/P2.
build_context_block <- function(module_obj, study_context) {

  scalar <- function(x, na = "NA") {
    if (is.null(x) || length(x) == 0) return(na)
    x <- trimws(as.character(x[1]))
    if (!nzchar(x) || is.na(x)) na else x
  }

  # Origin profile
  op <- module_obj$origin_profile
  if (!is.null(op) && !is.na(op$n_annotated) && as.integer(op$n_annotated) > 0) {
    n  <- as.integer(op$n_annotated)
    en <- as.integer(op$endogenous_n %||% 0L)
    ex <- as.integer(op$exogenous_n  %||% 0L)
    uk <- as.integer(op$unknown_n    %||% 0L)
    origin_text <- glue::glue(
      "Origin profile (HMDB-annotated, n={n}):
  Endogenous: {round(en/n*100,1)}% ({en}/{n})
  Exogenous (dietary/microbial/drug/env): {round(ex/n*100,1)}% ({ex}/{n})
  Unknown: {round(uk/n*100,1)}% ({uk}/{n})"
    )
  } else {
    origin_text <- "Origin profile: not available."
  }

  # Species and biospecimen only (NO phenotype)
  sp  <- scalar(study_context$species)
  bio <- scalar(study_context$biospecimen)

  glue::glue(
    "--- Biological context ---
{origin_text}
Species: {sp}
Biospecimen: {bio}"
  )
}


# ============================================================
# SECTION 3 – Hypothesis state helpers
# ============================================================

# Format hypotheses for display (title + full biochemical summary + anchor + tension).
# Used by all agents (S1/S2/P1/P2 and Arbiter).
# claim has been removed; the full summary is the primary mechanistic statement.
format_hypotheses_compact <- function(hypotheses) {
  format_hypotheses_full(hypotheses)
}

# Format hypotheses with full mechanism (title + biochemical summary + anchor + tension).
# Used by all agents. No claim field (removed from pipeline).
format_hypotheses_full <- function(hypotheses) {
  if (length(hypotheses) == 0) return("(no hypotheses)")
  blocks <- lapply(hypotheses, function(h) {
    anc <- if (length(h$anchor_metabolites) > 0)
      paste(h$anchor_metabolites, collapse = ", ") else "none listed"
    ten <- if (length(h$tension_metabolites) > 0)
      paste(h$tension_metabolites, collapse = ", ") else "none listed"
    typ <- h$type %||% "not specified"
    if (is.na(typ) || !nzchar(typ)) typ <- "not specified"
    glue::glue(
      '{h$hypothesis_id} — "{h$title}"
  Type: {typ}
  Biochemical summary: {h$summary}
  Anchor metabolites: {anc}
  Tension metabolites: {ten}'
    )
  })
  paste(blocks, collapse = "\n\n")
}

# Apply Proposer revision decisions to the current hypothesis list.
# Returns updated hypothesis list with dropped ones removed.
apply_revisions <- function(hypotheses, revisions) {
  if (is.null(revisions) || length(revisions) == 0) return(hypotheses)

  # Build lookup: hypothesis_id -> revision
  rev_map <- setNames(revisions, sapply(revisions, `[[`, "hypothesis_id"))

  result <- list()
  for (h in hypotheses) {
    rev <- rev_map[[h$hypothesis_id]]
    if (is.null(rev)) {
      result[[length(result) + 1]] <- h  # no revision decision → keep as-is
      next
    }
    action <- trimws(tolower(rev$action %||% "keep"))
    if (action == "drop") next  # remove from list

    if (action == "modify") {
      # Apply field-by-field updates
      if (!is.null(rev$revised_type) && !is.na(rev$revised_type) && nzchar(rev$revised_type))
        h$type <- rev$revised_type
      if (!is.null(rev$revised_title) && nzchar(rev$revised_title))
        h$title   <- rev$revised_title
      if (!is.null(rev$revised_summary) && nzchar(rev$revised_summary))
        h$summary <- rev$revised_summary
      if (!is.null(rev$revised_anchor_metabolites) &&
          length(rev$revised_anchor_metabolites) > 0)
        h$anchor_metabolites  <- rev$revised_anchor_metabolites
      if (!is.null(rev$revised_tension_metabolites))
        h$tension_metabolites <- rev$revised_tension_metabolites
    }
    result[[length(result) + 1]] <- h
  }
  result
}

# Parse Proposer output into a clean list of hypothesis objects.
adv_top_metabolite_names <- function(module_obj) {
  top_df <- module_obj$topN_metabolites
  if (is.null(top_df) || !is.data.frame(top_df) || nrow(top_df) == 0 ||
      !"Compound_name" %in% colnames(top_df)) {
    return(character(0))
  }
  unique(trimws(as.character(top_df$Compound_name)))
}

restrict_to_top_metabolites <- function(x, valid_metabolites = NULL) {
  if (is.null(x)) return(character(0))
  x <- unique(trimws(as.character(x)))
  x <- x[nzchar(x)]
  if (is.null(valid_metabolites) || length(valid_metabolites) == 0) return(x)

  valid_metabolites <- unique(trimws(as.character(valid_metabolites)))
  valid_metabolites <- valid_metabolites[nzchar(valid_metabolites)]
  valid_key <- stats::setNames(valid_metabolites, tolower(valid_metabolites))
  matched <- valid_key[tolower(x)]
  unique(unname(matched[!is.na(matched)]))
}

normalize_hypothesis_type <- function(x) {
  allowed <- c(
    "Host-Endogenous",
    "Microbiome-Derived",
    "Diet-Derived",
    "Xenobiotic / Exogenous",
    "Host–Microbiome"
  )
  if (is.null(x) || length(x) == 0) return(NA_character_)
  raw <- trimws(as.character(x[1]))
  if (!nzchar(raw) || is.na(raw)) return(NA_character_)

  key <- tolower(gsub("[[:space:]]+", " ", raw))
  key <- gsub("\u2013|\u2014", "-", key)
  key <- gsub("/", " / ", key)
  key <- gsub("[[:space:]]+", " ", key)

  if (key %in% tolower(allowed)) {
    return(allowed[match(key, tolower(allowed))])
  }
  if (grepl("host", key) && grepl("microbiome|microbial", key)) return("Host–Microbiome")
  if (grepl("microbiome|microbial", key)) return("Microbiome-Derived")
  if (grepl("diet|food|nutr", key)) return("Diet-Derived")
  if (grepl("xenobiotic|exogenous|drug|pharma|environment|exposure", key)) {
    return("Xenobiotic / Exogenous")
  }
  if (grepl("host|endogenous", key)) return("Host-Endogenous")
  NA_character_
}

parse_proposer_output <- function(parsed, valid_metabolites = NULL) {
  if (is.null(parsed)) return(NULL)
  hyps_raw <- parsed$hypotheses
  if (is.null(hyps_raw) || length(hyps_raw) == 0) return(NULL)

  ensure_char_vec <- function(x) {
    if (is.null(x)) return(character(0))
    if (is.list(x)) return(as.character(unlist(x)))
    as.character(x)
  }

  out <- lapply(hyps_raw, function(h) {
    list(
      hypothesis_id        = trimws(as.character(h$hypothesis_id %||% "H?")),
      type                 = normalize_hypothesis_type(h$type %||% NA_character_),
      title                = trimws(as.character(h$title         %||% "")),
      summary              = trimws(as.character(h$summary       %||% "")),
      anchor_metabolites   = restrict_to_top_metabolites(
        ensure_char_vec(h$anchor_metabolites), valid_metabolites
      ),
      tension_metabolites  = restrict_to_top_metabolites(
        ensure_char_vec(h$tension_metabolites), valid_metabolites
      )
    )
  })
  clean_hypotheses(out)
}

clean_hypotheses <- function(hypotheses, max_hypotheses = 3L) {
  if (is.null(hypotheses) || length(hypotheses) == 0) return(list())

  valid <- Filter(function(h) {
    nzchar(trimws(h$title %||% "")) &&
      nzchar(trimws(h$summary %||% "")) &&
      length(h$anchor_metabolites) > 0
  }, hypotheses)

  if (length(valid) == 0) return(list())
  valid <- valid[seq_len(min(length(valid), max_hypotheses))]
  lapply(seq_along(valid), function(i) {
    valid[[i]]$hypothesis_id <- paste0("H", i)
    valid[[i]]
  })
}

# Parse Skeptic output into a clean list of challenge objects.
parse_skeptic_output <- function(parsed, valid_hypothesis_ids = NULL) {
  if (is.null(parsed)) return(list())
  chs_raw <- parsed$challenges
  if (is.null(chs_raw) || length(chs_raw) == 0) return(list())

  allowed_types <- c(
    "weak_evidence", "tension_unresolved", "origin_inconsistency",
    "species_biospecimen", "inadequate_revision", "no_challenge"
  )

  out <- lapply(chs_raw, function(c) {
    ct <- trimws(tolower(as.character(c$challenge_type %||% "")))
    sv <- trimws(tolower(as.character(c$severity      %||% "minor")))
    # Bidirectional normalization: challenge_type and severity must agree on no_challenge.
    # If either side says "no_challenge", force both sides to "no_challenge".
    if (identical(ct, "no_challenge") || identical(sv, "no_challenge")) {
      ct <- "no_challenge"
      sv <- "no_challenge"
    }
    list(
      target_hypothesis_id = trimws(as.character(c$target_hypothesis_id %||% "")),
      challenge_type       = ct,
      challenge            = trimws(as.character(c$challenge             %||% "")),
      severity             = sv
    )
  })
  # "no_challenge" is a valid severity alongside "major" and "minor".
  # no_challenge entries are KEPT in the parsed list so debate history
  # builders can show S1/S2 explicitly found no issue with a hypothesis.
  valid_severities <- c("major", "minor", "no_challenge")
  Filter(function(ch) {
    if (!ch$challenge_type %in% allowed_types) return(FALSE)
    if (!ch$severity %in% valid_severities)    return(FALSE)
    if (!nzchar(ch$target_hypothesis_id))      return(FALSE)
    # Allow empty challenge text only for no_challenge entries
    if (!nzchar(ch$challenge) && !identical(ch$severity, "no_challenge")) return(FALSE)
    if (is.null(valid_hypothesis_ids)) return(TRUE)
    parts <- strsplit(ch$target_hypothesis_id, "_", fixed = TRUE)[[1]]
    all(parts %in% valid_hypothesis_ids)
  }, out)
}

# Parse Proposer revision output into a clean list of revision objects.
parse_revision_output <- function(parsed,
                                  valid_hypothesis_ids = NULL,
                                  valid_metabolites = NULL) {
  if (is.null(parsed)) return(list())
  revs_raw <- parsed$revisions
  if (is.null(revs_raw) || length(revs_raw) == 0) return(list())

  ensure_char_vec <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.list(x)) return(as.character(unlist(x)))
    as.character(x)
  }

  out <- lapply(revs_raw, function(r) {
    revised_type <- normalize_hypothesis_type(r$revised_type %||% NA_character_)
    if (is.na(revised_type)) revised_type <- ""
    list(
      hypothesis_id               = trimws(as.character(r$hypothesis_id %||% "")),
      action                      = trimws(tolower(as.character(r$action %||% "keep"))),
      justification               = trimws(as.character(r$justification %||% "")),
      revised_type                = revised_type,
      revised_title               = trimws(as.character(r$revised_title %||% "")),
      revised_summary             = trimws(as.character(r$revised_summary %||% "")),
      revised_anchor_metabolites  = restrict_to_top_metabolites(
        ensure_char_vec(r$revised_anchor_metabolites), valid_metabolites
      ),
      revised_tension_metabolites = restrict_to_top_metabolites(
        ensure_char_vec(r$revised_tension_metabolites), valid_metabolites
      )
    )
  })
  Filter(function(r) {
    r$action %in% c("keep", "modify", "drop") &&
      nzchar(r$hypothesis_id) &&
      (is.null(valid_hypothesis_ids) || r$hypothesis_id %in% valid_hypothesis_ids)
  }, out)
}


# ============================================================
# SECTION 3b – Challenge routing helpers
# ============================================================

# Check whether a challenge object targets a given hypothesis_id.
# Handles both exact targets ("H1") and compound targets ("H1_H2") defensively.
#
# Returns TRUE if hid is the sole target OR one of the joined targets.
challenge_targets_hid <- function(ch, hid) {
  tid   <- ch$target_hypothesis_id %||% ""
  parts <- strsplit(tid, "_", fixed = TRUE)[[1]]
  identical(tid, hid) || hid %in% parts
}

# Format one challenge for display in a Proposer or Arbiter prompt.
# Adds a "[cross-hypothesis: H1_H2]" annotation when the challenge is compound,
# so the recipient knows the same challenge also applies to another hypothesis.
format_challenge_line <- function(ch, current_hid) {
  tid       <- ch$target_hypothesis_id %||% ""
  cross_tag <- if (!identical(tid, current_hid))
    sprintf(" [cross-hypothesis target: %s]", tid) else ""
  sprintf("  [%s] %s%s: %s",
          toupper(ch$severity %||% ""), ch$challenge_type %||% "", cross_tag,
          ch$challenge %||% "")
}


# ============================================================
# SECTION 4 – Round 2 decision logic
# ============================================================

#' Decide Whether to Run Skeptic 2 + Proposer 3
#'
#' @param surviving List of hypothesis objects surviving after Proposer revision 1.
#' @param s1_challenges List of challenge objects from Skeptic 1.
#' @param p1_revisions  List of revision objects from Proposer 1.
#' @return Logical.
should_run_round2 <- function(surviving, s1_challenges, p1_revisions) {

  # n=0: all dropped, go straight to Arbiter
  if (length(surviving) == 0) return(FALSE)

  # Build lookup maps
  rev_map <- setNames(
    lapply(p1_revisions, `[[`, "action"),
    sapply(p1_revisions, `[[`, "hypothesis_id")
  )

  # For each surviving hypothesis, find its max severity in S1.
  # Uses challenge_targets_hid() defensively for older compound targets.
  # "no_challenge" severity is treated the same as "none" — does not trigger round 2.
  max_sev <- function(hid) {
    relevant <- Filter(function(c) challenge_targets_hid(c, hid), s1_challenges)
    relevant_real <- Filter(
      function(c) !identical(c$severity, "no_challenge"), relevant
    )
    if (length(relevant_real) == 0) return("none")
    if (any(sapply(relevant_real, function(c) identical(c$severity, "major")))) "major" else "minor"
  }

  actions      <- sapply(surviving, function(h) rev_map[[h$hypothesis_id]] %||% "keep")
  severities   <- sapply(surviving, function(h) max_sev(h$hypothesis_id))

  # n=1 special case
  if (length(surviving) == 1) {
    # modify: content changed, Skeptic should verify the revision is adequate
    if (identical(actions[1], "modify")) return(TRUE)
    # keep + major: not yet addressed, needs second round
    return(identical(actions[1], "keep") && identical(severities[1], "major"))
    # keep + minor/none: no material challenge, skip round 2
  }

  # n>=2: trigger if any revision OR any keep-with-major
  any(actions == "modify") || any(actions == "keep" & severities == "major")
}


# ============================================================
# SECTION 5 – Debate history formatters
# ============================================================

# Build the Round 1 history block shown to S2.
# S2 needs to see:
#   - The original P0 hypothesis summary
#   - All S1 challenges (full text, including no_challenge)
#   - P1's decision and the revised summary (if MODIFY)
# This prevents S2 from re-raising already-addressed issues and helps it
# focus on whether P1's revisions introduced new problems.
build_round1_history_for_s2 <- function(hypotheses_v0,
                                          s1_challenges,
                                          p1_revisions,
                                          hypotheses_v1) {

  rev_map <- setNames(p1_revisions, sapply(p1_revisions, `[[`, "hypothesis_id"))

  fmt_vec <- function(v) {
    v <- trimws(v %||% character(0))
    v <- v[nzchar(v)]
    if (length(v) == 0) "none" else paste(v, collapse = ", ")
  }

  blocks <- lapply(hypotheses_v0, function(h) {
    hid <- h$hypothesis_id
    lines <- c(
      sprintf('%s — "%s" (P0 original)', hid, h$title),
      sprintf("  Original type:     %s", h$type %||% "not specified"),
      sprintf("  Original summary:  %s", h$summary %||% ""),
      sprintf("  Original anchors:  %s", fmt_vec(h$anchor_metabolites)),
      sprintf("  Original tensions: %s", fmt_vec(h$tension_metabolites))
    )

    # S1 challenges for this hypothesis
    chs <- Filter(function(c) challenge_targets_hid(c, hid), s1_challenges)
    if (length(chs) > 0) {
      lines <- c(lines, "  S1 challenges:")
      for (c in chs)
        lines <- c(lines, sprintf("    %s", format_challenge_line(c, hid)))
    } else {
      lines <- c(lines, "  S1: (no challenge objects returned for this hypothesis)")
    }

    # P1 decision — include full revised state
    rev <- rev_map[[hid]]
    if (!is.null(rev)) {
      lines <- c(lines, sprintf("  P1 decision (%s): %s",
                                toupper(rev$action), rev$justification))
      if (identical(rev$action, "modify")) {
        if (nzchar(rev$revised_type %||% ""))
          lines <- c(lines, sprintf("    → Revised type:     %s", rev$revised_type))
        if (nzchar(rev$revised_title %||% ""))
          lines <- c(lines, sprintf("    → Revised title:    \"%s\"", rev$revised_title))
        h_v1 <- Filter(function(x) x$hypothesis_id == hid, hypotheses_v1)
        if (length(h_v1) > 0) {
          if (nzchar(h_v1[[1]]$summary %||% ""))
            lines <- c(lines, sprintf("    → Revised summary:  %s", h_v1[[1]]$summary))
          rv_anc <- fmt_vec(h_v1[[1]]$anchor_metabolites)
          rv_ten <- fmt_vec(h_v1[[1]]$tension_metabolites)
          lines <- c(lines,
                     sprintf("    → Revised anchors:  %s", rv_anc),
                     sprintf("    → Revised tensions: %s", rv_ten))
        }
      } else if (identical(rev$action, "drop")) {
        lines <- c(lines, "    → DROPPED in Round 1")
      }
    } else {
      lines <- c(lines, "  P1: (no revision decision recorded)")
    }

    paste(lines, collapse = "\n")
  })

  paste(
    "--- Round 1 full record (for S2 reference) ---",
    paste(blocks, collapse = "\n\n"),
    sep = "\n"
  )
}


# Build the comprehensive Round 1 context shown to P2.
# P2 needs to see P0 originals + S1 challenges + P1 decisions so it does not
# inadvertently revert to the P0 state when revising P1's work.
build_round1_history_for_p2 <- function(hypotheses_v0,
                                          s1_challenges,
                                          p1_revisions,
                                          hypotheses_v1) {

  rev_map <- setNames(p1_revisions, sapply(p1_revisions, `[[`, "hypothesis_id"))

  blocks <- lapply(hypotheses_v0, function(h) {
    hid   <- h$hypothesis_id
    anc   <- paste(h$anchor_metabolites %||% character(0), collapse = ", ")
    ten   <- paste(h$tension_metabolites %||% character(0), collapse = ", ")
    lines <- c(
      sprintf('%s — "%s" (P0 original)', hid, h$title),
      sprintf("  Original type: %s", h$type %||% "not specified"),
      sprintf("  Original summary: %s", h$summary %||% ""),
      sprintf("  Original anchors: %s", if (nzchar(anc)) anc else "none"),
      sprintf("  Original tensions: %s", if (nzchar(ten)) ten else "none")
    )

    chs <- Filter(function(c) challenge_targets_hid(c, hid), s1_challenges)
    if (length(chs) > 0) {
      lines <- c(lines, "  S1 challenges:")
      for (c in chs)
        lines <- c(lines, sprintf("    %s", format_challenge_line(c, hid)))
    } else {
      lines <- c(lines, "  S1: no challenges")
    }

    rev <- rev_map[[hid]]
    if (!is.null(rev)) {
      lines <- c(lines, sprintf("  P1 (%s): %s", toupper(rev$action), rev$justification))
      if (identical(rev$action, "modify")) {
        if (nzchar(rev$revised_type %||% ""))
          lines <- c(lines, sprintf("    → Revised type: %s", rev$revised_type))
        if (nzchar(rev$revised_title %||% ""))
          lines <- c(lines, sprintf("    → Revised title: \"%s\"", rev$revised_title))
        h_v1 <- Filter(function(x) x$hypothesis_id == hid, hypotheses_v1)
        if (length(h_v1) > 0 && nzchar(h_v1[[1]]$summary %||% ""))
          lines <- c(lines, sprintf("    → Revised summary: %s", h_v1[[1]]$summary))
        anc_v1 <- paste(h_v1[[1]]$anchor_metabolites %||% character(0), collapse = ", ")
        if (nzchar(anc_v1))
          lines <- c(lines, sprintf("    → Revised anchors: %s", anc_v1))
        ten_v1 <- paste(h_v1[[1]]$tension_metabolites %||% character(0), collapse = ", ")
        if (nzchar(ten_v1))
          lines <- c(lines, sprintf("    → Revised tensions: %s", ten_v1))
      } else if (identical(rev$action, "drop")) {
        lines <- c(lines, "    → DROPPED in Round 1")
      }
    }

    paste(lines, collapse = "\n")
  })

  paste(
    "--- Round 1 full record (P0 originals + S1 challenges + P1 decisions) ---",
    paste(blocks, collapse = "\n\n"),
    sep = "\n"
  )
}


build_debate_history_digest <- function(hypotheses_v0,
                                         s1_challenges,
                                         p1_revisions,
                                         hypotheses_v1,
                                         s2_challenges = NULL,
                                         p2_revisions  = NULL,
                                         hypotheses_final = NULL,
                                         include_ids = NULL) {

  rev_map_r1 <- setNames(p1_revisions, sapply(p1_revisions, `[[`, "hypothesis_id"))
  rev_map_r2 <- if (!is.null(p2_revisions))
    setNames(p2_revisions, sapply(p2_revisions, `[[`, "hypothesis_id")) else list()

  surviving_ids <- if (!is.null(hypotheses_final))
    sapply(hypotheses_final, `[[`, "hypothesis_id") else character(0)

  hypotheses_for_digest <- hypotheses_v0
  if (!is.null(include_ids) && length(include_ids) > 0) {
    hypotheses_for_digest <- Filter(
      function(h) h$hypothesis_id %in% include_ids,
      hypotheses_for_digest
    )
  }

  blocks <- lapply(hypotheses_for_digest, function(h) {
    hid <- h$hypothesis_id
    lines <- c(
      sprintf('%s — "%s" [P0 original]', hid, h$title),
      sprintf("  P0 type: %s", h$type %||% "not specified"),
      sprintf("  P0 summary: %s", h$summary %||% "")
    )

    # Round 1 challenges for this hypothesis (compound targets included)
    r1_ch <- Filter(function(c) challenge_targets_hid(c, hid), s1_challenges)
    if (length(r1_ch) > 0) {
      lines <- c(lines, "  Round 1 — Skeptic challenges:")
      for (c in r1_ch)
        lines <- c(lines, paste0("    ", format_challenge_line(c, hid)))
    } else {
      lines <- c(lines, "  Round 1 — Skeptic: no challenges raised")
    }

    # Round 1 Proposer revision decision
    r1_rev <- rev_map_r1[[hid]]
    if (!is.null(r1_rev)) {
      action_label <- toupper(r1_rev$action)
      lines <- c(lines, sprintf("  Round 1 — Proposer (%s): %s",
                                action_label, r1_rev$justification))
      if (r1_rev$action == "modify") {
        if (nzchar(r1_rev$revised_type %||% ""))
          lines <- c(lines, sprintf("    → Revised type: %s", r1_rev$revised_type))
        if (nzchar(r1_rev$revised_title %||% ""))
          lines <- c(lines, sprintf("    → Revised title: \"%s\"", r1_rev$revised_title))
        # Retrieve P1-revised summary/anchors/tensions from hypotheses_v1 if available
        h_v1 <- Filter(function(x) x$hypothesis_id == hid, hypotheses_v1)
        if (length(h_v1) > 0 && nzchar(h_v1[[1]]$summary %||% ""))
          lines <- c(lines, sprintf("    → Revised summary: %s", h_v1[[1]]$summary))
        if (length(h_v1) > 0) {
          anc_v1 <- paste(h_v1[[1]]$anchor_metabolites %||% character(0), collapse = ", ")
          if (nzchar(anc_v1))
            lines <- c(lines, sprintf("    → Revised anchors: %s", anc_v1))
          ten_v1 <- paste(h_v1[[1]]$tension_metabolites %||% character(0), collapse = ", ")
          if (nzchar(ten_v1))
            lines <- c(lines, sprintf("    → Revised tensions: %s", ten_v1))
        }
      }
    }

    if (!is.null(r1_rev) && identical(r1_rev$action, "drop")) {
      lines <- c(lines, "  STATUS: Dropped in Round 1")
      return(paste(lines, collapse = "\n"))
    }

    # Round 2 (if ran)
    if (!is.null(s2_challenges)) {
      r2_ch <- Filter(function(c) challenge_targets_hid(c, hid), s2_challenges)
      if (length(r2_ch) > 0) {
        lines <- c(lines, "  Round 2 — Skeptic challenges:")
        for (c in r2_ch)
          lines <- c(lines, paste0("    ", format_challenge_line(c, hid)))
      } else {
        lines <- c(lines, "  Round 2 — Skeptic: no further challenges")
      }

      r2_rev <- rev_map_r2[[hid]]
      if (!is.null(r2_rev)) {
        action_label <- toupper(r2_rev$action)
        lines <- c(lines, sprintf("  Round 2 — Proposer (%s): %s",
                                  action_label, r2_rev$justification))
        if (r2_rev$action == "modify") {
          if (nzchar(r2_rev$revised_type %||% ""))
            lines <- c(lines, sprintf("    → Revised type: %s", r2_rev$revised_type))
          if (nzchar(r2_rev$revised_title %||% ""))
            lines <- c(lines, sprintf("    → Revised title: \"%s\"", r2_rev$revised_title))
          # Retrieve P2-revised summary/anchors/tensions from hypotheses_final if available
          h_fin <- Filter(function(x) x$hypothesis_id == hid, hypotheses_final %||% list())
          if (length(h_fin) > 0 && nzchar(h_fin[[1]]$summary %||% ""))
            lines <- c(lines, sprintf("    → Revised summary: %s", h_fin[[1]]$summary))
          if (length(h_fin) > 0) {
            anc_fin <- paste(h_fin[[1]]$anchor_metabolites %||% character(0), collapse = ", ")
            if (nzchar(anc_fin))
              lines <- c(lines, sprintf("    → Revised anchors: %s", anc_fin))
            ten_fin <- paste(h_fin[[1]]$tension_metabolites %||% character(0), collapse = ", ")
            if (nzchar(ten_fin))
              lines <- c(lines, sprintf("    → Revised tensions: %s", ten_fin))
          }
        }
        if (identical(r2_rev$action, "drop")) {
          lines <- c(lines, "  STATUS: Dropped in Round 2")
          return(paste(lines, collapse = "\n"))
        }
      }
    }

    final_status <- if (hid %in% surviving_ids) "SURVIVED" else "DROPPED"
    lines <- c(lines, sprintf("  FINAL STATUS: %s", final_status))
    paste(lines, collapse = "\n")
  })

  paste(blocks, collapse = "\n\n")
}


# ============================================================
# SECTION 6 – Prompt builders
# ============================================================

# ---- P0: Proposer initial generation ----

build_proposer_system_prompt <- function() {
  glue::glue(
    "You are an expert in metabolomics and biochemical pathway interpretation.
Your task is to interpret one dysregulated metabolic module and, only when supported
by the data, split it into 1-3 candidate biochemical process hypotheses.The goal is NOT 
to invent competition. The goal is to identify specific biological processes supported 
by the module data.

You will receive:
- Signal_concentration (high | medium | low): whether importance is concentrated in a few metabolites
- Top metabolites: rank, name, chemical class, priority_level
- Class-level structure: class name, proportion of module, metabolite names, priority_level
- Pathway evidence from KEGG and HMDB

Number-of-hypotheses guidance:
- Generate 1 when early-ranked high-priority metabolites and pathway evidence point to one dominant process.
- Generate 2 when there are two well-supported biochemical processes anchored by different early-ranked metabolites.
- Generate 3 only when there are three clearly separable, well-supported processes.
- Do NOT create extra hypotheses just to fill a count.

For each hypothesis provide:
1. hypothesis_id : \"H1\", \"H2\", or \"H3\"
2. type          : one of the five origin/scope types below
3. title         : KEGG-style pathway/process name (see title rules below)
4. summary       : 4-6 sentences describing the biochemical mechanism (see summary rules below)
5. anchor_metabolites  : metabolite names that STRONGLY support this hypothesis
6. tension_metabolites : optional metabolite names that materially conflict with, or remain
                         difficult to explain under, this hypothesis

Type rules:
   - Choose exactly one:
     Host-Endogenous | Microbiome-Derived | Diet-Derived | Xenobiotic / Exogenous | Host–Microbiome
   - Type is an interpretation scope for the hypothesis, not a confidence score.
   - Assign type from the module evidence, especially top-ranked high-priority metabolites;
     use chemical class and KEGG/HMDB pathway evidence as support.
   - Assign type by the biological origin of the hypothesis's core substrate, dominant metabolite
     class, or central biochemical theme, not merely by the organism/enzyme performing a
     downstream step.
   - Host-Endogenous: intrinsic host metabolites/processes such as central metabolism, redox,
     lipid/amino acid/nucleotide metabolism, bile acids, steroids, or cofactors.
   - Microbiome-Derived: microbial products or transformations are the central signal.
   - Diet-Derived: food, nutrient, beverage, or plant-derived compounds, even if later modified
     by host enzymes.
   - Xenobiotic / Exogenous: drugs, pollutants, pesticides, industrial/consumer chemicals,
     heavy metals, or other non-dietary external exposures.
   - Host–Microbiome: the core theme requires both host and microbial metabolism, e.g. host
     conjugation of microbial products or microbial remodeling of host compounds.
   - If multiple labels fit, choose the one that best explains the top anchor metabolites.
   - The type, title, and summary must describe the same central biochemical theme.
   - Do not use or mention phenotype/disease when assigning type.

Title rules
   - Name the module in KEGG pathway style: concise, standardized, and anchored to a specific biochemical pathway or process.
   - Follow KEGG naming conventions:
     - Use \"X metabolism\" when the module spans multiple steps of a compound class (e.g., \"Tyrosine metabolism\", \"Sphingolipid metabolism\").
     - Use \"X biosynthesis\" or \"X degradation\" when directionality is clear from the metabolites.
     - Use a specific reaction or enzyme-level name when anchored by a single well-defined step.
     - Use a physiological or system-level name when appropriate.
   - Choose the most specific level of granularity supported by the early-ranked high-priority metabolites.
   - If top metabolites point to a specific sub-pathway, reflect that specificity rather than defaulting to the parent pathway name.
   - Use pathway evidence (KEGG/HMDB) to support or refine the title when it aligns with the metabolite evidence.
   - If no canonical pathway fits, use the most specific neutral description of the dominant chemical transformation.

Summary rules:
   - Write 4-6 sentences as a factual biochemical annotation of the module.
   - Describe what the metabolites do biochemically: enzymatic conversions, pathway steps, metabolite roles, and physiological functions — as if writing an entry in a biochemical database.
   - Base the explanation mainly on early-ranked high-priority metabolites; use class structure and pathway evidence only to support or qualify.
   - Keep origin/scope interpretation in the type field; in the summary, mention source or host/microbial
     processing only when it is necessary to explain the biochemical mechanism.
   - Do NOT comment on the study phenotype or disease relevance.
   - Do NOT reference module-level statistics (Signal_concentration).
   - Mention secondary signals only briefly at the end if needed.

Rules:
- Every hypothesis must be anchored by at least one high-priority metabolite
- anchor_metabolites and tension_metabolites MUST be selected only from the provided Top
  metabolites list. Class-level metabolites, pathway members, or other non-top metabolites
  may be used only as supporting context in the summary, not as anchors or tensions.
- anchor_metabolites should be as non-overlapping as the data allow; if two hypotheses rely on
  the same central anchor, they are probably not separate hypotheses
- tension_metabolites should be empty unless a metabolite creates a real biochemical conflict
  with the proposed mechanism. Do NOT list peripheral or merely non-central metabolites as
  tensions simply because they are not part of the main pathway.
- Do NOT reference study phenotype or disease
- Do NOT invent metabolites not present in the module data
- Do NOT invent unsupported pathways or mechanisms

Return STRICT JSON only:
{{
  \"n_hypotheses\": 1,
  \"hypotheses\": [
    {{
      \"hypothesis_id\": \"H1\",
      \"type\": \"Host-Endogenous | Microbiome-Derived | Diet-Derived | Xenobiotic / Exogenous | Host–Microbiome\",
      \"title\": \"...\",
      \"summary\": \"...\",
      \"anchor_metabolites\": [\"...\"],
      \"tension_metabolites\": [\"...\"]
    }}
  ]
}}"
  )
}

# P0 user prompt reuses build_module_prompt() from 45_gen_prompt.R (same input format)


# ---- S1/S2: Skeptic ----

build_skeptic_system_prompt <- function(round = 1L) {
  round2_addon <- if (round == 2L) glue::glue(
    "\n\nADDITIONAL INSTRUCTION — ROUND 2:
You have access to Round 1 challenges and Proposer revision decisions.
- Do NOT repeat a challenge that was adequately addressed by a modification in Round 1
- If a Round 1 modification introduced new content, examine that new content for issues
- Use challenge_type \"inadequate_revision\" if a Round 1 major challenge was NOT adequately resolved
- Focus on genuinely new or unresolved issues, not previously settled ones"
  ) else ""

  glue::glue(
    "You are an expert in metabolomics and biochemical pathway interpretation.
Your task is to critically evaluate candidate biochemical hypotheses for a metabolic module.

You will receive:
- Module evidence (top metabolites, class structure, pathway coverage)
- Basic biological context: origin profile, species, biospecimen
- 1-3 hypotheses, each shown as: type, title, full biochemical summary, anchor metabolites, tension metabolites

Use your general biochemical and metabolomics knowledge to evaluate pathway plausibility,
metabolite-to-pathway fit, mechanism specificity, and whether the proposed anchors really
support the claim. You are explicitly allowed to use general external biochemical knowledge
about pathways, metabolites, enzyme logic, pathway crosstalk, biospecimen detectability, and
species plausibility.
{round2_addon}

For each hypothesis generate 0-3 specific biochemical challenges.
Only raise a challenge when it materially affects interpretation. If a hypothesis is
biochemically well-supported and no meaningful issue exists, do not invent a weakness;
either omit that hypothesis from the challenges array or use challenge_type \"no_challenge\".

Challenge types (use exactly one per challenge):
- no_challenge          : no material biochemical concern was identified; use only when no
                          weakness needs Proposer response. This is not a criticism.
- weak_evidence         : anchor metabolites do not strongly support the claimed pathway
- tension_unresolved    : the hypothesis lists tension metabolites but the summary does not
                           explain why they are difficult to reconcile with the mechanism
- origin_inconsistency  : the current origin profile does not support the hypothesis and its type
- species_biospecimen   : the claimed pathway is implausible for the stated species or biospecimen
- inadequate_revision   : round 2 only; a Round 1 revision failed to adequately address a
                           major Skeptic challenge

Severity rules:
- major        : the challenge undermines the core mechanism; the hypothesis cannot stand without addressing this
- minor        : a notable gap or weakness that should be acknowledged, but does not invalidate the hypothesis
- no_challenge : no material biochemical concern was identified for this hypothesis; use ONLY with
                 challenge_type \"no_challenge\". This tells the next Proposer and the Arbiter that
                 you explicitly reviewed this hypothesis and found no issue.

Anchor each challenge to specific metabolites, pathway names, origin/type evidence, species,
biospecimen, or biochemical mechanisms.

Return STRICT JSON only. Do NOT force one challenge per hypothesis. If no substantive
challenge exists for a hypothesis, either omit it from the challenges array or return one
no_challenge entry for that hypothesis.
{{
  \"challenges\": [
    {{
      \"target_hypothesis_id\": \"H1\",
      \"challenge_type\": \"weak_evidence | tension_unresolved | origin_inconsistency | species_biospecimen | inadequate_revision | no_challenge\",
      \"challenge\": \"2-3 sentences with specific evidence, or 'No material biochemical concern identified.' for no_challenge\",
      \"severity\": \"major | minor | no_challenge\"
    }}
  ]
}}"
  )
}

build_skeptic_user_prompt <- function(module_evidence,
                                       context_block,
                                       hypotheses,
                                       history_text = NULL) {
  hyp_text <- format_hypotheses_compact(hypotheses)

  history_section <- if (!is.null(history_text) && nzchar(history_text))
    glue::glue("\n\n--- Round 1 history (for reference) ---\n{history_text}") else ""

  glue::glue(
    "{module_evidence}

{context_block}

--- Current hypotheses to evaluate ---
{hyp_text}{history_section}"
  )
}


# ---- P1/P2: Proposer revision ----

build_proposer_revision_system_prompt <- function() {
  glue::glue(
    "You are an expert in metabolomics and biochemical pathway interpretation.
Your task is to respond to critical challenges against proposed biological hypotheses.

You will receive:
- Module evidence (top metabolites, class structure, pathway coverage) — biochemical grounding
- Biological context: origin profile, species, biospecimen (no phenotype)
- Current hypotheses with their type, full biochemical summaries, and anchor/tension metabolites
- Skeptic challenges (with severity: major or minor)

For each hypothesis, decide: keep | modify | drop

Rules:
- keep   : hypothesis is valid; provide a counter-argument for every major challenge.
           You may cite the module evidence, general biochemical knowledge
           (pathway logic, metabolite roles, species or biospecimen biology), or point
           to specific weaknesses in the Skeptic's evidence.
- modify : revise P0-generated content to address reasonable Skeptic challenges; provide
           updated type, title, summary, anchor_metabolites, and tension_metabolites.
           A modification must change
           either the type, title, OR the mechanism text — renaming alone is not valid.
           revised_type must be exactly one of:
           Host-Endogenous | Microbiome-Derived | Diet-Derived | Xenobiotic / Exogenous | Host–Microbiome.
           Revised anchor_metabolites and revised_tension_metabolites MUST be selected only
           from the current hypotheses' top-metabolite-derived anchors/tensions or from the
           module's Top metabolites list. Do NOT promote class-level or other non-top
           metabolites into anchors or tensions.
           Remove tension_metabolites when they are merely peripheral rather than genuinely
           conflicting with the revised mechanism.
- drop   : hypothesis is fundamentally flawed; provide a clear reason.

You MUST respond to every major-severity challenge explicitly.
Silence on a major challenge is not acceptable.
You may drop all hypotheses if you judge them all to be fundamentally unsupported —
the Arbiter will handle the all-dropped case.

Return STRICT JSON only:
{{
  \"revisions\": [
    {{
      \"hypothesis_id\": \"H1\",
      \"action\": \"keep | modify | drop\",
      \"justification\": \"...\",
      \"revised_type\": \"Host-Endogenous | Microbiome-Derived | Diet-Derived | Xenobiotic / Exogenous | Host–Microbiome\",
      \"revised_title\": \"...\",
      \"revised_summary\": \"...\",
      \"revised_anchor_metabolites\": [\"...\"],
      \"revised_tension_metabolites\": [\"...\"]
    }}
  ]
}}"
  )
}

build_proposer_revision_user_prompt <- function(hypotheses,
                                                 challenges,
                                                 module_evidence = NULL,
                                                 context_block   = NULL,
                                                 history_text    = NULL) {
  hyp_text <- format_hypotheses_full(hypotheses)

  # Format challenges grouped by hypothesis.
  # challenge_targets_hid() ensures compound targets ("H1_H2") are distributed
  # to both H1 and H2 entries.  format_challenge_line() adds a cross-hypothesis
  # annotation so the Proposer can see which other hypothesis shares the issue.
  ch_blocks <- lapply(hypotheses, function(h) {
    hid  <- h$hypothesis_id
    chs  <- Filter(function(c) challenge_targets_hid(c, hid), challenges)
    # Strip no_challenge entries — Proposer need not respond to "no issues found"
    real_chs <- Filter(function(c) !identical(c$severity, "no_challenge"), chs)
    if (length(real_chs) == 0) return(sprintf("  %s: no challenges requiring response", hid))
    ch_lines <- sapply(real_chs, function(c) format_challenge_line(c, hid))
    paste(c(sprintf("  Challenges for %s:", hid), ch_lines), collapse = "\n")
  })
  ch_text <- paste(ch_blocks, collapse = "\n\n")

  history_section <- if (!is.null(history_text) && nzchar(history_text))
    glue::glue("\n\n--- Prior round history (your previous commitments) ---\n{history_text}") else ""

  # Module evidence + context block (raw biochemical grounding for rebuttal)
  evidence_section <- if (!is.null(module_evidence) && nzchar(module_evidence))
    glue::glue("{module_evidence}\n\n") else ""
  context_section  <- if (!is.null(context_block)   && nzchar(context_block))
    glue::glue("{context_block}\n\n") else ""

  glue::glue(
    "{evidence_section}{context_section}--- Current hypotheses (full details) ---
{hyp_text}

--- Skeptic challenges ---
{ch_text}{history_section}"
  )
}


# ---- A0: Arbiter ----

build_arbiter_system_prompt <- function() {
  glue::glue(
    "You are an expert in metabolomics and biochemical pathway interpretation.
Your task is to make final verdicts on candidate biological hypotheses for a metabolic module,
and to assess phenotypic relevance and confidence for surviving hypotheses.

You will receive:
- Module evidence (top metabolites, class structure, pathway coverage) — raw biochemical grounding
- Final surviving hypotheses in full (or all original hypotheses if all were dropped)
- A complete debate history digest that includes for each hypothesis:
    • P0 original type, title, and full biochemical summary
    • All Skeptic challenges with severity and full challenge text (including no_challenge entries
      which indicate the Skeptic explicitly found no issue)
    • Proposer revision decisions with justifications, revised types, titles, and summaries
    • Round 2 challenges and revisions (if the second round was triggered)
    • Final survival status
- Origin profile, species, biospecimen, and phenotype (full study context)
- Pre-computed phenotype relationship summaries based on PubMed evidence

Your role is to make a final verdict on the CURRENT (final) state of each hypothesis, informed
by the complete debate record. Specifically:
- Assess whether the Skeptic's challenges were biochemically valid
- Assess whether the Proposer's revisions adequately addressed the challenges
- Use the debate trajectory to calibrate confidence (well-contested and well-defended = higher confidence;
  major challenges that were only superficially addressed = lower confidence or rejection)

Note: global_direction is intentionally excluded — module-level direction is too complex to use
as a mechanistic constraint.

Verdict options:
- accepted  : biologically plausible and useful as the final module annotation, with
              evidence strength captured by confidence_label
- rejected  : biochemically implausible or contradicted by data within the hypotheses
              currently under your review, or too unsupported to be useful

For accepted hypotheses, also provide:
- confidence_label       : \"high\" | \"medium\" | \"low\"
- confidence_rationale   : 2-3 sentences explaining the strength of evidence supporting
  the accepted annotation.

Confidence reflects evidence strength for the accepted annotation. It is secondary to
status and should not be used as an independent verdict.

Confidence rules:
- high   : Strong module evidence AND adequate phenotype/literature support. Top metabolites,
           hypothesis type, class/pathway evidence, species, biospecimen, debate record, and
           PubMed phenotype evidence are coherent. High confidence requires at least indirect
           PubMed support; if literature_support is insufficient, confidence_label must not
           be high.
- medium : Reasonable module interpretation, but with partial evidence, indirect pathway support,
           incomplete explanation, limited literature support, or minor unresolved limitations.
- low    : Weak but still useful interpretation. Evidence is sparse, indirect, or partly
           unresolved, but the hypothesis remains more informative than rejecting the module
           annotation entirely.

Origin rule:
- Origin profile is an interpretation-scope check, not an automatic penalty.
- High exogenous fraction does NOT automatically lower confidence if the hypothesis type
  is Diet-Derived, Xenobiotic / Exogenous, Microbiome-Derived, or otherwise consistent with
  the origin profile.
- Lower confidence or reject only when the hypothesis type or summary conflicts with the
  origin profile.

Phenotype evidence rule:
- If phenotype is provided but PubMed phenotype summaries state insufficient evidence,
  you MAY infer a cautious phenotype relationship from general biochemical knowledge.
  When doing so, explicitly state that the relationship is inferred and lacks sufficient
  direct literature support.
- If phenotype is not provided, do NOT infer a phenotype relationship.

Recovery rule (when all hypotheses were dropped in debate):
- You MAY recover at most ONE hypothesis to \"accepted\" status with low confidence if the
  Proposer was overly aggressive and the hypothesis remains a useful, biologically plausible
  annotation
- If you agree all drops were justified, return all as \"rejected\"

Return STRICT JSON only:
{{
  \"verdicts\": [
    {{
      \"hypothesis_id\": \"H1\",
      \"status\": \"accepted | rejected\",
      \"verdict_rationale\": \"...\",
      \"confidence_label\": \"high | medium | low\",
      \"confidence_rationale\": \"...\"
    }}
  ],
  \"residual_ambiguity\": \"One sentence describing unresolved tension across hypotheses, or null.\"
}}"
  )
}

build_arbiter_user_prompt <- function(hypotheses_final,
                                       hypotheses_v0,
                                       debate_digest,
                                       module_obj,
                                       study_context,
                                       module_evidence         = NULL,
                                       phenotype_relationships = NULL) {

  scalar <- function(x, na = "NA") {
    if (is.null(x) || length(x) == 0) return(na)
    x <- trimws(as.character(x[1]))
    if (!nzchar(x) || is.na(x)) na else x
  }

  # Final hypotheses section (full format for Arbiter)
  if (length(hypotheses_final) == 0) {
    hyp_section <- glue::glue(
      "NOTE: All hypotheses were dropped during the debate.
The following are the ORIGINAL proposed hypotheses for your consideration:

{format_hypotheses_full(hypotheses_v0)}"
    )
  } else {
    hyp_section <- glue::glue(
      "Final surviving hypotheses:

{format_hypotheses_full(hypotheses_final)}"
    )
  }

  # Origin profile
  op <- module_obj$origin_profile
  if (!is.null(op) && !is.na(op$n_annotated) && as.integer(op$n_annotated) > 0) {
    n  <- as.integer(op$n_annotated)
    en <- as.integer(op$endogenous_n %||% 0L)
    ex <- as.integer(op$exogenous_n  %||% 0L)
    uk <- as.integer(op$unknown_n    %||% 0L)
    origin_text <- glue::glue(
      "Origin profile (n={n}): endogenous {round(en/n*100,1)}% | exogenous {round(ex/n*100,1)}% | unknown {round(uk/n*100,1)}%"
    )
  } else {
    origin_text <- "Origin profile: not available."
  }

  # Study context (full, including phenotype)
  ctx_text <- glue::glue(
    "Study context:
  Species:     {scalar(study_context$species)}
  Biospecimen: {scalar(study_context$biospecimen)}
  Phenotype:   {scalar(study_context$phenotype)}"
  )

  phenotype_text  <- format_phenotype_relationships_for_arbiter(phenotype_relationships)

  # Raw module evidence section (biochemical grounding, no global_direction)
  evidence_section <- if (!is.null(module_evidence) && nzchar(module_evidence))
    glue::glue("\n{module_evidence}\n") else ""

  glue::glue(
    "{evidence_section}
{hyp_section}

--- Debate history ---
{debate_digest}

--- Study context & biology ---
{origin_text}
{ctx_text}

--- PubMed-derived phenotype relationship summaries ---
{phenotype_text}"
  )
}


# ============================================================
# SECTION 7 – PubMed orchestration for advanced mode
# ============================================================

# Two-layer PubMed search strategy:
#   Layer 1 (shared): all unique anchor metabolites + phenotype — called ONCE.
#   Layer 2 (per-hypothesis): hypothesis title + phenotype — one call per hypothesis.
#   Anchor records are assigned only to hypotheses that use the matching anchors, then
#   merged with per-hypothesis title evidence and deduplicated by PMID.
#
# Searching and summarizing are kept separate:
#   search_pubmed_advanced()              → fetches raw PubMed evidence per hypothesis
#   summarize_pubmed_relationships_advanced() → LLM summarizes each hypothesis independently
search_pubmed_advanced <- function(hypotheses_final,
                                    phenotype            = NULL,
                                    NCBI_api_key         = NULL,
                                    retmax_per_query     = 3L,
                                    years                = 5) {

  if (length(hypotheses_final) == 0) return(list())

  clean_terms <- function(x) {
    if (is.null(x) || length(x) == 0) return(character(0))
    x <- unique(trimws(as.character(unlist(x, use.names = FALSE))))
    x[!is.na(x) & nzchar(x) & x != "NA"]
  }

  # Layer 1: search all unique anchors once, then assign matching source_terms
  # back to each hypothesis below.
  all_anchor_names <- clean_terms(lapply(hypotheses_final, function(h) h$anchor_metabolites))
  shared_anchors <- tryCatch(
    search_pubmed_for_annotation_evidence(
      module_title         = NULL,
      top_metabolite_names = all_anchor_names,
      phenotype            = phenotype,
      NCBI_api_key         = NCBI_api_key,
      retmax_per_query     = retmax_per_query,
      years                = years
    ),
    error = function(e) NULL
  )

  filter_pubmed_by_source_term <- function(evidence, terms) {
    terms <- clean_terms(terms)
    if (is.null(evidence) || length(evidence$pmids) == 0 || length(terms) == 0) {
      return(list(pmids = character(0), titles = character(0),
                  abstracts = character(0), source_type = character(0),
                  source_term = character(0)))
    }
    n <- min(length(evidence$pmids), length(evidence$titles), length(evidence$abstracts))
    if (n == 0) {
      return(list(pmids = character(0), titles = character(0),
                  abstracts = character(0), source_type = character(0),
                  source_term = character(0)))
    }
    source_term <- if (!is.null(evidence$source_term)) {
      as.character(evidence$source_term[seq_len(n)])
    } else {
      rep("", n)
    }
    keep <- tolower(trimws(source_term)) %in% tolower(terms)
    if (!any(keep, na.rm = TRUE)) {
      return(list(pmids = character(0), titles = character(0),
                  abstracts = character(0), source_type = character(0),
                  source_term = character(0)))
    }
    list(
      pmids       = as.character(evidence$pmids[seq_len(n)][keep]),
      titles      = as.character(evidence$titles[seq_len(n)][keep]),
      abstracts   = as.character(evidence$abstracts[seq_len(n)][keep]),
      source_type = if (!is.null(evidence$source_type)) {
        as.character(evidence$source_type[seq_len(n)][keep])
      } else {
        rep("unknown", sum(keep))
      },
      source_term = source_term[keep]
    )
  }

  # Layer 2: per-hypothesis title search + merge with that hypothesis's anchor records.
  per_hyp <- lapply(hypotheses_final, function(h) {
    title <- if (!is.null(h$title) && nzchar(h$title)) h$title else NULL
    hyp_specific <- if (!is.null(title)) {
      tryCatch(
        search_pubmed_for_annotation_evidence(
          module_title         = title,
          top_metabolite_names = NULL,
          phenotype            = phenotype,
          NCBI_api_key         = NCBI_api_key,
          retmax_per_query     = retmax_per_query,
          years                = years
        ),
        error = function(e) NULL
      )
    } else NULL

    hyp_anchor_evidence <- filter_pubmed_by_source_term(shared_anchors, h$anchor_metabolites)
    adv_merge_pubmed(hyp_specific, hyp_anchor_evidence)
  })
  names(per_hyp) <- sapply(hypotheses_final, `[[`, "hypothesis_id")
  per_hyp
}

# Merge two PubMed result lists, deduplicating by PMID.
# Hypothesis-specific results (a) take priority over shared results (b).
adv_merge_pubmed <- function(a, b) {
  to_df <- function(x) {
    if (is.null(x) || length(x$pmids) == 0) return(NULL)
    n <- min(length(x$pmids), length(x$titles), length(x$abstracts))
    if (n == 0) return(NULL)
    data.frame(
      pmid        = as.character(x$pmids[seq_len(n)]),
      title       = as.character(x$titles[seq_len(n)]),
      abstract    = as.character(x$abstracts[seq_len(n)]),
      source_type = if (!is.null(x$source_type)) as.character(x$source_type[seq_len(n)]) else rep("unknown", n),
      source_term = if (!is.null(x$source_term)) as.character(x$source_term[seq_len(n)]) else rep("", n),
      stringsAsFactors = FALSE
    )
  }
  da <- to_df(a); db <- to_df(b)
  combined <- rbind(da, db)
  if (is.null(combined) || nrow(combined) == 0)
    return(list(pmids = character(0), titles = character(0),
                abstracts = character(0), source_type = character(0), source_term = character(0)))
  combined <- combined[!duplicated(combined$pmid), , drop = FALSE]
  list(
    pmids       = combined$pmid,
    titles      = combined$title,
    abstracts   = combined$abstract,
    source_type = combined$source_type,
    source_term = combined$source_term
  )
}


build_pubmed_relationship_system_prompt <- function() {
  glue::glue(
    "You are an expert in metabolomics and biochemical pathway interpretation.
Your task is to summarize how one biochemical hypothesis relates to the study phenotype,
based ONLY on the PubMed evidence provided.

Return one field:

PHENOTYPE_RELATIONSHIP
- Write 2-4 sentences summarizing how the hypothesis metabolites or pathway relate to the
  study phenotype, based ONLY on the provided PubMed abstracts.
- If no PubMed evidence is available or relevant, write exactly:
  \"Insufficient literature evidence to establish phenotype relationship.\"
- Do NOT speculate beyond what is explicitly stated in the abstracts.
- Do NOT assess whether the hypothesis should be accepted or rejected.
- Do NOT reference origin (endogenous/exogenous), confidence, or module-level statistics.

Return STRICT JSON only:
{{
  \"phenotype_relationship\": \"...\",
  \"literature_support\": \"direct | indirect | insufficient\",
  \"evidence_note\": \"one concise sentence describing the strength or gap in the evidence\"
}}"
  )
}

build_pubmed_relationship_user_prompt <- function(hypothesis,
                                                   study_context,
                                                   pubmed_evidence = NULL) {
  scalar <- function(x, na = "NA") {
    if (is.null(x) || length(x) == 0) return(na)
    x <- trimws(as.character(x[1]))
    if (!nzchar(x) || is.na(x)) na else x
  }

  if (is.null(pubmed_evidence) || length(pubmed_evidence$pmids) == 0) {
    pubmed_text <- "No PubMed evidence retrieved."
  } else {
    n <- min(length(pubmed_evidence$pmids),
             length(pubmed_evidence$titles),
             length(pubmed_evidence$abstracts))
    papers <- vapply(seq_len(n), function(i) {
      abs <- trimws(as.character(pubmed_evidence$abstracts[i]))
      # Truncate abstract to 800 chars to keep prompt manageable
      if (nchar(abs) > 800) abs <- paste0(substr(abs, 1, 800), "...")
      sprintf("  [PMID:%s] %s\n    %s",
              pubmed_evidence$pmids[i],
              pubmed_evidence$titles[i],
              abs)
    }, character(1))
    pubmed_text <- paste(papers, collapse = "\n")
  }

  glue::glue(
    "--- Study phenotype ---
Species: {scalar(study_context$species)}
Biospecimen: {scalar(study_context$biospecimen)}
Phenotype: {scalar(study_context$phenotype)}

--- Hypothesis ---
ID: {hypothesis$hypothesis_id}
Type: {hypothesis$type %||% 'not specified'}
Title: {hypothesis$title}
Biochemical summary: {hypothesis$summary}

--- PubMed evidence ---
{pubmed_text}"
  )
}

parse_pubmed_relationship_output <- function(parsed) {
  if (is.null(parsed)) {
    return(list(
      phenotype_relationship = "Insufficient literature evidence to establish phenotype relationship.",
      literature_support = "insufficient",
      evidence_note = "No valid literature summary was returned."
    ))
  }
  support <- trimws(tolower(as.character(parsed$literature_support %||% "insufficient")))
  if (!support %in% c("direct", "indirect", "insufficient")) support <- "insufficient"
  rel <- trimws(as.character(parsed$phenotype_relationship %||% ""))
  if (!nzchar(rel)) rel <- "Insufficient literature evidence to establish phenotype relationship."
  list(
    phenotype_relationship = rel,
    literature_support = support,
    evidence_note = trimws(as.character(parsed$evidence_note %||% ""))
  )
}

adv_has_context_value <- function(x) {
  if (is.null(x) || length(x) == 0) return(FALSE)
  val <- trimws(as.character(x[1]))
  nzchar(val) && !is.na(val) && !identical(toupper(val), "NA")
}

adv_insufficient_phenotype_relationship <- function(reason = NULL) {
  note <- reason %||% "No phenotype-specific PubMed evidence was available."
  list(
    phenotype_relationship = "Insufficient literature evidence to establish phenotype relationship.",
    literature_support = "insufficient",
    evidence_note = note
  )
}

adv_inferable_phenotype_relationship <- function(reason = NULL) {
  note <- reason %||%
    "No phenotype-specific PubMed evidence was available; Arbiter may make a cautious general-knowledge inference."
  list(
    phenotype_relationship = paste(
      "No sufficient PubMed evidence was available for a direct phenotype relationship.",
      "A cautious phenotype relationship may be inferred from general biochemical knowledge if the study phenotype is provided;"
    ),
    literature_support = "insufficient",
    evidence_note = note
  )
}

adv_build_insufficient_relationships <- function(hypotheses, reason = NULL) {
  if (length(hypotheses) == 0) return(list())
  out <- lapply(hypotheses, function(h) adv_insufficient_phenotype_relationship(reason))
  names(out) <- vapply(hypotheses, `[[`, character(1), "hypothesis_id")
  out
}

adv_build_inferable_relationships <- function(hypotheses, reason = NULL) {
  if (length(hypotheses) == 0) return(list())
  out <- lapply(hypotheses, function(h) adv_inferable_phenotype_relationship(reason))
  names(out) <- vapply(hypotheses, `[[`, character(1), "hypothesis_id")
  out
}

summarize_pubmed_relationships_advanced <- function(hypotheses,
                                                     study_context,
                                                     pubmed_per_hyp = NULL,
                                                     provider,
                                                     llm_api_key,
                                                     model,
                                                     max_tokens = 1200L,
                                                     temperature = 0,
                                                     schema = NULL) {
  if (length(hypotheses) == 0) return(list(summaries = list(), raw_text = NA_character_))

  sys <- build_pubmed_relationship_system_prompt()
  raw_texts <- character(0)
  summaries <- lapply(hypotheses, function(h) {
    hid <- h$hypothesis_id
    pm <- if (!is.null(pubmed_per_hyp)) pubmed_per_hyp[[hid]] else NULL
    usr <- build_pubmed_relationship_user_prompt(h, study_context, pm)
    res <- adv_request_json(
      messages = list(list(role = "system", content = sys),
                      list(role = "user", content = usr)),
      required_keys = "phenotype_relationship",
      schema = schema,
      retry_instr = "Return only valid JSON with key \"phenotype_relationship\".",
      provider = provider,
      llm_api_key = llm_api_key,
      model = model,
      max_tokens = max_tokens,
      temperature = temperature
    )
    raw_texts <<- c(raw_texts, sprintf("[%s]\n%s", hid, res$raw_text %||% ""))
    parse_pubmed_relationship_output(res$parsed)
  })
  names(summaries) <- vapply(hypotheses, `[[`, character(1), "hypothesis_id")
  list(summaries = summaries, raw_text = paste(raw_texts, collapse = "\n\n"))
}

format_phenotype_relationships_for_arbiter <- function(phenotype_relationships) {
  if (is.null(phenotype_relationships) || length(phenotype_relationships) == 0) {
    return("Phenotype relationship summaries were not generated.")
  }
  blocks <- lapply(names(phenotype_relationships), function(hid) {
    x <- phenotype_relationships[[hid]]
    glue::glue(
      "[{hid}]
  phenotype_relationship: {x$phenotype_relationship %||% 'NA'}
  literature_support: {x$literature_support %||% 'insufficient'}
  evidence_note: {x$evidence_note %||% ''}"
    )
  })
  paste(blocks, collapse = "\n\n")
}

# Helper: collapse all PMIDs across ALL hypotheses (full-archive, same for every output row)
collapse_pubmed_pmids_advanced <- function(pph) {
  if (is.null(pph) || length(pph) == 0) return(NA_character_)
  all_pmids <- unique(unlist(lapply(pph, function(h) h$pmids %||% character(0))))
  all_pmids <- all_pmids[!is.na(all_pmids) & nzchar(trimws(all_pmids))]
  if (length(all_pmids) == 0) return(NA_character_)
  paste(all_pmids, collapse = "; ")
}

# Helper: serialize full pubmed_per_hyp to JSON (all hypotheses, full archive)
serialize_pubmed_evidence_advanced <- function(pph) {
  if (is.null(pph) || length(pph) == 0) return(NA_character_)
  tryCatch(
    as.character(jsonlite::toJSON(pph, auto_unbox = TRUE)),
    error = function(e) NA_character_
  )
}

# Helper: serialize ONE hypothesis's PubMed evidence to JSON (per-row column)
serialize_pubmed_evidence_one_advanced <- function(hyp_pph) {
  if (is.null(hyp_pph)) return(NA_character_)
  tryCatch(
    as.character(jsonlite::toJSON(hyp_pph, auto_unbox = TRUE)),
    error = function(e) NA_character_
  )
}


# ============================================================
# SECTION 8 – Gemini response schemas
# ============================================================

adv_gemini_schema_proposer <- function() {
  list(
    type = "OBJECT",
    properties = list(
      n_hypotheses = list(type = "INTEGER"),
      hypotheses = list(
        type = "ARRAY",
        items = list(
          type = "OBJECT",
          properties = list(
            hypothesis_id       = list(type = "STRING"),
            type                = list(
              type = "STRING",
              enum = c("Host-Endogenous", "Microbiome-Derived", "Diet-Derived",
                       "Xenobiotic / Exogenous", "Host–Microbiome")
            ),
            title               = list(type = "STRING"),
            summary             = list(type = "STRING"),
            anchor_metabolites  = list(type = "ARRAY", items = list(type = "STRING")),
            tension_metabolites = list(type = "ARRAY", items = list(type = "STRING"))
          ),
          required = c("hypothesis_id", "type", "title", "summary",
                       "anchor_metabolites", "tension_metabolites")
        )
      )
    ),
    required = c("n_hypotheses", "hypotheses")
  )
}

adv_gemini_schema_skeptic <- function() {
  list(
    type = "OBJECT",
    properties = list(
      challenges = list(
        type = "ARRAY",
        items = list(
          type = "OBJECT",
          properties = list(
            target_hypothesis_id = list(type = "STRING"),
            challenge_type       = list(type = "STRING"),
            challenge            = list(type = "STRING"),
            severity             = list(type = "STRING",
                                        enum = c("major", "minor", "no_challenge"))
          ),
          required = c("target_hypothesis_id", "challenge_type",
                       "challenge", "severity")
        )
      )
    ),
    required = "challenges"
  )
}

adv_gemini_schema_revision <- function() {
  list(
    type = "OBJECT",
    properties = list(
      revisions = list(
        type = "ARRAY",
        items = list(
          type = "OBJECT",
          properties = list(
            hypothesis_id               = list(type = "STRING"),
            action                      = list(type = "STRING",
                                               enum = c("keep", "modify", "drop")),
            justification               = list(type = "STRING"),
            revised_type                = list(
              type = "STRING",
              enum = c("Host-Endogenous", "Microbiome-Derived", "Diet-Derived",
                       "Xenobiotic / Exogenous", "Host–Microbiome")
            ),
            revised_title               = list(type = "STRING"),
            revised_summary             = list(type = "STRING"),
            revised_anchor_metabolites  = list(type = "ARRAY",
                                               items = list(type = "STRING")),
            revised_tension_metabolites = list(type = "ARRAY",
                                               items = list(type = "STRING"))
          ),
          required = c("hypothesis_id", "action", "justification")
        )
      )
    ),
    required = "revisions"
  )
}

adv_gemini_schema_arbiter <- function() {
  list(
    type = "OBJECT",
    properties = list(
      verdicts = list(
        type = "ARRAY",
        items = list(
          type = "OBJECT",
          properties = list(
            hypothesis_id        = list(type = "STRING"),
            status               = list(type = "STRING",
                                        enum = c("accepted", "rejected")),
            verdict_rationale    = list(type = "STRING"),
            confidence_label     = list(type = "STRING",
                                        enum = c("high", "medium", "low")),
            confidence_rationale = list(type = "STRING")
          ),
          required = c("hypothesis_id", "status", "verdict_rationale")
        )
      ),
      residual_ambiguity    = list(type = "STRING")
    ),
    required = c("verdicts")
  )
}

adv_gemini_schema_pubmed_relationship <- function() {
  list(
    type = "OBJECT",
    properties = list(
      phenotype_relationship = list(type = "STRING"),
      literature_support = list(type = "STRING",
                                enum = c("direct", "indirect", "insufficient")),
      evidence_note = list(type = "STRING")
    ),
    required = c("phenotype_relationship")
  )
}


# ============================================================
# SECTION 9 – Output assembly
# ============================================================

assemble_advanced_output_row <- function(module_obj,
                                          arbiter_result,
                                          hypotheses_v0,
                                          hypotheses_final,
                                          debate_digest,
                                          pubmed_per_hyp,
                                          phenotype_relationships = NULL,
                                          raw_texts) {

  mid <- module_obj$module_id %||% NA_character_

  verdicts     <- arbiter_result$verdicts     %||% list()
  residual_amb <- arbiter_result$residual_ambiguity %||% NA_character_

  # Status counts
  statuses   <- if (length(verdicts) > 0) sapply(verdicts, `[[`, "status") else character(0)
  n_accepted <- sum(statuses == "accepted",  na.rm = TRUE)
  n_rejected <- sum(statuses == "rejected",  na.rm = TRUE)

  # ---- Module-level fields (same across all output rows) ----

  # Global direction (computed, never from LLM)
  global_direction <- NA_character_
  if (!is.null(module_obj$global_direction))
    global_direction <- module_obj$global_direction$direction_label %||% NA_character_

  # Core drivers (high-priority metabolites)
  core_drivers <- NA_character_
  top_df <- module_obj$topN_metabolites
  if (!is.null(top_df) && nrow(top_df) > 0 && "Compound_name" %in% colnames(top_df)) {
    if ("priority_level" %in% colnames(top_df)) {
      hi <- top_df[top_df$priority_level == "high", , drop = FALSE]
      if (nrow(hi) == 0) hi <- head(top_df, 3)
    } else {
      hi <- head(top_df, 3)
    }
    core_drivers <- paste(hi$Compound_name, collapse = "; ")
  }

  # Origin fractions
  op <- module_obj$origin_profile
  origin_endo <- op$endogenous_frac %||% NA_real_
  origin_exo  <- op$exogenous_frac  %||% NA_real_

  # Pathway summaries
  kegg_pw <- NA_character_; hmdb_pw <- NA_character_
  pc <- module_obj$pathway_coverage
  if (!is.null(pc)) {
    fmt_pw <- function(df) {
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NA_character_)
      paste(apply(df, 1, function(r)
        sprintf("%s (%d/%d, %.0f%%)",
                r[["pathway_name"]], as.integer(r[["n_hit"]]),
                as.integer(r[["n_metabolites"]]),
                as.numeric(r[["coverage"]]) * 100)), collapse = "; ")
    }
    kegg_pw <- fmt_pw(pc$kegg)
    hmdb_pw <- fmt_pw(pc$hmdb)
  }

  # Full hypothesis JSON (all verdicts, for archival)
  all_hyp_json <- tryCatch(
    as.character(jsonlite::toJSON(lapply(verdicts, function(v) {
      h_final <- Filter(function(h) h$hypothesis_id == v$hypothesis_id, hypotheses_final)
      h_obj   <- if (length(h_final) > 0) h_final[[1]] else {
        h_v0 <- Filter(function(h) h$hypothesis_id == v$hypothesis_id, hypotheses_v0)
        if (length(h_v0) > 0) h_v0[[1]] else list(hypothesis_id = v$hypothesis_id)
      }
      rel <- if (!is.null(phenotype_relationships)) {
        phenotype_relationships[[v$hypothesis_id]]
      } else NULL
      c(h_obj, list(
        status                 = v$status,
        verdict_rationale      = v$verdict_rationale,
        phenotype_relationship = rel$phenotype_relationship %||%
          v$phenotype_relationship %||% NA_character_,
        phenotype_literature_support = rel$literature_support %||% NA_character_,
        phenotype_evidence_note = rel$evidence_note %||% NA_character_,
        confidence_label       = v$confidence_label       %||% NA_character_,
        confidence_rationale   = v$confidence_rationale   %||% NA_character_
      ))
    }), auto_unbox = TRUE)),
    error = function(e) NA_character_
  )

  # Hypothesis lookup (final version takes priority over v0)
  hyp_lookup <- c(
    setNames(hypotheses_final, sapply(hypotheses_final, `[[`, "hypothesis_id")),
    setNames(hypotheses_v0,    sapply(hypotheses_v0,    `[[`, "hypothesis_id"))
  )
  # Remove duplicate keys (final wins)
  hyp_lookup <- hyp_lookup[!duplicated(names(hyp_lookup))]

  # ---- Helper: build one output row for a verdict ----
  make_row <- function(v) {
    hid      <- v$hypothesis_id %||% NA_character_
    hobj     <- hyp_lookup[[hid]]
    hyp_pub  <- if (!is.null(pubmed_per_hyp)) pubmed_per_hyp[[hid]] else NULL
    hyp_rel  <- if (!is.null(phenotype_relationships)) phenotype_relationships[[hid]] else NULL
    data.frame(
      module_id                  = mid,
      hypothesis_id              = hid,
      hypothesis_type            = hobj$type %||% NA_character_,
      verdict_status             = v$status %||% NA_character_,
      title                      = hobj$title   %||% NA_character_,
      functional_summary         = hobj$summary %||% NA_character_,
      phenotype_relationship     = hyp_rel$phenotype_relationship %||%
        v$phenotype_relationship %||% NA_character_,
      phenotype_literature_support = hyp_rel$literature_support %||% NA_character_,
      phenotype_evidence_note    = hyp_rel$evidence_note %||% NA_character_,
      confidence_label           = v$confidence_label       %||% NA_character_,
      rationale                  = v$confidence_rationale   %||% NA_character_,
      global_direction           = global_direction,
      core_drivers               = core_drivers,
      origin_endogenous_frac     = origin_endo,
      origin_exogenous_frac      = origin_exo,
      top_kegg_pathways          = kegg_pw,
      top_hmdb_pathways          = hmdb_pw,
      # per-hypothesis PubMed (this row only)
      pubmed_pmids               = collapse_pubmed_pmids_advanced(list(hyp_pub)),
      pubmed_evidence_json       = serialize_pubmed_evidence_one_advanced(hyp_pub),
      # full archive: all hypotheses' PubMed evidence
      all_pubmed_json            = serialize_pubmed_evidence_advanced(pubmed_per_hyp),
      raw_response_call1         = raw_texts$p0      %||% NA_character_,
      raw_response_call2         = raw_texts$arbiter %||% NA_character_,
      raw_response_p0            = raw_texts$p0      %||% NA_character_,
      raw_response_s1            = raw_texts$s1      %||% NA_character_,
      raw_response_p1            = raw_texts$p1      %||% NA_character_,
      raw_response_s2            = raw_texts$s2      %||% NA_character_,
      raw_response_p2            = raw_texts$p2      %||% NA_character_,
      raw_response_pubmed_relationship = raw_texts$pubmed_relationship %||% NA_character_,
      raw_response_arbiter       = raw_texts$arbiter %||% NA_character_,
      debate_digest              = debate_digest    %||% NA_character_,
      annotation_mode            = "advanced",
      n_proposed                 = length(hypotheses_v0),
      n_survived_debate          = length(hypotheses_final),
      n_accepted                 = n_accepted,
      n_rejected                 = n_rejected,
      residual_ambiguity         = residual_amb,
      all_hypotheses_json        = all_hyp_json,
      stringsAsFactors           = FALSE
    )
  }

  # ---- One row per accepted hypothesis ----
  accepted <- Filter(function(v) identical(v$status, "accepted"), verdicts)

  if (length(accepted) > 0) {
    do.call(rbind, lapply(accepted, make_row))
  } else {
    # All hypotheses rejected → one summary row (no per-hypothesis PubMed)
    module_rel <- list(
      phenotype_relationship = NA_character_,
      literature_support     = NA_character_,
      evidence_note          = NA_character_
    )
    if (!is.null(phenotype_relationships) && length(phenotype_relationships) > 0) {
      rel_support <- vapply(
        phenotype_relationships,
        function(x) x$literature_support %||% NA_character_,
        character(1)
      )
      if (all(rel_support == "insufficient", na.rm = TRUE)) {
        rel_notes <- unique(vapply(
          phenotype_relationships,
          function(x) x$evidence_note %||% NA_character_,
          character(1)
        ))
        rel_notes <- rel_notes[!is.na(rel_notes) & nzchar(trimws(rel_notes))]
        module_rel <- adv_insufficient_phenotype_relationship(
          reason = if (length(rel_notes) > 0) paste(rel_notes, collapse = "; ") else NULL
        )
      }
    }
    data.frame(
      module_id                  = mid,
      hypothesis_id              = NA_character_,
      hypothesis_type            = NA_character_,
      verdict_status             = "rejected",
      title                      = NA_character_,
      functional_summary         = NA_character_,
      phenotype_relationship     = module_rel$phenotype_relationship %||% NA_character_,
      phenotype_literature_support = module_rel$literature_support %||% NA_character_,
      phenotype_evidence_note    = module_rel$evidence_note %||% NA_character_,
      confidence_label           = NA_character_,
      rationale                  = NA_character_,
      global_direction           = global_direction,
      core_drivers               = core_drivers,
      origin_endogenous_frac     = origin_endo,
      origin_exogenous_frac      = origin_exo,
      top_kegg_pathways          = kegg_pw,
      top_hmdb_pathways          = hmdb_pw,
      pubmed_pmids               = collapse_pubmed_pmids_advanced(pubmed_per_hyp),
      pubmed_evidence_json       = NA_character_,
      all_pubmed_json            = serialize_pubmed_evidence_advanced(pubmed_per_hyp),
      raw_response_call1         = raw_texts$p0      %||% NA_character_,
      raw_response_call2         = raw_texts$arbiter %||% NA_character_,
      raw_response_p0            = raw_texts$p0      %||% NA_character_,
      raw_response_s1            = raw_texts$s1      %||% NA_character_,
      raw_response_p1            = raw_texts$p1      %||% NA_character_,
      raw_response_s2            = raw_texts$s2      %||% NA_character_,
      raw_response_p2            = raw_texts$p2      %||% NA_character_,
      raw_response_pubmed_relationship = raw_texts$pubmed_relationship %||% NA_character_,
      raw_response_arbiter       = raw_texts$arbiter %||% NA_character_,
      debate_digest              = debate_digest    %||% NA_character_,
      annotation_mode            = "advanced",
      n_proposed                 = length(hypotheses_v0),
      n_survived_debate          = length(hypotheses_final),
      n_accepted                 = 0L,
      n_rejected                 = n_rejected,
      residual_ambiguity         = residual_amb,
      all_hypotheses_json        = all_hyp_json,
      stringsAsFactors           = FALSE
    )
  }
}


# ============================================================
# SECTION 10 – Main orchestrator: one module
# ============================================================

#' Annotate a Single Module Using the Advanced Multi-Agent Debate Pipeline
#'
#' @param module_obj  List from \code{prepare_module()}.
#' @param study_context List with \code{species}, \code{biospecimen}, \code{phenotype}.
#' @param provider,llm_api_key,model  LLM configuration.
#' @param max_tokens  Max tokens per LLM response (default 2000; Gemini auto-expanded).
#' @param proposer_temperature  Sampling temperature for Proposer P0 (default 0.3).
#' @param temperature  Temperature for Skeptic / Proposer revision / Arbiter (default 0).
#' @param pubmed  Logical. Whether to search PubMed before Arbiter.
#' @param NCBI_api_key  NCBI API key.
#' @param years  PubMed look-back years.
#' @param debate_html Logical. Whether to export a per-module HTML debate trace.
#' @param debate_html_dir Directory for debate HTML files.
#' @return Data frame compatible with fast mode output plus advanced columns.
#'   Accepted hypotheses are returned as separate rows; rejected
#'   modules return one rejected summary row. PubMed evidence is exported in
#'   \code{pubmed_pmids} and \code{pubmed_evidence_json}. Debate trace is
#'   exported in \code{raw_response_p0}, \code{raw_response_s1},
#'   \code{raw_response_p1}, \code{raw_response_s2},
#'   \code{raw_response_p2}, \code{raw_response_pubmed_relationship},
#'   \code{raw_response_arbiter}, and \code{debate_digest}.
annotate_one_module_advanced <- function(module_obj,
                                          study_context,
                                          provider,
                                          llm_api_key,
                                          model,
                                          max_tokens           = 2000L,
                                          proposer_temperature = 0.3,
                                          temperature          = 0,
                                          pubmed               = FALSE,
                                          NCBI_api_key         = NULL,
                                          years                = 5,
                                          debate_html          = FALSE,
                                          debate_html_dir      = getwd()) {

  mid <- module_obj$module_id %||% "?"
  is_gemini <- identical(provider, "gemini")

  # Gemini schemas (NULL for non-Gemini)
  sch_proposer  <- if (is_gemini) adv_gemini_schema_proposer()  else NULL
  sch_skeptic   <- if (is_gemini) adv_gemini_schema_skeptic()   else NULL
  sch_revision  <- if (is_gemini) adv_gemini_schema_revision()  else NULL
  sch_arbiter   <- if (is_gemini) adv_gemini_schema_arbiter()   else NULL
  sch_pubmed_rel <- if (is_gemini) adv_gemini_schema_pubmed_relationship() else NULL

  raw_texts <- list()

  # ----------------------------------------------------------------
  # Build reusable inputs (computed once, reused across all agents)
  # ----------------------------------------------------------------
  skeptic_evidence <- build_skeptic_module_evidence(module_obj)
  context_block    <- build_context_block(module_obj, study_context)
  valid_top_metabolites <- adv_top_metabolite_names(module_obj)

  # ----------------------------------------------------------------
  # P0 – Proposer: initial hypothesis generation
  # ----------------------------------------------------------------
  message("[advanced] Module ", mid, " | P0 Proposer...")

  p0_sys  <- build_proposer_system_prompt()
  p0_usr  <- build_module_prompt(module_obj)   # reuse fast-mode prompt builder

  p0_res  <- adv_request_json(
    messages    = list(list(role = "system", content = p0_sys),
                       list(role = "user",   content = p0_usr)),
    required_keys = "hypotheses",
    schema      = sch_proposer,
    retry_instr = "Return only valid JSON with key \"hypotheses\" (array of hypothesis objects).",
    provider    = provider, llm_api_key = llm_api_key, model = model,
    max_tokens  = max_tokens, temperature = proposer_temperature
  )
  raw_texts$p0 <- p0_res$raw_text

  hypotheses_v0 <- parse_proposer_output(
    p0_res$parsed,
    valid_metabolites = valid_top_metabolites
  )

  if (is.null(hypotheses_v0) || length(hypotheses_v0) == 0) {
    message("[advanced] Module ", mid, ": P0 failed, falling back to fast mode.")
    fast_row <- annotate_one_module(
      module_obj    = module_obj,
      study_context = study_context,
      provider      = provider,
      llm_api_key   = llm_api_key,
      model         = model,
      max_tokens    = max_tokens,
      temperature   = temperature,
      pubmed        = pubmed,
      NCBI_api_key  = NCBI_api_key,
      years         = years
    )
    fast_row$annotation_mode       <- "fast_fallback"
    fast_row$hypothesis_id         <- NA_character_
    fast_row$hypothesis_type       <- NA_character_
    fast_row$verdict_status        <- NA_character_
    fast_row$n_proposed            <- 0L
    fast_row$n_survived_debate     <- NA_integer_
    fast_row$n_accepted            <- NA_integer_
    fast_row$n_rejected            <- NA_integer_
    fast_row$residual_ambiguity    <- NA_character_
    fast_row$all_hypotheses_json   <- NA_character_
    return(fast_row)
  }

  # ----------------------------------------------------------------
  # S1 – Skeptic round 1
  # ----------------------------------------------------------------
  message("[advanced] Module ", mid, " | S1 Skeptic (round 1)...")

  s1_sys <- build_skeptic_system_prompt(round = 1L)
  s1_usr <- build_skeptic_user_prompt(skeptic_evidence, context_block,
                                       hypotheses_v0, history_text = NULL)

  s1_res <- adv_request_json(
    messages    = list(list(role = "system", content = s1_sys),
                       list(role = "user",   content = s1_usr)),
    required_keys = "challenges",
    schema      = sch_skeptic,
    retry_instr = "Return only valid JSON with key \"challenges\" (array of challenge objects).",
    provider    = provider, llm_api_key = llm_api_key, model = model,
    max_tokens  = max_tokens, temperature = temperature
  )
  raw_texts$s1 <- s1_res$raw_text
  hids_v0 <- vapply(hypotheses_v0, `[[`, character(1), "hypothesis_id")
  s1_challenges <- parse_skeptic_output(s1_res$parsed, valid_hypothesis_ids = hids_v0)

  # ----------------------------------------------------------------
  # P1 – Proposer revision 1
  # ----------------------------------------------------------------
  message("[advanced] Module ", mid, " | P1 Proposer revision (round 1)...")

  p1_sys <- build_proposer_revision_system_prompt()
  p1_usr <- build_proposer_revision_user_prompt(
    hypotheses_v0, s1_challenges,
    module_evidence = skeptic_evidence,
    context_block   = context_block,
    history_text    = NULL
  )

  p1_res <- adv_request_json(
    messages    = list(list(role = "system", content = p1_sys),
                       list(role = "user",   content = p1_usr)),
    required_keys = "revisions",
    schema      = sch_revision,
    retry_instr = "Return only valid JSON with key \"revisions\" (array of revision objects).",
    provider    = provider, llm_api_key = llm_api_key, model = model,
    max_tokens  = max_tokens, temperature = temperature
  )
  raw_texts$p1 <- p1_res$raw_text
  p1_revisions      <- parse_revision_output(
    p1_res$parsed,
    valid_hypothesis_ids = hids_v0,
    valid_metabolites = valid_top_metabolites
  )
  hypotheses_after_r1 <- apply_revisions(hypotheses_v0, p1_revisions)

  # ----------------------------------------------------------------
  # Round 2 decision
  # ----------------------------------------------------------------
  do_r2 <- should_run_round2(hypotheses_after_r1, s1_challenges, p1_revisions)

  hypotheses_final <- hypotheses_after_r1
  s2_challenges    <- NULL
  p2_revisions     <- NULL

  if (do_r2) {
    hids_r1 <- vapply(hypotheses_after_r1, `[[`, character(1), "hypothesis_id")

    # Build comprehensive Round 1 histories for S2 and P2.
    # S2 sees: P0 originals + S1 full challenges + P1 decisions (including revised summaries).
    # P2 sees: same as S2 plus all P0 anchor/tension details to prevent reverting.
    r1_history_for_s2 <- build_round1_history_for_s2(
      hypotheses_v0  = hypotheses_v0,
      s1_challenges  = s1_challenges,
      p1_revisions   = p1_revisions,
      hypotheses_v1  = hypotheses_after_r1
    )
    r1_history_for_p2 <- build_round1_history_for_p2(
      hypotheses_v0  = hypotheses_v0,
      s1_challenges  = s1_challenges,
      p1_revisions   = p1_revisions,
      hypotheses_v1  = hypotheses_after_r1
    )

    # S2 – Skeptic round 2
    message("[advanced] Module ", mid, " | S2 Skeptic (round 2)...")

    s2_sys <- build_skeptic_system_prompt(round = 2L)
    s2_usr <- build_skeptic_user_prompt(
      skeptic_evidence, context_block,
      hypotheses_after_r1, history_text = r1_history_for_s2
    )

    s2_res <- adv_request_json(
      messages    = list(list(role = "system", content = s2_sys),
                         list(role = "user",   content = s2_usr)),
      required_keys = "challenges",
      schema      = sch_skeptic,
      retry_instr = "Return only valid JSON with key \"challenges\".",
      provider    = provider, llm_api_key = llm_api_key, model = model,
      max_tokens  = max_tokens, temperature = temperature
    )
    raw_texts$s2 <- s2_res$raw_text
    s2_challenges <- parse_skeptic_output(s2_res$parsed, valid_hypothesis_ids = hids_r1)

    # P2 – Proposer revision 2
    message("[advanced] Module ", mid, " | P2 Proposer revision (round 2)...")

    p2_sys <- build_proposer_revision_system_prompt()
    p2_usr <- build_proposer_revision_user_prompt(
      hypotheses_after_r1, s2_challenges,
      module_evidence = skeptic_evidence,
      context_block   = context_block,
      history_text    = r1_history_for_p2
    )

    p2_res <- adv_request_json(
      messages    = list(list(role = "system", content = p2_sys),
                         list(role = "user",   content = p2_usr)),
      required_keys = "revisions",
      schema      = sch_revision,
      retry_instr = "Return only valid JSON with key \"revisions\".",
      provider    = provider, llm_api_key = llm_api_key, model = model,
      max_tokens  = max_tokens, temperature = temperature
    )
    raw_texts$p2 <- p2_res$raw_text
    p2_revisions     <- parse_revision_output(
      p2_res$parsed,
      valid_hypothesis_ids = hids_r1,
      valid_metabolites = valid_top_metabolites
    )
    hypotheses_final <- apply_revisions(hypotheses_after_r1, p2_revisions)
  }

  # ----------------------------------------------------------------
  # PubMed (optional, always before Arbiter regardless of debate outcome)
  # When all hypotheses were dropped, search the original v0 hypotheses so
  # the Arbiter has evidence to support or reject any potential recovery.
  # ----------------------------------------------------------------
  pubmed_per_hyp <- NULL
  phenotype_relationships <- NULL
  hyps_for_phenotype <- if (length(hypotheses_final) > 0) hypotheses_final else hypotheses_v0
  has_phenotype <- adv_has_context_value(study_context$phenotype)

  if (!has_phenotype) {
    phenotype_relationships <- adv_build_insufficient_relationships(
      hyps_for_phenotype,
      reason = "Study phenotype was not provided, so phenotype-specific literature evidence was not searched."
    )
  } else if (!isTRUE(pubmed)) {
    phenotype_relationships <- adv_build_inferable_relationships(
      hyps_for_phenotype,
      reason = paste(
        "PubMed phenotype evidence was not requested (pubmed = FALSE).",
        "General-knowledge phenotype inference is allowed but must be labeled as inferred and not sufficiently literature-supported."
      )
    )
  } else {
    hyps_for_pubmed <- hyps_for_phenotype
    message("[advanced] Module ", mid, " | PubMed search for ",
            length(hyps_for_pubmed), " hypothesis(es) (pre-Arbiter)...")

    pubmed_per_hyp <- search_pubmed_advanced(
      hypotheses_final     = hyps_for_pubmed,
      phenotype            = study_context$phenotype,
      NCBI_api_key         = NCBI_api_key,
      retmax_per_query     = 3L,
      years                = years
    )

    message("[advanced] Module ", mid, " | PubMed phenotype relationship summary...")
    pubmed_rel <- summarize_pubmed_relationships_advanced(
      hypotheses      = hyps_for_pubmed,
      study_context   = study_context,
      pubmed_per_hyp  = pubmed_per_hyp,
      provider        = provider,
      llm_api_key     = llm_api_key,
      model           = model,
      max_tokens      = 1200L,
      temperature     = temperature,
      schema          = sch_pubmed_rel
    )
    phenotype_relationships <- pubmed_rel$summaries
    raw_texts$pubmed_relationship <- pubmed_rel$raw_text
  }

  # ----------------------------------------------------------------
  # Build debate history digest for Arbiter
  # ----------------------------------------------------------------
  arbiter_scope_ids <- if (length(hypotheses_final) > 0) {
    vapply(hypotheses_final, `[[`, character(1), "hypothesis_id")
  } else {
    vapply(hypotheses_v0, `[[`, character(1), "hypothesis_id")
  }

  debate_digest <- build_debate_history_digest(
    hypotheses_v0    = hypotheses_v0,
    s1_challenges    = s1_challenges,
    p1_revisions     = p1_revisions,
    hypotheses_v1    = hypotheses_after_r1,
    s2_challenges    = s2_challenges,
    p2_revisions     = p2_revisions,
    hypotheses_final = hypotheses_final,
    include_ids      = arbiter_scope_ids
  )

  # ----------------------------------------------------------------
  # A0 – Arbiter
  # ----------------------------------------------------------------
  message("[advanced] Module ", mid, " | A0 Arbiter (final verdict)...")

  a0_sys <- build_arbiter_system_prompt()
  a0_usr <- build_arbiter_user_prompt(
    hypotheses_final        = hypotheses_final,
    hypotheses_v0           = hypotheses_v0,
    debate_digest           = debate_digest,
    module_obj              = module_obj,
    study_context           = study_context,
    module_evidence         = skeptic_evidence,
    phenotype_relationships = phenotype_relationships
  )

  a0_res <- adv_request_json(
    messages    = list(list(role = "system", content = a0_sys),
                       list(role = "user",   content = a0_usr)),
    required_keys = "verdicts",
    schema      = sch_arbiter,
    retry_instr = "Return only valid JSON with key \"verdicts\" (array) and \"residual_ambiguity\" (string or null).",
    provider    = provider, llm_api_key = llm_api_key, model = model,
    max_tokens  = max_tokens, temperature = temperature
  )
  raw_texts$arbiter <- a0_res$raw_text

  # Parse Arbiter verdicts
  arbiter_result <- list(
    verdicts           = list(),
    residual_ambiguity = NA_character_
  )
  if (!is.null(a0_res$parsed)) {
    parsed_verdicts <- a0_res$parsed$verdicts
    if (!is.null(parsed_verdicts) && length(parsed_verdicts) > 0) {
      arbiter_result$verdicts <- lapply(parsed_verdicts, function(v) {
        list(
          hypothesis_id          = adv_scalar(v$hypothesis_id, ""),
          status                 = tolower(adv_scalar(v$status, "rejected")),
          verdict_rationale      = adv_scalar(v$verdict_rationale, ""),
          phenotype_relationship = NA_character_,
          confidence_label       = tolower(adv_scalar(v$confidence_label, "")),
          confidence_rationale   = adv_scalar(v$confidence_rationale, "")
        )
      })
    }
    arbiter_result$verdicts <- Filter(function(v) {
      v$status %in% c("accepted", "rejected") &&
        (v$status == "rejected" || v$confidence_label %in% c("high", "medium", "low"))
    }, arbiter_result$verdicts)
    arbiter_result$verdicts <- Filter(
      function(v) v$hypothesis_id %in% arbiter_scope_ids,
      arbiter_result$verdicts
    )
    rab <- adv_scalar(a0_res$parsed$residual_ambiguity, "")
    if (nzchar(rab) && !identical(toupper(rab), "NA"))
      arbiter_result$residual_ambiguity <- rab
  }

  # If Arbiter failed or returned empty verdicts, build minimal rejected output
  if (length(arbiter_result$verdicts) == 0) {
    message("[advanced] Module ", mid, ": Arbiter returned no verdicts; marking all rejected.")
    arbiter_result$verdicts <- lapply(
      c(hypotheses_final, if (length(hypotheses_final) == 0) hypotheses_v0 else list()),
      function(h) list(
        hypothesis_id     = h$hypothesis_id,
        status            = "rejected",
        verdict_rationale = "Arbiter call failed or returned no output.",
        phenotype_relationship = NA_character_,
        confidence_label  = NA_character_,
        confidence_rationale = NA_character_
      )
    )
  }

  if (isTRUE(debate_html)) {
    # --------------------------------------------------------------
    # Generate debate HTML visualization
    # --------------------------------------------------------------
    tryCatch(
      generate_debate_html(
        module_id        = mid,
        module_obj       = module_obj,
        study_context    = study_context,
        hypotheses_v0    = hypotheses_v0,
        s1_challenges    = s1_challenges,
        p1_revisions     = p1_revisions,
        hypotheses_v1    = hypotheses_after_r1,
        s2_challenges    = s2_challenges,
        p2_revisions     = p2_revisions,
        hypotheses_final = hypotheses_final,
        arbiter_result   = arbiter_result,
        phenotype_relationships = phenotype_relationships,
        pubmed_per_hyp   = pubmed_per_hyp,
        output_dir       = debate_html_dir
      ),
      error = function(e) message("[HTML] Skipped (", e$message, ")")
    )
  }

  # ----------------------------------------------------------------
  # Assemble output row
  # ----------------------------------------------------------------
  assemble_advanced_output_row(
    module_obj       = module_obj,
    arbiter_result   = arbiter_result,
    hypotheses_v0    = hypotheses_v0,
    hypotheses_final = hypotheses_final,
    debate_digest    = debate_digest,
    pubmed_per_hyp   = pubmed_per_hyp,
    phenotype_relationships = phenotype_relationships,
    raw_texts        = raw_texts
  )
}


# ============================================================
# SECTION 11 – Loop over all modules
# ============================================================

#' Annotate All Modules Using the Advanced Multi-Agent Debate Pipeline
#'
#' @param all_modules_ready Named list of module objects from \code{prepare_module()}.
#' @param study_context  List with \code{species}, \code{biospecimen}, \code{phenotype}.
#' @param provider,llm_api_key,model  LLM configuration.
#' @param max_tokens  Default 2000.
#' @param proposer_temperature  Default 0.3.
#' @param temperature  Default 0.
#' @param sleep_sec  Pause between modules (default 0).
#' @param pubmed  Logical.
#' @param NCBI_api_key  NCBI API key.
#' @param years  PubMed look-back years.
#' @param debate_html Logical. Whether to export per-module HTML debate traces.
#' @param debate_html_dir Directory for debate HTML files.
#' @return Data frame with one row per module.
annotate_all_modules_advanced <- function(all_modules_ready,
                                           study_context,
                                           provider,
                                           llm_api_key,
                                           model,
                                           max_tokens           = 2000L,
                                           proposer_temperature = 0.3,
                                           temperature          = 0,
                                           sleep_sec            = 0,
                                           pubmed               = FALSE,
                                           NCBI_api_key         = NULL,
                                           years                = 5,
                                           debate_html          = FALSE,
                                           debate_html_dir      = getwd()) {

  n_mod <- length(all_modules_ready)
  results_list <- lapply(seq_len(n_mod), function(i) {
    mod <- all_modules_ready[[i]]
    message("[advanced] Processing module ", mod$module_id %||% i,
            " (", i, "/", n_mod, ")")

    out <- tryCatch(
      annotate_one_module_advanced(
        module_obj           = mod,
        study_context        = study_context,
        provider             = provider,
        llm_api_key          = llm_api_key,
        model                = model,
        max_tokens           = max_tokens,
        proposer_temperature = proposer_temperature,
        temperature          = temperature,
        pubmed               = pubmed,
        NCBI_api_key         = NCBI_api_key,
        years                = years,
        debate_html          = debate_html,
        debate_html_dir      = debate_html_dir
      ),
      error = function(e) {
        message("[advanced] Module ", mod$module_id %||% i,
                " failed: ", e$message)
        NULL
      }
    )

    if (sleep_sec > 0) Sys.sleep(sleep_sec)
    out
  })

  dplyr::bind_rows(Filter(Negate(is.null), results_list))
}
