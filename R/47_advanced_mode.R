# ============================================================
# Advanced Mode Pipeline
#
# Multi-agent debate for metabolic module annotation.
#
# Agents (per module):
#   P0  Proposer initial      – generates 2-3 competing hypotheses
#   S1  Skeptic round 1       – challenges each hypothesis
#   P1  Proposer revision 1   – keep / modify / drop
#   S2  Skeptic round 2       – (conditional) re-challenges revised set
#   P2  Proposer revision 2   – (conditional) final keep / modify / drop
#   A0  Arbiter               – accepted / uncertain / rejected + confidence
#
# Information visible per agent:
#   P0  : full module data (metabolites, class, pathway, Signal_concentration)
#   S1/S2: compressed module summary + origin_profile + species + biospecimen
#          + current hypothesis titles/claims/anchors/tensions (NO full mechanism)
#   P1/P2: same as S1/S2 input + Skeptic challenges + full hypothesis text
#   A0  : surviving hypotheses (full) + debate history digest
#          + origin_profile + species + biospecimen + phenotype
#          + per-hypothesis PubMed evidence
#
# global_direction: computed from module data, stored in output row ONLY.
# It is NOT passed to any agent — module-level direction is too complex
# to use as a mechanistic constraint.
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

# Build a ~150-token compressed reference summary used by S1/S2/P1/P2.
# Does NOT include phenotype or global_direction.
build_compressed_module_summary <- function(module_obj) {

  scalar <- function(x, na = "NA") {
    if (is.null(x) || length(x) == 0) return(na)
    x <- trimws(as.character(x[1]))
    if (!nzchar(x) || is.na(x)) na else x
  }

  # Signal concentration
  sig <- scalar(module_obj$module_stats$Signal_concentration)

  # Top metabolites (up to 8, compact)
  top_df <- module_obj$topN_metabolites
  if (!is.null(top_df) && nrow(top_df) > 0) {
    top_lines <- apply(head(top_df, 8), 1, function(r) {
      sprintf("  %s. %s [%s]",
              r[["top_rank"]],
              r[["Compound_name"]],
              r[["priority_level"]])
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
    "--- Module compressed summary ---
Signal_concentration: {sig}

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

# Format hypotheses compactly for Skeptic / Proposer revision.
# Returns a character string (not the full mechanism, just claim + anchors + tensions).
format_hypotheses_compact <- function(hypotheses) {
  if (length(hypotheses) == 0) return("(no hypotheses)")
  blocks <- lapply(hypotheses, function(h) {
    anc <- if (length(h$anchor_metabolites) > 0)
      paste(h$anchor_metabolites, collapse = ", ") else "none listed"
    ten <- if (length(h$tension_metabolites) > 0)
      paste(h$tension_metabolites, collapse = ", ") else "none listed"
    glue::glue(
      '{h$hypothesis_id} — "{h$title}"
  Core claim: {h$claim}
  Anchor metabolites: {anc}
  Tension metabolites: {ten}'
    )
  })
  paste(blocks, collapse = "\n\n")
}

# Format hypotheses with full mechanism (used only for Proposer revision prompt
# where the Proposer may need to rewrite the mechanism).
format_hypotheses_full <- function(hypotheses) {
  if (length(hypotheses) == 0) return("(no hypotheses)")
  blocks <- lapply(hypotheses, function(h) {
    anc <- if (length(h$anchor_metabolites) > 0)
      paste(h$anchor_metabolites, collapse = ", ") else "none listed"
    ten <- if (length(h$tension_metabolites) > 0)
      paste(h$tension_metabolites, collapse = ", ") else "none listed"
    glue::glue(
      '{h$hypothesis_id} — "{h$title}"
  Mechanism: {h$summary}
  Core claim: {h$claim}
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
      if (!is.null(rev$revised_title)   && nzchar(rev$revised_title))
        h$title   <- rev$revised_title
      if (!is.null(rev$revised_summary) && nzchar(rev$revised_summary))
        h$summary <- rev$revised_summary
      if (!is.null(rev$revised_claim)   && nzchar(rev$revised_claim))
        h$claim   <- rev$revised_claim
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
parse_proposer_output <- function(parsed) {
  if (is.null(parsed)) return(NULL)
  hyps_raw <- parsed$hypotheses
  if (is.null(hyps_raw) || length(hyps_raw) == 0) return(NULL)

  ensure_char_vec <- function(x) {
    if (is.null(x)) return(character(0))
    if (is.list(x)) return(as.character(unlist(x)))
    as.character(x)
  }

  lapply(hyps_raw, function(h) {
    list(
      hypothesis_id        = trimws(as.character(h$hypothesis_id %||% "H?")),
      title                = trimws(as.character(h$title         %||% "")),
      summary              = trimws(as.character(h$summary       %||% "")),
      claim                = trimws(as.character(h$claim         %||% "")),
      anchor_metabolites   = ensure_char_vec(h$anchor_metabolites),
      tension_metabolites  = ensure_char_vec(h$tension_metabolites)
    )
  })
}

# Parse Skeptic output into a clean list of challenge objects.
parse_skeptic_output <- function(parsed) {
  if (is.null(parsed)) return(list())
  chs_raw <- parsed$challenges
  if (is.null(chs_raw) || length(chs_raw) == 0) return(list())

  lapply(chs_raw, function(c) {
    list(
      target_hypothesis_id = trimws(as.character(c$target_hypothesis_id %||% "")),
      challenge_type       = trimws(as.character(c$challenge_type       %||% "")),
      challenge            = trimws(as.character(c$challenge            %||% "")),
      severity             = trimws(tolower(as.character(c$severity     %||% "minor")))
    )
  })
}

# Parse Proposer revision output into a clean list of revision objects.
parse_revision_output <- function(parsed) {
  if (is.null(parsed)) return(list())
  revs_raw <- parsed$revisions
  if (is.null(revs_raw) || length(revs_raw) == 0) return(list())

  ensure_char_vec <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.list(x)) return(as.character(unlist(x)))
    as.character(x)
  }

  lapply(revs_raw, function(r) {
    list(
      hypothesis_id               = trimws(as.character(r$hypothesis_id %||% "")),
      action                      = trimws(tolower(as.character(r$action %||% "keep"))),
      justification               = trimws(as.character(r$justification %||% "")),
      revised_title               = trimws(as.character(r$revised_title %||% "")),
      revised_summary             = trimws(as.character(r$revised_summary %||% "")),
      revised_claim               = trimws(as.character(r$revised_claim %||% "")),
      revised_anchor_metabolites  = ensure_char_vec(r$revised_anchor_metabolites),
      revised_tension_metabolites = ensure_char_vec(r$revised_tension_metabolites)
    )
  })
}


# ============================================================
# SECTION 3b – Challenge routing helpers
# ============================================================

# Check whether a challenge object targets a given hypothesis_id.
# Handles both exact targets ("H1") and compound targets ("H1_H2") produced
# by the "redundant" challenge type, where Skeptic flags two hypotheses at once.
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
  # Uses challenge_targets_hid() to handle compound "H1_H2" redundant targets.
  max_sev <- function(hid) {
    relevant <- Filter(function(c) challenge_targets_hid(c, hid), s1_challenges)
    if (length(relevant) == 0) return("none")
    if (any(sapply(relevant, function(c) identical(c$severity, "major")))) "major" else "minor"
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
# SECTION 5 – Debate history formatter
# ============================================================

build_debate_history_digest <- function(hypotheses_v0,
                                         s1_challenges,
                                         p1_revisions,
                                         hypotheses_v1,
                                         s2_challenges = NULL,
                                         p2_revisions  = NULL,
                                         hypotheses_final = NULL) {

  rev_map_r1 <- setNames(p1_revisions, sapply(p1_revisions, `[[`, "hypothesis_id"))
  rev_map_r2 <- if (!is.null(p2_revisions))
    setNames(p2_revisions, sapply(p2_revisions, `[[`, "hypothesis_id")) else list()

  surviving_ids <- if (!is.null(hypotheses_final))
    sapply(hypotheses_final, `[[`, "hypothesis_id") else character(0)

  blocks <- lapply(hypotheses_v0, function(h) {
    hid <- h$hypothesis_id
    lines <- c(sprintf('%s — "%s" (original claim: %s)', hid, h$title, h$claim))

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
      if (r1_rev$action == "modify" && nzchar(r1_rev$revised_title))
        lines <- c(lines, sprintf("    → Revised title: \"%s\"", r1_rev$revised_title))
      if (r1_rev$action == "modify" && nzchar(r1_rev$revised_claim))
        lines <- c(lines, sprintf("    → Revised claim: %s", r1_rev$revised_claim))
    }

    if (r1_rev$action == "drop") {
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
        if (r2_rev$action == "modify" && nzchar(r2_rev$revised_title))
          lines <- c(lines, sprintf("    → Revised title: \"%s\"", r2_rev$revised_title))
        if (r2_rev$action == "modify" && nzchar(r2_rev$revised_claim))
          lines <- c(lines, sprintf("    → Revised claim: %s", r2_rev$revised_claim))
        if (r2_rev$action == "drop") {
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
Your task is to generate 2-3 competing biological hypotheses for a dysregulated metabolic module.

Each hypothesis must represent a GENUINELY DISTINCT mechanistic explanation —
different biochemical pathways or processes, not just different names for the same mechanism.
Two hypotheses invoking the same central metabolite via the same enzymatic route are NOT distinct.

You will receive:
- Signal_concentration (high | medium | low): whether importance is concentrated in a few metabolites
- Top metabolites: rank, name, chemical class, priority_level
- Class-level structure: class name, proportion of module, metabolite names, priority_level
- Pathway evidence from KEGG and HMDB

Number-of-hypotheses guidance:
- Generate 2 when Signal_concentration is high OR metabolites clearly cluster around one or two related pathways
- Generate 3 when Signal_concentration is medium/low OR metabolites span distinctly different chemical classes

For each hypothesis provide:
1. hypothesis_id : \"H1\", \"H2\", or \"H3\"
2. title         : KEGG-style concise pathway name
3. summary       : 3-4 sentences describing the biochemical mechanism (factual, database-style)
4. claim         : ONE sentence capturing the core mechanistic claim
5. anchor_metabolites  : metabolite names that STRONGLY support this hypothesis
6. tension_metabolites : metabolite names that are DIFFICULT to explain under this hypothesis
                         (be honest — tension_metabolites are mandatory for every hypothesis)

Rules:
- Every hypothesis must be anchored by at least one high-priority metabolite
- anchor_metabolites must not overlap more than 50% between any two hypotheses
- Do NOT reference study phenotype, disease, or origin (endogenous/exogenous)
- Do NOT invent metabolites not present in the module data

Return STRICT JSON only:
{{
  \"n_hypotheses\": 2,
  \"hypotheses\": [
    {{
      \"hypothesis_id\": \"H1\",
      \"title\": \"...\",
      \"summary\": \"...\",
      \"claim\": \"...\",
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
- Focus on genuinely new issues, not previously settled ones"
  ) else ""

  glue::glue(
    "You are an expert in metabolomics and biochemical pathway interpretation.
Your task is to critically evaluate competing biological hypotheses for a metabolic module.

You will receive:
- Compressed module summary (top metabolites, class structure, pathway coverage, Signal_concentration)
- Biological context: origin profile, species, biospecimen
- 2-3 hypotheses, each shown as: title, core claim, anchor metabolites, tension metabolites

For each hypothesis generate 1-3 specific challenges grounded ONLY in the provided information.

Challenge types (use exactly one per challenge):
- weak_evidence         : anchor metabolites do not strongly support the claimed pathway
- metabolite_mismatch   : high-priority metabolites are better explained by a different mechanism
- pathway_gap           : KEGG/HMDB coverage is absent or very low for the claimed pathway
- origin_inconsistency  : the mechanism (endo/exo) conflicts with the origin profile proportions
- species_biospecimen   : the claimed pathway is implausible for the stated species or biospecimen
- redundant             : two hypotheses are mechanistically equivalent (use target_hypothesis_id \"H1_H2\")
- tension_unresolved    : the hypothesis lists metabolites as tensions but its claim ignores them entirely
- inadequate_revision   : (round 2 only) a Round 1 revision failed to address the original challenge{round2_addon}

Severity rules:
- major : the challenge undermines the core claim; the hypothesis cannot stand without addressing this
- minor : a notable gap or weakness that should be acknowledged, but does not invalidate the hypothesis

FORBIDDEN:
- Do not challenge based on phenotype relevance or disease association
- Do not cite literature, papers, or external studies
- Every challenge must reference specific metabolite names, pathway names, or numeric values from the prompt

Return STRICT JSON only:
{{
  \"challenges\": [
    {{
      \"target_hypothesis_id\": \"H1\",
      \"challenge_type\": \"...\",
      \"challenge\": \"2-3 sentences with specific evidence...\",
      \"severity\": \"major | minor\"
    }}
  ]
}}"
  )
}

build_skeptic_user_prompt <- function(compressed_summary,
                                       context_block,
                                       hypotheses,
                                       history_text = NULL) {
  hyp_text <- format_hypotheses_compact(hypotheses)

  history_section <- if (!is.null(history_text) && nzchar(history_text))
    glue::glue("\n\n--- Round 1 history (for reference) ---\n{history_text}") else ""

  glue::glue(
    "{compressed_summary}

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
- Compressed module summary and biological context (origin, species, biospecimen)
- Current hypotheses with their full mechanisms, claims, and anchor/tension metabolites
- Skeptic challenges (with severity: major or minor)

For each hypothesis, decide: keep | modify | drop

Rules:
- keep   : hypothesis is valid; provide a counter-argument for every major challenge,
           citing specific evidence from the module data
- modify : revise to address the challenges; provide updated title, summary, claim,
           anchor_metabolites, tension_metabolites; a modification must change either
           the title OR the mechanism text — renaming alone is not valid
- drop   : hypothesis is fundamentally flawed; provide a clear reason

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
      \"revised_title\": \"...\",
      \"revised_summary\": \"...\",
      \"revised_claim\": \"...\",
      \"revised_anchor_metabolites\": [\"...\"],
      \"revised_tension_metabolites\": [\"...\"]
    }}
  ]
}}"
  )
}

build_proposer_revision_user_prompt <- function(compressed_summary,
                                                 context_block,
                                                 hypotheses,
                                                 challenges,
                                                 history_text = NULL) {
  hyp_text <- format_hypotheses_full(hypotheses)

  # Format challenges grouped by hypothesis.
  # challenge_targets_hid() ensures compound targets ("H1_H2") are distributed
  # to both H1 and H2 entries.  format_challenge_line() adds a cross-hypothesis
  # annotation so the Proposer can see which other hypothesis shares the issue.
  ch_blocks <- lapply(hypotheses, function(h) {
    hid  <- h$hypothesis_id
    chs  <- Filter(function(c) challenge_targets_hid(c, hid), challenges)
    if (length(chs) == 0) return(sprintf("  %s: no challenges", hid))
    ch_lines <- sapply(chs, function(c) format_challenge_line(c, hid))
    paste(c(sprintf("  Challenges for %s:", hid), ch_lines), collapse = "\n")
  })
  ch_text <- paste(ch_blocks, collapse = "\n\n")

  history_section <- if (!is.null(history_text) && nzchar(history_text))
    glue::glue("\n\n--- Prior round history (your previous commitments) ---\n{history_text}") else ""

  glue::glue(
    "{compressed_summary}

{context_block}

--- Current hypotheses (full details) ---
{hyp_text}

--- Skeptic challenges ---
{ch_text}{history_section}"
  )
}


# ---- A0: Arbiter ----

build_arbiter_system_prompt <- function() {
  glue::glue(
    "You are an expert in metabolomics and biochemical pathway interpretation.
Your task is to make final verdicts on competing biological hypotheses for a metabolic module,
and to assess phenotypic relevance and confidence for surviving hypotheses.

You will receive:
- Final surviving hypotheses (or all original hypotheses if all were dropped during debate)
- A debate history digest showing how each hypothesis evolved through challenges and revisions
- Origin profile, global direction, species, biospecimen, and phenotype (full study context)
- PubMed evidence: shared (from top metabolite searches) and per-hypothesis (from title searches)

Verdict options:
- accepted  : biochemically sound, survived debate, plausible phenotype relevance
- uncertain : has merit but carries unresolved issues, weak literature support,
              or ambiguous evidence (including hypotheses recovered after all-drop)
- rejected  : biochemically implausible, contradicted by data, or dropped during debate
              for reasons you agree with

For accepted and uncertain hypotheses, also provide:
- phenotype_relationship : 2-4 sentences based ONLY on provided PubMed evidence.
  If no evidence, write exactly: \"Insufficient literature evidence to establish phenotype relationship.\"
- confidence_label       : \"high\" | \"medium\" | \"low\"
- confidence_rationale   : 2-3 sentences covering species validity, biospecimen appropriateness,
  phenotype relevance, and how the origin profile influenced the confidence decision.

Confidence rules:
- high   : ALL dimensions hold; early-ranked metabolites are recognized endogenous substrates/products;
           exogenous_frac > 0.5 CANNOT be high
- medium : pathway fits species and biospecimen but one dimension is indirect/speculative,
           OR exogenous metabolites dominate
- low    : pathway fails species validity OR module is overwhelmingly exogenous/synthetic

Recovery rule (when all hypotheses were dropped in debate):
- You MAY recover at most ONE hypothesis to \"uncertain\" status if the Proposer was overly aggressive
- If you agree all drops were justified, return all as \"rejected\"

Return STRICT JSON only:
{{
  \"verdicts\": [
    {{
      \"hypothesis_id\": \"H1\",
      \"status\": \"accepted | uncertain | rejected\",
      \"verdict_rationale\": \"...\",
      \"phenotype_relationship\": \"...\",
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
                                       pubmed_per_hyp = NULL) {

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

  # PubMed section
  pubmed_text <- format_pubmed_for_arbiter(pubmed_per_hyp)

  glue::glue(
    "{hyp_section}

--- Debate history ---
{debate_digest}

--- Study context & biology ---
{origin_text}
{ctx_text}

{pubmed_text}"
  )
}


# ============================================================
# SECTION 7 – PubMed orchestration for advanced mode
# ============================================================

# Layer 1: shared metabolite search (reuses existing function, called once)
# Layer 2: per-hypothesis title search (lightweight, retmax=3)
# Each hypothesis gets: its own title papers + shared metabolite papers, deduplicated by PMID
search_pubmed_advanced <- function(hypotheses_final,
                                    top_metabolite_names = NULL,
                                    phenotype            = NULL,
                                    NCBI_api_key         = NULL,
                                    retmax_per_query     = 3L,
                                    years                = 5) {

  if (length(hypotheses_final) == 0) return(list())

  # Layer 1: shared metabolite search
  shared <- tryCatch(
    search_pubmed_for_annotation_evidence(
      module_title         = NULL,
      top_metabolite_names = top_metabolite_names,
      phenotype            = phenotype,
      NCBI_api_key         = NCBI_api_key,
      retmax_per_query     = retmax_per_query,
      years                = years
    ),
    error = function(e) NULL
  )

  # Layer 2: per-hypothesis title search
  per_hyp <- lapply(hypotheses_final, function(h) {
    title <- if (nzchar(h$title)) h$title else NULL
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

    adv_merge_pubmed(hyp_specific, shared)
  })
  names(per_hyp) <- sapply(hypotheses_final, `[[`, "hypothesis_id")
  per_hyp
}

# Merge two pubmed result lists, deduplicating by PMID.
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
    return(list(pmids=character(0), titles=character(0),
                abstracts=character(0), source_type=character(0), source_term=character(0)))
  combined <- combined[!duplicated(combined$pmid), , drop = FALSE]
  list(
    pmids       = combined$pmid,
    titles      = combined$title,
    abstracts   = combined$abstract,
    source_type = combined$source_type,
    source_term = combined$source_term
  )
}

# Format per-hypothesis PubMed results for Arbiter prompt (~300 chars per abstract)
format_pubmed_for_arbiter <- function(pubmed_per_hyp) {
  if (is.null(pubmed_per_hyp) || length(pubmed_per_hyp) == 0)
    return("Literature evidence: not retrieved (pubmed = FALSE).")

  sections <- lapply(names(pubmed_per_hyp), function(hid) {
    pm <- pubmed_per_hyp[[hid]]
    if (is.null(pm) || length(pm$pmids) == 0)
      return(glue::glue("[{hid}] No PubMed evidence retrieved."))

    n <- min(length(pm$pmids), length(pm$titles), length(pm$abstracts))
    papers <- vapply(seq_len(n), function(i) {
      abs_trunc <- if (nchar(pm$abstracts[i]) > 300)
        paste0(substr(pm$abstracts[i], 1, 300), "...") else pm$abstracts[i]
      sprintf("  [PMID:%s] %s\n    %s", pm$pmids[i], pm$titles[i], abs_trunc)
    }, character(1))
    paste(c(glue::glue("[{hid}] Literature evidence:"), papers), collapse = "\n")
  })

  paste(c("--- PubMed evidence (per hypothesis) ---",
          paste(sections, collapse = "\n\n")), collapse = "\n")
}

# Helper: collapse all PMIDs across ALL hypotheses (full-archive, same for every output row)
collapse_pubmed_pmids <- function(pph) {
  if (is.null(pph) || length(pph) == 0) return(NA_character_)
  all_pmids <- unique(unlist(lapply(pph, function(h) h$pmids %||% character(0))))
  all_pmids <- all_pmids[!is.na(all_pmids) & nzchar(trimws(all_pmids))]
  if (length(all_pmids) == 0) return(NA_character_)
  paste(all_pmids, collapse = "; ")
}

# Helper: serialize full pubmed_per_hyp to JSON (all hypotheses, full archive)
serialize_pubmed_evidence <- function(pph) {
  if (is.null(pph) || length(pph) == 0) return(NA_character_)
  tryCatch(
    as.character(jsonlite::toJSON(pph, auto_unbox = TRUE)),
    error = function(e) NA_character_
  )
}

# Helper: serialize ONE hypothesis's PubMed evidence to JSON (per-row column)
serialize_pubmed_evidence_one <- function(hyp_pph) {
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
            title               = list(type = "STRING"),
            summary             = list(type = "STRING"),
            claim               = list(type = "STRING"),
            anchor_metabolites  = list(type = "ARRAY", items = list(type = "STRING")),
            tension_metabolites = list(type = "ARRAY", items = list(type = "STRING"))
          ),
          required = c("hypothesis_id", "title", "summary", "claim",
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
                                        enum = c("major", "minor"))
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
            revised_title               = list(type = "STRING"),
            revised_summary             = list(type = "STRING"),
            revised_claim               = list(type = "STRING"),
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
                                        enum = c("accepted", "uncertain", "rejected")),
            verdict_rationale    = list(type = "STRING"),
            phenotype_relationship = list(type = "STRING"),
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


# ============================================================
# SECTION 9 – Output assembly
# ============================================================

assemble_advanced_output_row <- function(module_obj,
                                          arbiter_result,
                                          hypotheses_v0,
                                          hypotheses_final,
                                          debate_digest,
                                          pubmed_per_hyp,
                                          raw_texts) {

  mid <- module_obj$module_id %||% NA_character_

  verdicts     <- arbiter_result$verdicts     %||% list()
  residual_amb <- arbiter_result$residual_ambiguity %||% NA_character_

  # Status counts
  statuses   <- if (length(verdicts) > 0) sapply(verdicts, `[[`, "status") else character(0)
  n_accepted <- sum(statuses == "accepted",  na.rm = TRUE)
  n_uncertain<- sum(statuses == "uncertain", na.rm = TRUE)
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
      c(h_obj, list(
        status                 = v$status,
        verdict_rationale      = v$verdict_rationale,
        phenotype_relationship = v$phenotype_relationship %||% NA_character_,
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
    data.frame(
      module_id                  = mid,
      hypothesis_id              = hid,
      verdict_status             = v$status %||% NA_character_,
      title                      = hobj$title   %||% NA_character_,
      functional_summary         = hobj$summary %||% NA_character_,
      phenotype_relationship     = v$phenotype_relationship %||% NA_character_,
      confidence_label           = v$confidence_label       %||% NA_character_,
      rationale                  = v$confidence_rationale   %||% NA_character_,
      global_direction           = global_direction,
      core_drivers               = core_drivers,
      origin_endogenous_frac     = origin_endo,
      origin_exogenous_frac      = origin_exo,
      top_kegg_pathways          = kegg_pw,
      top_hmdb_pathways          = hmdb_pw,
      # per-hypothesis PubMed (this row only)
      pubmed_pmids               = collapse_pubmed_pmids(list(hyp_pub)),
      pubmed_evidence_json       = serialize_pubmed_evidence_one(hyp_pub),
      # full archive: all hypotheses' PubMed evidence
      all_pubmed_json            = serialize_pubmed_evidence(pubmed_per_hyp),
      raw_response_call1         = raw_texts$p0      %||% NA_character_,
      raw_response_call2         = raw_texts$arbiter %||% NA_character_,
      annotation_mode            = "advanced",
      n_proposed                 = length(hypotheses_v0),
      n_survived_debate          = length(hypotheses_final),
      n_accepted                 = n_accepted,
      n_uncertain                = n_uncertain,
      n_rejected                 = n_rejected,
      residual_ambiguity         = residual_amb,
      all_hypotheses_json        = all_hyp_json,
      stringsAsFactors           = FALSE
    )
  }

  # ---- One row per accepted/uncertain hypothesis ----
  non_rejected <- Filter(function(v) v$status %in% c("accepted", "uncertain"), verdicts)

  if (length(non_rejected) > 0) {
    do.call(rbind, lapply(non_rejected, make_row))
  } else {
    # All hypotheses rejected → one summary row (no per-hypothesis PubMed)
    data.frame(
      module_id                  = mid,
      hypothesis_id              = NA_character_,
      verdict_status             = "rejected",
      title                      = NA_character_,
      functional_summary         = NA_character_,
      phenotype_relationship     = NA_character_,
      confidence_label           = "rejected",
      rationale                  = NA_character_,
      global_direction           = global_direction,
      core_drivers               = core_drivers,
      origin_endogenous_frac     = origin_endo,
      origin_exogenous_frac      = origin_exo,
      top_kegg_pathways          = kegg_pw,
      top_hmdb_pathways          = hmdb_pw,
      pubmed_pmids               = collapse_pubmed_pmids(pubmed_per_hyp),
      pubmed_evidence_json       = NA_character_,
      all_pubmed_json            = serialize_pubmed_evidence(pubmed_per_hyp),
      raw_response_call1         = raw_texts$p0      %||% NA_character_,
      raw_response_call2         = raw_texts$arbiter %||% NA_character_,
      annotation_mode            = "advanced",
      n_proposed                 = length(hypotheses_v0),
      n_survived_debate          = length(hypotheses_final),
      n_accepted                 = 0L,
      n_uncertain                = 0L,
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
#' @return Data frame compatible with fast mode output plus advanced columns.
#'   Accepted and uncertain hypotheses are returned as separate rows; rejected
#'   modules return one rejected summary row. PubMed evidence is exported in
#'   \code{pubmed_pmids} and \code{pubmed_evidence_json}.
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
                                          years                = 5) {

  mid <- module_obj$module_id %||% "?"
  is_gemini <- identical(provider, "gemini")

  # Gemini schemas (NULL for non-Gemini)
  sch_proposer  <- if (is_gemini) adv_gemini_schema_proposer()  else NULL
  sch_skeptic   <- if (is_gemini) adv_gemini_schema_skeptic()   else NULL
  sch_revision  <- if (is_gemini) adv_gemini_schema_revision()  else NULL
  sch_arbiter   <- if (is_gemini) adv_gemini_schema_arbiter()   else NULL

  raw_texts <- list()

  # ----------------------------------------------------------------
  # Build reusable inputs (computed once, reused across all agents)
  # ----------------------------------------------------------------
  compressed_summary <- build_compressed_module_summary(module_obj)
  context_block      <- build_context_block(module_obj, study_context)

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

  hypotheses_v0 <- parse_proposer_output(p0_res$parsed)

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
    fast_row$verdict_status        <- NA_character_
    fast_row$n_proposed            <- 0L
    fast_row$n_survived_debate     <- NA_integer_
    fast_row$n_accepted            <- NA_integer_
    fast_row$n_uncertain           <- NA_integer_
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
  s1_usr <- build_skeptic_user_prompt(compressed_summary, context_block,
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
  s1_challenges <- parse_skeptic_output(s1_res$parsed)

  # ----------------------------------------------------------------
  # P1 – Proposer revision 1
  # ----------------------------------------------------------------
  message("[advanced] Module ", mid, " | P1 Proposer revision (round 1)...")

  p1_sys <- build_proposer_revision_system_prompt()
  p1_usr <- build_proposer_revision_user_prompt(
    compressed_summary, context_block,
    hypotheses_v0, s1_challenges, history_text = NULL
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
  p1_revisions      <- parse_revision_output(p1_res$parsed)
  hypotheses_after_r1 <- apply_revisions(hypotheses_v0, p1_revisions)

  # ----------------------------------------------------------------
  # Round 2 decision
  # ----------------------------------------------------------------
  do_r2 <- should_run_round2(hypotheses_after_r1, s1_challenges, p1_revisions)

  hypotheses_final <- hypotheses_after_r1
  s2_challenges    <- NULL
  p2_revisions     <- NULL

  if (do_r2) {
    # Build compact history for round 2 prompts
    r1_history_text <- paste(
      sapply(p1_revisions, function(r)
        sprintf("%s: %s — %s", r$hypothesis_id, toupper(r$action), r$justification)),
      collapse = "\n"
    )

    # S2 – Skeptic round 2
    message("[advanced] Module ", mid, " | S2 Skeptic (round 2)...")

    s2_sys <- build_skeptic_system_prompt(round = 2L)
    s2_usr <- build_skeptic_user_prompt(
      compressed_summary, context_block,
      hypotheses_after_r1, history_text = r1_history_text
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
    s2_challenges <- parse_skeptic_output(s2_res$parsed)

    # P2 – Proposer revision 2
    message("[advanced] Module ", mid, " | P2 Proposer revision (round 2)...")

    p2_sys <- build_proposer_revision_system_prompt()
    p2_usr <- build_proposer_revision_user_prompt(
      compressed_summary, context_block,
      hypotheses_after_r1, s2_challenges,
      history_text = r1_history_text
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
    p2_revisions     <- parse_revision_output(p2_res$parsed)
    hypotheses_final <- apply_revisions(hypotheses_after_r1, p2_revisions)
  }

  # ----------------------------------------------------------------
  # PubMed (optional, always before Arbiter regardless of debate outcome)
  # When all hypotheses were dropped, search the original v0 hypotheses so
  # the Arbiter has evidence to support or reject any potential recovery.
  # ----------------------------------------------------------------
  pubmed_per_hyp <- NULL
  if (pubmed) {
    hyps_for_pubmed <- if (length(hypotheses_final) > 0) hypotheses_final else hypotheses_v0
    message("[advanced] Module ", mid, " | PubMed search for ",
            length(hyps_for_pubmed), " hypothesis(es) (pre-Arbiter)...")
    top_names <- NULL
    if (!is.null(module_obj$topN_metabolites) &&
        "Compound_name" %in% colnames(module_obj$topN_metabolites))
      top_names <- unique(module_obj$topN_metabolites$Compound_name)

    pubmed_per_hyp <- search_pubmed_advanced(
      hypotheses_final     = hyps_for_pubmed,
      top_metabolite_names = top_names,
      phenotype            = study_context$phenotype,
      NCBI_api_key         = NCBI_api_key,
      retmax_per_query     = 3L,
      years                = years
    )
  }

  # ----------------------------------------------------------------
  # Build debate history digest for Arbiter
  # ----------------------------------------------------------------
  debate_digest <- build_debate_history_digest(
    hypotheses_v0    = hypotheses_v0,
    s1_challenges    = s1_challenges,
    p1_revisions     = p1_revisions,
    hypotheses_v1    = hypotheses_after_r1,
    s2_challenges    = s2_challenges,
    p2_revisions     = p2_revisions,
    hypotheses_final = hypotheses_final
  )

  # ----------------------------------------------------------------
  # A0 – Arbiter
  # ----------------------------------------------------------------
  message("[advanced] Module ", mid, " | A0 Arbiter (final verdict)...")

  a0_sys <- build_arbiter_system_prompt()
  a0_usr <- build_arbiter_user_prompt(
    hypotheses_final = hypotheses_final,
    hypotheses_v0    = hypotheses_v0,
    debate_digest    = debate_digest,
    module_obj       = module_obj,
    study_context    = study_context,
    pubmed_per_hyp   = pubmed_per_hyp
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
          hypothesis_id          = trimws(as.character(v$hypothesis_id %||% "")),
          status                 = trimws(tolower(as.character(v$status %||% "rejected"))),
          verdict_rationale      = trimws(as.character(v$verdict_rationale %||% "")),
          phenotype_relationship = trimws(as.character(v$phenotype_relationship %||% "")),
          confidence_label       = trimws(tolower(as.character(v$confidence_label %||% ""))),
          confidence_rationale   = trimws(as.character(v$confidence_rationale %||% ""))
        )
      })
    }
    rab <- a0_res$parsed$residual_ambiguity
    if (!is.null(rab) && !is.na(rab) && nzchar(trimws(as.character(rab))))
      arbiter_result$residual_ambiguity <- trimws(as.character(rab))
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
                                           years                = 5) {

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
        years                = years
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
