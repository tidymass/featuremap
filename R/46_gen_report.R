#' Annotate a Single Module Using Two Sequential LLM Calls  (Fast Mode)
#'
#' Implements the Fast Mode two-call pipeline for one metabolic module:
#'
#' \enumerate{
#'   \item \strong{Call 1} — Biochemical annotation (no study context):
#'         generates \code{title} and \code{functional_summary} from metabolite
#'         structure, class composition, and pathway coverage.
#'   \item \strong{PubMed search} (optional) — retrieves phenotype-relevant
#'         abstracts using the module title and study phenotype.
#'   \item \strong{Call 2} — Phenotype and confidence assessment:
#'         generates \code{phenotype_relationship}, \code{confidence_label},
#'         and \code{rationale} using origin profile, study context, and the
#'         PubMed evidence collected in step 2.
#' }
#'
#' @param module_obj List produced by \code{prepare_module()}. Must contain:
#'   \code{module_id}, \code{module_stats}, \code{topN_metabolites},
#'   \code{other_metabolites}, \code{pathway_coverage}, \code{global_direction},
#'   \code{origin_profile}.
#' @param study_context List with \code{species}, \code{biospecimen},
#'   \code{phenotype}.  Used only in Call 2.
#' @param provider Character. LLM provider (\code{"openai"}, \code{"siliconflow"},
#'   \code{"gemini"}).
#' @param llm_api_key Character. API key for the LLM service.
#' @param model Character. Model identifier (e.g. \code{"gpt-4.1-mini"}).
#' @param max_tokens Integer. Maximum response tokens (default: 5000).
#' @param temperature Numeric. Sampling temperature, 0 for deterministic
#'   (default: 0).
#' @param pubmed Logical. Whether to search PubMed after Call 1 (default: FALSE).
#' @param NCBI_api_key Character or NULL. NCBI API key for PubMed.
#' @param years Numeric. PubMed look-back years (default: 5).
#'
#' @return A one-row data frame with columns:
#'   \code{module_id}, \code{title}, \code{functional_summary},
#'   \code{phenotype_relationship}, \code{confidence_label}, \code{rationale},
#'   \code{global_direction}, \code{core_drivers},
#'   \code{origin_endogenous_frac}, \code{origin_exogenous_frac},
#'   \code{top_kegg_pathways}, \code{top_hmdb_pathways},
#'   \code{raw_response_call1}, \code{raw_response_call2}.
annotate_one_module <- function(module_obj,
                                study_context,
                                provider,
                                llm_api_key,
                                model,
                                max_tokens   = 5000,
                                temperature  = 0,
                                pubmed       = FALSE,
                                NCBI_api_key = NULL,
                                years        = 5) {

  mid <- module_obj$module_id %||% NA_character_

  # ------------------------------------------------------------------
  # Internal: normalise a parsed LLM response to a named list
  # ------------------------------------------------------------------
  normalize_parsed <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.data.frame(x) && nrow(x) >= 1) return(as.list(x[1, , drop = FALSE]))
    if (is.list(x) && !is.null(names(x))) return(x)
    if (is.list(x) && length(x) >= 1 && is.list(x[[1]]) && !is.null(names(x[[1]])))
      return(x[[1]])
    NULL
  }

  # ------------------------------------------------------------------
  # Internal: extract the first JSON object from raw text
  # ------------------------------------------------------------------
  parse_json_response <- function(txt) {
    direct <- tryCatch(
      jsonlite::fromJSON(txt, simplifyVector = FALSE),
      error = function(e) NULL
    )
    direct <- normalize_parsed(direct)
    if (!is.null(direct)) return(direct)

    json_obj_txt <- tryCatch(extract_json_object(txt), error = function(e) NULL)
    if (!is.null(json_obj_txt)) {
      obj <- tryCatch(
        jsonlite::fromJSON(json_obj_txt, simplifyVector = FALSE),
        error = function(e) NULL
      )
      obj <- normalize_parsed(obj)
      if (!is.null(obj)) return(obj)
    }

    json_fence <- regmatches(txt, regexpr("```\\s*[\\s\\S]*?\\s*```", txt, perl = TRUE))
    if (length(json_fence) > 0 && nzchar(json_fence[1])) {
      inner <- sub("^```[a-zA-Z]*\\s*", "", json_fence[1], perl = TRUE)
      inner <- sub("\\s*```$", "", inner, perl = TRUE)
      fobj  <- tryCatch(
        jsonlite::fromJSON(inner, simplifyVector = FALSE),
        error = function(e) NULL
      )
      fobj <- normalize_parsed(fobj)
      if (!is.null(fobj)) return(fobj)
    }
    NULL
  }

  # ------------------------------------------------------------------
  # Internal: partial-JSON key extractor
  # ------------------------------------------------------------------
  parse_partial_json_fallback <- function(txt) {
    if (is.null(txt) || !nzchar(txt)) return(NULL)

    extract_str <- function(key) {
      pattern <- paste0('"', key, '"\\s*:\\s*"')
      m <- regexpr(pattern, txt, perl = TRUE)
      if (m[1] == -1) return(NA_character_)
      start <- m[1] + attr(m, "match.length")
      if (start > nchar(txt)) return(NA_character_)
      chars <- strsplit(substr(txt, start, nchar(txt)), "", fixed = TRUE)[[1]]
      out <- character(0); esc <- FALSE; closed <- FALSE
      for (ch in chars) {
        if (!esc && ch == "\\") { esc <- TRUE; out <- c(out, ch); next }
        if (!esc && ch == "\"") { closed <- TRUE; break }
        out <- c(out, ch); esc <- FALSE
      }
      val <- trimws(gsub("\\\\n", " ", paste(out, collapse = ""), perl = TRUE))
      if (!nzchar(val)) NA_character_ else val
    }

    keys <- c("title", "functional_summary",
              "phenotype_relationship", "confidence_label", "rationale")
    res <- setNames(lapply(keys, extract_str), keys)

    for (kk in c("confidence_label")) {
      v <- res[[kk]]
      if (!is.na(v)) {
        if      (grepl("\\bhigh\\b",   v, ignore.case = TRUE)) res[[kk]] <- "high"
        else if (grepl("\\bmedium\\b", v, ignore.case = TRUE)) res[[kk]] <- "medium"
        else if (grepl("\\blow\\b",    v, ignore.case = TRUE)) res[[kk]] <- "low"
        else res[[kk]] <- NA_character_
      }
    }

    found_n <- sum(!is.na(unlist(res)))
    if (found_n == 0) NULL else res
  }

  # ------------------------------------------------------------------
  # Internal: validate all required keys are present and non-empty
  # ------------------------------------------------------------------
  validate_required_keys <- function(x, required_keys) {
    if (is.null(x)) return(NULL)
    ok <- vapply(required_keys, function(k) {
      v <- x[[k]]
      !is.null(v) && length(v) >= 1 && !is.na(v[[1]]) && nzchar(trimws(as.character(v[[1]])))
    }, logical(1))
    if (all(ok)) x else NULL
  }

  # ------------------------------------------------------------------
  # Internal: send one request with 2-attempt retry
  # ------------------------------------------------------------------
  request_json <- function(messages, schema, required_keys, retry_instruction,
                            stage_label) {
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
          is.character(res_call$content) && length(res_call$content) == 1 &&
          nzchar(res_call$content)) {
        raw_text <- res_call$content
        parsed   <- validate_required_keys(parse_json_response(raw_text), required_keys)
        if (!is.null(parsed)) break

        message("[annotate_module] No valid JSON (module ", mid, ", ",
                stage_label, ", attempt ", attempt, ")")

        # Gemini truncation workaround
        if (identical(provider, "gemini") &&
            !is.null(res_call$finish_reason) &&
            identical(toupper(as.character(res_call$finish_reason)), "MAX_TOKENS")) {
          cur_max <- min(max(3000L, cur_max * 2L), 8192L)
          message("[annotate_module] Gemini MAX_TOKENS; retrying with max_tokens = ", cur_max)
        }
      } else {
        message("[annotate_module] Empty response (module ", mid, ", ",
                stage_label, ", attempt ", attempt, ")")
      }

      if (attempt == 1) {
        messages <- c(messages,
                      list(list(role = "user", content = retry_instruction)))
      }
    }

    # Fallback parsing
    if (is.null(parsed))
      parsed <- validate_required_keys(parse_partial_json_fallback(raw_text), required_keys)

    list(parsed = parsed, raw_text = raw_text)
  }

  # ------------------------------------------------------------------
  # Helper: format pathway coverage as compact text
  # ------------------------------------------------------------------
  pathway_to_text <- function(df) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NA_character_)
    rows <- apply(df, 1, function(r) {
      r <- as.list(r)
      cov_pct <- if (!is.null(r$coverage) && is.finite(as.numeric(r$coverage)))
        paste0(round(as.numeric(r$coverage) * 100, 1), "%")
      else "NA"
      sprintf("%s (%d/%d, %s)",
              trimws(as.character(r$pathway_name %||% "?")),
              as.integer(r$n_hit %||% 0L),
              as.integer(r$n_metabolites %||% 0L),
              cov_pct)
    })
    paste(rows, collapse = "; ")
  }

  # ------------------------------------------------------------------
  # Gemini response schema helpers
  # ------------------------------------------------------------------
  schema_call1 <- NULL
  schema_call2 <- NULL
  if (identical(provider, "gemini")) {
    schema_call1 <- list(
      type = "OBJECT",
      properties = list(
        title              = list(type = "STRING"),
        functional_summary = list(type = "STRING")
      ),
      required = c("title", "functional_summary")
    )
    schema_call2 <- list(
      type = "OBJECT",
      properties = list(
        phenotype_relationship = list(type = "STRING"),
        confidence_label       = list(type = "STRING",
                                      enum = c("high", "medium", "low")),
        rationale              = list(type = "STRING")
      ),
      required = c("phenotype_relationship", "confidence_label", "rationale")
    )
  }

  # ==================================================================
  # CALL 1 – Biochemical annotation (title + functional_summary)
  # No study context, no confidence judgment
  # ==================================================================
  sys1  <- build_system_prompt()
  usr1  <- build_module_prompt(module_obj)

  res1  <- request_json(
    messages = list(
      list(role = "system", content = sys1),
      list(role = "user",   content = usr1)
    ),
    schema           = schema_call1,
    required_keys    = c("title", "functional_summary"),
    retry_instruction = paste(
      "Return ONLY a valid JSON object with exactly two keys:",
      "\"title\" and \"functional_summary\"."
    ),
    stage_label = "call1 (title+summary)"
  )

  title              <- NA_character_
  functional_summary <- NA_character_
  if (!is.null(res1$parsed)) {
    title              <- as.character(res1$parsed[["title"]]              %||% NA_character_)
    functional_summary <- as.character(res1$parsed[["functional_summary"]] %||% NA_character_)
  } else {
    message("[annotate_module] Call 1 failed for module ", mid,
            ". Returning NA for title/summary; continuing to Call 2.")
  }

  # ==================================================================
  # PUBMED SEARCH (optional) – module title plus every top metabolite
  # ==================================================================
  pubmed_context <- NULL

  if (pubmed && !is.na(title) && nzchar(trimws(title))) {
    top_met_names <- NULL
    if (!is.null(module_obj$topN_metabolites) &&
        nrow(module_obj$topN_metabolites) > 0 &&
        "Compound_name" %in% colnames(module_obj$topN_metabolites)) {
      top_met_names <- unique(module_obj$topN_metabolites$Compound_name)
    }
    phenotype_val <- NULL
    if (!is.null(study_context$phenotype)) {
      ph <- trimws(as.character(study_context$phenotype[1]))
      if (nzchar(ph) && !identical(ph, "NA")) phenotype_val <- ph
    }

    pubmed_context <- tryCatch(
      search_pubmed_for_annotation_evidence(
        module_title         = title,
        top_metabolite_names = top_met_names,
        phenotype            = phenotype_val,
        NCBI_api_key         = NCBI_api_key,
        retmax_per_query     = 3,
        years                = years
      ),
      error = function(e) {
        message("[annotate_module] PubMed search failed for module ", mid, ": ", e$message)
        NULL
      }
    )
  }

  # ==================================================================
  # CALL 2 – Phenotype relationship + confidence + rationale
  # Uses origin profile, study context, PubMed evidence
  # ==================================================================
  sys2 <- build_phenotype_confidence_system_prompt()
  usr2 <- build_phenotype_confidence_prompt(
    module_obj     = module_obj,
    title          = title,
    study_context  = study_context,
    pubmed_context = pubmed_context
  )

  res2 <- request_json(
    messages = list(
      list(role = "system", content = sys2),
      list(role = "user",   content = usr2)
    ),
    schema            = schema_call2,
    required_keys     = c("phenotype_relationship", "confidence_label", "rationale"),
    retry_instruction = paste(
      "Return ONLY a valid JSON object with exactly three keys:",
      "\"phenotype_relationship\", \"confidence_label\", \"rationale\"."
    ),
    stage_label = "call2 (confidence)"
  )

  phenotype_relationship <- NA_character_
  confidence_label       <- NA_character_
  rationale              <- NA_character_
  if (!is.null(res2$parsed)) {
    phenotype_relationship <-
      as.character(res2$parsed[["phenotype_relationship"]] %||% NA_character_)
    confidence_label <-
      as.character(res2$parsed[["confidence_label"]]       %||% NA_character_)
    rationale <-
      as.character(res2$parsed[["rationale"]]              %||% NA_character_)
  }

  # ==================================================================
  # Computed fields (no LLM)
  # ==================================================================

  # Global direction
  global_direction <- NA_character_
  if (!is.null(module_obj$global_direction)) {
    global_direction <- module_obj$global_direction$direction_label %||% NA_character_
  }

  # Core drivers (high-priority top metabolites)
  core_drivers <- NA_character_
  if (!is.null(module_obj$topN_metabolites) &&
      is.data.frame(module_obj$topN_metabolites) &&
      nrow(module_obj$topN_metabolites) > 0 &&
      "Compound_name" %in% colnames(module_obj$topN_metabolites)) {

    tdf <- module_obj$topN_metabolites
    if ("priority_level" %in% colnames(tdf)) {
      drivers <- tdf[tdf$priority_level == "high", , drop = FALSE]
      if ("top_rank" %in% colnames(drivers))
        drivers <- drivers[order(as.integer(drivers$top_rank)), , drop = FALSE]
    } else {
      drivers <- tdf
    }
    if (nrow(drivers) == 0)
      drivers <- head(tdf[order(as.integer(tdf$top_rank %||% seq_len(nrow(tdf)))), ], 3)
    core_drivers <- paste(drivers$Compound_name, collapse = "; ")
  }

  # Origin fractions
  origin_endo_frac <- NA_real_
  origin_exo_frac  <- NA_real_
  if (!is.null(module_obj$origin_profile)) {
    origin_endo_frac <- module_obj$origin_profile$endogenous_frac %||% NA_real_
    origin_exo_frac  <- module_obj$origin_profile$exogenous_frac  %||% NA_real_
  }

  # Pathway text summaries
  kegg_pathways <- NA_character_
  hmdb_pathways <- NA_character_
  if (!is.null(module_obj$pathway_coverage)) {
    kegg_pathways <- pathway_to_text(module_obj$pathway_coverage$kegg)
    hmdb_pathways <- pathway_to_text(module_obj$pathway_coverage$hmdb)
  }

  # ==================================================================
  # Return
  # ==================================================================
  data.frame(
    module_id              = mid,
    title                  = title,
    functional_summary     = functional_summary,
    phenotype_relationship = phenotype_relationship,
    confidence_label       = confidence_label,
    rationale              = rationale,
    global_direction       = global_direction,
    core_drivers           = core_drivers,
    origin_endogenous_frac = origin_endo_frac,
    origin_exogenous_frac  = origin_exo_frac,
    top_kegg_pathways      = kegg_pathways,
    top_hmdb_pathways      = hmdb_pathways,
    raw_response_call1     = res1$raw_text  %||% NA_character_,
    raw_response_call2     = res2$raw_text  %||% NA_character_,
    stringsAsFactors = FALSE
  )
}


#' Annotate All Modules Using the Two-Call Fast Mode Pipeline
#'
#' Iterates over a list of prepared module objects and generates functional
#' annotations for each using the sequential LLM Call 1 → PubMed → Call 2
#' workflow.
#'
#' @param all_modules_ready Named list of module objects from
#'   \code{prepare_module()}.
#' @param study_context List with \code{species}, \code{biospecimen},
#'   \code{phenotype}.
#' @param provider Character. LLM provider.
#' @param llm_api_key Character. API key.
#' @param model Character. Model name.
#' @param max_tokens Integer. Max response tokens (default: 2000).
#' @param temperature Numeric. Sampling temperature (default: 0).
#' @param sleep_sec Numeric. Pause between modules in seconds (default: 0).
#' @param pubmed Logical. Whether to search PubMed after Call 1 (default: FALSE).
#' @param NCBI_api_key Character or NULL. NCBI API key.
#' @param years Numeric. PubMed look-back years (default: 5).
#'
#' @return A data frame with one row per module containing all annotation
#'   columns produced by \code{annotate_one_module()}.
annotate_all_modules <- function(all_modules_ready,
                                  study_context,
                                  provider,
                                  llm_api_key,
                                  model,
                                  max_tokens   = 2000,
                                  temperature  = 0,
                                  sleep_sec    = 0,
                                  pubmed       = FALSE,
                                  NCBI_api_key = NULL,
                                  years        = 5) {

  n_mod <- length(all_modules_ready)

  results_list <- lapply(seq_len(n_mod), function(i) {
    mod <- all_modules_ready[[i]]
    message("[annotate_module] Processing module ", mod$module_id %||% i,
            " (", i, "/", n_mod, ")")

    out <- annotate_one_module(
      module_obj    = mod,
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

    if (sleep_sec > 0) Sys.sleep(sleep_sec)
    out
  })

  dplyr::bind_rows(results_list)
}
