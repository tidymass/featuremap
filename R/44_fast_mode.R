
# ============================================================
# LLM Prompt Builders  (Fast Mode)
#
# Call 1:  build_system_prompt()  +  build_module_prompt()
#          → {type, title, functional_summary}
#          Inputs: Signal_concentration, top metabolites, class structure,
#                  pathway coverage (KEGG + HMDB)
#          No study context; no confidence assessment
#
# Call 2:  build_phenotype_confidence_system_prompt()
#          + build_phenotype_confidence_prompt()
#          → {phenotype_relationship, confidence_label, rationale}
#          Inputs: module type/title/summary, origin profile, study context,
#                  PubMed abstracts (optional)
# ============================================================


# ============================================================
# CALL 1 – system prompt
# ============================================================

build_system_prompt <- function() {
  glue::glue(
    "You are an expert in metabolomics and biochemical pathway interpretation.
Your task is to name and summarize one dysregulated metabolic module from an omics study.

For each module, you will receive:
- Module-level label: Signal_concentration (high | medium | low)
- Top metabolites: rank, metabolite name, chemical class, Priority_level
- Module classes: class name, member metabolite names, proportion, Priority_level
- Pathway evidence from KEGG and HMDB databases

Definitions:
- Signal_concentration: whether phenotype-associated Importance is concentrated in a few metabolites (high) or distributed across many metabolites (low).
- Priority_level: interpretation priority within the module; early-ranked high-priority metabolites should anchor the main interpretation.

Using only the provided information, return three fields:

1) TYPE
   - Choose exactly one:
     Host-Endogenous | Microbiome-Derived | Diet-Derived | Xenobiotic / Exogenous | Host–Microbiome
   - Type is an interpretation scope for the module annotation, not a confidence score.
   - Assign type from the module evidence, especially top-ranked high-priority metabolites;
     use chemical class and KEGG/HMDB pathway evidence as support.
   - Assign type by the biological origin of the module's core substrate, dominant metabolite
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
   - The type, title, and functional_summary must describe the same central biochemical theme.
   - Do not use or mention phenotype/disease when assigning type.

2) TITLE
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

3) FUNCTIONAL_SUMMARY
   - Write 4-6 sentences as a factual biochemical annotation of the module.
   - Describe what the metabolites do biochemically: enzymatic conversions, pathway steps, metabolite roles, and physiological functions — as if writing an entry in a biochemical database.
   - Base the explanation mainly on early-ranked high-priority metabolites; use class structure and pathway evidence only to support or qualify.
   - Keep origin/scope interpretation in the type field; in the summary, mention source or host/microbial
     processing only when it is necessary to explain the biochemical mechanism.
   - Do NOT comment on the study phenotype or disease relevance.
   - Do NOT reference module-level statistics (Signal_concentration).
   - Mention secondary signals only briefly at the end if needed.

4) ANCHOR_METABOLITES
   - List 2-5 metabolites that most strongly anchor the module interpretation.
   - Choose from the provided top metabolites, prioritising early-ranked high-priority ones.
   - These will be used for literature retrieval; prefer metabolites with established names
     (avoid generic class names).
   - Return as a JSON array of strings.

Return STRICT JSON only:
{{
  \"type\": \"Host-Endogenous | Microbiome-Derived | Diet-Derived | Xenobiotic / Exogenous | Host–Microbiome\",
  \"title\": \"...\",
  \"functional_summary\": \"...\",
  \"anchor_metabolites\": [\"metabolite1\", \"metabolite2\"]
}}

Do not invent metabolites, unsupported pathways, or unsupported mechanisms."
  )
}


# ============================================================
# CALL 1 – module user prompt
# ============================================================

build_module_prompt <- function(module_obj) {

  scalar_text <- function(x, na = "NA") {
    if (is.null(x) || length(x) == 0) return(na)
    x <- as.character(x[1])
    x <- trimws(x)
    if (!nzchar(x) || is.na(x) || toupper(x) %in% c("NA", "NULL", "UNKNOWN")) na else x
  }

  # ---- Module ID ----
  mid <- scalar_text(module_obj$module_id, na = "NA")

  # ---- Signal_concentration (only stat in fast mode) ----
  stats <- module_obj$module_stats %||% list()
  signal_level <- scalar_text(stats$Signal_concentration)

  # ---- Top metabolites ----
  top_df <- module_obj$topN_metabolites

  if (is.null(top_df) || !is.data.frame(top_df) || nrow(top_df) == 0) {
    top_text <- "NA (no top metabolites available)"
  } else {
    if (!"functional_summary" %in% colnames(top_df)) {
      top_df$functional_summary <- NA_character_
    }

    top_blocks <- apply(top_df, 1, function(row) {
      row <- as.list(row)
      glue::glue(
        "Metabolite (top):
- Top_rank: {scalar_text(row$top_rank)}
- Name: {scalar_text(row$Compound_name)}
- Class: {scalar_text(row$class, na = 'not annotated')}
- Priority_level: {scalar_text(row$priority_level)}",
        .na = "NA"
      )
    })
    top_text <- paste(top_blocks, collapse = "\n\n")
  }

  # ---- Class-level structure ----
  class_df <- module_obj$other_metabolites

  if (is.null(class_df) || !is.data.frame(class_df) || nrow(class_df) == 0) {
    class_text <- "No class-level structure available."
  } else {
    class_blocks <- apply(class_df, 1, function(row) {
      row <- as.list(row)
      met_vec <- tryCatch(
        clean_field(row$metabolites),
        error = function(e) row$metabolites
      )
      met_str <- if (is.null(met_vec) || length(met_vec) == 0) {
        "NA"
      } else {
        paste(as.character(met_vec), collapse = ", ")
      }
      glue::glue(
        "Metabolite class:
- Class_name: {scalar_text(row$class)}
- Proportion_of_module: {scalar_text(row$n_metabolites)}
- Metabolite_names: {met_str}
- Priority_level: {scalar_text(row$priority_level)}",
        .na = "NA"
      )
    })
    class_text <- paste(class_blocks, collapse = "\n\n")
  }

  # ---- Pathway evidence ----
  pathway_text <- build_pathway_text(module_obj$pathway_coverage)

  # ---- Assemble prompt ----
  glue::glue("
Module information:
- Module ID: {mid}
- Signal_concentration: {signal_level}

Top metabolites (strongest statistical signals in the module):
Listed in descending priority order (smaller Top_rank = stronger evidence).

{top_text}

Class-level structure of the module (all metabolites grouped by class):

{class_text}

{pathway_text}")
}


# Helper: format pathway coverage for prompt
build_pathway_text <- function(pathway_coverage) {
  if (is.null(pathway_coverage)) {
    return("Pathway evidence: not available.")
  }

  format_pathway_table <- function(df, label) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
      return(paste0(label, ": no pathways matched."))
    }
    rows <- apply(df, 1, function(r) {
      r <- as.list(r)
      n_hit  <- as.integer(r$n_hit %||% 0)
      n_tot  <- as.integer(r$n_metabolites %||% 0)
      cov_pct <- if (!is.null(r$coverage) && is.finite(as.numeric(r$coverage))) {
        paste0(round(as.numeric(r$coverage) * 100, 1), "%")
      } else "NA"
      pname <- trimws(as.character(r$pathway_name %||% "Unknown"))
      sprintf("  %s — %d of %d module metabolites (coverage: %s)",
              pname, n_hit, n_tot, cov_pct)
    })
    paste0(label, ":\n", paste(rows, collapse = "\n"))
  }

  kegg_text <- format_pathway_table(pathway_coverage$kegg,
                                     "Pathway evidence (KEGG)")
  hmdb_text <- format_pathway_table(pathway_coverage$hmdb,
                                     "Pathway evidence (HMDB, supplemental)")

  paste(kegg_text, hmdb_text, sep = "\n\n")
}


# ============================================================
# CALL 2 – system prompt
# ============================================================

build_phenotype_confidence_system_prompt <- function() {
  glue::glue(
    "You are an expert in metabolomics and biochemical pathway interpretation.
Your task is to assess phenotype relationship and evidence strength for a previously
generated fast-mode module annotation.

You will receive:
- Module type, title, and functional summary from prior biochemical analysis
- Origin profile: proportions of endogenous, exogenous, and unknown-origin metabolites
- Study context: species, biospecimen, phenotype
- Literature evidence: PubMed abstracts retrieved by searching the module title,
  anchor metabolites (key metabolites identified in Call 1), and phenotype
  (may be absent if pubmed = FALSE or no papers were found)

Definitions:
- Endogenous: metabolites documented as produced by the host organism.
- Exogenous: metabolites of dietary, microbial, pharmaceutical, plant, or environmental origin.
- Unknown: HMDB-annotated metabolites with no recorded origin.

Return three fields:

1) PHENOTYPE_RELATIONSHIP
   - Write 2-4 sentences summarizing how the module's metabolites relate to the study phenotype,
     based on the provided PubMed evidence when available.
   - If phenotype is provided but PubMed evidence is absent or insufficient, you MAY infer a
     cautious phenotype relationship from general biochemical knowledge. When doing so,
     explicitly state that the relationship is inferred and lacks sufficient direct literature support.
   - If phenotype is not provided, write exactly:
     \"Insufficient literature evidence to establish phenotype relationship.\"

2) CONFIDENCE
   Confidence reflects the strength of evidence supporting the fast-mode annotation. It is
   not an accepted/rejected verdict.

   Scoring rules:
   - High: Strong module evidence AND adequate phenotype/literature support. Top metabolites,
     annotation type, class/pathway evidence, species, biospecimen, and available
     phenotype/literature evidence are coherent. High confidence requires at least indirect
     PubMed support; if no PubMed evidence is retrieved or the evidence is insufficient,
     confidence_label must not be high.
   - Medium: Reasonable module interpretation, but with partial evidence, indirect pathway
     support, incomplete explanation, limited literature support, or minor limitations.
   - Low: Weak but still useful interpretation. Evidence is sparse, indirect, partly unresolved,
     or the annotation type/summary has notable but non-fatal limitations.

   Origin rule:
   - Origin profile is an interpretation-scope check, not an automatic penalty.
   - High exogenous fraction does NOT automatically lower confidence if the annotation type is
     Diet-Derived, Xenobiotic / Exogenous, Microbiome-Derived, or otherwise consistent with
     the origin profile.
   - Lower confidence only when the annotation type or summary conflicts with the origin profile.

3) RATIONALE
   - 2-3 sentences explaining why the annotation received this evidence-strength label,
     including species, biospecimen, phenotype/literature support, and origin/type consistency
     when relevant.

Return STRICT JSON only:
{{
  \"phenotype_relationship\": \"...\",
  \"confidence_label\": \"high\" | \"medium\" | \"low\",
  \"rationale\": \"2-3 sentences.\"
}}

Do not invent PubMed support or claim direct literature evidence when it is not provided."
  )
}


# ============================================================
# CALL 2 – module user prompt
# ============================================================

build_phenotype_confidence_prompt <- function(module_obj,
                                               title,
                                               hypothesis_type = NULL,
                                               functional_summary = NULL,
                                               study_context  = NULL,
                                               pubmed_context = NULL) {

  scalar_text <- function(x, na = "NA") {
    if (is.null(x) || length(x) == 0) return(na)
    x <- as.character(x[1])
    x <- trimws(x)
    if (!nzchar(x) || is.na(x)) na else x
  }

  title_str <- scalar_text(title, na = "Unknown module")
  type_str  <- scalar_text(hypothesis_type, na = "not specified")
  summary_str <- scalar_text(functional_summary, na = "Not available.")

  # ---- Origin profile ----
  op <- module_obj$origin_profile
  if (!is.null(op) && !is.na(op$n_annotated) && as.integer(op$n_annotated) > 0) {
    n_ann  <- as.integer(op$n_annotated)
    endo_n <- as.integer(op$endogenous_n %||% NA_integer_)
    exo_n  <- as.integer(op$exogenous_n  %||% NA_integer_)
    unk_n  <- as.integer(op$unknown_n    %||% NA_integer_)
    endo_pct <- if (!is.na(endo_n)) paste0(round(endo_n / n_ann * 100, 1), "%") else "NA"
    exo_pct  <- if (!is.na(exo_n))  paste0(round(exo_n  / n_ann * 100, 1), "%") else "NA"
    unk_pct  <- if (!is.na(unk_n))  paste0(round(unk_n  / n_ann * 100, 1), "%") else "NA"
    origin_text <- glue::glue(
      "Origin profile (HMDB-annotated metabolites, n = {n_ann}):
- Endogenous: {endo_pct} ({endo_n}/{n_ann})
- Exogenous (dietary/microbial/drug/environmental): {exo_pct} ({exo_n}/{n_ann})
- Unknown origin: {unk_pct} ({unk_n}/{n_ann})"
    )
  } else {
    origin_text <- "Origin profile: not available (no HMDB-annotated metabolites with origin data)."
  }

  # ---- Study context ----
  has_context <- !is.null(study_context) && any(vapply(
    study_context[c("species", "biospecimen", "phenotype")],
    function(x) {
      if (is.null(x) || length(x) == 0) return(FALSE)
      val <- trimws(as.character(x[1]))
      nzchar(val) && !identical(val, "NA")
    },
    logical(1)
  ))

  if (has_context) {
    context_text <- glue::glue(
      "Study context:
- Species: {scalar_text(study_context$species)}
- Biospecimen: {scalar_text(study_context$biospecimen)}
- Phenotype: {scalar_text(study_context$phenotype)}"
    )
  } else {
    context_text <- "Study context: not provided."
  }

  # ---- PubMed evidence ----
  if (!is.null(pubmed_context) &&
      length(pubmed_context$abstracts) > 0 &&
      length(pubmed_context$titles)    > 0) {

    n_papers   <- min(length(pubmed_context$titles), length(pubmed_context$abstracts))
    ctx_value <- function(field, i, default = "NA") {
      x <- pubmed_context[[field]]
      if (is.null(x) || length(x) < i) return(default)
      x[[i]] %||% default
    }
    paper_blocks <- vapply(seq_len(n_papers), function(i) {
      ttl <- trimws(as.character(ctx_value("titles", i, "")))
      abs <- trimws(as.character(ctx_value("abstracts", i, "")))
      pmid <- trimws(as.character(ctx_value("pmids", i, "NA")))
      src_type <- trimws(as.character(ctx_value("source_type", i, "unknown")))
      src_term <- trimws(as.character(ctx_value("source_term", i, "NA")))
      # Truncate abstract to 800 chars to keep prompt manageable
      if (nchar(abs) > 800) abs <- paste0(substr(abs, 1, 800), "...")
      if (!nzchar(ttl) && !nzchar(abs)) return(NA_character_)
      sprintf("[%d] PMID: %s\n    Source: %s (%s)\n    Title: %s\n    Abstract excerpt: %s",
              i, pmid, src_type, src_term, ttl, abs)
    }, character(1))
    paper_blocks <- paper_blocks[!is.na(paper_blocks)]

    if (length(paper_blocks) > 0) {
      lit_text <- paste0(
        "Literature evidence (PubMed, deduplicated by PMID after module title and top metabolite searches):\n",
        paste(paper_blocks, collapse = "\n\n")
      )
    } else {
      lit_text <- "Literature evidence: no relevant abstracts retrieved."
    }
  } else {
    lit_text <- "Literature evidence: not retrieved (pubmed = FALSE or no papers found)."
  }

  # ---- Assemble prompt ----
  as.character(glue::glue("
Module title: {title_str}
Module type: {type_str}
Functional summary: {summary_str}

{origin_text}

{context_text}

{lit_text}
"))
}


#' Annotate a Single Module Using Two Sequential LLM Calls  (Fast Mode)
#'
#' Implements the Fast Mode two-call pipeline for one metabolic module:
#'
#' \enumerate{
#'   \item \strong{Call 1} — Biochemical annotation (no study context):
#'         generates \code{hypothesis_type}, \code{title}, and
#'         \code{functional_summary} from metabolite structure, class
#'         composition, and pathway coverage.
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
#'   \code{module_id}, \code{hypothesis_type}, \code{title}, \code{functional_summary},
#'   \code{phenotype_relationship}, \code{confidence_label}, \code{rationale},
#'   \code{global_direction},
#'   \code{origin_endogenous_frac}, \code{origin_exogenous_frac},
#'   \code{top_kegg_pathways}, \code{top_hmdb_pathways},
#'   \code{pubmed_pmids}, \code{pubmed_evidence_json},
#'   \code{raw_response_call1}, \code{raw_response_call2}.
annotate_one_module <- function(module_obj,
                                study_context,
                                provider,
                                llm_api_key,
                                model,
                                base_url     = NULL,
                                max_tokens   = 5000,
                                temperature  = 0,
                                pubmed       = FALSE,
                                NCBI_api_key = NULL,
                                years        = 5) {
  
  mid <- module_obj$module_id %||% NA_character_

  normalize_fast_hypothesis_type <- function(x) {
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
        base_url        = base_url,
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
        type               = list(
          type = "STRING",
          enum = c("Host-Endogenous", "Microbiome-Derived", "Diet-Derived",
                   "Xenobiotic / Exogenous", "Host–Microbiome")
        ),
        title              = list(type = "STRING"),
        functional_summary = list(type = "STRING"),
        anchor_metabolites = list(type = "ARRAY",
                                  items = list(type = "STRING"))
      ),
      required = c("type", "title", "functional_summary", "anchor_metabolites")
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
  # CALL 1 – Biochemical annotation (type + title + functional_summary)
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
    required_keys    = c("type", "title", "functional_summary"),
    retry_instruction = paste(
      "Return ONLY a valid JSON object with exactly four keys:",
      "\"type\", \"title\", \"functional_summary\", and \"anchor_metabolites\"."
    ),
    stage_label = "call1 (type+title+summary)"
  )

  hypothesis_type    <- NA_character_
  title              <- NA_character_
  functional_summary <- NA_character_
  anchor_metabolites <- NA_character_
  if (!is.null(res1$parsed)) {
    hypothesis_type    <- normalize_fast_hypothesis_type(res1$parsed[["type"]] %||% NA_character_)
    title              <- as.character(res1$parsed[["title"]]              %||% NA_character_)
    functional_summary <- as.character(res1$parsed[["functional_summary"]] %||% NA_character_)
    anc_raw            <- res1$parsed[["anchor_metabolites"]]
    if (!is.null(anc_raw) && length(anc_raw) > 0) {
      anchor_metabolites <- paste(as.character(unlist(anc_raw)), collapse = "; ")
    }
  } else {
    message("[annotate_module] Call 1 failed for module ", mid,
            ". Returning NA for title/summary; continuing to Call 2.")
  }
  
  # ==================================================================
  # PUBMED SEARCH (optional) – module title plus anchor metabolites
  # ==================================================================
  pubmed_context <- NULL

  if (pubmed && !is.na(title) && nzchar(trimws(title))) {
    # Use anchor metabolites from Call 1 for focused literature retrieval
    anchor_met_names <- NULL
    if (!is.na(anchor_metabolites) && nzchar(anchor_metabolites)) {
      anchor_met_names <- trimws(strsplit(anchor_metabolites, ";\\s*")[[1]])
      anchor_met_names <- anchor_met_names[nzchar(anchor_met_names)]
    }
    phenotype_val <- NULL
    if (!is.null(study_context$phenotype)) {
      ph <- trimws(as.character(study_context$phenotype[1]))
      if (nzchar(ph) && !identical(ph, "NA")) phenotype_val <- ph
    }

    pubmed_context <- tryCatch(
      search_pubmed_for_annotation_evidence(
        module_title         = title,
        top_metabolite_names = anchor_met_names,
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
    hypothesis_type = hypothesis_type,
    functional_summary = functional_summary,
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

  # Class names only
  class_names <- NA_character_
  if (!is.null(module_obj$other_metabolites) &&
      is.data.frame(module_obj$other_metabolites) &&
      nrow(module_obj$other_metabolites) > 0 &&
      "class" %in% colnames(module_obj$other_metabolites)) {
    cls <- unique(trimws(as.character(module_obj$other_metabolites$class)))
    cls <- cls[!is.na(cls) & nzchar(cls) & toupper(cls) != "NA"]
    if (length(cls) > 0) class_names <- paste(cls, collapse = "; ")
  }
  
  # ==================================================================
  # Return
  # ==================================================================
  data.frame(
    module_id              = mid,
    hypothesis_type        = hypothesis_type,
    title                  = title,
    functional_summary     = functional_summary,
    anchor_metabolites     = anchor_metabolites,
    phenotype_relationship = phenotype_relationship,
    confidence_label       = confidence_label,
    rationale              = rationale,
    class_names            = class_names,
    global_direction       = global_direction,
    origin_endogenous_frac = origin_endo_frac,
    origin_exogenous_frac  = origin_exo_frac,
    top_kegg_pathways      = kegg_pathways,
    top_hmdb_pathways      = hmdb_pathways,
    pubmed_pmids           = collapse_pubmed_pmids(pubmed_context),
    pubmed_evidence_json   = serialize_pubmed_evidence(pubmed_context),
    raw_response_call1     = res1$raw_text  %||% NA_character_,
    raw_response_call2     = res2$raw_text  %||% NA_character_,
    annotation_mode        = "fast",
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
                                 base_url     = NULL,
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
      base_url      = base_url,
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
