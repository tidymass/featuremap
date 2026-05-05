
# ============================================================
# LLM Prompt Builders  (Fast Mode)
#
# Call 1:  build_system_prompt()  +  build_module_prompt()
#          → {title, functional_summary}
#          Inputs: Signal_concentration, top metabolites, class structure,
#                  pathway coverage (KEGG + HMDB)
#          No study context; no confidence assessment
#
# Call 2:  build_phenotype_confidence_system_prompt()
#          + build_phenotype_confidence_prompt()
#          → {phenotype_relationship, confidence_label, rationale}
#          Inputs: module title, origin profile, study context,
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

Using only the provided information, return two fields:

1) TITLE
   - Name the pathway in KEGG style: concise, standardized, and anchored to a specific biochemical pathway or process.
   - Follow KEGG naming conventions:
     - Use \"X metabolism\" when the module spans multiple steps of a compound class (e.g., \"Tyrosine metabolism\", \"Sphingolipid metabolism\").
     - Use \"X biosynthesis\" or \"X degradation\" when directionality is clear from the metabolites.
     - Use a specific reaction or enzyme-level name when anchored by a single well-defined step.
     - Use a physiological or system-level name when appropriate.
   - Choose the most specific level of granularity supported by the early-ranked high-priority metabolites.
   - If top metabolites point to a specific sub-pathway, reflect that specificity rather than defaulting to the parent pathway name.
   - Use pathway evidence (KEGG/HMDB) to support or refine the title when it aligns with the metabolite evidence.
   - If no canonical pathway fits, use the most specific neutral description of the dominant chemical transformation.

2) FUNCTIONAL_SUMMARY
   - Write 4-6 sentences as a factual biochemical annotation of the module.
   - Describe what the metabolites do biochemically: enzymatic conversions, pathway steps, metabolite roles, and physiological functions — as if writing an entry in a biochemical database.
   - Base the explanation mainly on early-ranked high-priority metabolites; use class structure and pathway evidence only to support or qualify.
   - STRICTLY FORBIDDEN: do NOT state or imply whether metabolites are endogenous or exogenous, dietary, microbial, plant-derived, or environmental.
   - Do NOT comment on the study phenotype or disease relevance.
   - Do NOT reference module-level statistics (Signal_concentration).
   - Mention secondary signals only briefly at the end.

Return STRICT JSON only:
{{
  \"title\": \"...\",
  \"functional_summary\": \"...\"
}}

Do not invent unsupported pathways or mechanisms."
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
    if (!nzchar(x) || is.na(x)) na else x
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
- Class: {scalar_text(row$class)}
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
Your task is to assess the phenotypic relevance and confidence of a previously named metabolic module.

You will receive:
- Module title (from prior biochemical analysis)
- Origin profile: proportions of endogenous, exogenous, and unknown-origin metabolites
  (based on HMDB-annotated module members only)
- Study context: species, biospecimen, phenotype
- Literature evidence: PubMed abstracts retrieved by searching the module title,
  all top metabolites, and phenotype
  (may be absent if pubmed = FALSE or no papers were found)

Definitions:
- Endogenous: metabolites documented as produced by the host organism.
- Exogenous: metabolites of dietary, microbial, pharmaceutical, plant, or environmental origin.
- Unknown: HMDB-annotated metabolites with no recorded origin.

Return three fields:

1) PHENOTYPE_RELATIONSHIP
   - Write 2-4 sentences summarizing how the module's metabolites relate to the study phenotype,
     based ONLY on the provided PubMed evidence.
   - If no PubMed evidence is available, write exactly:
     \"Insufficient literature evidence to establish phenotype relationship.\"
   - Do NOT speculate beyond what is stated in the abstracts.

2) CONFIDENCE
   Evaluate across three dimensions using the study context AND origin profile:
   (a) Species validity: Is this an established pathway in the study species?
   (b) Biospecimen appropriateness: Is this pathway known to be active/detectable in the study
       biospecimen?
   (c) Phenotype relevance: Is the pathway plausibly perturbed in the study phenotype?

   Scoring rules:
   - High: ALL three dimensions hold AND early-ranked metabolites are recognized endogenous
     substrates/products in that species. If exogenous_frac > 0.5, confidence CANNOT be high.
   - Medium: pathway exists in the species and fits the biospecimen, but at least one dimension
     is indirect/speculative, OR exogenous metabolites dominate.
   - Low: pathway fails species validity OR the module is overwhelmingly exogenous/synthetic.

   Origin profile is a primary signal: high exogenous fraction shifts confidence toward medium
   or low regardless of pathway plausibility.

3) RATIONALE
   - 2-3 sentences explicitly addressing species validity, biospecimen appropriateness,
     phenotype relevance, and how the origin profile influenced the confidence decision.

Return STRICT JSON only:
{{
  \"phenotype_relationship\": \"...\",
  \"confidence_label\": \"high\" | \"medium\" | \"low\",
  \"rationale\": \"2-3 sentences.\"
}}

Do not invent information not present in the provided evidence."
  )
}


# ============================================================
# CALL 2 – module user prompt
# ============================================================

build_phenotype_confidence_prompt <- function(module_obj,
                                               title,
                                               study_context  = NULL,
                                               pubmed_context = NULL) {

  scalar_text <- function(x, na = "NA") {
    if (is.null(x) || length(x) == 0) return(na)
    x <- as.character(x[1])
    x <- trimws(x)
    if (!nzchar(x) || is.na(x)) na else x
  }

  title_str <- scalar_text(title, na = "Unknown module")

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
      # Truncate abstract to 400 chars to keep prompt manageable
      if (nchar(abs) > 400) abs <- paste0(substr(abs, 1, 400), "...")
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

{origin_text}

{context_text}

{lit_text}
"))
}
