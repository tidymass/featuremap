# ============================================================
# Module Synthesis Agent  (advanced mode only)
#
# Requires annotate_module(mode = "advanced") output.
# Fast mode does not produce hypothesis-level debate verdicts and
# therefore cannot be synthesized at this level of granularity.
#
# Takes all accepted hypotheses from annotate_module() output
# and produces:
#   1. Thematic clusters  — grouping of hypotheses by shared biology
#   2. Sparse relationship table — only direct / indirect pairs;
#      pairs not listed are completed as unrelated (not_identified)
#   3. Overall summary narrative — dominant themes, cross-module
#      connections, and integrated phenotype interpretation
#
# Node unit: individual hypothesis (hypothesis_id within a module).
#   Fast-mode output  → node_id = module_id
#   Advanced-mode     → node_id = "module_id::hypothesis_id"
#   Same-module hypotheses share module_id; readers can see which
#   hypothesis is more central in the network.
#
# Relationship taxonomy:
#   direct   / shared_pathway      — same or tightly adjacent pathway
#   direct   / containment         — one is a sub-process of the other
#   direct   / cross_talk          — explicit substrate/product/regulatory link
#   indirect / parallel_under_theme — same higher-level theme, different routes
#   indirect / phenotype_convergent — different mechanisms, same phenotype target
#   (unrelated pairs are completed by R, not listed by the LLM)
# ============================================================


# ============================================================
# SECTION 1 – Input preparation
# ============================================================

#' Prepare node list from module_function data frame
#'
#' Filters to accepted rows and assigns stable node_ids.
#' Fast-mode rows (no hypothesis_id column, or hypothesis_id NA) use
#' module_id as node_id.  Advanced-mode rows with a real hypothesis_id
#' get a compound node_id "module_id::hypothesis_id".
#'
#' @param module_df  Data frame from \code{res$module_function}.
#' @return Data frame with columns: node_id, module_id, hypothesis_id,
#'   verdict_status, title, functional_summary, phenotype_relationship,
#'   confidence_label, rationale, core_drivers, pathway summaries, and
#'   phenotype literature support fields when available.
prepare_synthesis_nodes <- function(module_df) {
  if (is.null(module_df) || nrow(module_df) == 0)
    stop("module_df is empty — run annotate_module(mode = \"advanced\") first.")

  # Require advanced mode output
  if (!"verdict_status" %in% colnames(module_df))
    stop(paste(
      "module_df does not contain a 'verdict_status' column.",
      "synthesize_modules() requires advanced mode output",
      "(annotate_module(mode = \"advanced\"))."
    ))

  # Exclude fast_fallback rows (P0 failed, reverted to fast mode internally)
  if ("annotation_mode" %in% colnames(module_df)) {
    module_df <- module_df[module_df$annotation_mode == "advanced", , drop = FALSE]
    if (nrow(module_df) == 0)
      stop("No advanced-mode rows found (all modules were fast_fallback).")
  }

  # Keep accepted rows only
  df <- module_df[module_df$verdict_status == "accepted", , drop = FALSE]

  if (nrow(df) == 0)
    stop("No accepted hypotheses found. Nothing to synthesize.")

  # node_id = "module_id::hypothesis_id" (hypothesis_id is always present in advanced mode)
  df$node_id <- paste(df$module_id, df$hypothesis_id, sep = "::")

  # Deduplicate (guard against accidental duplicates)
  df <- df[!duplicated(df$node_id), , drop = FALSE]

  optional_cols <- c(
    "hypothesis_type", "confidence_label", "rationale", "core_drivers",
    "top_kegg_pathways", "top_hmdb_pathways",
    "phenotype_literature_support", "phenotype_evidence_note"
  )
  for (nm in optional_cols) {
    if (!nm %in% colnames(df)) df[[nm]] <- NA_character_
  }

  df[, c("node_id", "module_id",
         "hypothesis_id", "verdict_status",
         "title", "functional_summary", "phenotype_relationship",
         optional_cols)]
}


# ============================================================
# SECTION 2 – Prompt builders
# ============================================================

build_synthesis_system_prompt <- function() {
  glue::glue(
    "You are an expert in metabolomics and systems biology.
Your task is to synthesize a set of annotated metabolic module hypotheses from one omics study
into an integrated biological picture.

You will receive:
- A list of hypothesis nodes, each with: node_id, module_id, verdict_status,
  hypothesis type, title, functional_summary, phenotype_relationship, confidence/rationale,
  core driver metabolites, pathway evidence, and phenotype literature support
- Study context: species, biospecimen, phenotype

Your output has three parts:

────────────────────────────────────────
PART 1 — CLUSTERS
────────────────────────────────────────
Group nodes into thematic clusters based on shared or closely related biology.
A cluster must have at least 2 nodes; singleton nodes belong to their own cluster only
if they represent a clearly distinct biological theme.

For each cluster provide:
- cluster_id    : \"C1\", \"C2\", ...
- cluster_theme : concise biological theme name (e.g., \"Aromatic amino acid catabolism\")
- node_ids      : list of node_ids in this cluster

────────────────────────────────────────
PART 2 — RELATIONSHIPS (cross-module only, sparse)
────────────────────────────────────────
List node pairs that have a direct or indirect biological relationship.

Relationship taxonomy (use exactly one subtype):
  direct   | shared_pathway      : same or tightly adjacent pathway
  direct   | containment         : one node is a sub-process of the other
  direct   | cross_talk          : explicit substrate/product/regulatory biochemical link
  indirect | parallel_under_theme: different routes under the same biological theme
  indirect | phenotype_convergent: different mechanisms converging on the same phenotype process

For each relationship provide:
- node_1         : node_id
- node_2         : node_id
- relation_class : \"direct\" | \"indirect\"
- subtype        : one of the five subtypes above
- rationale      : 1-2 sentences citing specific metabolites, pathways, or phenotype evidence.
                   Do not claim metabolite overlap unless the exact same metabolite names are
                   visible in both node records.

────────────────────────────────────────
PART 3 — SUMMARY
────────────────────────────────────────
Write a 4-7 sentence narrative integrating the findings. Structure it as:
  Sentence 1-2 : Dominant biological themes (reference cluster names)
  Sentence 3-4 : Key cross-module connections (reference specific node titles and relations)
  Sentence 5-6 : Integrated phenotype interpretation (tie the modules to the study phenotype)
  Sentence 7   : (optional) Caveats or unresolved questions across modules

RULES:
- Base all statements on the provided node information only
- Do NOT invent metabolites, pathways, or literature not present in the inputs
- Treat modules as distinct community-detected modules. Do NOT describe two modules as
  having \"overlapping\" or \"shared\" metabolites unless the exact same metabolite
  name is explicitly present in both nodes' Core drivers or pathway evidence.
- For shared_pathway relationships, cite shared or adjacent pathway evidence,
  biochemical class, or process similarity. Same pathway evidence is not the same
  as overlapping metabolites.
- Use core drivers and pathway evidence when deciding shared pathways, containment, or cross-talk
- Use phenotype relationship and literature support when deciding phenotype_convergent relationships
- Do NOT rank accepted nodes by verdict status in the summary (treat all accepted nodes equally)

Return STRICT JSON only:
{{
  \"clusters\": [
    {{
      \"cluster_id\": \"C1\",
      \"cluster_theme\": \"...\",
      \"node_ids\": [\"...\"]
    }}
  ],
  \"relationships\": [
    {{
      \"node_1\": \"...\",
      \"node_2\": \"...\",
      \"relation_class\": \"direct | indirect\",
      \"subtype\": \"...\",
      \"rationale\": \"...\"
    }}
  ],
  \"summary\": \"...\"
}}"
  )
}


build_synthesis_user_prompt <- function(nodes_df, study_context = NULL) {

  scalar <- function(x, na = "NA") {
    if (is.null(x) || length(x) == 0) return(na)
    x <- trimws(as.character(x[1]))
    if (!nzchar(x) || is.na(x)) na else x
  }

  # Study context block
  if (!is.null(study_context) &&
      any(vapply(study_context[c("species", "biospecimen", "phenotype")],
                 function(x) !is.null(x) && nzchar(trimws(as.character(x[1] %||% ""))),
                 logical(1)))) {
    ctx_text <- glue::glue(
      "Study context:
  Species:     {scalar(study_context$species)}
  Biospecimen: {scalar(study_context$biospecimen)}
  Phenotype:   {scalar(study_context$phenotype)}"
    )
  } else {
    ctx_text <- "Study context: not provided."
  }

  # Node blocks — group same-module nodes together
  nodes_df <- nodes_df[order(nodes_df$module_id, nodes_df$node_id), ]
  modules  <- unique(nodes_df$module_id)

  node_blocks <- unlist(lapply(modules, function(mid) {
    rows <- nodes_df[nodes_df$module_id == mid, , drop = FALSE]
    lapply(seq_len(nrow(rows)), function(i) {
      r <- rows[i, ]
      fs <- scalar(r$functional_summary)
      pr <- scalar(r$phenotype_relationship)
      rat <- scalar(r$rationale)
      # Truncate to keep prompt manageable
      if (nchar(fs) > 500) fs <- paste0(substr(fs, 1, 500), "...")
      if (nchar(pr) > 300) pr <- paste0(substr(pr, 1, 300), "...")
      if (nchar(rat) > 300) rat <- paste0(substr(rat, 1, 300), "...")
      glue::glue(
        "Node: {r$node_id}
  Module: {r$module_id}  |  Status: {scalar(r$verdict_status)}
  Type: {scalar(r$hypothesis_type)}
  Title: {scalar(r$title)}
  Confidence: {scalar(r$confidence_label)}
  Confidence rationale: {rat}
  Core drivers: {scalar(r$core_drivers)}
  Top KEGG pathways: {scalar(r$top_kegg_pathways)}
  Top HMDB pathways: {scalar(r$top_hmdb_pathways)}
  Functional summary: {fs}
  Phenotype relationship: {pr}
  Phenotype literature support: {scalar(r$phenotype_literature_support)}
  Phenotype evidence note: {scalar(r$phenotype_evidence_note)}"
      )
    })
  }), recursive = FALSE)

  node_text <- paste(node_blocks, collapse = "\n\n")

  n_nodes   <- nrow(nodes_df)
  n_modules <- length(unique(nodes_df$module_id))

  glue::glue(
    "{ctx_text}

Total nodes: {n_nodes} hypothesis/module entries across {n_modules} modules.

--- Hypothesis nodes ---

{node_text}"
  )
}


# ============================================================
# SECTION 3 – Gemini response schema
# ============================================================

synth_gemini_schema <- function() {
  list(
    type = "OBJECT",
    properties = list(
      clusters = list(
        type = "ARRAY",
        items = list(
          type = "OBJECT",
          properties = list(
            cluster_id    = list(type = "STRING"),
            cluster_theme = list(type = "STRING"),
            node_ids      = list(type = "ARRAY", items = list(type = "STRING"))
          ),
          required = c("cluster_id", "cluster_theme", "node_ids")
        )
      ),
      relationships = list(
        type = "ARRAY",
        items = list(
          type = "OBJECT",
          properties = list(
            node_1         = list(type = "STRING"),
            node_2         = list(type = "STRING"),
            relation_class = list(type = "STRING",
                                  enum = c("direct", "indirect")),
            subtype        = list(type = "STRING",
                                  enum = c("shared_pathway", "containment",
                                           "cross_talk", "parallel_under_theme",
                                           "phenotype_convergent")),
            rationale      = list(type = "STRING")
          ),
          required = c("node_1", "node_2", "relation_class", "subtype", "rationale")
        )
      ),
      summary = list(type = "STRING")
    ),
    required = c("clusters", "relationships", "summary")
  )
}


# ============================================================
# SECTION 4a – Relationship table completion
# ============================================================

# After the LLM returns sparse cross-module direct/indirect pairs, this
# function completes the full N*(N-1)/2 pairwise table by:
#
#   1. same_module pairs  → relation_class = "same_module"
#                           subtype        = "same_module_hypotheses"
#                           rationale      = (auto-generated note)
#
#   2. LLM direct/indirect pairs → kept as-is (cross-module only)
#
#   3. Remaining cross-module pairs → relation_class = "unrelated"
#                                     subtype        = "not_identified"
#                                     rationale      = NA
#
# Returns a complete data frame with one row per unique pair.
complete_relationship_table <- function(llm_relationships, nodes_df) {

  node_ids <- nodes_df$node_id
  n        <- length(node_ids)

  empty_complete <- data.frame(
    node_1 = character(), node_2 = character(),
    module_id_1 = character(), module_id_2 = character(),
    same_module = logical(), relation_class = character(),
    subtype = character(), rationale = character(),
    stringsAsFactors = FALSE
  )
  if (n < 2) return(empty_complete)

  # All unique pairs (lexicographic order)
  pairs <- combn(sort(node_ids), 2, simplify = FALSE)

  # Module lookup
  mid_lookup <- setNames(nodes_df$module_id, nodes_df$node_id)

  # LLM relationship lookup: "node_1|||node_2" → row index
  if (!is.null(llm_relationships) && nrow(llm_relationships) > 0) {
    llm_key <- paste(llm_relationships$node_1, llm_relationships$node_2, sep = "|||")
  } else {
    llm_key <- character(0)
  }

  rows <- lapply(pairs, function(p) {
    n1 <- p[1]; n2 <- p[2]
    key <- paste(n1, n2, sep = "|||")

    m1 <- unname(mid_lookup[n1] %||% NA_character_)
    m2 <- unname(mid_lookup[n2] %||% NA_character_)

    if (!is.na(m1) && !is.na(m2) && identical(m1, m2)) {
      # Same module — hypotheses from the same debate; may be complementary or alternative.
      data.frame(
        node_1         = n1,
        node_2         = n2,
        module_id_1    = m1,
        module_id_2    = m2,
        same_module    = TRUE,
        relation_class = "same_module",
        subtype        = "same_module_hypotheses",
        rationale      = paste0("Both nodes are hypotheses for module ", m1,
                                " generated in the same module-level debate."),
        stringsAsFactors = FALSE
      )
    } else if (key %in% llm_key) {
      # LLM-labeled direct or indirect relationship
      r <- llm_relationships[llm_key == key, , drop = FALSE][1, ]
      data.frame(
        node_1         = n1,
        node_2         = n2,
        module_id_1    = m1,
        module_id_2    = m2,
        same_module    = FALSE,
        relation_class = r$relation_class,
        subtype        = r$subtype,
        rationale      = r$rationale,
        stringsAsFactors = FALSE
      )
    } else {
      # Cross-module, not mentioned by LLM → unrelated
      data.frame(
        node_1         = n1,
        node_2         = n2,
        module_id_1    = m1,
        module_id_2    = m2,
        same_module    = FALSE,
        relation_class = "unrelated",
        subtype        = "not_identified",
        rationale      = "No direct or indirect relationship was identified by the synthesis agent.",
        stringsAsFactors = FALSE
      )
    }
  })

  do.call(rbind, rows)
}


# ============================================================
# SECTION 4b – Output parsing & assembly
# ============================================================

#' Parse synthesis LLM output into tidy data frames
#'
#' @param parsed   Parsed list from \code{adv_parse_json_response()}.
#' @param nodes_df Node data frame from \code{prepare_synthesis_nodes()}.
#' @return List with \code{$clusters}, \code{$relationships_sparse},
#'   \code{$relationships}, \code{$summary}, \code{$nodes}.
parse_synthesis_output <- function(parsed, nodes_df) {

  valid_nodes <- nodes_df$node_id

  # ---- clusters ----
  clusters_df <- tryCatch({
    cl <- parsed$clusters
    if (is.null(cl) || length(cl) == 0) stop("empty")
    do.call(rbind, lapply(cl, function(c) {
      nids <- trimws(as.character(unlist(c$node_ids %||% list())))
      nids <- nids[nids %in% valid_nodes]
      data.frame(
        cluster_id    = trimws(as.character(c$cluster_id    %||% "")),
        cluster_theme = trimws(as.character(c$cluster_theme %||% "")),
        node_ids      = paste(nids, collapse = "; "),
        stringsAsFactors = FALSE
      )
    }))
  }, error = function(e) {
    data.frame(cluster_id = NA_character_, cluster_theme = NA_character_,
               node_ids = NA_character_, stringsAsFactors = FALSE)[0, ]
  })
  if (nrow(clusters_df) > 0) {
    clusters_df <- clusters_df[nzchar(clusters_df$cluster_id) & nzchar(clusters_df$node_ids), , drop = FALSE]
  }

  # ---- relationships ----
  relationships_df <- tryCatch({
    rl <- parsed$relationships
    if (is.null(rl) || length(rl) == 0) {
      data.frame(node_1 = character(), node_2 = character(),
                 relation_class = character(), subtype = character(),
                 rationale = character(), stringsAsFactors = FALSE)
    } else {
      do.call(rbind, lapply(rl, function(r) {
        n1 <- trimws(as.character(r$node_1 %||% ""))
        n2 <- trimws(as.character(r$node_2 %||% ""))
        # Enforce lexicographic ordering to avoid duplicate pairs
        if (n1 > n2) { tmp <- n1; n1 <- n2; n2 <- tmp }
        data.frame(
          node_1         = n1,
          node_2         = n2,
          relation_class = trimws(tolower(as.character(r$relation_class %||% ""))),
          subtype        = trimws(as.character(r$subtype        %||% "")),
          rationale      = trimws(as.character(r$rationale      %||% "")),
          stringsAsFactors = FALSE
        )
      }))
    }
  }, error = function(e) {
    data.frame(node_1 = character(), node_2 = character(),
               relation_class = character(), subtype = character(),
               rationale = character(), stringsAsFactors = FALSE)
  })

  # Deduplicate pairs (same node_1/node_2)
  if (nrow(relationships_df) > 0) {
    allowed_subtypes <- c(
      "shared_pathway", "containment", "cross_talk",
      "parallel_under_theme", "phenotype_convergent"
    )
    relationships_df <- relationships_df[
      relationships_df$node_1 %in% valid_nodes &
        relationships_df$node_2 %in% valid_nodes &
        relationships_df$node_1 != relationships_df$node_2 &
        relationships_df$relation_class %in% c("direct", "indirect") &
        relationships_df$subtype %in% allowed_subtypes,
      ,
      drop = FALSE
    ]

    mid_lookup <- setNames(nodes_df$module_id, nodes_df$node_id)
    module_1 <- unname(mid_lookup[relationships_df$node_1])
    module_2 <- unname(mid_lookup[relationships_df$node_2])
    relationships_df <- relationships_df[
      !is.na(module_1) & !is.na(module_2) & module_1 != module_2,
      ,
      drop = FALSE
    ]
  }

  if (nrow(relationships_df) > 0) {
    pair_key <- paste(relationships_df$node_1, relationships_df$node_2, sep = "|||")
    relationships_df <- relationships_df[!duplicated(pair_key), , drop = FALSE]
  }

  # Enrich sparse LLM relationships with module_id info for each node
  mid_lookup <- setNames(nodes_df$module_id, nodes_df$node_id)
  if (nrow(relationships_df) > 0) {
    relationships_df$module_id_1 <- unname(mid_lookup[relationships_df$node_1])
    relationships_df$module_id_2 <- unname(mid_lookup[relationships_df$node_2])
    relationships_df$same_module  <- relationships_df$module_id_1 == relationships_df$module_id_2 &
      !is.na(relationships_df$module_id_1)
  } else {
    relationships_df$module_id_1  <- character()
    relationships_df$module_id_2  <- character()
    relationships_df$same_module   <- logical()
  }

  relationships_sparse <- relationships_df
  relationships_complete <- complete_relationship_table(relationships_sparse, nodes_df)

  # Enrich nodes with cluster assignment
  nodes_out <- nodes_df
  nodes_out$cluster_id <- NA_character_
  for (i in seq_len(nrow(clusters_df))) {
    nids <- trimws(strsplit(clusters_df$node_ids[i], "; ", fixed = TRUE)[[1]])
    nodes_out$cluster_id[nodes_out$node_id %in% nids] <- clusters_df$cluster_id[i]
  }

  list(
    nodes                = nodes_out,
    clusters             = clusters_df,
    relationships_sparse = relationships_sparse,
    relationships        = relationships_complete,
    summary              = trimws(as.character(parsed$summary %||% ""))
  )
}


# ============================================================
# SECTION 5 – Main orchestrator
# ============================================================

#' Synthesize All Annotated Module Hypotheses
#'
#' Takes the \code{module_function} data frame from an annotated featuremap
#' result and runs a single LLM call to produce:
#' \itemize{
#'   \item Thematic clusters of hypotheses
#'   \item A sparse pairwise relationship table (direct/indirect only)
#'   \item An overall narrative summary
#' }
#'
#' @param res            featuremap result object after
#'   \code{annotate_module(mode = "advanced")}.
#'   \strong{Not compatible with fast mode output.}
#' @param study_context  Named list with \code{species}, \code{biospecimen},
#'   \code{phenotype}.  If \code{NULL}, extracted from \code{res} if available.
#' @param provider       LLM provider: \code{"openai"} or \code{"gemini"}.
#' @param llm_api_key    API key string.
#' @param llm_model      Model name (e.g. \code{"gpt-4.1-mini"}).
#' @param max_tokens     Max tokens for the LLM response (default 4000).
#' @param temperature    Sampling temperature (default 0).
#' @param sleep_sec      Seconds to pause after the API call (default 0).
#' @return The input \code{res} object with a \code{module_synthesis} list
#'   appended, containing \code{$nodes}, \code{$clusters},
#'   \code{$relationships_sparse}, \code{$relationships},
#'   \code{$summary}, and \code{$raw_response}.
#' @export
synthesize_modules <- function(res,
                                study_context  = NULL,
                                provider       = "openai",
                                llm_api_key    = "",
                                llm_model      = "gpt-4.1-mini",
                                base_url       = NULL,
                                max_tokens     = 4000L,
                                temperature    = 0,
                                sleep_sec      = 0) {

  # ---- Validate input ----
  if (is.null(res$module_function) || nrow(res$module_function) == 0)
    stop("res$module_function is empty. Run annotate_module() first.")

  # synthesize_modules() is only defined for advanced mode output.
  # Fast mode produces one flat annotation per module with no hypothesis debate,
  # so cross-module synthesis at the hypothesis level is not applicable.
  if (!"verdict_status" %in% colnames(res$module_function)) {
    stop(paste(
      "synthesize_modules() requires mode = \"advanced\" output.",
      "Fast mode (annotate_module(mode = \"fast\")) does not produce",
      "hypothesis-level verdicts and cannot be synthesized this way.",
      "Re-run annotate_module() with mode = \"advanced\" first."
    ))
  }

  # ---- Extract study context from res if not supplied ----
  if (is.null(study_context) && !is.null(res$study_context))
    study_context <- res$study_context

  # ---- Prepare nodes ----
  module_function_for_synthesis <- res$module_function
  if (!is.null(res$module_pubmed_evidence) &&
      is.data.frame(res$module_pubmed_evidence) &&
      nrow(res$module_pubmed_evidence) > 0) {
    pub_cols <- c(
      "module_id", "hypothesis_id",
      "phenotype_literature_support", "phenotype_evidence_note"
    )
    pub_cols <- pub_cols[pub_cols %in% colnames(res$module_pubmed_evidence)]
    join_by <- intersect(c("module_id", "hypothesis_id"), pub_cols)
    if (length(join_by) > 0 &&
        any(c("phenotype_literature_support", "phenotype_evidence_note") %in% pub_cols)) {
      module_function_for_synthesis <- dplyr::left_join(
        module_function_for_synthesis,
        dplyr::distinct(res$module_pubmed_evidence[, pub_cols, drop = FALSE]),
        by = join_by
      )
    }
  }

  nodes_df <- tryCatch(
    prepare_synthesis_nodes(module_function_for_synthesis),
    error = function(e) {
      cli::cli_abort("synthesize_modules: failed to prepare nodes — {e$message}")
    }
  )

  n_nodes   <- nrow(nodes_df)
  n_modules <- length(unique(nodes_df$module_id))

  if (n_modules <= 1) {
    cli::cli_alert_info(
      "Only one module has accepted hypotheses; skipping cross-module synthesis."
    )
    res$module_synthesis <- list(
      nodes                = nodes_df,
      clusters             = data.frame(
        cluster_id = character(),
        cluster_theme = character(),
        node_ids = character(),
        stringsAsFactors = FALSE
      ),
      relationships_sparse = data.frame(
        node_1 = character(),
        node_2 = character(),
        relation_class = character(),
        subtype = character(),
        rationale = character(),
        module_id_1 = character(),
        module_id_2 = character(),
        same_module = logical(),
        stringsAsFactors = FALSE
      ),
      relationships        = complete_relationship_table(NULL, nodes_df),
      summary              = NA_character_,
      raw_response         = NA_character_
    )
    return(res)
  }

  cli::cli_alert_info(
    "Synthesizing {.val {n_nodes}} hypothesis node(s) across {.val {n_modules}} module(s)."
  )

  # ---- Build prompts ----
  sys_prompt  <- build_synthesis_system_prompt()
  user_prompt <- build_synthesis_user_prompt(nodes_df, study_context)

  is_gemini <- identical(provider, "gemini")
  schema    <- if (is_gemini) synth_gemini_schema() else NULL

  # ---- LLM call ----
  result <- adv_request_json(
    messages      = list(list(role = "system", content = sys_prompt),
                         list(role = "user",   content = user_prompt)),
    required_keys = c("clusters", "relationships", "summary"),
    schema        = schema,
    retry_instr   = paste(
      "Return only valid JSON with keys \"clusters\" (array),",
      "\"relationships\" (array), and \"summary\" (string)."
    ),
    provider    = provider,
    llm_api_key = llm_api_key,
    model       = llm_model,
    base_url    = base_url,
    max_tokens  = max_tokens,
    temperature = temperature
  )

  if (sleep_sec > 0) Sys.sleep(sleep_sec)

  # ---- Parse output ----
  synthesis <- if (!is.null(result$parsed)) {
    tryCatch(
      parse_synthesis_output(result$parsed, nodes_df),
      error = function(e) {
        cli::cli_warn("synthesize_modules: output parsing failed — {e$message}")
        list(nodes = nodes_df, clusters = NULL, relationships_sparse = NULL,
             relationships = complete_relationship_table(NULL, nodes_df),
             summary = NA_character_)
      }
    )
  } else {
    cli::cli_warn("synthesize_modules: LLM returned no parseable output.")
    list(nodes = nodes_df, clusters = NULL, relationships_sparse = NULL,
         relationships = complete_relationship_table(NULL, nodes_df),
         summary = NA_character_)
  }

  synthesis$raw_response <- result$raw_text %||% NA_character_

  n_rel_sparse <- if (!is.null(synthesis$relationships_sparse)) nrow(synthesis$relationships_sparse) else 0L
  n_rel_full   <- if (!is.null(synthesis$relationships))        nrow(synthesis$relationships)        else 0L
  n_cl         <- if (!is.null(synthesis$clusters))             nrow(synthesis$clusters)             else 0L
  cli::cli_alert_success(
    "Synthesis complete: {.val {n_cl}} cluster(s), {.val {n_rel_sparse}} LLM relationship(s), {.val {n_rel_full}} completed pair(s)."
  )

  res$module_synthesis <- synthesis
  return(res)
}
