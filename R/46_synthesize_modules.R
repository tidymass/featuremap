# ============================================================
# Module Synthesis Agent  (advanced mode only)
#
# Requires annotate_module(mode = "advanced") output.
# Fast mode does not produce hypothesis-level debate verdicts and
# therefore cannot be synthesized at this level of granularity.
#
# Takes all accepted hypotheses from annotate_module() output
# and produces:
#   1. Thematic clusters  \u2014 grouping of hypotheses by shared biology
#   2. Evidence-backed relationship tables:
#      - biochemical links are computed deterministically from direct MRN
#        top-metabolite edges plus metabolic pathway ID overlap
#      - biology links are inferred by the synthesis LLM from title/summary
#      - final relation_type is biochemical / biology / both / none
#   3. Overall summary narrative \u2014 dominant themes, cross-module
#      connections, and integrated phenotype interpretation
#
# Node unit: individual hypothesis (hypothesis_id within a module).
#   Fast-mode output  \u2192 node_id = module_id
#   Advanced-mode     \u2192 node_id = "module_id::hypothesis_id"
#   Same-module hypotheses share module_id; readers can see which
#   hypothesis is more central in the network.
#
# Relationship taxonomy:
#   biochemical \u2014 direct top-metabolite edge in MRN and metabolic pathway ID overlap
#   biology     \u2014 LLM-inferred functional relationship from hypothesis titles/summaries
#   both        \u2014 both evidence sources support the pair
#   none        \u2014 no biochemical or biology link identified
# ============================================================


# ============================================================
# SECTION 1 \u2013 Input preparation
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
#'   confidence_label, rationale, anchor/top metabolites, pathway summaries,
#'   and phenotype literature support fields when available.
prepare_synthesis_nodes <- function(module_df) {
  if (is.null(module_df) || nrow(module_df) == 0)
    stop("module_df is empty \u2014 run annotate_module(mode = \"advanced\") first.")

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
    "hypothesis_type", "confidence_label", "rationale", "anchor_metabolites",
    "anchor_metabolite_ids", "top_metabolites", "top_metabolite_ids",
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
# SECTION 2 \u2013 Prompt builders
# ============================================================

build_synthesis_system_prompt <- function() {
  glue::glue(
    "You are an expert in metabolomics and systems biology.
Your task is to synthesize a set of annotated metabolic module hypotheses from one omics study
into an integrated biological picture.

You will receive:
- A list of accepted hypothesis nodes, each with: node_id, module_id,
  hypothesis type, title, and functional_summary
- Study context: species, biospecimen, phenotype

Your output has three parts:

\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
PART 1 \u2014 CLUSTERS
\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
Group nodes into thematic clusters based on shared or closely related biology.
A cluster must have at least 2 nodes; singleton nodes belong to their own cluster only
if they represent a clearly distinct biological theme.
Each node may appear in at most one cluster. Choose the single most likely cluster for
ambiguous nodes; do not duplicate a node across clusters.

For each cluster provide:
- cluster_id    : \"C1\", \"C2\", ...
- cluster_theme : concise biological theme name (e.g., \"Aromatic amino acid catabolism\")
- node_ids      : list of node_ids in this cluster

\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
PART 2 \u2014 BIOLOGY LINKS (cross-module only, sparse)
\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
List cross-module node pairs that have a clear functional or biological relationship
based only on their titles and functional summaries.

For each relationship provide:
- node_1         : node_id
- node_2         : node_id
- rationale      : 1-2 sentences explaining the shared functional theme or biological
                   dependency implied by the two titles/summaries.

\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
PART 3 \u2014 SUMMARY
\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
Write a 4-7 sentence narrative integrating the findings. Structure it as:
  Sentence 1-2 : Dominant biological themes (reference cluster names)
  Sentence 3-4 : Key cross-module biology links (reference specific node titles)
  Sentence 5-6 : Integrated phenotype interpretation (tie the modules to the study phenotype)
  Sentence 7   : (optional) Caveats or unresolved questions across modules

RULES:
- Base all statements on the provided node information only
- Do NOT infer biochemical adjacency, pathway overlap, substrate/product edges, or MRN links;
  those are computed separately by R code.
- Do NOT invent metabolites, pathways, or literature.
- Do NOT rank accepted nodes by verdict status; treat all accepted nodes equally.
- Only list a biology link when the titles/summaries imply a meaningful functional connection.
- Be conservative. Do not link nodes merely because they share a broad umbrella label
  such as lipid metabolism, amino acid metabolism, xenobiotic metabolism, or host-endogenous
  metabolism.
- Prefer no biology link unless the pair shares a specific process, substrate family,
  pathway-level function, detoxification mechanism, or clear upstream/downstream biological
  dependency stated in the titles/summaries.

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

  # Node blocks \u2014 group same-module nodes together.  Biology-link inference
  # intentionally sees only title/summary-level semantics; deterministic
  # biochemical evidence is computed separately in R.
  nodes_df <- nodes_df[order(nodes_df$module_id, nodes_df$node_id), ]
  modules  <- unique(nodes_df$module_id)

  node_blocks <- unlist(lapply(modules, function(mid) {
    rows <- nodes_df[nodes_df$module_id == mid, , drop = FALSE]
    lapply(seq_len(nrow(rows)), function(i) {
      r <- rows[i, ]
      fs <- scalar(r$functional_summary)
      # Truncate to keep prompt manageable
      if (nchar(fs) > 500) fs <- paste0(substr(fs, 1, 500), "...")
      glue::glue(
        "Node: {r$node_id}
  Module: {r$module_id}
  Type: {scalar(r$hypothesis_type)}
  Title: {scalar(r$title)}
  Functional summary: {fs}"
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
# SECTION 3 \u2013 Gemini response schema
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
            rationale      = list(type = "STRING")
          ),
          required = c("node_1", "node_2", "rationale")
        )
      ),
      summary = list(type = "STRING")
    ),
    required = c("clusters", "relationships", "summary")
  )
}


# ============================================================
# SECTION 4a \u2013 Deterministic biochemical evidence
# ============================================================

.syn_split_field <- function(x) {
  if (is.null(x) || length(x) == 0L) return(character(0))
  x <- as.character(x[1])
  if (is.na(x) || !nzchar(trimws(x)) || toupper(trimws(x)) %in% c("NA", "NULL")) {
    return(character(0))
  }
  # InChI strings may contain semicolons in charge/protonation layers, so split
  # only on the package's human-readable delimiter ("; ").
  out <- trimws(strsplit(x, ";\\s+|\\{\\}")[[1]])
  out[!is.na(out) & nzchar(out) & toupper(out) != "NA"]
}

.syn_pair_key <- function(a, b) {
  paste(pmin(a, b), pmax(a, b), sep = "|||")
}

.syn_empty_biochemical_links <- function() {
  data.frame(
    node_1 = character(), node_2 = character(),
    module_id_1 = character(), module_id_2 = character(),
    biochemical_mrn_edge = logical(),
    biochemical_mrn_pairs = character(),
    biochemical_pathway_overlap = logical(),
    biochemical_pathway_ids = character(),
    biochemical_pathway_names = character(),
    biochemical_evidence_type = character(),
    biochemical_evidence_detail = character(),
    stringsAsFactors = FALSE
  )
}

.syn_empty_biology_links <- function() {
  data.frame(
    node_1 = character(), node_2 = character(),
    module_id_1 = character(), module_id_2 = character(),
    same_module = logical(),
    has_biochemical_link = logical(),
    biochemical_evidence_type = character(),
    biochemical_evidence_detail = character(),
    biochemical_mrn_edge = logical(),
    biochemical_mrn_pairs = character(),
    biochemical_pathway_overlap = logical(),
    biochemical_pathway_ids = character(),
    biochemical_pathway_names = character(),
    has_biology_link = logical(),
    biology_rationale = character(),
    relation_type = character(),
    relation_class = character(),
    subtype = character(),
    rationale = character(),
    stringsAsFactors = FALSE
  )
}

.syn_extract_mrn_edge_keys <- function(mrn) {
  if (is.null(mrn)) return(character(0))
  edge_df <- tryCatch(igraph::as_data_frame(mrn, what = "edges"), error = function(e) NULL)
  if (is.null(edge_df) || !is.data.frame(edge_df) || nrow(edge_df) == 0L) {
    return(character(0))
  }

  keys <- character(0)
  add_keys <- function(from, to) {
    from <- trimws(as.character(from))
    to   <- trimws(as.character(to))
    ok <- !is.na(from) & !is.na(to) & nzchar(from) & nzchar(to)
    .syn_pair_key(from[ok], to[ok])
  }

  if (all(c("INCHI_From", "INCHI_To") %in% colnames(edge_df))) {
    keys <- c(keys, add_keys(edge_df$INCHI_From, edge_df$INCHI_To))
  }
  if (all(c("from", "to") %in% colnames(edge_df))) {
    keys <- c(keys, add_keys(edge_df$from, edge_df$to))
  }

  unique(keys)
}

.syn_pathway_records <- function(pathway_coverage) {
  if (is.null(pathway_coverage) || !is.list(pathway_coverage)) {
    return(data.frame(
      module_id = character(), db = character(),
      pathway_id = character(), pathway_name = character(),
      pathway_key = character(), stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(c("kegg", "hmdb"), function(db) {
    df <- pathway_coverage[[db]]
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L ||
        !"module_id" %in% colnames(df) || !"pathway_id" %in% colnames(df)) {
      return(NULL)
    }
    pathway_name <- if ("pathway_name" %in% colnames(df)) df$pathway_name else NA_character_
    out <- data.frame(
      module_id = trimws(as.character(df$module_id)),
      db = db,
      pathway_id = trimws(as.character(df$pathway_id)),
      pathway_name = trimws(as.character(pathway_name)),
      stringsAsFactors = FALSE
    )
    ok <- !is.na(out$module_id) & nzchar(out$module_id) &
      !is.na(out$pathway_id) & nzchar(out$pathway_id) &
      toupper(out$pathway_id) != "NA"
    out <- out[ok, , drop = FALSE]
    if (nrow(out) == 0L) return(NULL)
    out$pathway_key <- paste(out$db, out$pathway_id, sep = ":")
    out
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame(
      module_id = character(), db = character(),
      pathway_id = character(), pathway_name = character(),
      pathway_key = character(), stringsAsFactors = FALSE
    ))
  }
  unique(do.call(rbind, rows))
}

#' Build deterministic biochemical links for synthesis nodes
#'
#' A pair has a biochemical link if the two nodes come from different modules
#' and their top metabolites have a direct MRN edge while their modules also
#' share at least one metabolic pathway ID in the pathway-coverage evidence.
#' @noRd
build_synthesis_biochemical_links <- function(nodes_df,
                                              mrn = NULL,
                                              pathway_coverage = NULL) {
  out_empty <- .syn_empty_biochemical_links()
  if (is.null(nodes_df) || !is.data.frame(nodes_df) || nrow(nodes_df) < 2L) {
    return(out_empty)
  }

  node_ids <- nodes_df$node_id
  edge_keys <- .syn_extract_mrn_edge_keys(mrn)
  pathway_df <- .syn_pathway_records(pathway_coverage)
  pathways_by_module <- split(pathway_df, pathway_df$module_id)

  top_metabolites_by_node <- setNames(lapply(seq_len(nrow(nodes_df)), function(i) {
    ids <- .syn_split_field(nodes_df$top_metabolite_ids[i])
    nms <- .syn_split_field(nodes_df$top_metabolites[i])
    n <- max(length(ids), length(nms), 0L)
    if (n == 0L) {
      return(data.frame(id = character(), name = character(), stringsAsFactors = FALSE))
    }
    metabolite_df <- data.frame(
      id = vapply(seq_len(n), function(j) if (j <= length(ids)) ids[[j]] else NA_character_, character(1)),
      name = vapply(seq_len(n), function(j) if (j <= length(nms)) nms[[j]] else NA_character_, character(1)),
      stringsAsFactors = FALSE
    )
    metabolite_df[!is.na(metabolite_df$id) & nzchar(metabolite_df$id) &
                    toupper(metabolite_df$id) != "NA", , drop = FALSE]
  }), node_ids)

  mid_lookup <- setNames(nodes_df$module_id, nodes_df$node_id)
  pairs <- combn(sort(node_ids), 2, simplify = FALSE)

  rows <- lapply(pairs, function(p) {
    n1 <- p[[1]]; n2 <- p[[2]]
    m1 <- unname(mid_lookup[n1] %||% NA_character_)
    m2 <- unname(mid_lookup[n2] %||% NA_character_)
    if (is.na(m1) || is.na(m2) || identical(m1, m2)) return(NULL)

    a1 <- top_metabolites_by_node[[n1]]
    a2 <- top_metabolites_by_node[[n2]]
    mrn_pairs <- character(0)
    mrn_name_pairs <- character(0)
    if (nrow(a1) > 0L && nrow(a2) > 0L && length(edge_keys) > 0L) {
      for (i in seq_len(nrow(a1))) {
        for (j in seq_len(nrow(a2))) {
          key <- .syn_pair_key(a1$id[i], a2$id[j])
          if (key %in% edge_keys) {
            left <- if (!is.na(a1$name[i]) && nzchar(a1$name[i])) {
              paste0(a1$name[i], " [", a1$id[i], "]")
            } else a1$id[i]
            right <- if (!is.na(a2$name[j]) && nzchar(a2$name[j])) {
              paste0(a2$name[j], " [", a2$id[j], "]")
            } else a2$id[j]
            mrn_pairs <- c(mrn_pairs, paste(left, right, sep = " -- "))

            left_name <- if (!is.na(a1$name[i]) && nzchar(a1$name[i])) a1$name[i] else "Unnamed metabolite"
            right_name <- if (!is.na(a2$name[j]) && nzchar(a2$name[j])) a2$name[j] else "Unnamed metabolite"
            mrn_name_pairs <- c(mrn_name_pairs, paste(left_name, right_name, sep = " <-> "))
          }
        }
      }
    }

    p1 <- pathways_by_module[[as.character(m1)]] %||% pathway_df[0, , drop = FALSE]
    p2 <- pathways_by_module[[as.character(m2)]] %||% pathway_df[0, , drop = FALSE]
    overlap_keys <- intersect(p1$pathway_key, p2$pathway_key)
    overlap <- unique(rbind(
      p1[p1$pathway_key %in% overlap_keys, , drop = FALSE],
      p2[p2$pathway_key %in% overlap_keys, , drop = FALSE]
    ))
    pathway_ids <- unique(overlap$pathway_key)
    pathway_names <- unique(vapply(seq_len(nrow(overlap)), function(i) {
      db_id <- paste0(toupper(overlap$db[i]), ": ", overlap$pathway_id[i])
      if (is.na(overlap$pathway_name[i]) || !nzchar(overlap$pathway_name[i])) {
        db_id
      } else {
        paste0(overlap$pathway_name[i], " (", db_id, ")")
      }
    }, character(1)))

    has_mrn <- length(mrn_pairs) > 0L
    has_pw  <- length(pathway_ids) > 0L
    if (!has_mrn || !has_pw) return(NULL)

    evidence_type <- "mrn_top_metabolite_edge+pathway_overlap"
    unique_name_pairs <- unique(mrn_name_pairs)
    details <- c(
      if (has_mrn) paste0(
        if (length(unique_name_pairs) == 1L) "Module key-metabolite link: " else "Module key-metabolite links: ",
        paste(unique_name_pairs, collapse = "; ")
      ),
      if (has_pw) paste0(
        if (length(pathway_names) == 1L) "Shared metabolic pathway: " else "Shared metabolic pathways: ",
        paste(pathway_names, collapse = "; ")
      )
    )

    data.frame(
      node_1 = n1,
      node_2 = n2,
      module_id_1 = m1,
      module_id_2 = m2,
      biochemical_mrn_edge = has_mrn,
      biochemical_mrn_pairs = if (has_mrn) paste(unique(mrn_pairs), collapse = "; ") else NA_character_,
      biochemical_pathway_overlap = has_pw,
      biochemical_pathway_ids = if (has_pw) paste(pathway_ids, collapse = "; ") else NA_character_,
      biochemical_pathway_names = if (has_pw) paste(pathway_names, collapse = "; ") else NA_character_,
      biochemical_evidence_type = evidence_type,
      biochemical_evidence_detail = paste(details, collapse = " | "),
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(out_empty)
  do.call(rbind, rows)
}


# ============================================================
# SECTION 4b \u2013 Relationship table completion
# ============================================================

# After the LLM returns sparse cross-module biology links and R computes
# biochemical links, this function completes the full N*(N-1)/2 pairwise table.
# Final relation_type is derived, never directly assigned by the LLM.
complete_relationship_table <- function(biology_relationships,
                                        nodes_df,
                                        biochemical_links = NULL) {

  node_ids <- nodes_df$node_id
  n        <- length(node_ids)

  empty_complete <- .syn_empty_biology_links()
  if (n < 2) return(empty_complete)

  # All unique pairs (lexicographic order)
  pairs <- combn(sort(node_ids), 2, simplify = FALSE)

  # Module lookup
  mid_lookup <- setNames(nodes_df$module_id, nodes_df$node_id)

  # Relationship lookups: "node_1|||node_2" -> row index
  if (!is.null(biology_relationships) && nrow(biology_relationships) > 0) {
    bio_key <- paste(biology_relationships$node_1, biology_relationships$node_2, sep = "|||")
  } else {
    bio_key <- character(0)
  }
  if (!is.null(biochemical_links) && nrow(biochemical_links) > 0) {
    biochem_key <- paste(biochemical_links$node_1, biochemical_links$node_2, sep = "|||")
  } else {
    biochem_key <- character(0)
  }

  rows <- lapply(pairs, function(p) {
    n1 <- p[1]; n2 <- p[2]
    key <- paste(n1, n2, sep = "|||")

    m1 <- unname(mid_lookup[n1] %||% NA_character_)
    m2 <- unname(mid_lookup[n2] %||% NA_character_)

    if (!is.na(m1) && !is.na(m2) && identical(m1, m2)) {
      # Same module \u2014 hypotheses from the same debate; may be complementary or alternative.
      data.frame(
        node_1         = n1,
        node_2         = n2,
        module_id_1    = m1,
        module_id_2    = m2,
        same_module    = TRUE,
        has_biochemical_link = FALSE,
        biochemical_evidence_type = NA_character_,
        biochemical_evidence_detail = NA_character_,
        biochemical_mrn_edge = FALSE,
        biochemical_mrn_pairs = NA_character_,
        biochemical_pathway_overlap = FALSE,
        biochemical_pathway_ids = NA_character_,
        biochemical_pathway_names = NA_character_,
        has_biology_link = FALSE,
        biology_rationale = NA_character_,
        relation_type = "same_module",
        relation_class = "same_module",
        subtype = "same_module_hypotheses",
        rationale = paste0("Both nodes are hypotheses for module ", m1,
                           " generated in the same module-level debate."),
        stringsAsFactors = FALSE
      )
    } else {
      bio <- if (key %in% bio_key) {
        biology_relationships[bio_key == key, , drop = FALSE][1, ]
      } else NULL
      biochem <- if (key %in% biochem_key) {
        biochemical_links[biochem_key == key, , drop = FALSE][1, ]
      } else NULL

      has_bio <- !is.null(bio)
      has_biochem <- !is.null(biochem)
      relation_type <- if (has_bio && has_biochem) {
        "both"
      } else if (has_biochem) {
        "biochemical"
      } else if (has_bio) {
        "biology"
      } else {
        "none"
      }

      biochem_detail <- if (has_biochem) biochem$biochemical_evidence_detail else NA_character_
      bio_rationale <- if (has_bio) bio$biology_rationale %||% bio$rationale else NA_character_
      rationale <- paste(c(
        if (has_biochem) biochem_detail,
        if (has_bio) bio_rationale
      ), collapse = " | ")
      if (!nzchar(rationale)) {
        rationale <- "No biochemical or biology link was identified."
      }

      data.frame(
        node_1         = n1,
        node_2         = n2,
        module_id_1    = m1,
        module_id_2    = m2,
        same_module    = FALSE,
        has_biochemical_link = has_biochem,
        biochemical_evidence_type = if (has_biochem) biochem$biochemical_evidence_type else NA_character_,
        biochemical_evidence_detail = biochem_detail,
        biochemical_mrn_edge = if (has_biochem) biochem$biochemical_mrn_edge else FALSE,
        biochemical_mrn_pairs = if (has_biochem) biochem$biochemical_mrn_pairs else NA_character_,
        biochemical_pathway_overlap = if (has_biochem) biochem$biochemical_pathway_overlap else FALSE,
        biochemical_pathway_ids = if (has_biochem) biochem$biochemical_pathway_ids else NA_character_,
        biochemical_pathway_names = if (has_biochem) biochem$biochemical_pathway_names else NA_character_,
        has_biology_link = has_bio,
        biology_rationale = bio_rationale,
        relation_type = relation_type,
        relation_class = relation_type,
        subtype = relation_type,
        rationale = rationale,
        stringsAsFactors = FALSE
      )
    }
  })

  do.call(rbind, rows)
}

.syn_empty_module_relationships <- function() {
  data.frame(
    module_id_1 = character(), module_id_2 = character(),
    has_biochemical_link = logical(),
    has_biology_link = logical(),
    relation_type = character(),
    biochemical_mrn_edge = logical(),
    biochemical_pathway_overlap = logical(),
    biochemical_evidence_type = character(),
    biochemical_evidence_detail = character(),
    biology_rationale = character(),
    contributing_node_pairs = character(),
    stringsAsFactors = FALSE
  )
}

#' Summarize hypothesis-pair relationships to module-pair relationships
#'
#' The synthesis graph remains hypothesis-level, but this table gives one row
#' per module pair for downstream module-level interpretation.
#' @noRd
summarize_module_relationships <- function(relationships) {
  out_empty <- .syn_empty_module_relationships()
  if (is.null(relationships) || !is.data.frame(relationships) ||
      nrow(relationships) == 0L) {
    return(out_empty)
  }

  rel <- relationships[
    !is.na(relationships$module_id_1) &
      !is.na(relationships$module_id_2) &
      relationships$module_id_1 != relationships$module_id_2,
    ,
    drop = FALSE
  ]
  if (nrow(rel) == 0L) return(out_empty)

  m1 <- pmin(as.character(rel$module_id_1), as.character(rel$module_id_2))
  m2 <- pmax(as.character(rel$module_id_1), as.character(rel$module_id_2))
  pair_key <- paste(m1, m2, sep = "|||")

  rows <- lapply(unique(pair_key), function(key) {
    idx <- pair_key == key
    df <- rel[idx, , drop = FALSE]
    mods <- strsplit(key, "|||", fixed = TRUE)[[1]]
    has_biochem <- any(df$has_biochemical_link, na.rm = TRUE)
    has_bio <- any(df$has_biology_link, na.rm = TRUE)
    relation_type <- if (has_biochem && has_bio) {
      "both"
    } else if (has_biochem) {
      "biochemical"
    } else if (has_bio) {
      "biology"
    } else {
      "none"
    }

    mrn_edge <- any(df$biochemical_mrn_edge, na.rm = TRUE)
    pathway_overlap <- any(df$biochemical_pathway_overlap, na.rm = TRUE)
    clean_join <- function(x) {
      x <- unique(trimws(as.character(x)))
      x <- x[!is.na(x) & nzchar(x) & toupper(x) != "NA"]
      if (length(x) == 0L) NA_character_ else paste(x, collapse = " | ")
    }

    contributing <- paste0(df$node_1, " -- ", df$node_2, " [", df$relation_type, "]")

    data.frame(
      module_id_1 = mods[[1]],
      module_id_2 = mods[[2]],
      has_biochemical_link = has_biochem,
      has_biology_link = has_bio,
      relation_type = relation_type,
      biochemical_mrn_edge = mrn_edge,
      biochemical_pathway_overlap = pathway_overlap,
      biochemical_evidence_type = clean_join(df$biochemical_evidence_type),
      biochemical_evidence_detail = clean_join(df$biochemical_evidence_detail),
      biology_rationale = clean_join(df$biology_rationale),
      contributing_node_pairs = clean_join(contributing),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


# ============================================================
# SECTION 4c \u2013 Output parsing & assembly
# ============================================================

#' Parse synthesis LLM output into tidy data frames
#'
#' @param parsed   Parsed list from \code{adv_parse_json_response()}.
#' @param nodes_df Node data frame from \code{prepare_synthesis_nodes()}.
#' @param biochemical_links Data frame from \code{build_synthesis_biochemical_links()}.
#' @return List with \code{$clusters}, \code{$biology_links},
#'   \code{$biochemical_links}, \code{$relationships},
#'   \code{$module_relationships}, \code{$summary}, \code{$nodes}.
parse_synthesis_output <- function(parsed, nodes_df, biochemical_links = NULL) {

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
  if (nrow(clusters_df) > 0) {
    assigned_nodes <- character(0)
    for (i in seq_len(nrow(clusters_df))) {
      nids <- .syn_split_field(clusters_df$node_ids[i])
      nids <- unique(nids[nids %in% valid_nodes])
      keep <- setdiff(nids, assigned_nodes)
      clusters_df$node_ids[i] <- paste(keep, collapse = "; ")
      assigned_nodes <- unique(c(assigned_nodes, keep))
    }
    clusters_df <- clusters_df[nzchar(clusters_df$node_ids), , drop = FALSE]
  }

  # ---- relationships ----
  relationships_df <- tryCatch({
    rl <- parsed$relationships
    if (is.null(rl) || length(rl) == 0) {
      data.frame(
        node_1 = character(), node_2 = character(),
        biology_rationale = character(),
        stringsAsFactors = FALSE
      )
    } else {
      do.call(rbind, lapply(rl, function(r) {
        n1 <- trimws(as.character(r$node_1 %||% ""))
        n2 <- trimws(as.character(r$node_2 %||% ""))
        # Enforce lexicographic ordering to avoid duplicate pairs
        if (n1 > n2) { tmp <- n1; n1 <- n2; n2 <- tmp }
        data.frame(
          node_1 = n1,
          node_2 = n2,
          biology_rationale = trimws(as.character(r$rationale %||% "")),
          stringsAsFactors = FALSE
        )
      }))
    }
  }, error = function(e) {
    data.frame(
      node_1 = character(), node_2 = character(),
      biology_rationale = character(),
      stringsAsFactors = FALSE
    )
  })

  # Deduplicate pairs (same node_1/node_2)
  if (nrow(relationships_df) > 0) {
    relationships_df <- relationships_df[
      relationships_df$node_1 %in% valid_nodes &
        relationships_df$node_2 %in% valid_nodes &
        relationships_df$node_1 != relationships_df$node_2 &
        nzchar(relationships_df$biology_rationale),
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
    relationships_df$has_biology_link <- TRUE
    relationships_df$relation_type <- "biology"
    relationships_df$relation_class <- "biology"
    relationships_df$subtype <- "biology_link"
    relationships_df$rationale <- relationships_df$biology_rationale
  } else {
    relationships_df$module_id_1  <- character()
    relationships_df$module_id_2  <- character()
    relationships_df$same_module   <- logical()
    relationships_df$has_biology_link <- logical()
    relationships_df$relation_type <- character()
    relationships_df$relation_class <- character()
    relationships_df$subtype <- character()
    relationships_df$rationale <- character()
  }

  relationships_sparse <- relationships_df
  if (is.null(biochemical_links)) biochemical_links <- .syn_empty_biochemical_links()
  relationships_complete <- complete_relationship_table(
    biology_relationships = relationships_sparse,
    nodes_df = nodes_df,
    biochemical_links = biochemical_links
  )
  module_relationships <- summarize_module_relationships(relationships_complete)

  # Enrich nodes with cluster assignment
  nodes_out <- nodes_df
  nodes_out$cluster_id <- NA_character_
  for (i in seq_len(nrow(clusters_df))) {
    nids <- .syn_split_field(clusters_df$node_ids[i])
    needs_cluster <- nodes_out$node_id %in% nids & is.na(nodes_out$cluster_id)
    nodes_out$cluster_id[needs_cluster] <- clusters_df$cluster_id[i]
  }

  list(
    nodes                = nodes_out,
    clusters             = clusters_df,
    biology_links        = relationships_sparse,
    biochemical_links    = biochemical_links,
    relationships_sparse = relationships_sparse,
    relationships        = relationships_complete,
    module_relationships = module_relationships,
    summary              = trimws(as.character(parsed$summary %||% ""))
  )
}


# ============================================================
# SECTION 5 \u2013 Main orchestrator
# ============================================================

#' Synthesize All Annotated Module Hypotheses
#'
#' Takes the \code{module_function} data frame from an annotated featuremap
#' result and runs a single LLM call to produce:
#' \itemize{
#'   \item Thematic clusters of hypotheses
#'   \item LLM-inferred biology links from hypothesis titles/summaries
#'   \item Deterministic biochemical links from direct MRN top-metabolite edges
#'         and metabolic pathway ID overlap
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
#' @param base_url       Character or NULL. Optional API base URL for an
#'   OpenAI-compatible provider.
#' @param max_tokens     Max tokens for the LLM response (default 4000).
#' @param temperature    Sampling temperature (default 0).
#' @param sleep_sec      Seconds to pause after the API call (default 0).
#' @return The input \code{res} object with a \code{module_synthesis} list
#'   appended, containing \code{$nodes}, \code{$clusters},
#'   \code{$biology_links}, \code{$biochemical_links},
#'   \code{$relationships_sparse}, \code{$relationships},
#'   \code{$module_relationships},
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
      cli::cli_abort("synthesize_modules: failed to prepare nodes \u2014 {e$message}")
    }
  )

  n_nodes   <- nrow(nodes_df)
  n_modules <- length(unique(nodes_df$module_id))
  pathway_coverage <- NULL
  if (!is.null(res$module_evidence) && !is.null(res$module_evidence$pathway_coverage)) {
    pathway_coverage <- res$module_evidence$pathway_coverage
  }
  biochemical_links <- build_synthesis_biochemical_links(
    nodes_df = nodes_df,
    mrn = res$mrn %||% NULL,
    pathway_coverage = pathway_coverage
  )

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
        node_1 = character(), node_2 = character(),
        biology_rationale = character(), module_id_1 = character(),
        module_id_2 = character(), same_module = logical(),
        has_biology_link = logical(), relation_type = character(),
        relation_class = character(), subtype = character(),
        rationale = character(), stringsAsFactors = FALSE
      ),
      biology_links        = .syn_empty_biology_links()[0, , drop = FALSE],
      biochemical_links    = biochemical_links,
      relationships        = complete_relationship_table(NULL, nodes_df, biochemical_links),
      module_relationships = summarize_module_relationships(
        complete_relationship_table(NULL, nodes_df, biochemical_links)
      ),
      summary              = NA_character_,
      raw_response         = NA_character_
    )
    return(res)
  }

  cli::cli_alert_info(
    "Synthesizing {.val {n_nodes}} hypothesis node(s) across {.val {n_modules}} module(s); precomputed {.val {nrow(biochemical_links)}} biochemical link(s)."
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
        "\"relationships\" (array of biology links), and \"summary\" (string)."
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
        parse_synthesis_output(result$parsed, nodes_df,
                               biochemical_links = biochemical_links),
      error = function(e) {
        cli::cli_warn("synthesize_modules: output parsing failed \u2014 {e$message}")
        list(nodes = nodes_df, clusters = NULL,
             biology_links = .syn_empty_biology_links()[0, , drop = FALSE],
             biochemical_links = biochemical_links,
             relationships_sparse = NULL,
             relationships = complete_relationship_table(NULL, nodes_df, biochemical_links),
             module_relationships = summarize_module_relationships(
               complete_relationship_table(NULL, nodes_df, biochemical_links)
             ),
             summary = NA_character_)
      }
    )
  } else {
    cli::cli_warn("synthesize_modules: LLM returned no parseable output.")
    list(nodes = nodes_df, clusters = NULL,
         biology_links = .syn_empty_biology_links()[0, , drop = FALSE],
         biochemical_links = biochemical_links,
         relationships_sparse = NULL,
         relationships = complete_relationship_table(NULL, nodes_df, biochemical_links),
         module_relationships = summarize_module_relationships(
           complete_relationship_table(NULL, nodes_df, biochemical_links)
         ),
         summary = NA_character_)
  }

  synthesis$raw_response <- result$raw_text %||% NA_character_

  n_rel_sparse <- if (!is.null(synthesis$relationships_sparse)) nrow(synthesis$relationships_sparse) else 0L
  n_rel_full   <- if (!is.null(synthesis$relationships))        nrow(synthesis$relationships)        else 0L
  n_rel_biochem <- if (!is.null(synthesis$biochemical_links))    nrow(synthesis$biochemical_links)    else 0L
  n_cl         <- if (!is.null(synthesis$clusters))             nrow(synthesis$clusters)             else 0L
  cli::cli_alert_success(
    "Synthesis complete: {.val {n_cl}} cluster(s), {.val {n_rel_sparse}} biology link(s), {.val {n_rel_biochem}} biochemical link(s), {.val {n_rel_full}} completed pair(s)."
  )

  res$module_synthesis <- synthesis
  return(res)
}
