
#' Flatten PubMed Data for a Single Module
#'
#' Converts a nested list of PubMed entries (one per HMDB ID) into a flat tibble
#' containing HMDB IDs and their functional summaries.
#'
#' @param pubmed_list_for_module A named list where names are HMDB IDs and values
#'   are lists containing metabolite information including functional_summary
#'
#' @return A tibble with columns:
#'   \item{HMDB_ID}{Character. HMDB identifier for the metabolite}
#'   \item{functional_summary}{Character. Functional summary from PubMed, or NA if missing}
#'
#' @details Uses purrr::imap_dfr to iterate over the list and bind rows.
#'   The %||% operator provides NA_character_ as default when functional_summary is NULL.
#'
#' @examples
#' \dontrun{
#' pubmed_data <- list(
#'   "HMDB0000001" = list(functional_summary = "Amino acid metabolism"),
#'   "HMDB0000002" = list(functional_summary = NULL)
#' )
#' flatten_pubmed_for_one_module(pubmed_data)
#' }
flatten_pubmed_for_one_module <- function(pubmed_list_for_module) {
  imap_dfr(pubmed_list_for_module, function(entry, id) {
    tibble(
      Compound_ID            = id,
      functional_summary = entry$functional_summary %||% NA_character_
    )
  })
}

#' Combine Top 10 Metabolites with PubMed Functional Summaries
#'
#' Performs a left join to merge the top metabolites data frame with
#' PubMed functional summaries, then selects and orders relevant columns.
#'
#' @param module_top_df A data frame containing the top metabolites with columns:
#'   top_id, top_rank, Compound_ID, Compound_name, class, priority_level, and
#'   optional INCHI_ID, HMDB_ID, KEGG_ID
#' @param pubmed_top_df A data frame containing PubMed data with columns:
#'   HMDB_ID, functional_summary
#'
#' @return A tibble with columns ordered as: top_id, top_rank, Compound_ID,
#'   Compound_name, class, priority_level, optional identifier columns, and
#'   functional_summary
#'
#' @details The left join ensures all metabolites from module_top_df are retained,
#'   even if they don't have corresponding PubMed data.
combine_topN_with_pubmed <- function(module_top_df, pubmed_top_df) {
  module_top_df %>%
    left_join(pubmed_top_df, by = "Compound_ID") %>%
    dplyr::select(
      dplyr::any_of(c(
        "top_id", "top_rank", "Compound_ID", "INCHI_ID", "HMDB_ID", "KEGG_ID",
        "Compound_name", "class", "priority_level", "functional_summary"
      ))
    )
}


# ============================================================
# Top metabolite reference helpers
# ============================================================

add_top_ids <- function(top_df) {
  if (is.null(top_df) || !is.data.frame(top_df) || nrow(top_df) == 0) {
    return(top_df)
  }
  top_df <- tibble::as_tibble(top_df)
  if (!"top_id" %in% colnames(top_df)) {
    if ("top_rank" %in% colnames(top_df)) {
      top_df <- dplyr::mutate(top_df, top_id = paste0("T", as.integer(top_rank)))
    } else {
      top_df <- dplyr::mutate(top_df, top_id = paste0("T", dplyr::row_number()))
    }
  }
  top_df
}

normalize_top_id_tokens <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  x <- unique(trimws(as.character(unlist(x))))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(character(0))

  vapply(x, function(tok) {
    up <- toupper(tok)
    up <- gsub("[[:space:]_-]+", "", up)
    if (grepl("^T[0-9]+$", up)) {
      return(paste0("T", as.integer(sub("^T", "", up))))
    }
    digits <- regmatches(up, regexpr("[0-9]+", up))
    if (length(digits) == 1 && nzchar(digits) &&
        grepl("^(TOP|RANK)?[0-9]+$", up)) {
      return(paste0("T", as.integer(digits)))
    }
    tok
  }, character(1), USE.NAMES = FALSE)
}

resolve_top_metabolite_refs <- function(top_ids = NULL,
                                        top_df,
                                        ref_type = "anchor") {
  empty <- list(
    top_ids      = character(0),
    names        = character(0),
    compound_ids = character(0),
    unmatched    = character(0),
    table        = tibble::tibble()
  )
  if (is.null(top_df) || !is.data.frame(top_df) || nrow(top_df) == 0) {
    unmatched <- unique(trimws(as.character(unlist(top_ids))))
    unmatched <- unmatched[!is.na(unmatched) & nzchar(unmatched)]
    empty$unmatched <- unmatched
    return(empty)
  }

  top_df <- add_top_ids(top_df)
  required <- c("top_id", "Compound_name", "Compound_ID")
  for (col in setdiff(required, colnames(top_df))) top_df[[col]] <- NA_character_

  valid_ids <- unique(trimws(as.character(top_df$top_id)))
  valid_ids <- valid_ids[!is.na(valid_ids) & nzchar(valid_ids)]
  requested_ids_raw <- unique(trimws(as.character(unlist(top_ids))))
  requested_ids_raw <- requested_ids_raw[!is.na(requested_ids_raw) & nzchar(requested_ids_raw)]
  requested_ids <- normalize_top_id_tokens(requested_ids_raw)
  matched_ids <- requested_ids[requested_ids %in% valid_ids]

  unmatched <- requested_ids_raw[!(normalize_top_id_tokens(requested_ids_raw) %in% matched_ids)]

  matched_ids <- unique(matched_ids[!is.na(matched_ids) & nzchar(matched_ids)])
  matched_rows <- top_df[top_df$top_id %in% matched_ids, , drop = FALSE]
  if (nrow(matched_rows) > 0) {
    matched_rows <- matched_rows[match(matched_ids, matched_rows$top_id), , drop = FALSE]
    matched_rows <- matched_rows[!is.na(matched_rows$top_id), , drop = FALSE]
  }

  out_table <- if (nrow(matched_rows) == 0) {
    tibble::tibble()
  } else {
    tibble::tibble(
      top_id = as.character(matched_rows$top_id),
      metabolite_name = as.character(matched_rows$Compound_name),
      Compound_ID = as.character(matched_rows$Compound_ID),
      ref_type = ref_type
    )
  }

  list(
    top_ids      = out_table$top_id %||% character(0),
    names        = out_table$metabolite_name %||% character(0),
    compound_ids = out_table$Compound_ID %||% character(0),
    unmatched    = unique(unmatched[!is.na(unmatched) & nzchar(unmatched)]),
    table        = out_table
  )
}



# ============================================================
# New helper: compute_global_direction
# ============================================================

#' Compute Global Directionality for a Metabolic Module
#'
#' Matches module metabolites to features via anno_feature_table2 and computes
#' a weighted mean of condition values (log2FC or correlation coefficients).
#' The result is stored as an output label only — it does NOT enter LLM prompts.
#'
#' @param module_df Data frame of module metabolites. Must contain HMDB_ID and
#'   KEGG_ID columns.
#' @param feature_table Data frame with at least variable_id (character) and
#'   condition (numeric: log2FC or correlation/regression coefficient).
#' @param anno_feature_table2 Data frame linking features to compounds. Must
#'   contain variable_id, and at least one of HMDB.ID / HMDB_ID, KEGG.ID /
#'   KEGG_ID, plus \code{weighted_anno_score}.
#'
#' @return A named list:
#'   \item{direction_label}{Character: "up", "down", or "mixed". NA if input
#'     is unavailable.}
#'   \item{weighted_mean_condition}{Numeric weighted mean of condition, or NA.}
#'   \item{n_features}{Integer. Number of matched features used.}
compute_global_direction <- function(module_df, feature_table, anno_feature_table2) {

  na_result <- list(direction_label = NA_character_,
                    weighted_mean_condition = NA_real_,
                    n_features = 0L)

  if (is.null(feature_table) || !is.data.frame(feature_table) ||
      !all(c("variable_id", "condition") %in% colnames(feature_table))) {
    return(na_result)
  }
  if (is.null(anno_feature_table2) || !is.data.frame(anno_feature_table2) ||
      !"variable_id" %in% colnames(anno_feature_table2)) {
    return(na_result)
  }

  if (!"weighted_anno_score" %in% colnames(anno_feature_table2)) {
    return(na_result)
  }

  # Detect ID columns in anno_feature_table2
  hmdb_col <- NULL
  for (hc in c("HMDB.ID", "HMDB_ID")) {
    if (hc %in% colnames(anno_feature_table2)) { hmdb_col <- hc; break }
  }
  kegg_col <- NULL
  for (kc in c("KEGG.ID", "KEGG_ID")) {
    if (kc %in% colnames(anno_feature_table2)) { kegg_col <- kc; break }
  }

  hmdb_ids <- unique(module_df$HMDB_ID[
    !is.na(module_df$HMDB_ID) & nzchar(as.character(module_df$HMDB_ID))])
  kegg_ids <- unique(module_df$KEGG_ID[
    !is.na(module_df$KEGG_ID) & nzchar(as.character(module_df$KEGG_ID))])

  # Build row-selection filter
  keep <- rep(FALSE, nrow(anno_feature_table2))
  if (!is.null(hmdb_col) && length(hmdb_ids) > 0) {
    v <- anno_feature_table2[[hmdb_col]]
    keep <- keep | (!is.na(v) & v %in% hmdb_ids)
  }
  if (!is.null(kegg_col) && length(kegg_ids) > 0) {
    v <- anno_feature_table2[[kegg_col]]
    keep <- keep | (!is.na(v) & v %in% kegg_ids)
  }
  if (!any(keep)) return(na_result)

  anno_sub <- anno_feature_table2[keep, , drop = FALSE]
  anno_sub <- anno_sub[, c("variable_id", "weighted_anno_score"), drop = FALSE]
  colnames(anno_sub)[2] <- "score"
  anno_sub <- anno_sub[!duplicated(anno_sub$variable_id), , drop = FALSE]

  combined <- merge(
    anno_sub,
    feature_table[, c("variable_id", "condition"), drop = FALSE],
    by = "variable_id"
  )
  combined <- combined[
    !is.na(combined$condition) & is.finite(combined$condition) &
    !is.na(combined$score)     & is.finite(combined$score) &
    combined$score > 0, , drop = FALSE]

  if (nrow(combined) == 0) return(na_result)

  wmean <- sum(combined$condition * combined$score) / sum(combined$score)

  direction_label <- if (wmean > 0) "up" else if (wmean < 0) "down" else "mixed"

  list(
    direction_label         = direction_label,
    weighted_mean_condition = round(wmean, 4),
    n_features              = nrow(combined)
  )
}


# ============================================================
# New helper: compute_pathway_coverage
# ============================================================

#' Compute Pathway Coverage for a Module
#'
#' For each pathway database (HMDB and KEGG), counts how many module
#' metabolites appear in each pathway and computes the hit/total coverage
#' fraction.  Returns the top \code{top_n} pathways sorted by hit count then
#' coverage.
#'
#' @param module_df Data frame of module metabolites. Must contain HMDB_ID and
#'   KEGG_ID columns.
#' @param pathway_hmdb Data frame loaded via
#'   \code{featuremap_load_resource("pathway_hmdb")} with columns
#'   pathway_id, pathway_name, HMDB_ID, n_metabolites. Pass \code{NULL} to
#'   skip.
#' @param pathway_kegg Data frame loaded via
#'   \code{featuremap_load_resource("pathway_kegg")} with columns
#'   pathway_id, pathway_name, KEGG_ID, n_metabolites. Pass \code{NULL} to
#'   skip.
#' @param top_n Integer. Maximum number of pathways to return per database
#'   (default: 3).
#'
#' @return A named list with components \code{hmdb} and \code{kegg}, each a
#'   tibble with columns pathway_id, pathway_name, n_hit, n_metabolites,
#'   coverage.
#'   Empty tibble when the database is unavailable or no hits are found.
compute_pathway_coverage <- function(module_df,
                                     pathway_hmdb = NULL,
                                     pathway_kegg = NULL,
                                     top_n = 3) {

  compute_one <- function(module_ids, pathway_db, id_col) {
    if (is.null(pathway_db) || !is.data.frame(pathway_db)) return(tibble::tibble())
    if (!id_col %in% colnames(pathway_db))                 return(tibble::tibble())
    if (!all(c("pathway_id", "pathway_name", "n_metabolites") %in% colnames(pathway_db)))
      return(tibble::tibble())
    if (length(module_ids) == 0)                           return(tibble::tibble())

    # Per-pathway total size (one canonical value per pathway_id)
    pathway_sizes <- pathway_db[!duplicated(pathway_db$pathway_id),
                                c("pathway_id", "n_metabolites"), drop = FALSE]

    hits <- pathway_db[pathway_db[[id_col]] %in% module_ids, , drop = FALSE]
    if (nrow(hits) == 0) return(tibble::tibble())

    agg <- aggregate(
      as.formula(paste(id_col, "~ pathway_id + pathway_name")),
      data = hits,
      FUN  = function(x) length(unique(x))
    )
    colnames(agg)[colnames(agg) == id_col] <- "n_hit"

    agg <- merge(agg, pathway_sizes, by = "pathway_id", all.x = TRUE)
    agg$coverage <- round(agg$n_hit / pmax(agg$n_metabolites, 1L), 3)
    agg <- agg[order(-agg$coverage, -agg$n_hit), , drop = FALSE]
    agg <- agg[seq_len(min(top_n, nrow(agg))), , drop = FALSE]
    tibble::as_tibble(agg[, c("pathway_id", "pathway_name", "n_hit",
                              "n_metabolites", "coverage"), drop = FALSE])
  }

  hmdb_ids <- unique(module_df$HMDB_ID[
    !is.na(module_df$HMDB_ID) & nzchar(as.character(module_df$HMDB_ID))])
  kegg_ids <- unique(module_df$KEGG_ID[
    !is.na(module_df$KEGG_ID) & nzchar(as.character(module_df$KEGG_ID))])

  list(
    hmdb = compute_one(hmdb_ids, pathway_hmdb, "HMDB_ID"),
    kegg = compute_one(kegg_ids, pathway_kegg, "KEGG_ID")
  )
}


# ============================================================
# New helper: compute_origin_profile
# ============================================================

#' Compute Origin Profile for a Module
#'
#' Looks up the compound origin (endogenous vs. exogenous) for module
#' metabolites in the HMDB functional reference and returns proportions.
#' Only HMDB-annotated metabolites are counted.
#'
#' @param module_df Data frame of module metabolites with an HMDB_ID column.
#' @param hmdb_function_ref Data frame. HMDB functional reference loaded via
#'   \code{featuremap_load_resource("functional_reference_hmdb")}. Must
#'   contain an \code{id} column (HMDB IDs) and an origin / source column
#'   (detected automatically by name pattern).
#'
#' @return A named list with numeric slots endogenous_frac, exogenous_frac,
#'   and integer slots endogenous_n, exogenous_n, unknown_n, n_annotated.
#'   All slots are NA when the reference or origin column cannot be found.
compute_origin_profile <- function(module_df, hmdb_function_ref = NULL) {

  na_result <- list(
    endogenous_frac = NA_real_,
    exogenous_frac  = NA_real_,
    endogenous_n    = NA_integer_,
    exogenous_n     = NA_integer_,
    unknown_n       = NA_integer_,
    n_annotated     = NA_integer_
  )

  if (is.null(hmdb_function_ref) || !is.data.frame(hmdb_function_ref)) {
    return(na_result)
  }

  # Normalise id column
  ref <- hmdb_function_ref
  if (!"id" %in% colnames(ref)) {
    id_cand <- grep("^hmdb_?id$", colnames(ref), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(id_cand)) {
      colnames(ref)[colnames(ref) == id_cand] <- "id"
    } else {
      warning("[compute_origin_profile] Cannot find id column. Available: ",
              paste(colnames(ref), collapse = ", "))
      return(na_result)
    }
  }

  # Detect origin column
  origin_col <- grep("origin|source", colnames(ref), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(origin_col)) {
    warning("[compute_origin_profile] No origin/source column found. Available: ",
            paste(colnames(ref), collapse = ", "))
    return(na_result)
  }

  hmdb_ids <- unique(module_df$HMDB_ID[
    !is.na(module_df$HMDB_ID) & nzchar(as.character(module_df$HMDB_ID))])

  zero_result <- list(endogenous_frac = NA_real_, exogenous_frac = NA_real_,
                      endogenous_n = 0L, exogenous_n = 0L,
                      unknown_n = 0L, n_annotated = 0L)

  if (length(hmdb_ids) == 0) return(zero_result)

  module_ref <- ref[ref$id %in% hmdb_ids,
                    c("id", origin_col), drop = FALSE]
  module_ref <- module_ref[!duplicated(module_ref$id), , drop = FALSE]

  n_annotated <- nrow(module_ref)
  if (n_annotated == 0) return(zero_result)

  origins_lc <- tolower(as.character(module_ref[[origin_col]]))

  endogenous_n <- sum(grepl("endogenous", origins_lc), na.rm = TRUE)
  exogenous_n  <- sum(
    grepl("exogenous|dietary|microbial|drug|xenobiotic|environmental|food|plant",
          origins_lc),
    na.rm = TRUE)
  unknown_n <- as.integer(max(0L, n_annotated - endogenous_n - exogenous_n))

  list(
    endogenous_frac = round(endogenous_n / n_annotated, 3),
    exogenous_frac  = round(exogenous_n  / n_annotated, 3),
    endogenous_n    = as.integer(endogenous_n),
    exogenous_n     = as.integer(exogenous_n),
    unknown_n       = unknown_n,
    n_annotated     = as.integer(n_annotated)
  )
}


# ============================================================
# Generic anchor metabolite ID constants
# Ubiquitous energy/redox cofactors that appear in virtually
# every metabolic pathway and must not dominate annotations.
# ============================================================

.FEATUREMAP_GENERIC_ANCHOR_HMDB <- c(
  "HMDB0000538",  # ATP
  "HMDB0001341",  # ADP
  "HMDB0000045",  # AMP
  "HMDB0001273",  # GTP
  "HMDB0000244",  # GDP
  "HMDB0001429",  # Inorganic phosphate (Pi)
  "HMDB0000250",  # Pyrophosphate (PPi)
  "HMDB0000902",  # NAD+
  "HMDB0001487",  # NADH
  "HMDB0000217",  # NADP+
  "HMDB0000221",  # NADPH
  "HMDB0001248",  # FAD
  "HMDB0001197",  # FADH2
  "HMDB0001301",  # Coenzyme A (CoA)
  "HMDB0001353",  # CO2
  "HMDB0002111"   # H2O
)

.FEATUREMAP_GENERIC_ANCHOR_KEGG <- c(
  "C00002",  # ATP
  "C00008",  # ADP
  "C00020",  # AMP
  "C00044",  # GTP
  "C00035",  # GDP
  "C00009",  # Inorganic phosphate (Pi)
  "C00013",  # Pyrophosphate (PPi)
  "C00003",  # NAD+
  "C00004",  # NADH
  "C00006",  # NADP+
  "C00005",  # NADPH
  "C00016",  # FAD
  "C01352",  # FADH2
  "C00010",  # Coenzyme A (CoA)
  "C00011",  # CO2
  "C00001"   # H2O
)


# ============================================================
# New helper: compute_generic_anchor_flags
# ============================================================

#' Check Whether Top Metabolites Are Dominated by Generic Cofactors
#'
#' Matches the topN metabolite HMDB/KEGG IDs against a predefined list of
#' ubiquitous energy/redox cofactors (ATP, ADP, NAD+, CoA, etc.) that should
#' not serve as primary annotation anchors.
#'
#' @param topN_df Data frame of top metabolites from prepare_module(). Must
#'   contain HMDB_ID/KEGG_ID columns for database-ID matching.
#' @param hmdb_ids Character vector of HMDB IDs considered generic anchors.
#'   Defaults to the package-level .FEATUREMAP_GENERIC_ANCHOR_HMDB list.
#' @param kegg_ids Character vector of KEGG IDs considered generic anchors.
#'   Defaults to the package-level .FEATUREMAP_GENERIC_ANCHOR_KEGG list.
#'
#' @return A named list:
#'   \item{has_generic}{Logical. TRUE if any topN metabolite is generic.}
#'   \item{generic_names}{Character vector of generic metabolite names found.}
#'   \item{generic_ids}{Character vector of matched generic HMDB/KEGG IDs.}
#'   \item{n_generic}{Integer. Count of generic metabolites in topN.}
#'   \item{n_top}{Integer. Total size of topN.}
compute_generic_anchor_flags <- function(
    topN_df,
    hmdb_ids            = .FEATUREMAP_GENERIC_ANCHOR_HMDB,
    kegg_ids            = .FEATUREMAP_GENERIC_ANCHOR_KEGG
) {
  na_result <- list(
    has_generic   = FALSE,
    generic_names = character(0),
    generic_ids   = character(0),
    n_generic     = 0L,
    n_top         = 0L
  )

  if (is.null(topN_df) || !is.data.frame(topN_df) || nrow(topN_df) == 0) return(na_result)
  if (!any(c("Compound_ID", "HMDB_ID", "KEGG_ID") %in% colnames(topN_df))) return(na_result)

  clean_ids <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x)] <- ""
    x
  }

  n_top <- nrow(topN_df)
  hmdb <- if ("HMDB_ID" %in% colnames(topN_df)) {
    clean_ids(topN_df$HMDB_ID)
  } else {
    rep("", n_top)
  }
  kegg <- if ("KEGG_ID" %in% colnames(topN_df)) {
    clean_ids(topN_df$KEGG_ID)
  } else {
    rep("", n_top)
  }
  is_gen <- hmdb %in% hmdb_ids | kegg %in% kegg_ids
  n_gen  <- sum(is_gen)

  generic_names <- character(0)
  if ("Compound_name" %in% colnames(topN_df) && n_gen > 0L) {
    generic_names <- unique(trimws(as.character(topN_df$Compound_name[is_gen])))
    generic_names <- generic_names[!is.na(generic_names) & nzchar(generic_names)]
  }

  list(
    has_generic   = n_gen > 0L,
    generic_names = generic_names,
    generic_ids   = unique(c(hmdb[hmdb %in% hmdb_ids], kegg[kegg %in% kegg_ids])),
    n_generic     = as.integer(n_gen),
    n_top         = as.integer(n_top)
  )
}


# ============================================================
# New helper: compute_isolated_outliers
# ============================================================

#' Identify Isolated Outlier Metabolites in the Top Metabolite List
#'
#' A top metabolite is classified as an isolated outlier when it satisfies
#' ALL three conditions simultaneously:
#'   (1) Its HMDB_ID is absent from the member set of the top
#'       \code{top_n_pathways} HMDB pathways (by hit count) for this module.
#'   (2) It cannot be mapped via \code{module_df} to a KEGG ID present in
#'       the top \code{top_n_pathways} KEGG pathways.
#'   (3) Its chemical class is not the dominant (highest-priority) class
#'       reported in \code{other_metabolites}.
#'
#' Such metabolites have limited support from the current pathway/class
#' context. They may be noise, exposure-specific signals, or cross-pathway
#' biology, and should not dominate LLM annotation without corroboration.
#'
#' @param topN_df Data frame of top metabolites with at least Compound_ID,
#'   HMDB_ID, KEGG_ID, Compound_name, and class columns when available.
#' @param module_df Data frame of the full module containing HMDB_ID and
#'   KEGG_ID columns.  Used to look up KEGG IDs for topN members and to
#'   supply the full KEGG ID set for pathway ranking.  Pass NULL to skip
#'   KEGG pathway membership checks.
#' @param pathway_hmdb HMDB pathway database data frame (pathway_id,
#'   pathway_name, HMDB_ID, n_metabolites).  Pass NULL to skip.
#' @param pathway_kegg KEGG pathway database data frame (pathway_id,
#'   pathway_name, KEGG_ID, n_metabolites).  Pass NULL to skip.
#' @param other_metabolites Class-level summary data frame from
#'   prepare_module().  Used to determine the dominant chemical class.
#' @param top_n_pathways Integer. Number of top-hit pathways whose member
#'   sets define "covered" metabolites (default 3).
#'
#' @return A named list:
#'   \item{outlier_names}{Character vector of isolated metabolite names.}
#'   \item{outlier_ids}{Character vector of their Compound_IDs.}
#'   \item{n_outliers}{Integer.}
#'   \item{n_top}{Integer. Total size of topN.}
#'   \item{top_pathway_hmdb_ids}{All HMDB IDs covered by the top HMDB pathways.}
#'   \item{top_pathway_kegg_ids}{All KEGG IDs covered by the top KEGG pathways.}
compute_isolated_outliers <- function(
    topN_df,
    module_df         = NULL,
    pathway_hmdb      = NULL,
    pathway_kegg      = NULL,
    other_metabolites = NULL,
    top_n_pathways    = 3L
) {
  na_result <- list(
    outlier_names        = character(0),
    outlier_ids          = character(0),
    n_outliers           = 0L,
    n_top                = 0L,
    top_pathway_hmdb_ids = character(0),
    top_pathway_kegg_ids = character(0)
  )

  if (is.null(topN_df) || !is.data.frame(topN_df) || nrow(topN_df) == 0) return(na_result)
  if (!any(c("Compound_ID", "HMDB_ID", "KEGG_ID") %in% colnames(topN_df))) return(na_result)

  clean_ids <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x)] <- ""
    x
  }

  n_top <- nrow(topN_df)
  top_ids <- if ("Compound_ID" %in% colnames(topN_df)) {
    clean_ids(topN_df$Compound_ID)
  } else {
    rep("", n_top)
  }
  top_hmdb_ids <- if ("HMDB_ID" %in% colnames(topN_df)) {
    clean_ids(topN_df$HMDB_ID)
  } else {
    rep("", n_top)
  }
  top_kegg_ids <- if ("KEGG_ID" %in% colnames(topN_df)) {
    clean_ids(topN_df$KEGG_ID)
  } else {
    rep("", n_top)
  }

  # ---- Inner helper: IDs that appear in the top-N hit pathways ----
  get_top_pathway_members <- function(pathway_db, query_ids, id_col) {
    if (is.null(pathway_db) || !is.data.frame(pathway_db)) return(character(0))
    if (!id_col %in% colnames(pathway_db))                 return(character(0))
    if (!"pathway_id" %in% colnames(pathway_db))           return(character(0))
    if (length(query_ids) == 0)                            return(character(0))

    hits <- pathway_db[pathway_db[[id_col]] %in% query_ids, , drop = FALSE]
    if (nrow(hits) == 0) return(character(0))

    agg <- aggregate(
      as.formula(paste(id_col, "~ pathway_id")),
      data = hits,
      FUN  = function(x) length(unique(x))
    )
    colnames(agg)[2] <- "n_hit"
    agg <- agg[order(-agg$n_hit), , drop = FALSE]
    top_pids <- agg$pathway_id[seq_len(min(top_n_pathways, nrow(agg)))]

    mem <- pathway_db[[id_col]][pathway_db$pathway_id %in% top_pids]
    unique(as.character(mem[!is.na(mem) & nzchar(as.character(mem))]))
  }

  # ---- HMDB membership ----
  top_hmdb_members <- get_top_pathway_members(pathway_hmdb, top_hmdb_ids, "HMDB_ID")

  # ---- KEGG membership ----
  # Map topN HMDB IDs → KEGG IDs via module_df for per-metabolite check when
  # topN_df does not already carry KEGG_ID.
  # Use the full module KEGG ID set to rank which KEGG pathways are "top".
  top_kegg_for_top <- top_kegg_ids
  kegg_module_ids  <- character(0)

  if (!is.null(module_df) && is.data.frame(module_df) &&
      all(c("HMDB_ID", "KEGG_ID") %in% colnames(module_df))) {

    hmdb2kegg <- stats::setNames(
      as.character(module_df$KEGG_ID),
      as.character(module_df$HMDB_ID)
    )
    if (!any(nzchar(top_kegg_for_top))) {
      mapped <- hmdb2kegg[top_hmdb_ids]
      mapped[is.na(mapped)] <- ""
      top_kegg_for_top <- as.character(mapped)
    }

    kegg_module_ids <- unique(as.character(module_df$KEGG_ID[
      !is.na(module_df$KEGG_ID) & nzchar(as.character(module_df$KEGG_ID))]))
  }

  top_kegg_members <- get_top_pathway_members(pathway_kegg, kegg_module_ids, "KEGG_ID")

  # ---- Dominant class determination ----
  dominant_classes <- character(0)
  if (!is.null(other_metabolites) && is.data.frame(other_metabolites) &&
      "class" %in% colnames(other_metabolites)) {
    cls_df <- other_metabolites
    if ("priority_level" %in% colnames(cls_df)) {
      hi <- cls_df[cls_df$priority_level == "high", , drop = FALSE]
      if (nrow(hi) == 0) hi <- head(cls_df, 1L)
    } else {
      hi <- head(cls_df, 1L)
    }
    dominant_classes <- unique(trimws(as.character(hi$class)))
    dominant_classes <- dominant_classes[!is.na(dominant_classes) & nzchar(dominant_classes)]
  }

  # ---- Classify each topN metabolite ----
  in_hmdb_top  <- top_ids %in% top_hmdb_members
  in_kegg_top  <- nzchar(top_kegg_for_top) & (top_kegg_for_top %in% top_kegg_members)

  in_dom_class <- rep(FALSE, n_top)
  if (length(dominant_classes) > 0 && "class" %in% colnames(topN_df)) {
    met_class    <- trimws(as.character(topN_df$class))
    in_dom_class <- met_class %in% dominant_classes
  }

  is_isolated <- !in_hmdb_top & !in_kegg_top & !in_dom_class

  outlier_names <- character(0)
  if ("Compound_name" %in% colnames(topN_df) && any(is_isolated)) {
    outlier_names <- unique(trimws(as.character(topN_df$Compound_name[is_isolated])))
    outlier_names <- outlier_names[!is.na(outlier_names) & nzchar(outlier_names)]
  }

  list(
    outlier_names        = outlier_names,
    outlier_ids          = top_ids[is_isolated],
    n_outliers           = as.integer(sum(is_isolated)),
    n_top                = as.integer(n_top),
    top_pathway_hmdb_ids = top_hmdb_members,
    top_pathway_kegg_ids = top_kegg_members
  )
}


# ============================================================
# Updated: prepare_module  (adds pathway_coverage, global_direction, origin_profile)
# ============================================================

prepare_module <- function(
    k,
    module_list_top,
    module_list,
    module_pubmed_summaries,
    module_statistics,
    module_class,
    pathway_coverage_list = NULL,
    direction_list        = NULL,
    origin_profile_list   = NULL,
    pubmed       = TRUE,
    pathway_hmdb = NULL,
    pathway_kegg = NULL
) {
  make_lookup_keys <- function(module_id, fallback_key) {
    unique(as.character(stats::na.omit(c(
      module_id,
      fallback_key
    ))))
  }

  pick_named_entry <- function(x, keys, fallback = NULL) {
    if (is.null(x)) {
      return(NULL)
    }
    if (is.null(names(x))) {
      if (!is.null(fallback) && is.numeric(fallback) && fallback <= length(x)) {
        return(x[[fallback]])
      }
      return(NULL)
    }

    hit <- keys[keys %in% names(x)][1]
    if (is.na(hit) || !nzchar(hit)) {
      if (!is.null(fallback) && is.numeric(fallback) && fallback <= length(x)) {
        return(x[[fallback]])
      }
      return(NULL)
    }

    x[[hit]]
  }

  module_topN <- module_list_top[[k]]
  module_topN <- add_top_ids(module_topN)
  if (!"Compound_ID" %in% colnames(module_topN)) {
    if ("INCHI_ID" %in% colnames(module_topN)) {
      module_topN <- module_topN %>% dplyr::mutate(Compound_ID = INCHI_ID)
    } else {
      module_topN$Compound_ID <- NA_character_
    }
  }
  module_id <- unique(module_topN$Module)
  module_id <- as.character(module_id[1])
  lookup_keys <- make_lookup_keys(module_id, as.character(k))

  if (pubmed) {

    if (!is.null(module_pubmed_summaries[[k]]) &&
        length(module_pubmed_summaries[[k]]) > 0) {

      pubmed_df <- flatten_pubmed_for_one_module(module_pubmed_summaries[[k]])

      module_topN_merged <- combine_topN_with_pubmed(
        module_top_df = module_topN,
        pubmed_top_df = pubmed_df)
      
    } else {

      module_topN_merged <- module_topN %>%
        dplyr::select(
          dplyr::any_of(c(
            "top_id", "top_rank", "Compound_ID", "INCHI_ID", "HMDB_ID", "KEGG_ID",
            "Compound_name", "class", "priority_level"
          ))
        ) %>%
        mutate(functional_summary = NA_character_)
    }
    
  } else {
    module_topN_merged <- module_topN %>%
      dplyr::select(
        dplyr::any_of(c(
          "top_id", "top_rank", "Compound_ID", "INCHI_ID", "HMDB_ID", "KEGG_ID",
          "Compound_name", "class", "priority_level"
        ))
      ) %>%
      mutate(functional_summary = NA_character_)
  }

  # for other metabolites ----
  module_other <- NULL
  if (!is.null(module_class)) {
    module_other <- pick_named_entry(module_class, lookup_keys, fallback = k)
  }
  
  stats <- pick_named_entry(module_statistics, lookup_keys, fallback = k)

  # ---- New fields: pathway coverage, direction, origin ----
  pathway_cv <- pick_named_entry(pathway_coverage_list, lookup_keys, fallback = k)
  direction   <- pick_named_entry(direction_list,        lookup_keys, fallback = k)
  origin      <- pick_named_entry(origin_profile_list,   lookup_keys, fallback = k)

  # ---- Advanced-mode data-layer flags ----
  # generic_anchor_flags: detects cofactor/currency metabolite dominance in topN
  generic_flags <- compute_generic_anchor_flags(module_topN_merged)

  # isolation_flags: detects topN metabolites absent from top pathway hit sets
  # and not in the dominant chemical class (limited current context support)
  module_full_df <- pick_named_entry(module_list, lookup_keys, fallback = k)
  isolation_flags <- compute_isolated_outliers(
    topN_df           = module_topN_merged,
    module_df         = module_full_df,
    pathway_hmdb      = pathway_hmdb,
    pathway_kegg      = pathway_kegg,
    other_metabolites = module_other
  )

  list(
    module_id            = module_id,
    module_stats         = stats,
    topN_metabolites     = module_topN_merged,
    other_metabolites    = module_other,
    pathway_coverage     = pathway_cv,
    global_direction     = direction,
    origin_profile       = origin,
    generic_anchor_flags = generic_flags,
    isolation_flags      = isolation_flags
  )
}
