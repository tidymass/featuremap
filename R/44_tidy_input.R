
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
#'   top_rank, Compound_ID, Compound_name, class, priority_level
#' @param pubmed_top_df A data frame containing PubMed data with columns:
#'   HMDB_ID, functional_summary
#'
#' @return A tibble with columns ordered as: top_rank, Compound_ID, Compound_name,
#'   class, priority_level, functional_summary
#'
#' @details The left join ensures all metabolites from module_top_df are retained,
#'   even if they don't have corresponding PubMed data.
combine_topN_with_pubmed <- function(module_top_df, pubmed_top_df) {
  module_top_df %>%
    left_join(pubmed_top_df, by = "Compound_ID") %>%
    dplyr::select(
      top_rank,
      Compound_ID,
      Compound_name,
      class,
      priority_level,
      functional_summary
    )
}



#' Prepare Comprehensive Module Data for Analysis or Reporting
#'
#' Combines top metabolites, PubMed functional summaries, other metabolites,
#' and module statistics into a structured list for a single module.
#'
#' @param k Integer or character. Index or identifier for the module to prepare
#' @param module_list_top List of data frames, each containing top 10 metabolites
#'   for a module. Each data frame must have columns including Module, Compound_ID
#'   (or HMDB_ID / KEGG_ID), Compound_name, class, top_rank, and priority_level
#' @param module_list Named list of data frames containing all metabolites for each module.
#'   Names should match the module IDs returned by \code{\link{identify_module}}
#'   (for example \code{"M1"}, \code{"M2"}, ...)
#' @param module_pubmed_embed_refined List of PubMed data for each module. Each element
#'   is a named list where names are HMDB IDs
#' @param module_statistics List of statistics for each module. Names should match
#'   the character module IDs returned by \code{\link{identify_module}}
#' @param module_class Optional list of class-level summaries for each module.
#' @param pubmed Logical. If TRUE (default), attempts to merge PubMed functional summaries.
#'   If FALSE, functional_summary column will be NA
#'
#' @return A named list with four elements:
#'   \item{module_id}{Character. The module identifier (for example \code{"M1"})}
#'   \item{module_stats}{Statistics data for the module}
#'   \item{topN_metabolites}{Tibble of top 10 metabolites with functional summaries (if available)}
#'   \item{other_metabolites}{Tibble of remaining metabolites not in top 10}
#'
#' @details This function orchestrates the preparation of module data by:
#'   \enumerate{
#'     \item Extracting the top 10 metabolites for the specified module
#'     \item Optionally merging PubMed functional summaries (with warning if data missing)
#'     \item Extracting metabolites not in the top 10
#'     \item Retrieving module statistics
#'     \item Packaging everything into a structured list
#'   }
#'
#' @section Warning:
#'   If pubmed = TRUE but PubMed data is missing or empty for the module,
#'   a warning is issued and functional_summary will be set to NA for all metabolites.
#'
#' @examples
#' \dontrun{
#' # Prepare the first module entry with PubMed data
#' module_data <- prepare_module(
#'   k = 1,
#'   module_list_top = topN_list,
#'   module_list = all_modules,
#'   module_pubmed_embed_refined = pubmed_data,
#'   module_statistics = stats_list,
#'   pubmed = TRUE
#' )
#'
#' # Access components
#' module_data$module_id
#' module_data$topN_metabolites
#' module_data$other_metabolites
#' }


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
#'   KEGG_ID, plus a score column (weighted_anno_score preferred, else
#'   score_norm).
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

  # Detect the annotation score column
  score_col <- NULL
  for (sc in c("weighted_anno_score", "score_norm")) {
    if (sc %in% colnames(anno_feature_table2)) { score_col <- sc; break }
  }
  if (is.null(score_col)) return(na_result)

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
  anno_sub <- anno_sub[, c("variable_id", score_col), drop = FALSE]
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
#'   tibble with columns pathway_name, n_hit, n_metabolites, coverage.
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
      formula = as.formula(paste(id_col, "~ pathway_id + pathway_name")),
      data    = hits,
      FUN     = function(x) length(unique(x))
    )
    colnames(agg)[colnames(agg) == id_col] <- "n_hit"

    agg <- merge(agg, pathway_sizes, by = "pathway_id", all.x = TRUE)
    agg$coverage <- round(agg$n_hit / pmax(agg$n_metabolites, 1L), 3)
    agg <- agg[order(-agg$n_hit, -agg$coverage), , drop = FALSE]
    agg <- agg[seq_len(min(top_n, nrow(agg))), , drop = FALSE]
    tibble::as_tibble(agg[, c("pathway_name", "n_hit", "n_metabolites", "coverage"),
                          drop = FALSE])
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
# Updated: prepare_module  (adds pathway_coverage, global_direction, origin_profile)
# ============================================================

prepare_module <- function(
    k,
    module_list_top,
    module_list,
    module_pubmed_embed_refined,
    module_statistics,
    module_class,
    pathway_coverage_list = NULL,
    direction_list        = NULL,
    origin_profile_list   = NULL,
    pubmed = TRUE
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
  if (!"Compound_ID" %in% colnames(module_topN)) {
    if ("HMDB_ID" %in% colnames(module_topN)) {
      module_topN <- module_topN %>% dplyr::rename(Compound_ID = HMDB_ID)
    } else if ("KEGG_ID" %in% colnames(module_topN)) {
      module_topN <- module_topN %>% dplyr::rename(Compound_ID = KEGG_ID)
    } else {
      colnames(module_topN)[3] = "Compound_ID"
    }
  }
  module_id <- unique(module_topN$Module)
  module_id <- as.character(module_id[1])
  lookup_keys <- make_lookup_keys(module_id, as.character(k))

  if (pubmed) {

    if (!is.null(module_pubmed_embed_refined[[k]]) &&
        length(module_pubmed_embed_refined[[k]]) > 0) {

      pubmed_df <- flatten_pubmed_for_one_module(module_pubmed_embed_refined[[k]])

      module_topN_merged <- combine_topN_with_pubmed(
        module_top_df = module_topN,
        pubmed_top_df = pubmed_df)
      
    } else {

      module_topN_merged <- module_topN %>%
        dplyr::select(
          top_rank,
          Compound_ID,
          Compound_name,
          class,
          priority_level
        ) %>%
        mutate(functional_summary = NA_character_)
    }
    
  } else {
    module_topN_merged <- module_topN %>%
      dplyr::select(
        top_rank,
        Compound_ID,
        Compound_name,
        class,
        priority_level
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

  list(
    module_id          = module_id,
    module_stats       = stats,
    topN_metabolites   = module_topN_merged,
    other_metabolites  = module_other,
    pathway_coverage   = pathway_cv,
    global_direction   = direction,
    origin_profile     = origin
  )
}
