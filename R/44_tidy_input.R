
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
#'   HMDB_ID, Compound_name, Importance, Robustness, class
#' @param pubmed_top_df A data frame containing PubMed data with columns:
#'   HMDB_ID, functional_summary
#'
#' @return A tibble with columns ordered as: HMDB_ID, Compound_name, Importance,
#'   Robustness, class, functional_summary
#'
#' @details The left join ensures all metabolites from module_top_df are retained,
#'   even if they don't have corresponding PubMed data.
combine_topN_with_pubmed <- function(module_top_df, pubmed_top_df) {
  module_top_df %>%
    left_join(pubmed_top_df, by = "Compound_ID") %>%
    dplyr::select(
      Compound_ID,
      Compound_name,
      class,
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
#'   for a module. Each data frame must have columns: Module, HMDB_ID, Compound_name,
#'   Importance, Robustness, class
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
prepare_module <- function(
    k,
    module_list_top,
    module_list,
    module_pubmed_embed_refined,
    module_statistics,
    module_class,
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
          Compound_ID,
          Compound_name,
          class
        ) %>%
        mutate(functional_summary = NA_character_)
    }
    
  } else {
    module_topN_merged <- module_topN %>%
      dplyr::select(
        Compound_ID,
        Compound_name,
        class
      ) %>%
      mutate(functional_summary = NA_character_)
  }

  # for other metabolites ----
  module_other <- NULL
  if (!is.null(module_class)) {
    module_other <- pick_named_entry(module_class, lookup_keys, fallback = k)
  }
  
  stats <- pick_named_entry(module_statistics, lookup_keys, fallback = k)

  list(
    module_id          = module_id,
    module_stats       = stats,
    topN_metabolites   = module_topN_merged,
    other_metabolites  = module_other
  )
}
