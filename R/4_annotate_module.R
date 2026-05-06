#' Annotate functional modules using database, literature, and LLM models
#'
#' Perform functional annotation of statistically significant metabolic modules
#' by integrating database-derived annotations, optional PubMed literature
#' evidence, and large language model (LLM) functional interpretation.
#'
#' @param res A \code{featuremap_result} object returned by
#'   \code{\link{identify_module}} (or already module-annotated). Must contain
#'   \code{module_table} and \code{module_score}.
#' @param p Character. Column name in \code{res$module_score} used for module
#'   significance filtering (e.g. \code{"p"} or \code{"fdr"}).
#' @param pcutoff Numeric. Threshold applied to \code{res$module_score[[p]]}.
#' @param ref_data Functional reference metadata. Can be:
#'   \itemize{
#'     \item \code{NULL}: cached HMDB and KEGG references are loaded automatically.
#'     \item A single data frame with an \code{id} column (used for both lookups).
#'     \item A list with components \code{hmdb} and \code{kegg}.
#'   }
#' @param mode Character. \code{"fast"} (default) runs the two-call pipeline
#'   with pathway coverage, origin profiling, and directed confidence scoring.
#'   \code{"advanced"} runs the multi-agent debate pipeline: Proposer generates
#'   2-3 competing hypotheses, Skeptic challenges them (up to 2 rounds), and an
#'   Arbiter delivers final \code{accepted / uncertain / rejected} verdicts with
#'   confidence scores.
#' @param proposer_temperature Numeric. Sampling temperature for the Proposer's
#'   initial hypothesis generation in advanced mode (default 0.3). Skeptic,
#'   Proposer revision, and Arbiter always run at temperature 0.
#' @param sleep_sec Numeric. Pause in seconds between modules (default 0).
#' @param species Character. Species context for LLM confidence scoring.
#' @param biospecimen Character. Biospecimen type (e.g. plasma, urine, tissue).
#' @param phenotype Character. Phenotype/disease string used in PubMed retrieval
#'   (when \code{pubmed = TRUE}) and confidence assessment.
#' @param pubmed Logical. Whether to retrieve PubMed evidence after LLM Call 1.
#'   When \code{TRUE}, the module title and every top metabolite are searched
#'   separately; up to three papers are retained per query and deduplicated by
#'   PMID before LLM Call 2.
#' @param years Integer. PubMed look-back years when \code{pubmed = TRUE}.
#' @param provider Character. LLM provider (\code{"openai"}, \code{"siliconflow"},
#'   \code{"gemini"}).
#' @param llm_model Character. LLM model identifier.
#' @param NCBI_api_key Character. NCBI Entrez API key when \code{pubmed = TRUE}.
#' @param llm_api_key Character. API key for \code{provider}.
#'
#' @return The input \code{res} object updated with:
#' \itemize{
#'   \item \code{module_function}: data frame with one row per annotated module
#'         and columns covering all seven annotation dimensions (title,
#'         functional_summary, phenotype_relationship, confidence_label,
#'         rationale, global_direction, core_drivers, origin fractions,
#'         top pathway names, and PubMed evidence columns
#'         \code{pubmed_pmids} and \code{pubmed_evidence_json}.
#'   \item \code{module_top_metabolite_pmids}: compact PubMed evidence export
#'         with \code{module_id}, \code{pubmed_pmids}, and
#'         \code{pubmed_evidence_json}.
#'   \item \code{stage}: updated to \code{"module_annotated"}.
#' }
#'
#' @details
#' \strong{Fast Mode pipeline (default):}
#' \enumerate{
#'   \item Load functional reference databases (HMDB + KEGG).
#'   \item Annotate significant module metabolites via database lookup.
#'   \item Compute per-module: Signal_concentration (HHI), pathway coverage
#'         (KEGG + HMDB), global directionality, and origin profile.
#'   \item LLM Call 1 — title + functional_summary (no study context).
#'   \item Optional PubMed search using module title, all top metabolites,
#'         phenotype, and PMID deduplication.
#'   \item LLM Call 2 — phenotype_relationship + confidence_label + rationale.
#' }
#'
#' \strong{Advanced Mode pipeline:}
#' \enumerate{
#'   \item Same data preparation as fast mode (steps 1-3).
#'   \item P0 Proposer: generates 2-3 competing hypotheses from metabolite
#'         data only (no study context).
#'   \item S1 Skeptic: challenges each hypothesis on biochemical, origin,
#'         species, and biospecimen grounds (no phenotype, no PubMed).
#'   \item P1 Proposer revision: keep / modify / drop each hypothesis.
#'   \item (Conditional) S2 + P2: second debate round if revisions were made
#'         or major challenges were kept without revision.
#'   \item Optional PubMed: per-hypothesis title search + shared metabolite
#'         search, results held for Arbiter only.
#'   \item A0 Arbiter: final accepted / uncertain / rejected verdicts with
#'         full study context and PubMed evidence.
#' }
#'
#' @examples
#' \dontrun{
#' res <- annotate_module(
#'   res        = res,
#'   p          = "fdr",
#'   pcutoff    = 0.05,
#'   mode       = "fast",
#'   species    = "human",
#'   biospecimen = "plasma",
#'   phenotype  = "type 2 diabetes",
#'   pubmed     = FALSE,
#'   provider   = "openai",
#'   llm_model  = "gpt-4.1-mini",
#'   llm_api_key = Sys.getenv("OPENAI_API_KEY")
#' )
#' }
#'
#' @export
annotate_module <- function(res                  = NULL,
                             p                    = "fdr",
                             pcutoff              = 0.05,
                             ref_data             = NULL,
                             mode                 = c("fast", "advanced"),
                             species              = "human",
                             biospecimen          = NULL,
                             phenotype            = NULL,
                             pubmed               = TRUE,
                             years                = 5,
                             provider             = "openai",
                             llm_model            = "gpt-4.1-mini",
                             NCBI_api_key         = "",
                             llm_api_key          = "",
                             proposer_temperature = 0.3,
                             sleep_sec            = 0) {

  mode <- match.arg(mode)

  # ---- Inner helpers -------------------------------------------------------

  prepare_function_ref <- function(ref_df, label) {
    if (!is.data.frame(ref_df))
      stop("[Error] ", label, " functional reference must be a data.frame.")
    if (!"id" %in% colnames(ref_df))
      stop("[Error] ", label, " functional reference must contain an 'id' column.")
    ref_df <- tibble::as_tibble(ref_df)
    for (nm in setdiff(c("name", "class", "synonyms"), colnames(ref_df)))
      ref_df[[nm]] <- NA_character_
    ref_df
  }

  combine_function_refs <- function(module_df, hmdb_ref, kegg_ref) {
    hmdb_cols <- setdiff(colnames(hmdb_ref), "id")
    kegg_cols <- setdiff(colnames(kegg_ref), "id")
    out <- module_df %>%
      dplyr::left_join(
        hmdb_ref %>% dplyr::rename_with(~ paste0(.x, ".hmdb"), -id),
        by = c("HMDB_ID" = "id")
      ) %>%
      dplyr::left_join(
        kegg_ref %>% dplyr::rename_with(~ paste0(.x, ".kegg"), -id),
        by = c("KEGG_ID" = "id")
      )
    for (nm in intersect(hmdb_cols, kegg_cols))
      out[[nm]] <- dplyr::coalesce(out[[paste0(nm, ".hmdb")]], out[[paste0(nm, ".kegg")]])
    for (nm in setdiff(hmdb_cols, kegg_cols))
      out[[nm]] <- out[[paste0(nm, ".hmdb")]]
    for (nm in setdiff(kegg_cols, hmdb_cols))
      out[[nm]] <- out[[paste0(nm, ".kegg")]]
    out %>%
      dplyr::mutate(
        Compound_ID   = dplyr::coalesce(HMDB_ID, KEGG_ID),
        Compound_name = dplyr::coalesce(Compound_name, name)
      ) %>%
      dplyr::select(-dplyr::matches("\\.(hmdb|kegg)$"))
  }

  normalize_module_ids <- function(df) {
    if ("Module" %in% colnames(df)) df %>% dplyr::mutate(Module = as.character(Module)) else df
  }

  split_modules <- function(df) {
    split(df, factor(df$Module, levels = unique(df$Module)))
  }

  # ---- Validate input ------------------------------------------------------

  if (is.null(res)) stop("[Error] res cannot be NULL.")
  res <- ensure_featuremap_object(
    res, expected_stages = c("module_identified", "module_annotated")
  )
  res$module_table <- normalize_module_ids(res$module_table)
  res$module_score <- normalize_module_ids(res$module_score)

  cli::cli_h2("Step 4: Module Annotation  [mode = {mode}]")
  has_user_ref_data <- !is.null(ref_data)

  # ---- Step 1: Load reference databases ------------------------------------
  total_steps <- if (pubmed) 4L else 3L
  cli::cli_alert_info("[1/{total_steps}] Loading functional and pathway references.")

  if (is.null(ref_data)) {
    ref_data <- list(
      hmdb = featuremap_load_resource("functional_reference_hmdb"),
      kegg = featuremap_load_resource("functional_reference_kegg")
    )
  }
  if (is.data.frame(ref_data)) ref_data <- list(hmdb = ref_data, kegg = ref_data)
  if (!is.list(ref_data) || !all(c("hmdb", "kegg") %in% names(ref_data)))
    stop("[Error] ref_data must be NULL, a data.frame, or a list with 'hmdb' and 'kegg'.")
  ref_data$hmdb <- prepare_function_ref(ref_data$hmdb, "HMDB")
  ref_data$kegg <- prepare_function_ref(ref_data$kegg, "KEGG")

  # Load pathway databases (used for coverage, not for API calls)
  pathway_hmdb <- tryCatch(
    featuremap_load_resource("pathway_hmdb"),
    error = function(e) {
      cli::cli_alert_warning("pathway_hmdb resource not available: {e$message}")
      NULL
    }
  )
  pathway_kegg <- tryCatch(
    featuremap_load_resource("pathway_kegg"),
    error = function(e) {
      cli::cli_alert_warning("pathway_kegg resource not available: {e$message}")
      NULL
    }
  )

  if (has_user_ref_data) {
    cli::cli_alert_success("Using user-supplied functional references.")
  } else {
    cli::cli_alert_success("Loaded cached/remote functional and pathway references.")
  }

  start_time <- Sys.time()

  # ---- Step 2: Database annotation -----------------------------------------
  cli::cli_alert_info("[2/{total_steps}] Annotating significant modules.")

  module_levels <- unique(res$module_table$Module)
  sig_modules   <- res$module_score %>%
    dplyr::filter(.data[[p]] < pcutoff) %>%
    dplyr::pull(Module)

  module_table_for_annotation <- res$module_table
  if (!"Use_for_annotation" %in% colnames(module_table_for_annotation))
    module_table_for_annotation$Use_for_annotation <- TRUE
  module_table_for_annotation <- module_table_for_annotation %>%
    dplyr::mutate(Use_for_annotation = dplyr::coalesce(Use_for_annotation, TRUE))

  sig_module <- module_table_for_annotation %>%
    dplyr::mutate(
      .module_order = match(Module, module_levels),
      Robustness    = Node_freq,
      Compound_ID   = dplyr::coalesce(HMDB_ID, KEGG_ID)
    ) %>%
    dplyr::filter(Module %in% sig_modules, Use_for_annotation) %>%
    dplyr::arrange(.module_order, desc(Importance), desc(Robustness), Type) %>%
    dplyr::select(-.module_order) %>%
    combine_function_refs(hmdb_ref = ref_data$hmdb, kegg_ref = ref_data$kegg)

  if (is.null(sig_module) || nrow(sig_module) == 0) {
    cli::cli_alert_warning(
      "No modules passed the {.val {p}} < {.val {pcutoff}} threshold."
    )
    res$module_function            <- tibble::tibble()
    res$module_top_metabolite_pmids <- tibble::tibble()
    res <- set_featuremap_stage(res, "module_annotated")
    return(res)
  }

  module_list <- split_modules(sig_module)

  # ---- Top-N metabolite selection (GMM-based) ------------------------------
  module_list_top <- lapply(module_list, function(module) {
    module2 <- module %>%
      dplyr::mutate(weighted_score = 0.7 * Importance + 0.3 * Robustness) %>%
      dplyr::filter(is.finite(weighted_score)) %>%
      dplyr::arrange(dplyr::desc(weighted_score)) %>%
      dplyr::mutate(
        top_rank       = dplyr::row_number(),
        weighted_rank  = dplyr::percent_rank(weighted_score),
        priority_level = dplyr::case_when(
          weighted_rank >= 0.7 ~ "high",
          weighted_rank >= 0.3 ~ "medium",
          TRUE                 ~ "low"
        )
      )

    s <- module2$weighted_score
    n <- length(s)
    if (n == 0) return(module2)
    total <- sum(s)
    if (!is.finite(total) || total <= 0)
      return(dplyr::slice_head(module2, n = min(3L, n)))

    share <- s / total
    if (requireNamespace("mclust", quietly = TRUE) && n >= 8L) {
      fit <- tryCatch(
        suppressWarnings(mclust::Mclust(share, G = 1:3)),
        error = function(e) NULL
      )
      cls <- if (!is.null(fit)) fit$classification else NULL
      if (is.null(cls) || length(cls) != n || any(is.na(cls))) {
        k_cov <- which(cumsum(s) / total >= 0.20)[1]
        if (is.na(k_cov)) k_cov <- n
        topN  <- min(max(6L, k_cov), min(20L, n))
      } else if (max(cls) == 1L) {
        k_cov <- which(cumsum(s) / total >= 0.20)[1]
        if (is.na(k_cov)) k_cov <- n
        topN  <- min(max(5L, k_cov), min(20L, n))
      } else {
        means    <- tapply(share, cls, mean)
        head_comp <- which.max(means)
        topN     <- max(5L, min(20L, sum(cls == head_comp), n))
      }
    } else {
      k_cov <- which(cumsum(s) / total >= 0.20)[1]
      if (is.na(k_cov)) k_cov <- n
      topN  <- min(max(6L, k_cov), min(20L, n))
    }

    dplyr::slice_head(module2, n = topN)
  })

  # ---- Class-level structure (all metabolites, for prompt) -----------------
  module_class <- purrr::imap(module_list, ~ {
    if (nrow(.x) == 0) return(tibble::tibble())
    .x %>%
      dplyr::group_by(class) %>%
      dplyr::summarise(
        n_metabolites  = round(n() / nrow(.x), 2),
        metabolites    = paste(unique(Compound_name), collapse = "{}"),
        Compound_IDs   = paste(unique(Compound_ID),   collapse = "{}"),
        mean_priority  = mean(0.7 * Importance + 0.3 * Robustness, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        priority_rank  = dplyr::percent_rank(mean_priority),
        priority_level = dplyr::case_when(
          priority_rank >= 0.7 ~ "high",
          priority_rank >= 0.3 ~ "medium",
          TRUE                 ~ "low"
        )
      ) %>%
      dplyr::select(class, n_metabolites, metabolites, Compound_IDs, priority_level)
  })

  # ---- Module-level statistics: Signal_concentration (HHI) only -----------
  # Computed across ALL modules for relative calibration, then filtered to sig.
  module_statistics <- {
    full_mod_tbl <- module_table_for_annotation %>%
      dplyr::filter(Use_for_annotation) %>%
      dplyr::mutate(Robustness = Node_freq)

    tmp_all <- purrr::imap_dfr(split_modules(full_mod_tbl), function(module, name) {
      imp <- module$Importance[is.finite(module$Importance)]
      hhi <- if (length(imp) == 0) NA_real_ else {
        shr <- imp / sum(imp); sum(shr^2)
      }
      tibble::tibble(Module = name, HHI = hhi)
    })

    q_hhi <- quantile(tmp_all$HHI, probs = c(0.3, 0.7), na.rm = TRUE)
    level3 <- function(x, qs) dplyr::case_when(
      is.na(x)    ~ "NA",
      x >= qs[2]  ~ "high",
      x >= qs[1]  ~ "medium",
      TRUE        ~ "low"
    )

    tmp_sig <- tmp_all %>%
      dplyr::filter(Module %in% names(module_list)) %>%
      dplyr::mutate(Signal_concentration = level3(HHI, q_hhi))

    split_modules(tmp_sig) %>%
      purrr::map(~ list(Signal_concentration = .x$Signal_concentration))
  }

  # ---- Step 3: Per-module computed features --------------------------------
  cli::cli_alert_info("[3/{total_steps}] Computing pathway coverage, direction, and origin profiles.")

  # Check if feature_table and anno_feature_table2 are available for direction
  has_direction_data <-
    !is.null(res$feature_table)       && is.data.frame(res$feature_table) &&
    "condition" %in% colnames(res$feature_table) &&
    !is.null(res$anno_feature_table2) && is.data.frame(res$anno_feature_table2)

  pathway_coverage_list <- vector("list", length(module_list))
  direction_list        <- vector("list", length(module_list))
  origin_profile_list   <- vector("list", length(module_list))
  names(pathway_coverage_list) <- names(module_list)
  names(direction_list)        <- names(module_list)
  names(origin_profile_list)   <- names(module_list)

  for (mod_name in names(module_list)) {
    mod_df <- module_list[[mod_name]]

    pathway_coverage_list[[mod_name]] <- compute_pathway_coverage(
      module_df    = mod_df,
      pathway_hmdb = pathway_hmdb,
      pathway_kegg = pathway_kegg,
      top_n        = 3L
    )

    direction_list[[mod_name]] <- if (has_direction_data) {
      compute_global_direction(
        module_df          = mod_df,
        feature_table      = res$feature_table,
        anno_feature_table2 = res$anno_feature_table2
      )
    } else {
      list(direction_label = NA_character_,
           weighted_mean_condition = NA_real_, n_features = 0L)
    }

    origin_profile_list[[mod_name]] <- compute_origin_profile(
      module_df         = mod_df,
      hmdb_function_ref = ref_data$hmdb
    )
  }

  n_pathways_found <- sum(vapply(pathway_coverage_list, function(pc) {
    nrow(pc$kegg %||% tibble::tibble()) + nrow(pc$hmdb %||% tibble::tibble())
  }, integer(1)))
  cli::cli_alert_success(
    "Pathway coverage computed; {.val {n_pathways_found}} pathway hits across all modules."
  )

  # ---- PubMed info (per-module search is deferred to LLM annotation) -------
  res$module_top_metabolite_pmids <- tibble::tibble()

  # ---- Prepare module objects for LLM annotation ---------------------------
  all_modules_ready <- lapply(
    seq_along(module_list_top),
    function(i) {
      prepare_module(
        k                     = i,
        module_list_top       = module_list_top,
        module_list           = module_list,
        module_class          = module_class,
        module_pubmed_summaries = NULL,
        module_statistics     = module_statistics,
        pathway_coverage_list = pathway_coverage_list,
        direction_list        = direction_list,
        origin_profile_list   = origin_profile_list,
        pubmed                = FALSE         # pubmed handled inside annotate_one_module
      )
    }
  )
  names(all_modules_ready) <- names(module_list_top)

  study_context <- list(
    species     = species,
    biospecimen = biospecimen,
    phenotype   = phenotype
  )

  cli::cli_alert_success(
    "Prepared annotation inputs for {.val {length(all_modules_ready)}} significant modules."
  )

  # ---- Step 3/4: LLM annotation ----------------------------------------
  llm_step <- if (pubmed) "4" else "3"
  cli::cli_alert_info(
    "[{llm_step}/{total_steps}] Running LLM annotation (mode = {mode})."
  )

  if (length(all_modules_ready) == 0) {
    cli::cli_alert_warning("No modules available for LLM annotation.")
    res$module_function <- tibble::tibble()
    res <- set_featuremap_stage(res, "module_annotated")
    return(res)
  }

  if (pubmed && nzchar(NCBI_api_key)) {
    rentrez::set_entrez_key(NCBI_api_key)
  }

  all_results <- if (identical(mode, "advanced")) {
    annotate_all_modules_advanced(
      all_modules_ready    = all_modules_ready,
      study_context        = study_context,
      provider             = provider,
      llm_api_key          = llm_api_key,
      model                = llm_model,
      max_tokens           = if (identical(provider, "gemini")) 4000L else 2000L,
      proposer_temperature = proposer_temperature,
      temperature          = 0,
      sleep_sec            = sleep_sec,
      pubmed               = pubmed,
      NCBI_api_key         = NCBI_api_key,
      years                = years
    )
  } else {
    annotate_all_modules(
      all_modules_ready = all_modules_ready,
      study_context     = study_context,
      provider          = provider,
      llm_api_key       = llm_api_key,
      model             = llm_model,
      max_tokens        = if (identical(provider, "gemini")) 4000L else 2000L,
      temperature       = 0,
      sleep_sec         = sleep_sec,
      pubmed            = pubmed,
      NCBI_api_key      = NCBI_api_key,
      years             = years
    )
  }

  res$module_function <- all_results
  if ("pubmed_evidence_json" %in% colnames(all_results)) {
    res$module_top_metabolite_pmids <-
      all_results[, c("module_id", "pubmed_pmids", "pubmed_evidence_json"), drop = FALSE]
  }

  elapsed_min <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2)
  cli::cli_alert_success(
    "Annotated {.val {nrow(res$module_function)}} modules in {.val {elapsed_min}} minutes."
  )
  res <- set_featuremap_stage(res, "module_annotated")
  return(res)
}
