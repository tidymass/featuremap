#' Extract a sub-network using the PCSF algorithm
#'
#' Construct a prize-collecting Steiner forest (PCSF)-based metabolic
#' sub-network from annotated metabolomics features and a reference
#' metabolic reaction network (MRN).
#'
#' This function integrates feature-level statistical scores with a
#' metabolic reaction network and identifies a robust sub-network
#' through repeated PCSF optimization. Optional grid tuning of
#' parameters `b` and `w` is supported.
#'
#' @param res A \code{featuremap_result} object returned by
#'   \code{\link{annotate_feature}} (or a later stage). It must contain at
#'   least \code{feature_table}, \code{anno_feature_table1}, and \code{refdata},
#'   and should be at stage \code{"feature_annotated"} or later.
#' @param feature Character. Feature selection mode:
#'   \itemize{
#'   \item \code{"all"}: use all rows in \code{res$feature_table}
#'   \item \code{"sig"}: use only rows with \code{pvalue < pcutoff}
#'   }
#' @param pcutoff Numeric. Significance threshold used only when
#'   \code{feature = "sig"} to filter \code{res$feature_table$pvalue}.
#' @param mrn A \code{tidygraph::tbl_graph} object for the metabolic reaction
#'   network. If \code{NULL}, a cached/default MRN resource will be loaded
#'   based on \code{res$refdata}.
#' @param workers Integer. Number of parallel workers used in repeated PCSF runs.
#' @param times Integer. Number of PCSF repeats for final subnetwork construction.
#'   Larger values usually improve robustness but increase runtime.
#' @param b Numeric. Prize scaling factor for terminal nodes in PCSF.
#' @param mu Numeric. Hub penalty coefficient. Non-terminal nodes receive a
#'   negative prize proportional to their network degree; larger values penalize
#'   high-degree nodes more strongly.
#' @param w Numeric. Dummy edge cost in PCSF.
#' @param r Numeric. Relative Gaussian edge-cost noise applied across repeated
#'   PCSF runs. For example, \code{r = 0.20} perturbs each edge cost by a
#'   multiplicative factor sampled from \code{N(1, 0.20)}.
#' @param grid_tune Logical. Whether to tune \code{b} and \code{w} on a grid
#'   before the final run.
#' @param b_grid Numeric vector of candidate \code{b} values when
#'   \code{grid_tune = TRUE}.
#' @param w_grid Numeric vector of candidate \code{w} values when
#'   \code{grid_tune = TRUE}.
#' @param times_tune Integer. Number of repeats per \code{(b, w)} pair during
#'   grid tuning.
#'
#' @return The input \code{res} object updated in-place with:
#' \itemize{
#'   \item \code{feature_table}: with derived \code{importance} score columns
#'   \item \code{pcsf_subnet}: final PCSF-derived subnetwork (\code{igraph})
#'   \item \code{node_frequency}: node appearance frequency across repeats
#'   \item \code{edge_frequency}: edge appearance frequency across repeats
#'   \item \code{anno_feature_table2}: network-supported feature-metabolite
#'         annotation table including \code{INCHI_ID}, \code{edge_freq},
#'         \code{weighted_anno_score}, and feature-level \code{importance}
#'   \item \code{tuning_grid} and \code{tuning_best}: (if \code{grid_tune = TRUE})
#'         tuning summary
#'   \item \code{stage}: updated to \code{"subnetwork_extracted"}
#' }
#'
#' @details
#' Annotation rows are projected to the MRN at the metabolite-node level using
#' \code{HMDB.ID} first and \code{KEGG.ID} as fallback when no HMDB-based match
#' is available. Only annotation rows supported by the final PCSF subnetwork are
#' retained in \code{anno_feature_table2}.
#'
#' Feature importance scores are derived from effects value and p-value.
#' The PCSF algorithm is executed multiple times with relative Gaussian
#' perturbation of edge costs to estimate robustness. Final subnetworks are
#' constructed based on node and edge frequency across runs.
#' Network-supported annotation confidence combines normalized annotation score
#' and edge frequency with a global adaptive edge-frequency weight bounded
#' between 0.1 and 0.3.
#'
#' @examples
#' \dontrun{
#' # res should be generated from annotate_feature() first
#' res <- extract_subnetwork(
#'   res = res,
#'   feature = "sig",
#'   pcutoff = 0.05,
#'   mrn = NULL,
#'   workers = 4,
#'   times = 100,
#'   b = 20,
#'   mu = 1e-4,
#'   w = 25,
#'   r = 0.20,
#'   grid_tune = FALSE
#' )
#' }
#'
#' @export

extract_subnetwork = function(res = NULL,
                              feature = c("all", "sig"),
                              pcutoff = 0.05,
                              mrn = NULL,
                              workers = 4,
                              times = 100,
                              b = 20,
                              mu = 1e-4,
                              w = 25,
                              r = 0.20,
                              grid_tune = FALSE,
                              b_grid = seq(20, 40, by = 5),
                              w_grid = seq(20, 40, by = 5),
                              times_tune = 10){
  feature <- match.arg(feature)

  if (is.null(res)){
    stop("[Error] res cannot be NULL.")
  }
  res <- ensure_featuremap_object(
    res,
    expected_stages = c("feature_annotated", "subnetwork_extracted", "module_identified", "module_annotated")
  )
  if (!"feature_table" %in% names(res)){
    stop("[Error] res$feature_table is missing.")
  }
  cli::cli_h2("Step 2: Subnetwork Extraction")
  has_user_mrn <- !is.null(mrn)
  start_time <- Sys.time()
  total_steps <- if (isTRUE(grid_tune)) 4 else 3
  step_label <- function(step_id) sprintf("[%d/%d]", step_id, total_steps)

  cli::cli_alert_info("{step_label(1)} Loading the metabolic reaction network.")
  if (is.null(mrn)){
    ref_tag <- if (!is.null(res$refdata)) as.character(res$refdata) else NULL
    resource_name <- featuremap_default_resource_name("mrn", refdata = ref_tag)
    mrn <- featuremap_load_resource(resource_name)
  }
  if (has_user_mrn) {
    cli::cli_alert_success("Using a user-supplied metabolic reaction network.")
  } else {
    cli::cli_alert_success("No user-supplied metabolic reaction network detected; loaded cached or remote resources.")
  }

  res$mrn = mrn
  res$feature_table %>%
    dplyr::mutate(
      score_fc = {
        fc <- abs(condition)
        (fc - min(fc, na.rm = TRUE)) /
          (max(fc, na.rm = TRUE) - min(fc, na.rm = TRUE))
      },
      score_p = {
        sp <- -log10(pvalue)
        (sp - min(sp, na.rm = TRUE)) /
          (max(sp, na.rm = TRUE) - min(sp, na.rm = TRUE))
      },
      importance = sqrt(score_fc * score_p)
    ) %>%
    as.data.frame() -> res$feature_table


  ## prepare data for pcsf analysis ---
  if(feature == "all"){
    feature_sig <- res$feature_table
  }else if(feature == "sig"){
    feature_sig <- res$feature_table %>%
      dplyr::filter(pvalue < pcutoff)
    if (nrow(feature_sig) == 0)
      stop("[Error] No significant features found under pcutoff = ", pcutoff)
  }
  cli::cli_alert_info("{step_label(2)} Preparing subnetwork extraction.")

  ## build FMpair + mrnref ----
  anno_feature_table1_sub <- res$anno_feature_table1 %>%
    dplyr::filter(variable_id %in% feature_sig$variable_id) %>%
    dplyr::mutate(annotation_row_id = dplyr::row_number()) %>%
    as.data.frame()

  anno_feature_table1_sub %>%
    dplyr::select(
      annotation_row_id,
      variable_id,
      Compound.name,
      HMDB.ID,
      KEGG.ID,
      score_norm
    ) -> anno_FM

  mrn_nodes <- mrn %>%
    tidygraph::activate(nodes) %>%
    as.data.frame()

  mrnref_hmdb <- mrn_nodes %>%
    dplyr::filter(!is.na(HMDB_ID), HMDB_ID != "", HMDB_ID != "<NA>") %>%
    dplyr::select(HMDB_ID, INCHI_ID, Compound_name) %>%
    dplyr::distinct(HMDB_ID, .keep_all = TRUE)

  mrnref_kegg <- mrn_nodes %>%
    dplyr::filter(!is.na(KEGG_ID), KEGG_ID != "", KEGG_ID != "<NA>") %>%
    dplyr::select(KEGG_ID, INCHI_ID, Compound_name) %>%
    dplyr::distinct(KEGG_ID, .keep_all = TRUE)

  suppressWarnings({
    anno_hmdb <- anno_FM %>%
      dplyr::filter(!is.na(HMDB.ID), HMDB.ID != "") %>%
      dplyr::left_join(mrnref_hmdb, by = c("HMDB.ID" = "HMDB_ID")) %>%
      dplyr::filter(!is.na(INCHI_ID)) %>%
      dplyr::arrange(annotation_row_id, INCHI_ID) %>%
      dplyr::distinct(annotation_row_id, .keep_all = TRUE) %>%
      dplyr::mutate(mapping_source = "HMDB")

    hmdb_mapped_rows <- anno_hmdb$annotation_row_id

    anno_kegg <- anno_FM %>%
      dplyr::filter(!annotation_row_id %in% hmdb_mapped_rows) %>%
      dplyr::filter(!is.na(KEGG.ID), KEGG.ID != "") %>%
      dplyr::left_join(mrnref_kegg, by = c("KEGG.ID" = "KEGG_ID")) %>%
      dplyr::filter(!is.na(INCHI_ID)) %>%
      dplyr::arrange(annotation_row_id, INCHI_ID) %>%
      dplyr::distinct(annotation_row_id, .keep_all = TRUE) %>%
      dplyr::mutate(mapping_source = "KEGG")
  })

  anno_mapped <- dplyr::bind_rows(anno_hmdb, anno_kegg) %>%
    dplyr::arrange(variable_id, dplyr::desc(score_norm), mapping_source, INCHI_ID)

  mapped_feature_n <- dplyr::n_distinct(anno_mapped$variable_id)
  cli::cli_alert_success(
    "Mapped {.val {mapped_feature_n}} features with {.val {nrow(anno_mapped)}} annotation rows to the metabolic reaction network."
  )

  c_min = 0.02
  FMpair <- anno_mapped %>%
    dplyr::mutate(cost = 1 - score_norm + c_min) %>%
    dplyr::select(variable_id, INCHI_ID, cost) %>%
    dplyr::distinct(variable_id, INCHI_ID, .keep_all = TRUE) %>%
    purrr::set_names(c("from", "to", "weight"))

  if (nrow(FMpair) == 0){
    stop("[Error] No matched feature-metabolite pairs (FMpair); check annotation table or ID format.")
  }

  mrn %>%
    tidygraph::activate(edges) %>%
    as.data.frame() %>%
    dplyr::select(INCHI_From, INCHI_To) %>%
    purrr::set_names(c("from", "to")) %>%
    dplyr::mutate(weight = min(max(FMpair$weight), 0.99)) -> MMpair

  imrn = igraph::graph_from_data_frame(dplyr::bind_rows(FMpair, MMpair), directed = FALSE)

  ## prepare node weights (parts independent of b) ---
  terminals <- feature_sig %>%
    dplyr::filter(variable_id %in% FMpair$from) %>%
    dplyr::select(variable_id, importance) %>%
    tibble::deframe()

  node_names = igraph::V(imrn)$name
  node_prz <- tibble::tibble(name = node_names) %>%
    dplyr::left_join(
      tibble::enframe(terminals, name = "name", value = "prize"),
      by = "name"
    ) %>%
    tidyr::replace_na(list(prize = 0)) %>%
    dplyr::pull(prize)

  hub_penalty <- tibble::tibble(
    name = node_names,
    degree = igraph::degree(imrn)
  ) %>%
    dplyr::mutate(penalty = -mu * degree) %>%
    dplyr::pull(penalty)

  ## set dummy edge endpoints (independent of b, w) ---
  dummies = names(terminals)
  edges_df <- igraph::ends(imrn, es = igraph::E(imrn)) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rlang::set_names(c("from", "to"))
  final_edges <- dplyr::bind_rows(
    tibble::tibble(from = "DUMMY", to = dummies),
    edges_df
  )
  from <- final_edges$from
  to <- final_edges$to

  # prepare parallel computation ----
  if (workers > 1) {
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)
  }

  # define run pcsf once function ---
  run_pcsf_once <- function(b_val, w_val, times_val, show_progress = TRUE){

    node_prizes <- tibble::tibble(
      name = node_names,
      base_prize = b_val * node_prz,
      penalty = hub_penalty
    ) %>%
      dplyr::mutate(final_prize = dplyr::if_else(base_prize == 0, penalty, base_prize)) %>%
      dplyr::pull(final_prize)

    if(workers > 1){
      results_list <- furrr::future_map(
        1:times_val,
        function(i) {
          cost <- c(
            rep(w_val, length(dummies)),
            pmax(
              igraph::E(imrn)$weight *
                stats::rnorm(length(igraph::E(imrn)), mean = 1, sd = r),
              1e-4
            )
          )
          output <- pcsf2featuremap::call_sr(from, to, cost, node_names, node_prizes)
          edge <- data.frame(
            source = as.character(output[[1]]),
            target = as.character(output[[2]]),
            stringsAsFactors = FALSE
          )
          edge <- edge[edge$source != "DUMMY" & edge$target != "DUMMY", ]

          if (nrow(edge) > 0) {
            g <- igraph::graph_from_data_frame(edge, directed = FALSE)
            nodes <- igraph::V(g)$name
            return(list(graph = g, nodes = nodes, iter = i))
          } else {
            return(list(graph = NULL, nodes = NULL, iter = i))
          }
        },
        .options = furrr::furrr_options(
          seed = NULL,
          packages = c("igraph", "dplyr", "tibble", "tidyr"),
          globals = list(
            imrn = imrn, dummies = dummies, from = from, to = to, r = r,
            node_names = node_names, node_prz = node_prz, hub_penalty = hub_penalty
          )
        ),
        .progress = show_progress
      )
    }else{
      results_list <- purrr::map(
        1:times_val,
        function(i) {
          cost <- c(
            rep(w_val, length(dummies)),
            pmax(
              igraph::E(imrn)$weight *
                stats::rnorm(length(igraph::E(imrn)), mean = 1, sd = r),
              1e-4
            )
          )
          output <- pcsf2featuremap::call_sr(from, to, cost, node_names, node_prizes)
          edge <- data.frame(
            source = as.character(output[[1]]),
            target = as.character(output[[2]]),
            stringsAsFactors = FALSE
          )
          edge <- edge[edge$source != "DUMMY" & edge$target != "DUMMY", ]

          if (nrow(edge) > 0) {
            g <- igraph::graph_from_data_frame(edge, directed = FALSE)
            nodes <- igraph::V(g)$name
            return(list(graph = g, nodes = nodes, iter = i))
          } else {
            return(list(graph = NULL, nodes = NULL, iter = i))
          }
        }
      )
    }

    all_nodes <- unlist(lapply(results_list, function(x) x$nodes))
    graphs_list <- lapply(results_list, function(x) x$graph)
    names(graphs_list) <- paste0("graph_", seq_along(graphs_list))

    if (length(all_nodes) == 0) {
      stop("No nodes found across runs (all PCSF outputs empty).")
    }

    node_freq <- table(all_nodes) %>% sort(decreasing = TRUE)
    node_frequency <- tibble::tibble(
      name = names(node_freq),
      freq = as.numeric(node_freq)/times_val
    ) %>%
      dplyr::mutate(type = dplyr::case_when(
        stringr::str_detect(name, "InChI") ~ "metabolite",
        TRUE ~ "feature"
      ))

    all_edges <- purrr::map_dfr(seq_along(graphs_list), function(i) {
      g <- graphs_list[[i]]
      if (is.null(g) || igraph::gsize(g) == 0) {
        return(tibble::tibble(from = character(), to = character(), run = integer()))
      }
      igraph::ends(g, es = igraph::E(g)) %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        rlang::set_names(c("from", "to")) %>%
        dplyr::mutate(run = i)
    })

    edge_frequency <- all_edges %>%
      dplyr::mutate(edge = dplyr::if_else(from < to,
                                          paste(from, to, sep = "|"),
                                          paste(to, from, sep = "|"))) %>%
      dplyr::group_by(edge) %>%
      dplyr::tally(name = "freq") %>%
      tidyr::separate(edge, into = c("from", "to"), sep = "\\|") %>%
      dplyr::mutate(freq = freq/times_val)

    if (nrow(edge_frequency) == 0) {
      stop("No PCSF-based sub-network could be identified under current (b,w).")
    }

    pcsf_subnet <- igraph::graph_from_data_frame(edge_frequency, directed = FALSE)
    igraph::E(pcsf_subnet)$weight <- edge_frequency$freq

    pcsf_subnet <- pcsf_subnet %>%
      igraph::set_vertex_attr("prize", value = node_frequency$freq[match(igraph::V(.)$name, node_frequency$name)]) %>%
      igraph::set_vertex_attr("type", value = "Steiner")

    term_idx <- match(names(terminals), igraph::V(pcsf_subnet)$name)
    term_idx <- term_idx[!is.na(term_idx)]
    igraph::V(pcsf_subnet)$type[term_idx] <- "Terminal"
    class(pcsf_subnet) <- c("PCSF", "igraph")

    ## tidy pcsf results -> anno_feature_table2 ---
    feature_nodes <- igraph::V(pcsf_subnet)[igraph::V(pcsf_subnet)$type == "Terminal"]$name

    feature2metabolite <- list()
    for (f in feature_nodes) {
      neigh <- igraph::neighbors(pcsf_subnet, f, mode = "all")$name
      feature2metabolite[[f]] <- neigh
    }

    feature_metabolite_df <- data.frame(
      feature = rep(names(feature2metabolite), lengths(feature2metabolite)),
      metabolite = unlist(feature2metabolite)
    ) %>%
      dplyr::mutate(network_supported = TRUE) %>%
      dplyr::left_join(
        mrn_nodes %>%
          dplyr::select(INCHI_ID, Compound_name, HMDB_ID, KEGG_ID) %>%
          dplyr::distinct(INCHI_ID, .keep_all = TRUE),
        by = c("metabolite" = "INCHI_ID")
      )
    rownames(feature_metabolite_df) = NULL

    anno_feature_table2_extra <- anno_mapped %>%
      dplyr::rename(HMDB_ID = HMDB.ID, KEGG_ID = KEGG.ID) %>%
      dplyr::select(
        annotation_row_id,
        variable_id,
        INCHI_ID,
        score_norm,
        mapping_source
      ) %>%
      dplyr::left_join(
        feature_metabolite_df,
        by = c("variable_id" = "feature", "INCHI_ID" = "metabolite")
      ) %>%
      dplyr::filter(!is.na(network_supported) & network_supported) %>%
      dplyr::mutate(
        mapped_compound_name = Compound_name,
        mapped_hmdb_id = HMDB_ID,
        mapped_kegg_id = KEGG_ID
      ) %>%
      dplyr::select(
        annotation_row_id,
        variable_id,
        INCHI_ID,
        score_norm,
        mapping_source,
        network_supported,
        mapped_compound_name,
        mapped_hmdb_id,
        mapped_kegg_id
      ) %>%
      dplyr::arrange(variable_id, dplyr::desc(score_norm))

    # Feature-metabolite edge frequencies, resolved per edge: edge endpoints are
    # stored in string order, so the metabolite (InChI) may sit in either column
    # depending on the feature ID format.
    fm_edge_frequency <- edge_frequency %>%
      dplyr::mutate(
        from_is_inchi = grepl("^InChI", from),
        to_is_inchi   = grepl("^InChI", to)
      ) %>%
      dplyr::filter(xor(from_is_inchi, to_is_inchi)) %>%
      dplyr::transmute(
        variable_id = dplyr::if_else(from_is_inchi, to, from),
        INCHI_ID    = dplyr::if_else(from_is_inchi, from, to),
        freq        = freq
      )

    anno_feature_table2_extra <- anno_feature_table2_extra %>%
      dplyr::left_join(fm_edge_frequency, by = c("variable_id", "INCHI_ID"))

    anno_feature_table2_extra <- anno_feature_table2_extra %>%
      dplyr::select(
        annotation_row_id,
        variable_id,
        INCHI_ID,
        score_norm,
        mapping_source,
        network_supported,
        mapped_compound_name,
        mapped_hmdb_id,
        mapped_kegg_id,
        freq
      ) %>%
      dplyr::rename(edge_freq = freq) %>%
      dplyr::arrange(variable_id, dplyr::desc(score_norm), dplyr::desc(edge_freq))

    edge_info <- stats::sd(anno_feature_table2_extra$edge_freq, na.rm = TRUE)
    if (is.na(edge_info) || !is.finite(edge_info)) edge_info <- 0
    edge_weight <- 0.1 + 0.2 * pmin(edge_info / 0.25, 1)
    score_weight <- 1 - edge_weight

    anno_feature_table2_extra <- anno_feature_table2_extra %>%
      dplyr::mutate(
        weighted_anno_score = score_weight * score_norm + edge_weight * edge_freq
      ) %>%
      dplyr::group_by(variable_id) %>%
      dplyr::mutate(confidence_rank = dplyr::dense_rank(dplyr::desc(weighted_anno_score))) %>%
      dplyr::ungroup()

    anno_feature_table2 <- anno_feature_table1_sub %>%
      dplyr::filter(annotation_row_id %in% anno_feature_table2_extra$annotation_row_id) %>%
      dplyr::left_join(anno_feature_table2_extra, by = c("annotation_row_id", "variable_id", "score_norm")) %>%
      dplyr::left_join(
        res$feature_table[, c("variable_id", "importance")] %>%
          dplyr::distinct(variable_id, .keep_all = TRUE),
        by = "variable_id"
      ) %>%
      dplyr::mutate(
        Compound.name = dplyr::coalesce(Compound.name, mapped_compound_name),
        HMDB.ID = dplyr::coalesce(HMDB.ID, mapped_hmdb_id),
        KEGG.ID = dplyr::coalesce(KEGG.ID, mapped_kegg_id)
      ) %>%
      dplyr::select(
        -annotation_row_id,
        -network_supported,
        -mapped_compound_name,
        -mapped_hmdb_id,
        -mapped_kegg_id,
        -mapping_source
      ) %>%
      dplyr::arrange(variable_id, dplyr::desc(weighted_anno_score), dplyr::desc(edge_freq))

    list(
      pcsf_subnet = pcsf_subnet,
      node_frequency = node_frequency,
      edge_frequency = edge_frequency,
      anno_feature_table2 = anno_feature_table2
    )
  }

  # ---- optional grid tuning ----
  best_b <- b
  best_w <- w

  if (isTRUE(grid_tune)) {
    cli::cli_alert_info("{step_label(3)} Tuning PCSF parameters on the grid.")

    grid <- tidyr::expand_grid(b = b_grid, w = w_grid)
    grid$unique_variable_id <- NA_integer_
    grid$mean_weighted_score <- NA_real_
    grid$tuning_success <- FALSE
    grid$tuning_error <- NA_character_

    safe_cli_progress_update <- function(...) {
      try(cli::cli_progress_update(...), silent = TRUE)
      invisible(NULL)
    }
    safe_cli_progress_done <- function(id) {
      try(cli::cli_progress_done(id = id), silent = TRUE)
      invisible(NULL)
    }

    tuning_progress_id <- cli::cli_progress_bar(
      name = "Grid tuning",
      total = nrow(grid),
      format = "{cli::pb_name} [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} {cli::pb_percent} | {cli::pb_status}",
      current = FALSE,
      auto_terminate = FALSE,
      .auto_close = FALSE
    )
    on.exit(safe_cli_progress_done(tuning_progress_id), add = TRUE)

    for (i in seq_len(nrow(grid))) {
      bb <- grid$b[i]
      ww <- grid$w[i]
      safe_cli_progress_update(
        id = tuning_progress_id,
        status = paste("Evaluating b =", bb, ", w =", ww)
      )

      out <- tryCatch(
        run_pcsf_once(
          b_val = bb,
          w_val = ww,
          times_val = times_tune,
          show_progress = FALSE
        ),
        error = function(e) e
      )

      if (!inherits(out, "error")) {
        grid$unique_variable_id[i] <- dplyr::n_distinct(out$anno_feature_table2$variable_id)
        grid$mean_weighted_score[i] <- if (nrow(out$anno_feature_table2) > 0) {
          mean(out$anno_feature_table2$weighted_anno_score, na.rm = TRUE)
        } else {
          NA_real_
        }
        grid$tuning_success[i] <- TRUE
      } else {
        grid$unique_variable_id[i] <- 0L
        grid$mean_weighted_score[i] <- NA_real_
        grid$tuning_error[i] <- conditionMessage(out)
      }

      safe_cli_progress_update(id = tuning_progress_id, set = i)
    }

    safe_cli_progress_done(tuning_progress_id)

    grid2 <- grid %>%
      dplyr::arrange(
        dplyr::desc(tuning_success),
        dplyr::desc(unique_variable_id),
        dplyr::desc(mean_weighted_score),
        b,
        w
      )

    if (!any(grid2$tuning_success)) {
      stop("[Error] Grid tuning failed for all candidate (b, w) pairs.")
    }

    best_b <- grid2$b[1]
    best_w <- grid2$w[1]

    res$tuning_grid <- grid2
    res$tuning_best <- list(b = best_b, w = best_w,
                            unique_variable_id = grid2$unique_variable_id[1],
                            mean_weighted_score = grid2$mean_weighted_score[1],
                            times_tune = times_tune)

    cli::cli_alert_info(
      paste(
        "Grid tuning evaluated", nrow(grid2), "candidate (b, w) pairs;",
        sum(grid2$tuning_success), "succeeded.",
        "Selected b =", best_b, "and w =", best_w,
        "for the final run based on retained feature count and mean weighted annotation score."
      )
    )
  }


  # ---- final run with (best_b, best_w) and full times/workers ----
  cli::cli_alert_info("{step_label(total_steps)} Finalizing network-supported annotations.")
  final_out <- run_pcsf_once(b_val = best_b, w_val = best_w, times_val = times)

  res$pcsf_subnet <- final_out$pcsf_subnet
  res$node_frequency <- final_out$node_frequency
  res$edge_frequency <- final_out$edge_frequency
  res$anno_feature_table2 <- final_out$anno_feature_table2

  res <- set_featuremap_stage(res, "subnetwork_extracted")

  end_time <- Sys.time()
  elapsed_min <- round(as.numeric(end_time - start_time, units = "mins"), 2)
  cli::cli_alert_success(
    "Constructed a subnetwork with {.val {igraph::gorder(res$pcsf_subnet)}} nodes and {.val {igraph::gsize(res$pcsf_subnet)}} edges."
  )
  cli::cli_alert_success(
    "Retained {.val {nrow(res$anno_feature_table2)}} network-supported annotation rows for {.val {dplyr::n_distinct(res$anno_feature_table2$variable_id)}} features."
  )
  cli::cli_alert_info("Elapsed time: {.val {elapsed_min}} minutes.")

  return(res)
}
