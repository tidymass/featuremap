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
#' @param mrn An \code{igraph} object for the metabolic reaction network.
#'   If \code{NULL}, a cached/default MRN resource will be loaded based on
#'   \code{res$refdata}.
#' @param workers Integer. Number of parallel workers used in repeated PCSF runs.
#' @param times Integer. Number of PCSF repeats for final subnetwork construction.
#'   Larger values usually improve robustness but increase runtime.
#' @param b Numeric. Prize scaling factor for terminal nodes in PCSF.
#' @param mu Numeric. Hub penalty coefficient. Larger values penalize high-degree
#'   nodes more strongly.
#' @param w Numeric. Dummy edge cost in PCSF.
#' @param r Numeric. Random perturbation magnitude applied to edge weights across
#'   repeats (robustness resampling factor).
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
#'   \item \code{anno_feature_table2}: network-refined feature-metabolite table
#'   \item \code{tuning_grid} and \code{tuning_best}: (if \code{grid_tune = TRUE})
#'         tuning summary
#'   \item \code{stage}: updated to \code{"subnetwork_extracted"}
#' }
#'
#' @details
#' Feature importance scores are derived from effects value and p-value.
#' The PCSF algorithm is executed multiple times with randomized edge
#' perturbation to ensure robustness. Final subnetworks are constructed
#' based on node and edge frequency across runs.
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
#'   b = 40,
#'   mu = 1e-4,
#'   w = 35,
#'   r = 0.30,
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
                              b = 40,
                              mu = 1e-4,
                              w = 30,
                              r = 0.30,
                              grid_tune = FALSE,
                              b_grid = seq(20, 40, by = 5),
                              w_grid = seq(25, 40, by = 5),
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
  if (is.null(mrn)){
    ref_tag <- if (!is.null(res$refdata)) as.character(res$refdata) else NULL
    resource_name <- featuremap_default_resource_name("mrn", refdata = ref_tag)
    message("mrn is NULL. Loading resource from cache/GitHub: ", resource_name)
    mrn <- featuremap_load_resource(resource_name)
  }

  start_time <- Sys.time()

  res$mrn = mrn
  res$feature_table %>%
    dplyr::mutate(
      score_fc = {
        fc <- abs(condition)
        (fc - min(fc, na.rm = TRUE)) /
          (max(fc, na.rm = TRUE) - min(fc, na.rm = TRUE))
      },
      score_p = -log10(pvalue),
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

  ## build FMpair + mrnref ----
  if(res$refdata == "kegg"){
    res$anno_feature_table1 %>%
      dplyr::filter(variable_id %in% feature_sig$variable_id) %>%
      dplyr::select(variable_id, Compound.name, KEGG.ID, score_norm) -> anno_FM

    mrn %>%
      tidygraph::activate(nodes) %>%
      as.data.frame() %>%
      dplyr::filter(KEGG_ID != "<NA>") %>%
      dplyr::select(KEGG_ID, INCHI_ID, Compound_name) -> mrnref

    c_min = 0.02
    suppressWarnings({
      anno_FM %>%
        dplyr::left_join(mrnref, by = c("KEGG.ID" = "KEGG_ID")) %>%
        dplyr::filter(!is.na(INCHI_ID)) %>%
        dplyr::arrange(INCHI_ID) %>%
        dplyr::distinct(variable_id, KEGG.ID, .keep_all = TRUE) %>%
        dplyr::mutate(cost = 1-score_norm+c_min) %>%
        dplyr::select(variable_id, INCHI_ID, cost) %>%
        purrr::set_names(c("from", "to", "weight")) -> FMpair
    })

  }else if(res$refdata == "hmdb"){
    res$anno_feature_table1 %>%
      dplyr::filter(variable_id %in% feature_sig$variable_id) %>%
      dplyr::select(variable_id, Compound.name, HMDB.ID, score_norm) -> anno_FM

    mrn %>%
      tidygraph::activate(nodes) %>%
      as.data.frame() %>%
      dplyr::filter(HMDB_ID != "<NA>") %>%
      dplyr::select(HMDB_ID, INCHI_ID, Compound_name) -> mrnref

    c_min = 0.02
    suppressWarnings({
      anno_FM %>%
        dplyr::left_join(mrnref, by = c("HMDB.ID" = "HMDB_ID")) %>%
        dplyr::filter(!is.na(INCHI_ID)) %>%
        dplyr::arrange(INCHI_ID) %>%
        dplyr::distinct(variable_id, HMDB.ID, .keep_all = TRUE) %>%
        dplyr::mutate(cost = 1-score_norm+c_min) %>%
        dplyr::select(variable_id, INCHI_ID, cost) %>%
        purrr::set_names(c("from", "to", "weight")) -> FMpair
    })
  }

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
  message("Setting node and edge weights for pcsf ...")

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
    dplyr::mutate(penalty = -mu * sqrt(degree)) %>%
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
    message("Initializing parallel backend with ", workers, " workers...")
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)
  }

  # define run pcsf once function ---
  run_pcsf_once <- function(b_val, w_val, times_val){

    node_prizes <- tibble::tibble(
      name = node_names,
      base_prize = b_val * node_prz,
      penalty = hub_penalty
    ) %>%
      dplyr::mutate(final_prize = dplyr::if_else(base_prize == 0, penalty, base_prize)) %>%
      dplyr::pull(final_prize)

    message("Running PCSF (", times_val, " iterations; b=", b_val, ", w=", w_val,
            if (workers > 1) paste0(", parallel with ", workers, " workers") else ", sequential",
            ") ...")

    if(workers > 1){
      results_list <- furrr::future_map(
        1:times_val,
        function(i) {
          cost <- c(
            rep(w_val, length(dummies)),
            igraph::E(imrn)$weight + igraph::E(imrn)$weight *
              stats::runif(length(igraph::E(imrn)), 0, max(r, stats::sd(igraph::E(imrn)$weight)))
          )
          output <- PCSF2::call_sr(from, to, cost, node_names, node_prizes)
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
        .progress = TRUE
      )
    }else{
      results_list <- purrr::map(
        1:times_val,
        function(i) {
          cost <- c(
            rep(w_val, length(dummies)),
            igraph::E(imrn)$weight + igraph::E(imrn)$weight *
              stats::runif(length(igraph::E(imrn)), 0, max(r, stats::sd(igraph::E(imrn)$weight)))
          )
          output <- PCSF2::call_sr(from, to, cost, node_names, node_prizes)
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
      dplyr::left_join(mrnref, by = c("metabolite" = "INCHI_ID"))
    rownames(feature_metabolite_df) = NULL

    if(res$refdata == "hmdb"){
      anno_feature_table2 <- res$anno_feature_table1 %>%
        dplyr::select(variable_id, HMDB.ID, score_norm) %>%
        dplyr::left_join(feature_metabolite_df, by = c("variable_id" = "feature", "HMDB.ID" = "HMDB_ID")) %>%
        dplyr::filter(!is.na(metabolite)) %>%
        dplyr::rename(INCHI_ID = metabolite) %>%
        dplyr::rename(HMDB_ID = HMDB.ID) %>%
        dplyr::arrange(variable_id)
    }else if(res$refdata == "kegg"){
      anno_feature_table2 <- res$anno_feature_table1 %>%
        dplyr::select(variable_id, KEGG.ID, score_norm) %>%
        dplyr::left_join(feature_metabolite_df, by = c("variable_id" = "feature", "KEGG.ID" = "KEGG_ID")) %>%
        dplyr::filter(!is.na(metabolite)) %>%
        dplyr::rename(INCHI_ID = metabolite) %>%
        dplyr::rename(KEGG_ID = KEGG.ID) %>%
        dplyr::arrange(variable_id)
    }

    from_is_inchi <- all(grepl("^InChI", edge_frequency$from))
    to_is_inchi   <- all(grepl("^InChI", edge_frequency$to))

    if(from_is_inchi){
      anno_feature_table2 <- anno_feature_table2 %>%
        dplyr::left_join(edge_frequency %>%
                           dplyr::filter(to %in% anno_feature_table2$variable_id),
                         by = c("variable_id" = "to", "INCHI_ID" = "from"))
    }else if(to_is_inchi){
      anno_feature_table2 <- anno_feature_table2 %>%
        dplyr::left_join(edge_frequency %>%
                           dplyr::filter(to %in% anno_feature_table2$variable_id),
                         by = c("variable_id" = "from", "INCHI_ID" = "to"))
    }

    if (res$refdata == "hmdb") {
      anno_feature_table2 <- anno_feature_table2 %>%
        dplyr::select(variable_id, Compound_name, HMDB_ID, INCHI_ID, score_norm, freq) %>%
        dplyr::rename(edge_freq = freq) %>%
        dplyr::arrange(variable_id, dplyr::desc(score_norm), dplyr::desc(edge_freq))
    } else {
      anno_feature_table2 <- anno_feature_table2 %>%
        dplyr::select(variable_id, Compound_name, KEGG_ID, INCHI_ID, score_norm, freq) %>%
        dplyr::rename(edge_freq = freq) %>%
        dplyr::arrange(variable_id, dplyr::desc(score_norm), dplyr::desc(edge_freq))
    }

    anno_feature_table2 <- anno_feature_table2 %>%
      dplyr::mutate(weighted_anno_score = 0.7*score_norm + 0.3*edge_freq) %>%
      dplyr::group_by(variable_id) %>%
      dplyr::mutate(confidence_rank = dplyr::dense_rank(dplyr::desc(weighted_anno_score))) %>%
      dplyr::ungroup()

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
    message("\n[Grid tuning] Start tuning b (", min(b_grid), "-", max(b_grid),
            ") and w (", min(w_grid), "-", max(w_grid), ") ...")

    grid <- tidyr::expand_grid(b = b_grid, w = w_grid)
    grid$unique_variable_id <- NA_integer_
    grid$mean_weighted_score <- NA_real_

    for (i in seq_len(nrow(grid))) {
      bb <- grid$b[i]
      ww <- grid$w[i]

      out <- tryCatch(
        run_pcsf_once(b_val = bb, w_val = ww, times_val = times_tune),
        error = function(e) NULL
      )

      if (!is.null(out)) {
        grid$unique_variable_id[i] <- dplyr::n_distinct(out$anno_feature_table2$variable_id)
        grid$mean_weighted_score[i] <- mean(out$anno_feature_table2$weighted_anno_score, na.rm = TRUE)
      } else {
        grid$unique_variable_id[i] <- 0L
        grid$mean_weighted_score[i] <- NA_real_
      }
    }

    grid2 <- grid %>%
      dplyr::arrange(dplyr::desc(unique_variable_id), dplyr::desc(mean_weighted_score))

    best_b <- grid2$b[1]
    best_w <- grid2$w[1]

    res$tuning_grid <- grid2
    res$tuning_best <- list(b = best_b, w = best_w,
                            unique_variable_id = grid2$unique_variable_id[1],
                            mean_weighted_score = grid2$mean_weighted_score[1],
                            times_tune = times_tune)

    message("[Grid tuning] Best parameters: b = ", best_b, ", w = ", best_w,
            " (unique variable_id = ", grid2$unique_variable_id[1], ").\n")
  }


  # ---- final run with (best_b, best_w) and full times/workers ----
  final_out <- run_pcsf_once(b_val = best_b, w_val = best_w, times_val = times)

  res$pcsf_subnet <- final_out$pcsf_subnet
  res$node_frequency <- final_out$node_frequency
  res$edge_frequency <- final_out$edge_frequency
  res$anno_feature_table2 <- final_out$anno_feature_table2

  message("\nStep 2: Union PCSF-based subnetwork has been constructed: ",
          igraph::gorder(res$pcsf_subnet), " nodes and ", igraph::gsize(res$pcsf_subnet), " edges.")
  res <- set_featuremap_stage(res, "subnetwork_extracted")

  end_time <- Sys.time()
  cat("Execution time:", round(as.numeric(end_time - start_time, units = "mins"), 2), "minutes\n")

  return(res)
}
