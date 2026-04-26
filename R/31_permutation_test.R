

#' Permutation test with Gamma distribution
#' 
#' @param subnet Steiner subnetwork graph with vertex `type` labels.
#' @param mrn Full metabolic reaction network
#' @param module_table Module membership/annotation table containing
#'   `INCHI_ID`, `Module`, `Type`, and `Importance`.
#' @param module_score Observed module scores to be tested.
#' @param method Community detection method
#' @param min_size Minimum module size
#' @param n_perm Number of permutations
#' @param workers Number of parallel workers
#' @param seed Random seed
#' 
#' @return Data frame with p-values
permutation_test_gamma <- function(subnet,
                                   mrn,
                                   module_table,
                                   module_score,
                                   method = "louvain",
                                   min_size = 3,
                                   n_perm = 100, 
                                   workers = 4,
                                   seed = 123) {
  
  set.seed(seed)
  

  all_mrn_nodes <- V(mrn)$name
  n_annotated <- sum(V(subnet)$type == "annotated")
  n_hidden <- sum(V(subnet)$type == "hidden")
  
  # prepare the degree distribution ---------
  
  mrn_df <- tibble::tibble(
    name   = V(mrn)$name,
    degree = V(mrn)$degree
  )
  
  breaks <- c(-Inf, 1, 3, 5, 10, 50, Inf)
  mrn_df$bin <- cut(
    mrn_df$degree,
    breaks = breaks,
    labels = FALSE,
    right  = TRUE,        
    include.lowest = TRUE
  )
  
  subnet_df <- mrn_df %>%
    dplyr::filter(name %in% V(subnet)$name) %>%
    dplyr::mutate(
      type = V(subnet)$type[match(name, V(subnet)$name)]
    )
  
  type_bin_counts <- subnet_df %>%
    dplyr::count(type, bin, name = "n_real")
  
  # ---- importance pool (from real data) ----
  importance_pool <- module_table %>%
    dplyr::filter(Type == "annotated", is.finite(Importance)) %>%
    dplyr::distinct(INCHI_ID, .keep_all = TRUE) %>%
    dplyr::pull(Importance)
  

  # perform permutation test ------
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  
  perm_results <- future_map_dfr(1:n_perm, function(i) {
    
    # Extract fake Steiner subnet via stratified sampling ---
    avail_df <- mrn_df
    fake_annotated <- character(0)
    fake_hidden    <- character(0)
    
    for (b in sort(unique(mrn_df$bin))) {
      nA <- type_bin_counts %>%
        dplyr::filter(type == "annotated", bin == b) %>%
        dplyr::pull(n_real)
      nH <- type_bin_counts %>%
        dplyr::filter(type == "hidden", bin == b) %>%
        dplyr::pull(n_real)
      
      if (length(nA) == 0) nA <- 0
      if (length(nH) == 0) nH <- 0
      
      cand_bin <- avail_df %>% dplyr::filter(bin == b)
      
      if (nA > 0) {
        if (nrow(cand_bin) == 0) next
        idxA <- sample(seq_len(nrow(cand_bin)),
                       size    = nA,
                       replace = nA > nrow(cand_bin))
        selA <- cand_bin$name[idxA]
        fake_annotated <- c(fake_annotated, selA)
        avail_df <- avail_df %>%
          dplyr::filter(!name %in% selA)
        cand_bin <- avail_df %>% dplyr::filter(bin == b)
      }
      
      if (nH > 0) {
        if (nrow(cand_bin) == 0) next
        idxH <- sample(seq_len(nrow(cand_bin)),
                       size    = nH,
                       replace = nH > nrow(cand_bin))
        selH <- cand_bin$name[idxH]
        fake_hidden <- c(fake_hidden, selH)
        avail_df <- avail_df %>%
          dplyr::filter(!name %in% selH)
      }
    }
    
    # ---- assign importance to fake annotated nodes (shuffle / bootstrap) ----
    fake_importance_tbl <- tibble::tibble(
      INCHI_ID   = fake_annotated,
      Importance = sample(
        importance_pool,
        size    = length(fake_annotated),
        replace = TRUE))

    fake_steiner <- c(fake_annotated, fake_hidden)
    fake_steiner_subnet <- induced_subgraph(mrn, vids = fake_steiner)
    
    # Community detection for fake sub-network
    fake_comm <- tryCatch({
      if (method == "louvain") {
        comm <- igraph::cluster_louvain(fake_steiner_subnet)
      } else if (method == "fast_greedy") {
        comm <- igraph::cluster_fast_greedy(igraph::as_undirected(fake_steiner_subnet))
      } else if (method == "walktrap") {
        comm <- igraph::cluster_walktrap(fake_steiner_subnet, steps = 4)
      } else if (method == "leading_eigen") {
        comm <- igraph::cluster_leading_eigen(fake_steiner_subnet)
      } else if (method == "newman_spectral") {
        comm <- rSpectral::spectral_igraph_communities(
          fake_steiner_subnet,
          Cn_min = as.integer(min_size),
          names = 1L
        )
      }
    }, error = function(e) NULL)
    
    if (is.null(fake_comm)){
      return(data.frame(Perm = i, Module_score = 0))
    }
    
    fake_member <- igraph::membership(fake_comm)
    
    # Build fake module table
    fake_module_table <- data.frame(
      INCHI_ID = names(fake_member),
      Module = as.integer(fake_member)) %>% 
      group_by(Module) %>% 
      filter(n() >= min_size) %>% 
      ungroup()
    
    if (nrow(fake_module_table) == 0){
      return(data.frame(
        Perm = i,
        Module_score = 0
      ))
    }
    
    # Calculate module scores-----
    fake_scores <- lapply(unique(fake_module_table$Module), function(mod) {
      
      module_nodes <- fake_module_table %>% 
        filter(Module == mod) %>% 
        pull(INCHI_ID)
      
      
      activity_result <- calc_activity_score(
        graph = mrn,
        module_nodes = module_nodes,
        seed_nodes = fake_annotated,
        importance_tbl = fake_importance_tbl)
      
      data.frame(
        Perm = i,
        Module_score = activity_result$Module_score
      )
      
    }) %>% bind_rows()
    
    return(fake_scores)
    
  }, 
  .progress = TRUE, 
  .options = furrr_options(
    seed = TRUE,
    packages = c("igraph","dplyr","tibble","tidyr","MASS","rSpectral")
  ))
  
  future::plan(future::sequential)
  
  # Remove invalid results
  perm_results <- perm_results %>% 
    filter((!is.na(Module_score)) & is.finite(Module_score))
  
  # message(paste0("\n",sprintf("Generated %d null moduel scores", nrow(perm_results))))
  
  # Fit Gamma distribution
  null_scores <- perm_results$Module_score
  min_score <- min(null_scores)
  
  if (min_score <= 0) {
    shift <- abs(min_score) + 0.01
    null_scores_shifted <- null_scores + shift
  } else {
    shift <- 0
    null_scores_shifted <- null_scores
  }
  
  fit <- tryCatch({
    suppressWarnings(MASS::fitdistr(null_scores_shifted, "gamma"))
  }, error = function(e) {
    mean_score <- mean(null_scores_shifted)
    var_score <- var(null_scores_shifted)
    list(estimate = c(shape = mean_score^2 / var_score, rate = mean_score / var_score))
  })
  
  shape_param <- fit$estimate["shape"]
  rate_param <- fit$estimate["rate"]
  
  # message(sprintf(
  #   "Fitted Gamma distribution: shape = %.4f, rate = %.4f",
  #   shape_param, rate_param
  # ))
  
  # Calculate p-values---
  pvalue_results <- module_score %>%
    dplyr::rowwise() %>%
    mutate(
      p = pgamma(Module_score + shift, shape = shape_param, rate  = rate_param, lower.tail = FALSE)
    ) %>%
    dplyr::ungroup() %>%
    mutate(
      fdr = p.adjust(p, method = "BH")
    ) %>%
    dplyr::arrange(p)
  
  
  return(pvalue_results)
  
  
}
