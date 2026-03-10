#' Calculate module activity score
#' 
#' @param graph Full metabolic reaction network (mrn)
#' @param module_nodes Vector of node IDs in the module
#' @param seed_nodes Vector of annotated (seed) node IDs
#' @param importance_tbl Optional data frame with columns `INCHI_ID` and
#'   `Importance` used to compute module-level mean importance.
#' 
#' @return Data frame with module metrics
calc_activity_score <- function(graph, 
                                module_nodes, 
                                seed_nodes,
                                importance_tbl = NULL) {
  
  m   <- gsize(graph)  # Total edges in FULL network (mrn)
  N_I <- length(seed_nodes)  # Total annotated nodes
  N_G <- length(module_nodes)  # Nodes in module
  N_I_G <- sum(module_nodes %in% seed_nodes)  # Annotated nodes in module
  
  # Handle edge cases
  if (N_G < 2 || m == 0) {
    return(data.frame(
      Module_size = N_G, 
      Annotated_nodes = 0, 
      Modularity = NA_real_, 
      Activity_score = NA_real_
    ))
  }
  
  # Extract module subgraph from FULL network
  subg <- induced_subgraph(graph, vids = module_nodes)
  E_G  <- gsize(subg)  # Edges within module
  
  # Get degrees of module nodes in FULL network 
  degs <- degree(graph, v = module_nodes)
  
  # Calculate expected edges under null model
  expected <- 0
  for (i in seq_along(module_nodes)) {
    for (j in seq_along(module_nodes)) {
      if (i != j) {
        expected <- expected + as.numeric(degs[i] * degs[j])
      }
    }
  }
  expected <- expected / (4 * m)
  
  # Calculate modularity Q (adjusted by seed node proportion)
  Q <- sqrt(N_I / N_G) * (E_G - expected) / m
  
  # Calculate the mean importance for each module
  if (!is.null(importance_tbl)) {
    I <- importance_tbl %>%
      dplyr::filter(INCHI_ID %in% module_nodes) %>%
      dplyr::pull(Importance) %>%
      mean(na.rm = TRUE)
    if (is.nan(I)) I <- NA_real_
  } else {
    I <- NA_real_
  }
  
  # Calculate activity score (scaled by annotated node enrichment)
  A <- Q * (N_I_G / N_G) * I * 100
  
  return(data.frame(
    Module_size = N_G, 
    Annotated_nodes = N_I_G, 
    Module_score = A
  ))
}


#' Calculate activity scores for all modules
#' 
#' @param res Result object from pcsf_network
#' @param annotated_ID Vector of annotated node IDs
#' @param mrn Full metabolic reaction network
#' @param importance_tbl Data frame with columns `INCHI_ID` and `Importance`.
#' 
#' @return Data frame with module activity scores
calculate_module_activity <- function(res, annotated_ID, mrn, importance_tbl) {
  
  # Get unique modules
  modules <- unique(res$module_table$Module)
  
  # Calculate activity score for each module
  module_activity <- lapply(modules, function(mod) {
    
    # Get nodes in this module
    module_nodes <- res$module_table %>% 
      filter(Module == mod) %>% 
      pull(INCHI_ID)
    
    # Calculate activity score using FULL network (mrn)
    activity_result <- calc_activity_score(
      graph = mrn,  # Use full network here!
      module_nodes = module_nodes,
      seed_nodes = annotated_ID,
      importance_tbl = importance_tbl
    )
    
    # Add module ID
    activity_result$Module <- mod
    
    return(activity_result)
    
  }) %>% bind_rows()
  
  # Reorder columns
  module_activity <- module_activity %>% 
    dplyr::select(Module, Module_size, Annotated_nodes, Module_score) %>% 
    arrange(desc(Module_score))
  
  return(module_activity)
}
