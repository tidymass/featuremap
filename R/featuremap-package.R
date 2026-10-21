#' Network-Based Functional Module Identification in Metabolomics
#'
#' The featuremap package provides a staged workflow for annotating metabolite
#' features, extracting a network-supported subnetwork, identifying functional
#' modules, and assigning biological interpretations to those modules.
#'
#' The main workflow is:
#'
#' * [annotate_feature()] annotates metabolite features.
#' * [extract_subnetwork()] builds a robust metabolic subnetwork.
#' * [identify_module()] detects and tests functional modules.
#' * [annotate_module()] generates module-level interpretations.
#'
#' Optional literature retrieval and large language model services are used
#' only when requested by the annotation workflow.
#'
#' @section Startup message:
#' The startup message is shown once per interactive R session. Set
#' `options(featuremap.quiet = TRUE)` before calling [library()] to suppress it.
#'
#' @seealso [annotate_feature()], [extract_subnetwork()], [identify_module()],
#'   [annotate_module()]
#' @keywords internal
"_PACKAGE"

#' @importFrom magrittr %>%
#' @importFrom R6 R6Class
#' @importFrom dplyr coalesce if_else n
#' @importFrom furrr furrr_options future_map_dfr
#' @importFrom igraph E E<- V V<- degree gsize induced_subgraph vcount
#' @importFrom methods is
#' @importFrom purrr imap imap_dfr
#' @importFrom stats aggregate as.formula p.adjust pgamma quantile rt setNames var weighted.mean
#' @importFrom stringr str_detect
#' @importFrom tibble tibble
#' @importFrom utils combn head
NULL
