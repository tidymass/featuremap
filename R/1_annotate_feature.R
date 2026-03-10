
#' Annotate raw features via MS1 or MS2 information
#'
#' Perform preliminary metabolite annotation for raw MS1 features
#' by matching against a reference metabolite database using mass
#' accuracy, retention time tolerance, isotope patterns, and optional
#' MS2 confirmation.
#'
#' @param feature_table A feature-level data frame to annotate. Must satisfy
#'   \code{\link{check_feature_table}} requirements, including columns:
#'   \code{variable_id}, \code{mz}, \code{rt}, \code{polarity},
#'   \code{mean_intensity}, \code{condition}, and \code{pvalue}.
#' @param ms2_table Optional MS2 confirmation table. Can be \code{NULL}.
#'   If provided, it should contain \code{variable_id} and the identifier
#'   column matching \code{refdata} (that is, \code{HMDB.ID} for
#'   \code{refdata = "hmdb"} or \code{KEGG.ID} for
#'   \code{refdata = "kegg"}).
#' @param metabolite_database A \code{metid} database object (typically
#'   \code{databaseClass}) used for MS1 annotation. If \code{NULL},
#'   the function auto-loads a cached/reference database from the
#'   featuremap resource repository according to \code{refdata}.
#' @param refdata Character. Reference namespace for downstream ID mapping.
#'   Must be one of \code{"kegg"} or \code{"hmdb"}.
#' @param column Character. LC column mode used by \code{metid} adduct/matching
#'   rules. Must be one of \code{"rp"} or \code{"hilic"}.
#' @param ms1_match_ppm Numeric. MS1 mass tolerance (ppm) used in candidate
#'   matching.
#' @param mfc_rt_tol Numeric. Retention-time tolerance (same unit as
#'   \code{feature_table$rt}) used in metabolite-feature-cluster scoring
#'   and isotope consistency checks.
#' @param isotope_number Integer. Maximum isotope depth used in isotope
#'   annotation.
#' @param workers Integer. Number of workers for parallel steps (isotope
#'   annotation) and threads passed to \code{metid::annotate_metabolites()}.
#'
#' @return A \code{featuremap_result} R6 object containing at least:
#'   \itemize{
#'     \item \code{feature_table}: input feature table (with derived fields added later in pipeline)
#'     \item \code{refdata}: selected reference tag (\code{"kegg"} or \code{"hmdb"})
#'     \item \code{anno_feature_table1}: feature-metabolite annotation table after redundancy removal
#'     \item \code{anno_redundancy}: redundancy summary before/after filtering
#'     \item \code{object_type}: \code{"featuremap_result"}
#'     \item \code{stage}: \code{"feature_annotated"}
#'   }
#'
#' @examples
#' \dontrun{
#' data("demo_feature_table", package = "featuremap")
#'
#' res <- annotate_feature(
#'   feature_table = demo_feature_table,
#'   ms2_table = NULL,
#'   metabolite_database = NULL,
#'   refdata = "kegg",
#'   column = "rp",
#'   ms1_match_ppm = 15,
#'   mfc_rt_tol = 5,
#'   isotope_number = 3,
#'   workers = 4
#' )
#' }
#'
#' @export

annotate_feature <- function(feature_table,
                             ms2_table = NULL,
                             metabolite_database = NULL,
                             refdata = c("kegg","hmdb"),
                             column = c("rp", "hilic"),
                             ms1_match_ppm = 15,
                             mfc_rt_tol = 5,
                             isotope_number = 3,
                             workers = 4
){
  refdata <- match.arg(refdata)
  column <- match.arg(column)
  if (is.null(metabolite_database)) {
    resource_name <- featuremap_default_resource_name("metabolite_database", refdata = refdata)
    message("metabolite_database is NULL. Loading resource from cache/GitHub: ", resource_name)
    metabolite_database <- featuremap_load_resource(resource_name)
  }

  start_time <- Sys.time()

  res <- R6Class(
    "featuremap_result",
    public = list(initialize = function(...) {}),
    lock_objects = FALSE
  )$new()
  res <- set_featuremap_stage(res, "initialized")

  res$feature_table = feature_table
  res$refdata = refdata


  res$anno_feature_table1 <-
    annotate_feature_table(
      feature_table = feature_table,
      ms2_table = ms2_table,
      column = column,
      refdata = refdata,
      metabolite_database = metabolite_database,
      ms1_match_ppm = ms1_match_ppm,
      mfc_rt_tol = mfc_rt_tol,
      isotope_number = isotope_number,
      workers = workers
    )

  message("\nStarting to remove annotation redundancy ...")
  rr1 = calculate_redundancy(annotation_table = res$anno_feature_table1)

  res$anno_feature_table1 =
    remove_redundancy(annotation_table = res$anno_feature_table1)

  rr2 = calculate_redundancy(annotation_table = res$anno_feature_table1)

  res$anno_feature_table1 %>%
    dplyr::mutate(score_norm = score / max(res$anno_feature_table1$score, na.rm = TRUE)) %>%
    dplyr::arrange(variable_id, dplyr::desc(score)) -> res$anno_feature_table1

  res$anno_redundancy = data.frame(
    mfc_per_metabolite = as.numeric(c(rr1[1], rr2[1])),
    metabolite_per_feature = as.numeric(c(rr1[2], rr2[2]))
  )

  rownames(res$anno_redundancy) = c("before", "after")
  res <- set_featuremap_stage(res, "feature_annotated")

  message("\nStep 1: Feature annotation has finished!\n")

  end_time <- Sys.time()
  cat("Execution time:", round(as.numeric(end_time - start_time, units = "mins"), 2), "minutes\n")

  return(res)

}
