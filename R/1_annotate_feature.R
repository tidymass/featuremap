
#' Annotate raw MS1 features using m/z and retention time information
#'
#' Perform preliminary metabolite annotation for raw MS1 features
#' by matching against a reference metabolite database using mass
#' accuracy, retention time tolerance, isotope patterns, and optional
#' MS2 confirmation.
#'
#' @param feature_table A data frame containing raw MS1 feature information,
#'   including m/z and retention time columns.
#' @param ms2_table Optional data frame containing MS2 spectra information
#'   for feature-level confirmation. Can be NULL if MS2 is not available.
#' @param metabolite_database A data frame of candidate metabolites with
#'   theoretical m/z, formula, and other annotation fields.
#' @param refdata Character. Name of the reference database
#'   (e.g., "hmdb" or "kegg").
#' @param column Character. Name of the m/z column in `feature_table`.
#' @param ms1_match_ppm Numeric. Mass tolerance in parts-per-million (ppm)
#'   used for MS1 matching.
#' @param mfc_rt_tol Numeric. Retention time tolerance used for filtering
#'   candidate matches.
#' @param isotope_number Integer. Number of isotopic peaks considered
#'   during isotope pattern validation.
#' @param workers Integer. Number of parallel workers used for computation.
#'
#' @return A data frame containing annotated features with candidate
#'   metabolite matches and associated scoring metrics.
#'
#' @export

annotate_feature <- function(feature_table,
                             ms2_table = NULL,
                             metabolite_database,
                             refdata = c("kegg","hmdb"),
                             column = c("rp", "hilic"),
                             ms1_match_ppm = 15,
                             mfc_rt_tol = 5,
                             isotope_number = 3,
                             workers = 4
){

  start_time <- Sys.time()

  res <- R6Class("featureFMA_result", public = list(initialize = function(...) {}))

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

  message("\nStep 1: Feature annotation has finished!\n")

  end_time <- Sys.time()
  cat("Execution time:", round(as.numeric(end_time - start_time, units = "mins"), 2), "minutes\n")

  return(res)

}
