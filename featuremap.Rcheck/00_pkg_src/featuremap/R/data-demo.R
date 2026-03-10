#' Demo Feature Table for featuremap
#'
#' A small example feature table used to demonstrate the featuremap workflow.
#'
#' @format A data frame with 7 variables:
#' \describe{
#'   \item{variable_id}{Feature ID (character).}
#'   \item{mz}{Mass-to-charge ratio (numeric).}
#'   \item{rt}{Retention time (numeric).}
#'   \item{polarity}{Ion mode, either "positive" or "negative" (character).}
#'   \item{mean_intensity}{Mean feature intensity (numeric).}
#'   \item{condition}{Effect-size-like value used for network scoring (numeric).}
#'   \item{pvalue}{Statistical p-value (numeric).}
#' }
#' @source Internal demo dataset for package examples.
#' @examples
#' data("demo_feature_table", package = "featuremap")
#' head(demo_feature_table)
#' colnames(demo_feature_table)
#' @docType data
#' @name demo_feature_table
#' @usage data(demo_feature_table)
#' @keywords datasets
NULL
