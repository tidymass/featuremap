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
#' data("demo_ms2", package = "featuremap")
#' head(demo_ms2)
#' @docType data
#' @name demo_feature_table
#' @usage data(demo_feature_table)
#' @keywords datasets
NULL

#' Demo MS2 Confirmation Table for featuremap
#'
#' A small example MS2 confirmation table used together with
#' \code{demo_feature_table} to demonstrate workflows that include optional
#' MS2 support.
#'
#' @format A data frame with 6 variables:
#' \describe{
#'   \item{variable_id}{Feature ID (character).}
#'   \item{HMDB.ID}{HMDB identifier used as the primary MS2 matching key (character).}
#'   \item{KEGG.ID}{KEGG identifier used only when \code{HMDB.ID} is absent (character).}
#'   \item{Compound.name}{Annotated metabolite name (character).}
#'   \item{SS}{Similarity score from the MS2 annotation workflow (numeric).}
#'   \item{Level}{Annotation confidence level (numeric).}
#' }
#' @source Internal demo dataset for package examples.
#' @examples
#' data("demo_ms2", package = "featuremap")
#' head(demo_ms2)
#' colnames(demo_ms2)
#' @docType data
#' @name demo_ms2
#' @usage data(demo_ms2)
#' @keywords datasets
NULL
