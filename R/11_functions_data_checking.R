#' Check the Format and Integrity of a Feature Table
#'
#' Performs a series of checks on a feature table to ensure it meets
#' the expected structure and content required for downstream analysis.
#' The table must contain exactly six columns with specific names and types.
#'
#' @param feature_table A `data.frame` with exactly six columns named:
#'   `variable_id`, `mz`, `rt`, `polarity`, `condition`, and `pvalue`.
#'   Each column must conform to specific data types and constraints.
#'
#' @return Returns `TRUE` (invisibly) if all checks pass. Otherwise, the function
#'   stops with an informative error message.
#'
#' @details
#' The function checks:
#' \itemize{
#'   \item That `feature_table` is a `data.frame`
#'   \item That it has exactly six columns
#'   \item That column names and order match the expected names
#'   \item That data types are correct for each column
#'   \item That `variable_id` values are unique
#'   \item That `polarity` only contains "positive" or "negative"
#' }
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' \dontrun{
#' # Example feature table
#' feature_table <- data.frame(
#'   variable_id = c("M1", "M2"),
#'   mz = c(100.1, 200.2),
#'   rt = c(300, 400),
#'   polarity = c("positive", "negative"),
#'   condition = c(1, 2),
#'   pvalue = c(0.01, 0.2)
#' )
#'
#' check_feature_table(feature_table)
#' }
#'


check_feature_table <-
  function(feature_table) {
    # Check class
    if (!is.data.frame(feature_table)) {
      stop("feature_table must be a data.frame.")
    }

    # # Check number of columns
    # if (ncol(feature_table) != 6) {
    #   stop("feature_table must have exactly 5 columns.")
    # }

    # Expected column names
    expected_names <- c("variable_id",
                        "mz",
                        "rt",
                        "polarity",
                        "condition",
                        "pvalue")
    if (!all(expected_names == colnames(feature_table))) {
      stop(paste0(
        "feature_table columns must be named: ",
        paste(expected_names, collapse = ",")
      ))
    }

    # Check variable_id
    if (!is.character(feature_table$variable_id)) {
      stop("variable_id column must be character.")
    }
    if (anyDuplicated(feature_table$variable_id)) {
      stop("variable_id column contains duplicated values.")
    }

    # Check mz
    if (!is.numeric(feature_table$mz)) {
      stop("mz column must be numeric.")
    }

    # Check rt
    if (!is.numeric(feature_table$rt)) {
      stop("rt column must be numeric.")
    }

    # Check condition
    if (!is.numeric(feature_table$condition)) {
      stop("condition column must be numeric.")
    }

    # Check polarity
    if (!is.character(feature_table$polarity)) {
      stop("polarity column must be character.")
    }

    # Check pvalue
    if (!is.numeric(feature_table$pvalue)) {
      stop("pvalue column must be numeric.")
    }

    allowed_polarities <- c("positive", "negative")
    if (!all(feature_table$polarity %in% allowed_polarities)) {
      stop("polarity column can only contain 'positive' or 'negative'.")
    }

    message("[OK] feature_table passed all checks.")
    invisible(TRUE)
  }



#' Check the Format of an MS2 Annotation Table
#'
#' Validates that the MS2 annotation table meets required structural and type specifications.
#' The function checks for required columns and ensures proper data types for each.
#'
#' @param check_annotation_table_ms2 A `data.frame` containing MS2 annotation results.
#'   It must include the following columns: `variable_id`, `mz`, `rt`, `HMDB.ID`, and `KEGG.ID`.
#'
#' @return Returns `TRUE` (invisibly) if all checks pass. Otherwise, stops with an informative error message.
#'
#' @details The function performs the following checks:
#' \itemize{
#'   \item `check_annotation_table_ms2` must be a data frame
#'   \item Required columns: `variable_id`, `mz`, `rt`, `HMDB.ID`, `KEGG.ID`
#'   \item `variable_id`: character, unique
#'   \item `mz` and `rt`: numeric
#'   \item `HMDB.ID` and `KEGG.ID`: character (can include `NA`)
#' }
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' \dontrun{
#' # Example annotation table
#' annotation_table <- data.frame(
#'   variable_id = c("M1", "M2"),
#'   mz = c(100.1, 200.2),
#'   rt = c(300, 400),
#'   HMDB.ID = c("HMDB00001", NA),
#'   KEGG.ID = c("C00031", NA),
#'   stringsAsFactors = FALSE
#' )
#'
#' check_annotation_table_ms2(annotation_table)
#' }
#'


check_annotation_table_ms2 <-
  function(check_annotation_table_ms2) {
    # Check that it is a data.frame
    if (!is.data.frame(check_annotation_table_ms2)) {
      stop("check_annotation_table_ms2 must be a data.frame.")
    }

    # Required columns
    required_cols <- c("variable_id", "mz", "rt", "HMDB.ID", "KEGG.ID")
    missing_cols <- setdiff(required_cols, colnames(check_annotation_table_ms2))
    if (length(missing_cols) > 0) {
      stop(paste(
        "check_annotation_table_ms2 is missing required columns:",
        paste(missing_cols, collapse = ", ")
      ))
    }

    # Check variable_id
    if (!is.character(check_annotation_table_ms2$variable_id)) {
      stop("variable_id column must be character.")
    }
    if (anyDuplicated(check_annotation_table_ms2$variable_id)) {
      stop("variable_id column contains duplicated values.")
    }

    # Check mz
    if (!is.numeric(check_annotation_table_ms2$mz)) {
      stop("mz column must be numeric.")
    }

    # Check rt
    if (!is.numeric(check_annotation_table_ms2$rt)) {
      stop("rt column must be numeric.")
    }

    # Check HMDB.ID
    if (!is.character(check_annotation_table_ms2$HMDB.ID)) {
      stop(
        "HMDB.ID column must be character. It can contain NA but the column itself must be character type."
      )
    }

    # Check KEGG.ID
    if (!is.character(check_annotation_table_ms2$KEGG.ID)) {
      stop(
        "KEGG.ID column must be character. It can contain NA but the column itself must be character type."
      )
    }

    message("[OK] check_annotation_table_ms2 passed all checks.")
    invisible(TRUE)
  }
