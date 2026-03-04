#' @title calculate_redundancy
#' @description calculate_redundancy
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @param annotation_table annotation table
#' Calculate Annotation Redundancy Metrics
#'
#' Computes two types of annotation redundancy from an annotation table:
#' \itemize{
#'   \item \strong{Redundancy 1:} The average number of unique compound clusters (e.g., metabolite feature clusters)
#'     mapped to each compound (`Lab.ID`).
#'   \item \strong{Redundancy 2:} The average number of compound matches per metabolite peak (`variable_id`).
#' }
#'
#' @param annotation_table A data.frame or list of data.frames containing annotation results.
#'   Must include at least the columns `Lab.ID`, `metabolite_feature_cluster`, and `variable_id`.
#'
#' @return A named numeric vector with two elements:
#' \describe{
#'   \item{`redundancy1`}{Average number of unique compound classes per `Lab.ID`}
#'   \item{`redundancy2`}{Average number of compound matches per `variable_id`}
#' }
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' \dontrun{
#' # Example annotation table
#' annotation_table <- data.frame(
#'   Lab.ID = c("C001", "C001", "C002", "C002", "C002"),
#'   metabolite_feature_cluster = c("F1", "F2", "F1", "F1", "F3"),
#'   variable_id = c("M1", "M1", "M2", "M2", "M3"),
#'   stringsAsFactors = FALSE
#' )
#'
#' calculate_redundancy(annotation_table)
#' }
#'


calculate_redundancy <-
  function(annotation_table) {
    if (is(annotation_table, "list")) {
      annotation_table <-
        annotation_table %>%
        data.table::rbindlist()
    }

    ##redundancy1 means one compound contains how many compound class
    annotation_table <-
      data.table::as.data.table(annotation_table)

    redundancy1 <-
      annotation_table[, list(number = length(unique(metabolite_feature_cluster))), by = Lab.ID]

    redundancy1 <-
      mean(redundancy1$number)

    ##redundancy2 means one peak matches how many compound
    redundancy2 <-
      annotation_table[, .N, by = variable_id] %>%
      pull(N) %>%
      mean()

    c("mfc per metabolite" = redundancy1,
      "metabolite per feature" = redundancy2)
  }




#' Remove Redundant Annotations from an Annotation Table
#'
#' Iteratively filters out low-confidence annotations to reduce redundancy in metabolite identification results.
#' Redundancy is defined both at the compound level (`Lab.ID`) and the feature level (`name`).
#'
#' @param annotation_table A data frame containing metabolite annotations. Must include at least:
#'   `Lab.ID`, `name`, `score`, and `metabolite_feature_cluster`.
#'
#' @return A filtered annotation table with reduced redundancy and updated scores. The structure of the original
#'   table is preserved, but low-scoring or duplicate annotations are removed.
#'
#' @details
#' The function performs the following steps iteratively:
#' \itemize{
#'   \item For each compound (`Lab.ID`), if any annotation has a score > 100, annotations with score ≤ 20 are removed.
#'   \item For each peak (`name`), the same rule is applied.
#'   \item After filtering, scores for each metabolite feature cluster are recalculated using `score_mfc()`.
#'   \item Redundancy is recalculated using `calculate_redundancy()`. The process repeats until no improvement.
#' }
#'
#' @note This function depends on `calculate_redundancy()` and `score_mfc()` being defined and accessible in the environment.
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' \dontrun{
#' data("example_annotation_table")
#' filtered_table <- remove_redundancy(annotation_table = example_annotation_table)
#' }
#'
#' @importFrom dplyr group_by filter ungroup bind_rows
#' @importFrom purrr map
#' @importFrom plyr dlply
#'



remove_redundancy <-
  function(annotation_table) {
    ##for one compound, if one compound class with score > 100, then remove the
    ##compound class with score <= 20

    redundancy_diff = c(-1, -1)

    while (any(redundancy_diff < 0)) {
      # cat("i", " ")
      before_redundancy =
        calculate_redundancy(annotation_table = annotation_table)

      annotation_table <-
        annotation_table %>%
        dplyr::group_by(Lab.ID) %>%
        dplyr::filter(if (any(score > 100)) {
          score > 20
        } else{
          score > 0
        }) %>%
        dplyr::ungroup()

      ##for one peak, if it has a annotation with score > 100, then remove other
      ##annotations with score <= 20
      annotation_table <-
        annotation_table %>%
        dplyr::group_by(variable_id) %>%
        dplyr::filter(if (any(score > 100)) {
          score > 20
        } else{
          score > 0
        }) %>%
        dplyr::ungroup()


      ###re-calculated confidence score for each compound class.
      annotation_table =
        split(annotation_table, annotation_table$metabolite_feature_cluster) %>%
        purrr::map(function(x) {
          score <- score_mfc(x)
          x$score = score
          x
        }) %>%
        dplyr::bind_rows()

      after_redundancy =
        calculate_redundancy(annotation_table = annotation_table)

      redundancy_diff = after_redundancy - before_redundancy
    }

    annotation_table
  }
