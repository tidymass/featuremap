#' Group Peaks Based on Retention Time Proximity
#'
#' Clusters peaks into retention time (RT)-based groups. Peaks that fall within a specified RT window (`rt.tol`)
#' from one another are grouped into the same class. Each group is assigned a unique class label.
#'
#' @param rt A numeric vector of retention times.
#' @param rt.tol Numeric. RT tolerance (in seconds) for grouping nearby peaks. Peaks within `rt.tol` seconds are grouped together. Default is 10.
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{`rt`}{The original retention times}
#'   \item{`class`}{The assigned class/group ID for each RT}
#' }
#' The rows are sorted in the order of grouping.
#'
#' @details
#' Peaks are sequentially grouped by checking whether other peaks fall within `rt.tol` seconds of each peak.
#' Later groups exclude RTs already assigned to earlier clusters.
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' \dontrun{
#' rt_values <- c(100, 102, 109, 120, 123, 180)
#' group_peaks_rt(rt = rt_values, rt.tol = 10)
#' }
#'
#' @importFrom purrr map2

group_peaks_rt <- function(rt, rt.tol = 10) {
  rt_class <-
    lapply(rt, function(x) {
      which(rt >= x & rt < x + rt.tol)
    })

  rt_class <-
    lapply(seq_along(rt_class)[-1], function(i) {
      setdiff(rt_class[[i]], unlist(rt_class[1:(i - 1)]) %>% unique())
    }) %>%
    `c`(rt_class[1], .)

  rt_class <-
    rt_class[which(lapply(rt_class, length) != 0)]

  names(rt_class) <- seq_along(rt_class)

  rt_class <-
    purrr::map2(
      .x = rt_class,
      .y = names(rt_class),
      .f = function(x, y) {
        data.frame(rt = rt[x],
                   class = y,
                   stringsAsFactors = FALSE)
      }
    ) %>%
    do.call(rbind, .)
  rownames(rt_class) <- NULL
  return(rt_class)
}






###MS2 matching score 150
###positive
### [M+H] 50
### [M+H] isotope 20
### other adduct 20
### other adduct isotope 10

###negative
### [M-H] 50
### [M-H] isotope 20
### other adduct 20
### other adduct isotope 10

####The max score should be 500
#####The min score is 20

#' Score a Metabolite Feature Cluster (MFC)
#'
#' Computes a heuristic score for a metabolite feature cluster (`mfc`) based on MS annotation features such as
#' adduct type, polarity, isotopic composition, and MS2-supported annotation level.
#'
#' @param mfc A data frame representing a single metabolite feature cluster. Must contain at least the columns:
#'   `Adduct`, `polarity`, `isotope`, and optionally `ms2_confirmed` / `Level`.
#'
#' @return A numeric score representing the confidence or quality of the MFC annotation.
#'
#' @details
#' The scoring logic is based on heuristic rules:
#' \itemize{
#'   \item +150 points if any feature has `ms2_confirmed == TRUE`, indicating
#'         user-supplied MS2-backed confirmation
#'   \item +50 for canonical adducts: `(M+H)+` in positive mode or `(M-H)-` in negative mode
#'   \item +20 if isotopic peaks are observed for canonical adducts
#'   \item +20–10 for secondary adducts or isotopic forms
#' }
#' These rules help prioritize well-supported MFCs and can be refined in future implementations.
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' \dontrun{
#' mfc <- data.frame(
#'   Adduct = c("(M+H)+", "(M+H)+"),
#'   polarity = c("positive", "positive"),
#'   isotope = c("[M]", "[M+1]"),
#'   ms2_confirmed = c(TRUE, TRUE),
#'   stringsAsFactors = FALSE
#' )
#'
#' score_mfc(mfc)
#' }
#'


score_mfc <- function(mfc) {
  score <- 0

  if ("ms2_confirmed" %in% colnames(mfc) &&
      any(!is.na(mfc$ms2_confirmed) & mfc$ms2_confirmed)) {
    score <-
      score + 150
  }

  ##This should be optimized in the future
  ###if positive adduct are +H, score add 50
  if (any(mfc$Adduct == "(M+H)+")) {
    score <- score + 50
  }

  if (any(mfc$Adduct == "(M+H)+" &
          mfc$isotope != "[M]")) {
    score <- score + 20
  }

  if (any(mfc$polarity == "positive" &
          mfc$Adduct != "(M+H)+")) {
    score <- score + 20
  }

  if (any(mfc$polarity == "positive" & mfc$Adduct != "(M+H)+" &
          mfc$isotope != "[M]")) {
    score <- score + 10
  }

  if (any(mfc$Adduct == "(M-H)-")) {
    score <- score + 50
  }

  if (any(mfc$Adduct == "(M-H)-" &
          mfc$isotope != "[M]")) {
    score <- score + 20
  }

  if (any(mfc$polarity == "negative" &
          mfc$Adduct != "(M-H)-")) {
    score <- score + 20
  }

  if (any(mfc$polarity == "negative" & mfc$Adduct != "(M-H)-" &
          mfc$isotope != "[M]")) {
    score <- score + 10
  }

  return(score)
}
