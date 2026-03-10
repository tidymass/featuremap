ensure_featuremap_object <- function(res,
                                     expected_stages = NULL,
                                     allow_missing_stage = TRUE) {
  if (!is.null(res$object_type) && !identical(res$object_type, "featuremap_result")) {
    stop("[Error] Input object type is invalid. Expected object_type = 'featuremap_result'.")
  }

  if (is.null(res$object_type)) {
    res$object_type <- "featuremap_result"
  }

  if (!is.null(expected_stages)) {
    if (is.null(res$stage)) {
      if (!allow_missing_stage) {
        stop("[Error] Input object has no stage marker.")
      }
    } else if (!res$stage %in% expected_stages) {
      stop(
        "[Error] Input object stage is '", res$stage,
        "', but expected one of: ",
        paste(expected_stages, collapse = ", "), "."
      )
    }
  }

  res
}

set_featuremap_stage <- function(res, stage) {
  res$object_type <- "featuremap_result"
  res$stage <- stage
  res
}

