.onAttach <- function(libname, pkgname) {
  version <- as.character(utils::packageVersion(pkgname))

  # Silence startup banner in non-interactive sessions (e.g. future workers)
  # or when users explicitly disable it.
  if (!interactive() || isTRUE(getOption("featuremap.quiet", FALSE))) {
    return(invisible())
  }

  # Avoid repeated banner printing in the same R session.
  if (isTRUE(getOption("featuremap.banner_shown", FALSE))) {
    return(invisible())
  }
  options(featuremap.banner_shown = TRUE)

  cli::cli_h1("featuremap {version}")
  cli::cli_verbatim(paste(
    "Metabolite annotation",
    "Network extraction",
    "Module identification",
    "Function interpretation",
    sep = " \u2192 "
  ))
  cli::cli_alert_info("Run {.code ?featuremap} to get started.")
}
