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

  cli::cli_h1("featuremap v{version}")

  cli::cli_alert_success("Decode metabolic signals, module by module")
  cli::cli_alert_success("{.strong Follow pipeline below:}")
  cli::cli_ul(c(
    "annotate_feature()",
    "extract_subnetwork()",
    "identify_module()",
    "annotate_module()"
  ))
  
  cli::cli_alert_info("Type {.code ?featuremap} to get started!")
  cli::cli_alert_info("Please cite us: Yuting Wang,..., Xiaotao Shen.")
}
