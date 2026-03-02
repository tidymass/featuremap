.onAttach <- function(libname, pkgname) {
  version <- as.character(utils::packageVersion(pkgname))
  
  cli::cli_h1("featureMAP v{0.99}")
  
  cli::cli_alert_success("Decode metabolic signals, module by module")
  cli::cli_alert_success("{.strong Follow pipeline below:}")
  cli::cli_ul(c(
    "annotate_feature()",
    "extract_subnetwork()",
    "identify_module()",
    "annotate_module()"
  ))
  
  cli::cli_alert_info("Type {.code ?featureMAP} to get started!")
  cli::cli_alert_info("Please cite us: Yuting Wang,..., Xiaotao Shen.")
}