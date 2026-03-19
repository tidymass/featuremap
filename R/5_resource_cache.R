featuremap_resource_base_url <- function() {
  getOption(
    "featuremap.resources.base_url",
    "https://raw.githubusercontent.com/tidymass/featuremap_resources/main/reference"
  )
}


featuremap_cache_dir <- function() {
  dir <- getOption("featuremap.resources.cache_dir")
  if (is.null(dir) || !nzchar(dir)) {
    dir <- tools::R_user_dir("featuremap", which = "cache")
  }
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(dir)) {
    stop(
      "[Error] Cannot create cache directory: ", dir,
      ". You can set options(featuremap.resources.cache_dir = '/path/to/cache')."
    )
  }
  dir
}


featuremap_manifest_local_path <- function() {
  file.path(featuremap_cache_dir(), "resource_manifest.csv")
}


featuremap_manifest_url <- function() {
  paste0(featuremap_resource_base_url(), "/manifest.csv")
}


featuremap_download <- function(url, destfile) {
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(destfile, ".tmp")
  on.exit(if (file.exists(tmp)) file.remove(tmp), add = TRUE)

  if (grepl("^file://", url, ignore.case = TRUE)) {
    src <- sub("^file://", "", url, ignore.case = TRUE)
    ok <- file.copy(src, tmp, overwrite = TRUE)
  } else {
    ok <- tryCatch(
      utils::download.file(url, destfile = tmp, mode = "wb", quiet = TRUE) == 0,
      error = function(e) FALSE
    )
  }

  if (!ok || !file.exists(tmp)) {
    stop("[Error] Failed to download resource from URL: ", url, ".")
  }

  file.rename(tmp, destfile)
  invisible(destfile)
}


featuremap_read_manifest <- function(refresh = FALSE) {
  local_path <- featuremap_manifest_local_path()
  remote_url <- featuremap_manifest_url()

  if (refresh || !file.exists(local_path)) {
    tryCatch(
      featuremap_download(remote_url, local_path),
      error = function(e) {
        if (!file.exists(local_path)) stop(e$message)
        warning("Failed to refresh manifest; using cached manifest: ", e$message)
      }
    )
  }

  manifest <- utils::read.csv(local_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("name", "file")
  missing_cols <- setdiff(required_cols, colnames(manifest))
  if (length(missing_cols) > 0) {
    stop(
      "[Error] Invalid manifest. Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  manifest
}


featuremap_default_resource_name <- function(kind, refdata = NULL) {
  defaults <- list(
    metabolite_database = c(
      default = "featuremap_compound",
      all = "featuremap_compound",
      hmdb = "featuremap_compound",
      kegg = "featuremap_compound"
    ),
    mrn = c(
      default = "mrn",
      all = "mrn",
      hmdb = "mrn",
      kegg = "mrn"
    ),
    functional_reference = c(hmdb = "functional_reference_hmdb", kegg = "functional_reference_kegg")
  )

  if (!kind %in% names(defaults)) {
    stop("[Error] Unknown resource kind: ", kind)
  }

  if (is.null(refdata) || !refdata %in% names(defaults[[kind]])) {
    if ("default" %in% names(defaults[[kind]])) {
      return(unname(defaults[[kind]][["default"]]))
    }
    stop(
      "[Error] No default resource mapping for kind '", kind,
      "' and refdata '", refdata, "'."
    )
  }

  unname(defaults[[kind]][[refdata]])
}


featuremap_manifest_row <- function(name, version = NULL, refresh = FALSE) {
  manifest <- featuremap_read_manifest(refresh = refresh)
  cand <- manifest[manifest$name == name, , drop = FALSE]
  if (nrow(cand) == 0) {
    stop(
      "[Error] Resource '", name, "' not found in manifest: ",
      featuremap_manifest_url()
    )
  }

  if (!is.null(version)) {
    if (!"version" %in% colnames(cand)) {
      stop("[Error] Manifest has no 'version' column, cannot filter by version.")
    }
    cand <- cand[cand$version == version, , drop = FALSE]
    if (nrow(cand) == 0) {
      stop("[Error] Resource '", name, "' with version '", version, "' not found.")
    }
  }

  if ("version" %in% colnames(cand)) {
    cand <- cand[order(cand$version, decreasing = TRUE), , drop = FALSE]
  }

  cand[1, , drop = FALSE]
}


featuremap_resource_cache_path <- function(name, file_name, version = NULL) {
  version_tag <- if (is.null(version) || is.na(version) || !nzchar(version)) "latest" else version
  file.path(featuremap_cache_dir(), "resources", name, version_tag, basename(file_name))
}


featuremap_optional_hash_check <- function(path, row) {
  if ("md5" %in% colnames(row) && !is.na(row$md5) && nzchar(row$md5)) {
    md5_local <- unname(tools::md5sum(path))
    if (!identical(tolower(md5_local), tolower(row$md5))) {
      stop("[Error] MD5 mismatch for cached resource: ", path)
    }
  }
  invisible(TRUE)
}


featuremap_fetch_resource <- function(name, version = NULL, force = FALSE, refresh_manifest = FALSE) {
  row <- featuremap_manifest_row(name = name, version = version, refresh = refresh_manifest)
  version_tag <- if ("version" %in% colnames(row)) as.character(row$version[[1]]) else NULL
  remote_file <- as.character(row$file[[1]])

  url <- if ("url" %in% colnames(row) && !is.na(row$url[[1]]) && nzchar(row$url[[1]])) {
    as.character(row$url[[1]])
  } else {
    paste0(featuremap_resource_base_url(), "/", remote_file)
  }

  local_path <- featuremap_resource_cache_path(name, remote_file, version_tag)
  if (!file.exists(local_path) || isTRUE(force)) {
    message("Downloading resource: ", name, " -> ", local_path)
    featuremap_download(url, local_path)
  }

  featuremap_optional_hash_check(local_path, row)
  local_path
}


featuremap_load_resource <- function(name, version = NULL, force = FALSE, refresh_manifest = FALSE) {
  path <- featuremap_fetch_resource(
    name = name,
    version = version,
    force = force,
    refresh_manifest = refresh_manifest
  )

  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    return(readRDS(path))
  }

  if (ext %in% c("rda", "rdata")) {
    e <- new.env(parent = emptyenv())
    nm <- load(path, envir = e)
    if (length(nm) == 1) return(e[[nm]])
    return(mget(nm, envir = e))
  }

  path
}
