# ============================================================
# 1. 辅助工具函数
# ============================================================

#' Clean Field Values
#'
#' Splits and cleans multi-value character fields, removing blanks,
#' whitespace, duplicates, and single-character entries.
#'
#' @param x Character vector or scalar. Raw field value, possibly
#'   containing multiple entries separated by \code{\{\}} delimiters,
#'   \code{NA} values, or extra whitespace.
#'
#' @return Character vector of cleaned, unique, non-empty values with
#'   more than one character. Returns \code{character(0)} if input is
#'   \code{NULL}, empty, or all \code{NA}.
#'

clean_field <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    return(character(0))
  }
  x <- as.character(x)
  if (length(x) > 1) {
    parts <- x
  } else {
    parts <- unlist(strsplit(x, "\\{\\}"))
  }
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  parts <- unique(parts)
  parts <- parts[nchar(parts) > 1]
  return(parts)
}


#' Escape Special Regex Characters
#'
#' Escapes all non-word characters in a string so it can be safely
#' used as a literal pattern inside \code{grepl} or \code{stringr}
#' regex calls.
#'
#' @param x Character string. The raw text to escape.
#'
#' @return Character string with all non-word characters preceded by
#'   a double backslash.
#'

escape_regex <- function(x) {
  stringr::str_replace_all(x, "([\\W])", "\\\\\\1")
}


# ============================================================
# 2. 同义词过滤
# ============================================================

#' Filter Safe Synonyms
#'
#' Removes synonyms that are too short, purely uppercase abbreviations,
#' or very short mixed-case tokens that would produce noisy PubMed
#' queries (e.g. "ATP", "Glc", "H2O").
#'
#' @param synonyms Character vector. Candidate synonym strings to filter.
#'   Typically comes from a metabolite database field such as HMDB synonyms.
#'
#' @return Character vector of filtered synonyms. Returns \code{character(0)}
#'   if \code{synonyms} is \code{NULL}, empty, or all entries are removed
#'   by the filtering rules.
#'
#' @details
#' The following entries are removed:
#' \itemize{
#'   \item Duplicates and leading/trailing whitespace (via \code{unique} +
#'     \code{trimws}).
#'   \item Strings with 3 or fewer characters.
#'   \item Strings matching purely uppercase abbreviations of 1–4 letters
#'     (e.g. \code{"ATP"}, \code{"NAD"}).
#'   \item Strings matching a single uppercase letter followed by exactly
#'     two lowercase letters (e.g. \code{"Glc"}, \code{"Ser"}).
#'   \item Strings composed entirely of 1–4 letters in any case.
#' }
#'

#'
#' @examples
#' \dontrun{
#' filter_safe_synonyms(c("ATP", "D-glucose", "Glc", "dextrose", "H2O"))
#' }
filter_safe_synonyms <- function(synonyms) {

  if (is.null(synonyms) || length(synonyms) == 0)
    return(character(0))

  synonyms <- unique(trimws(synonyms))
  synonyms <- synonyms[nchar(synonyms) > 3]
  synonyms <- synonyms[!(grepl("^[A-Z]{1,4}$", synonyms))]
  synonyms <- synonyms[!(grepl("^[A-Z][a-z]{2}$", synonyms))]
  synonyms <- synonyms[!(grepl("^[A-Za-z]{1,4}$", synonyms))]

  return(synonyms)
}


# ============================================================
# 3. 数据结构转换
# ============================================================

#' Convert Module Data Frame to Metabolite List
#'
#' Transforms a module-level data frame into a named list with one
#' entry per unique metabolite, ready for downstream PubMed annotation.
#'
#' @param module_df Data frame. Metabolite annotation table. Required
#'   columns depend on \code{refdata}:
#'   \itemize{
#'     \item \strong{hmdb}: \code{HMDB_ID}, \code{name}, \code{synonyms}.
#'     \item \strong{kegg}: \code{KEGG_ID}, \code{name}, \code{synonyms}.
#'   }
#' @param refdata Character string. Reference database used to identify
#'   metabolites. Must be one of \code{"hmdb"} (default) or \code{"kegg"}.
#'   Controls which ID column is used as the list key and for deduplication.
#'
#' @return Named list where each element is named by its database ID
#'   (\code{HMDB_ID} or \code{KEGG_ID}) and contains:
#'   \itemize{
#'     \item \code{HMDB_ID} / \code{KEGG_ID} — the metabolite identifier.
#'     \item \code{name} — primary metabolite name.
#'     \item \code{synonyms} — cleaned synonym vector (via \code{clean_field}).
#'   }
#'
#' @details
#' When duplicate IDs are present in \code{module_df}, only the first row
#' for each ID is retained (\code{dplyr::slice(1)}).
#'

#'
#' @examples
#' \dontrun{
#' module_df <- data.frame(
#'   HMDB_ID = c("HMDB00001", "HMDB00002"),
#'   KEGG_ID = c("C00031", "C00022"),
#'   name = c("Glucose", "Pyruvate"),
#'   synonyms = c("D-Glucose;Dextrose", "2-Oxopropanoate"),
#'   stringsAsFactors = FALSE
#' )
#' convert_module_to_met_list(module_df, refdata = "hmdb")
#' convert_module_to_met_list(module_df, refdata = "kegg")
#' }
convert_module_to_met_list <- function(module_df, refdata = c("hmdb", "kegg")) {

  refdata <- match.arg(refdata)

  if (refdata == "hmdb") {
    hmdb_ids <- unique(module_df$HMDB_ID)
    met_list <- lapply(hmdb_ids, function(id) {
      df <- module_df %>% dplyr::filter(HMDB_ID == id) %>% dplyr::slice(1)
      list(
        HMDB_ID  = df$HMDB_ID,
        name     = df$name,
        synonyms = clean_field(df$synonyms)
      )
    })
    names(met_list) <- hmdb_ids

  } else if (refdata == "kegg") {
    kegg_ids <- unique(module_df$KEGG_ID)
    met_list <- lapply(kegg_ids, function(id) {
      df <- module_df %>% dplyr::filter(KEGG_ID == id) %>% dplyr::slice(1)
      list(
        KEGG_ID  = df$KEGG_ID,
        name     = df$name,
        synonyms = clean_field(df$synonyms)
      )
    })
    names(met_list) <- kegg_ids
  }

  return(met_list)
}


# ============================================================
# 4. PubMed 查询构建
# ============================================================

#' Build PubMed Query String
#'
#' Constructs a structured PubMed Boolean query for a given metabolite
#' by combining its primary name, optional synonyms, metabolism-related
#' topic terms, an optional phenotype constraint, and a publication date
#' range.
#'
#' @param met_name Character string. Primary name of the metabolite
#'   (e.g. \code{"glucose"}). Always included in the name search block.
#' @param synonyms Character vector or \code{NULL}. Alternative names or
#'   synonyms for the metabolite. Passed through \code{filter_safe_synonyms}
#'   before use; \code{NULL} or empty vectors are silently ignored.
#' @param phenotype Character string or \code{NULL}. Optional disease or
#'   phenotype term to restrict results (e.g. \code{"type 2 diabetes"}).
#'   When provided and longer than one character, added as an
#'   \code{AND "[Title/Abstract]"} constraint. Ignored when \code{NULL}.
#' @param years Numeric scalar. Number of years to look back from the
#'   current calendar year when setting the \code{pdat} (publication date)
#'   filter. For example, \code{years = 5} retrieves articles published
#'   in the last five years.
#'
#' @return Character string containing the fully formatted PubMed query,
#'   ready to pass to \code{rentrez::entrez_search}.
#'
#' @details
#' The query is assembled from up to four AND-joined blocks:
#' \enumerate{
#'   \item \strong{Name block}: metabolite name and synonyms joined with OR,
#'     each wrapped in \code{"..."[Title/Abstract]}.
#'   \item \strong{Topic block}: fixed metabolism-related terms
#'     (\code{metabolism}, \code{biochemical}, \code{pathway},
#'     \code{reaction}) joined with OR.
#'   \item \strong{Phenotype block} (optional): single
#'     \code{"phenotype"[Title/Abstract]} term.
#'   \item \strong{Date block}: \code{(start_year:current_year[pdat])}.
#' }
#'

#'
#' @examples
#' \dontrun{
#' build_pubmed_query("glucose", synonyms = c("dextrose", "D-glucose"), years = 3)
#' build_pubmed_query("lactate", phenotype = "sepsis", years = 5)
#' }
build_pubmed_query <- function(met_name,
                               synonyms  = NULL,
                               phenotype = NULL,
                               years     = NULL) {

  ## --- metabolite name block ---
  synonyms  <- filter_safe_synonyms(synonyms)
  name_part <- c(met_name, synonyms) %>%
    unique() %>%
    purrr::map_chr(~ sprintf('"%s"[Title/Abstract]', .x)) %>%
    paste(collapse = " OR ")

  ## --- metabolism-related topic block ---
  topic_part <- paste(
    c("metabolism", "biochemical", "pathway", "reaction") %>%
      sprintf('%s[Title/Abstract]', .),
    collapse = " OR "
  )

  ## --- optional phenotype constraint ---
  phenotype_part <- NULL
  if (!is.null(phenotype) && nchar(phenotype) > 1) {
    phenotype_part <- sprintf('"%s"[Title/Abstract]', phenotype)
  }

  ## --- publication date range ---
  start_year <- as.numeric(format(Sys.Date(), "%Y")) - years
  year_part  <- sprintf("(%d:%d[pdat])",
                        start_year,
                        as.numeric(format(Sys.Date(), "%Y")))

  ## --- assemble final query ---
  query_parts <- c(
    sprintf("(%s)", name_part),
    sprintf("(%s)", topic_part),
    phenotype_part,
    year_part
  )
  query_parts <- query_parts[!is.null(query_parts)]
  paste(query_parts, collapse = " AND ")
}


# ============================================================
# 5. PubMed 检索与记录获取
# ============================================================

#' Search PubMed Database
#'
#' Executes a search query against the PubMed database with automatic
#' retry logic and rate-limiting delays to comply with NCBI API policies.
#'
#' @param query Character string. A valid PubMed Boolean query string,
#'   typically produced by \code{build_pubmed_query}.
#' @param retmax Numeric scalar. Maximum number of PMIDs to return per
#'   search request (default: \code{10}). Passed directly to
#'   \code{rentrez::entrez_search}.
#' @param retries Numeric scalar. Maximum number of retry attempts when
#'   the NCBI request fails (default: \code{3}). Each failed attempt
#'   increments the counter before retrying.
#' @param NCBI_api_key Character string or \code{NULL}. Personal NCBI API
#'   key registered at \url{https://www.ncbi.nlm.nih.gov/account/}.
#'   When provided, raises the NCBI rate limit from 3 to 10 requests per
#'   second. Pass \code{NULL} (default) to use the anonymous limit.
#'
#' @return Character vector of PubMed IDs (PMIDs) matching the query.
#'   Returns \code{character(0)} if all retry attempts fail or no results
#'   are found.
#'
#' @details
#' A 0.5-second sleep is inserted before every request to stay within
#' the NCBI rate limit of 3 requests/second for unauthenticated users.
#' If an API key is supplied the delay can be safely reduced, but is kept
#' conservative here for robustness.
search_pubmed <- function(query,
                          retmax       = 10,
                          retries      = 3,
                          NCBI_api_key = NULL) {

  if (!is.null(NCBI_api_key)) rentrez::set_entrez_key(NCBI_api_key)

  attempt <- 1
  while (attempt <= retries) {
    Sys.sleep(0.5)
    res_pub <- try(
      rentrez::entrez_search(db = "pubmed", term = query, retmax = retmax),
      silent = TRUE
    )
    if (!inherits(res_pub, "try-error")) {
      return(res_pub$ids)
    }
    attempt <- attempt + 1
  }

  return(character(0))
}


#' Fetch PubMed Records
#'
#' Retrieves full article titles and abstracts from PubMed for a set of
#' PMIDs by downloading and parsing the corresponding XML records.
#'
#' @param pmids Character vector. One or more PubMed IDs (PMIDs) to fetch,
#'   as returned by \code{search_pubmed}. Returns empty lists immediately
#'   when \code{length(pmids) == 0}.
#' @param NCBI_api_key Character string or \code{NULL}. NCBI API key for
#'   elevated rate limits (see \code{search_pubmed} for details).
#'   Pass \code{NULL} (default) to use the anonymous limit.
#'
#' @return Named list with two character vectors:
#'   \itemize{
#'     \item \code{titles} — article titles, one element per PMID; empty
#'       string when the title node is absent.
#'     \item \code{abstracts} — full abstract text (structured abstract
#'       sections are concatenated with a space); empty string when no
#'       abstract is present.
#'   }
#'
#' @details
#' Records are fetched in a single \code{rentrez::entrez_fetch} call using
#' \code{rettype = "xml"}. XML parsing is performed with \pkg{xml2}.
#' When a \code{PubmedArticle} node lacks an \code{AbstractText} child,
#' the corresponding element is set to \code{""} to maintain positional
#' alignment between \code{titles} and \code{abstracts}.
fetch_pubmed_records <- function(pmids, NCBI_api_key = NULL) {

  if (length(pmids) == 0) {
    return(list(titles = character(0), abstracts = character(0)))
  }

  if (!is.null(NCBI_api_key)) rentrez::set_entrez_key(NCBI_api_key)

  rec <- rentrez::entrez_fetch(
    db      = "pubmed",
    id      = pmids,
    rettype = "xml",
    parsed  = FALSE
  )

  doc      <- xml2::read_xml(rec)
  articles <- xml2::xml_find_all(doc, ".//PubmedArticle")

  titles <- sapply(articles, function(a) {
    t <- xml2::xml_find_first(a, ".//ArticleTitle")
    if (is.na(t)) return("")
    xml2::xml_text(t)
  })

  abstracts <- sapply(articles, function(a) {
    abs_nodes <- xml2::xml_find_all(a, ".//Abstract/AbstractText")
    if (length(abs_nodes) == 0) return("")
    paste(xml2::xml_text(abs_nodes), collapse = " ")
  })

  return(list(titles = titles, abstracts = abstracts))
}


# ============================================================
# 6. 摘要相关性过滤
# ============================================================

#' Check Whether an Abstract Mentions a Metabolite
#'
#' Tests whether a PubMed abstract contains a whole-word match for a
#' metabolite name or any of its synonyms, using case-insensitive
#' boundary-anchored regex. Used to discard off-topic articles after
#' a PubMed search.
#'
#' @param abstract Character string. Full abstract text to search within.
#'   Returns \code{FALSE} immediately when \code{abstract} is \code{NA}
#'   or fewer than 30 characters long.
#' @param met_name Character string. Primary metabolite name to search for
#'   (e.g. \code{"glucose"}).
#' @param synonyms Character vector or \code{NULL}. Additional synonym
#'   strings to search for alongside \code{met_name}. \code{NA} entries
#'   and strings of length ≤ 1 are silently dropped before matching.
#'
#' @return Logical scalar. \code{TRUE} if at least one term produces a
#'   whole-word match in the abstract; \code{FALSE} otherwise.
#'
#' @details
#' All comparisons are performed in lower case. Special regex characters
#' in each term are escaped via \code{escape_regex} before pattern
#' construction, so metabolite names containing hyphens, parentheses,
#' or Greek letters are handled correctly.
abstract_mentions_metabolite <- function(abstract,
                                         met_name,
                                         synonyms = NULL) {

  if (is.na(abstract) || nchar(abstract) < 30) {
    return(FALSE)
  }

  text  <- tolower(abstract)
  terms <- tolower(c(met_name, synonyms))
  terms <- terms[!is.na(terms) & nchar(terms) > 1]

  any(sapply(terms, function(t) {
    pattern <- paste0("\\b", escape_regex(t), "\\b")
    grepl(pattern, text, perl = TRUE)
  }))
}


# ============================================================
# 7. 单个代谢物注释
# ============================================================

#' Annotate a Single Metabolite with PubMed Literature
#'
#' Runs the complete PubMed annotation pipeline for one metabolite:
#' query construction → search → record fetch → abstract filtering.
#'
#' @param met_info Named list. Metabolite metadata, typically one element
#'   from the list produced by \code{convert_module_to_met_list}. Must
#'   contain at minimum:
#'   \itemize{
#'     \item \code{name} — primary metabolite name (character string).
#'     \item \code{synonyms} — synonym vector (character or \code{NULL});
#'       processed through \code{clean_field} internally.
#'   }
#' @param phenotype Character string or \code{NULL}. Optional phenotype /
#'   disease term appended to the query (e.g. \code{"insulin resistance"}).
#'   Passed unchanged to \code{build_pubmed_query}. Use \code{NULL}
#'   (default) for a general metabolism query.
#' @param NCBI_api_key Character string or \code{NULL}. NCBI API key for
#'   elevated request rate limits. Forwarded to \code{search_pubmed} and
#'   \code{fetch_pubmed_records}. Use \code{NULL} (default) for anonymous
#'   access.
#' @param retmax Numeric scalar. Maximum number of PubMed results to
#'   retrieve for this metabolite (default: \code{10}). Passed to
#'   \code{search_pubmed}.
#' @param years Numeric scalar. Number of years to look back when building
#'   the publication date filter. Passed to \code{build_pubmed_query}.
#'
#' @return Named list with three elements:
#'   \itemize{
#'     \item \code{pmids} — character vector of PMIDs that passed the
#'       abstract relevance filter.
#'     \item \code{titles} — character vector of corresponding article
#'       titles.
#'     \item \code{abstracts} — character vector of corresponding
#'       abstracts.
#'   }
#'   All three vectors are \code{character(0)} when no relevant articles
#'   are found or when the search or fetch step fails.
#'
#' @details
#' After retrieving raw PubMed records, each abstract is passed to
#' \code{abstract_mentions_metabolite} to confirm that the metabolite
#' name or a synonym actually appears in the text, removing articles that
#' matched only on topic terms.
annotate_one_metabolite_pubmed <- function(met_info,
                                           phenotype    = NULL,
                                           NCBI_api_key = NULL,
                                           retmax       = 10,
                                           years) {

  query <- build_pubmed_query(
    met_name  = met_info$name,
    synonyms  = clean_field(met_info$synonyms),
    phenotype = phenotype,
    years     = years
  )

  pmids <- search_pubmed(query, retmax = retmax, NCBI_api_key = NCBI_api_key)

  recs <- tryCatch(
    fetch_pubmed_records(pmids),
    error = function(e) list(titles = character(0), abstracts = character(0))
  )

  if (length(pmids) == 0 || length(recs$abstracts) == 0) {
    return(list(pmids     = character(0),
                titles    = character(0),
                abstracts = character(0)))
  }

  keep <- mapply(
    abstract_mentions_metabolite,
    abstract  = recs$abstracts,
    MoreArgs  = list(met_name = met_info$name,
                     synonyms = clean_field(met_info$synonyms))
  )

  filtered_abstracts <- recs$abstracts[keep]
  filtered_titles    <- recs$titles[keep]
  filtered_pmids     <- pmids[keep]

  contexts <- mapply(
    extract_metabolite_context,
    abstract = filtered_abstracts,
    MoreArgs = list(
      met_name = met_info$name,
      synonyms = clean_field(met_info$synonyms),
      window   = 1
    )
  )

  contexts <- as.character(contexts)

  list(
    pmids        = filtered_pmids,
    titles       = filtered_titles,
    abstracts    = filtered_abstracts,
    contexts     = contexts
  )
}


extract_metabolite_context <- function(abstract,
                                       met_name,
                                       synonyms = NULL,
                                       window = 1) {

  if (is.na(abstract) || nchar(abstract) < 30) {
    return(NA_character_)
  }

  text <- tolower(abstract)
  safe_syn <- filter_safe_synonyms(synonyms)

  terms <- tolower(c(met_name, safe_syn))
  terms <- terms[!is.na(terms) & nchar(terms) > 1]

  sentences <- unlist(strsplit(abstract, "(?<=[.!?])\\s+", perl = TRUE))

  hit_idx <- which(sapply(sentences, function(s) {
    any(sapply(terms, function(t) {
      pattern <- paste0("\\b", escape_regex(t), "\\b")
      grepl(pattern, tolower(s), perl = TRUE)
    }))
  }))

  if (length(hit_idx) == 0) {
    return(NA_character_)
  }

  keep_idx <- unique(unlist(lapply(hit_idx, function(i) {
    seq(max(1, i - window), min(length(sentences), i + window))
  })))

  paste(sentences[keep_idx], collapse = "; ")
}


# ============================================================
# 8. 模块级批量注释（入口函数）
# ============================================================

#' Annotate All Metabolites in a Module with PubMed Literature
#'
#' Iterates over every metabolite in a module annotation list and runs
#' the full PubMed annotation pipeline for each one, collecting results
#' into a single named list.
#'
#' @param module_anno Named list. Each element represents one metabolite
#'   and must be structured as described in \code{annotate_one_metabolite_pubmed}
#'   (i.e. contains at minimum \code{name} and optionally \code{synonyms}).
#'   Typically produced by \code{convert_module_to_met_list}.
#' @param phenotype Character string or \code{NULL}. Optional phenotype
#'   constraint applied uniformly to every metabolite query
#'   (e.g. \code{"obesity"}). Pass \code{NULL} (default) to omit.
#' @param NCBI_api_key Character string or \code{NULL}. NCBI API key
#'   shared across all requests. Pass \code{NULL} (default) for anonymous
#'   access (≤ 3 requests/second).
#' @param retmax Numeric scalar. Maximum number of PubMed hits to retrieve
#'   per metabolite (default: \code{10}).
#' @param years Numeric scalar. Number of years to look back for the
#'   publication date filter applied to every query.
#'
#' @return Named list with one element per metabolite (names match
#'   \code{names(module_anno)}). Each element is the list returned by
#'   \code{annotate_one_metabolite_pubmed}: \code{pmids}, \code{titles},
#'   and \code{abstracts}.
#'
#' @details
#' Processing is intentionally sequential (via \code{lapply}) rather than
#' parallel to respect NCBI rate limits. If parallel execution is required,
#' supply a valid \code{NCBI_api_key} and implement appropriate throttling
#' (≤ 10 requests/second per key).
#'

#'
#' @examples
#' \dontrun{
#' module_df <- data.frame(
#'   HMDB_ID = c("HMDB00001", "HMDB00002"),
#'   name = c("Glucose", "Pyruvate"),
#'   synonyms = c("D-Glucose;Dextrose", "2-Oxopropanoate"),
#'   stringsAsFactors = FALSE
#' )
#' met_list <- convert_module_to_met_list(module_df, refdata = "hmdb")
#' results  <- annotate_pubmed_for_module(met_list, phenotype = "diabetes", years = 5)
#' }
annotate_pubmed_for_module <- function(module_anno,
                                       phenotype    = NULL,
                                       NCBI_api_key = NULL,
                                       retmax       = 10,
                                       years) {

  mets <- names(module_anno)

  res_mets <- lapply(mets, function(m) {
    annotate_one_metabolite_pubmed(
      met_info     = module_anno[[m]],
      phenotype    = phenotype,
      NCBI_api_key = NCBI_api_key,
      retmax       = retmax,
      years        = years
    )
  })

  names(res_mets) <- mets
  res_mets
}
