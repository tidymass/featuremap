#' Annotate a Single Module Using LLM
#'
#' Sends a single module's data to an LLM API (e.g., GPT) to generate a biological
#' interpretation including a title, summary, confidence label, and rationale.
#'
#' @param module_obj List containing module data with fields:
#'   - `module_id`: module identifier (for example `"M1"`).
#'   - `module_stats`: module-level statistics (Importance, Robustness, Annotated_fraction).
#'   - `top_N_metabolites`: data frame of top metabolites with
#'     `HMDB_ID`, `Compound_name`, `Importance`, `Robustness`, `class`,
#'     and optional `functional_summary`.
#'   - `other_metabolites`: data frame of additional metabolites aggregated at class level.
#' @param study_context List containing study metadata:
#'   - `biospecimen`: biological specimen type.
#'   - `phenotype`: phenotype of interest.
#'   - `design`: study design description (optional).
#' @param provider Character. LLM provider passed to `chat_api_call()`
#'   (e.g. `"openai"`, `"siliconflow"`, `"gemini"`).
#' @param llm_api_key Character. API key for the LLM service
#' @param model Character. Model name (e.g., "gpt-4", "gpt-3.5-turbo")
#' @param max_tokens Integer. Maximum number of tokens in the response (default: 800)
#' @param temperature Numeric. Sampling temperature for the LLM (default: 0 for deterministic output)
#'
#' @return A data frame with one row containing `module_id`, `title`,
#'   `summary`, `confidence_label`, `rationale`, and `raw_response`.
#'
#' @details This function:
#'   \enumerate{
#'     \item Constructs system and user prompts using \code{build_system_prompt()}
#'           and \code{build_module_prompt()}
#'     \item Sends the prompts to the LLM API via \code{gpt_api_call()}
#'     \item Parses the JSON response into structured fields
#'     \item Handles errors gracefully by returning NA values when parsing fails
#'   }
#'
#' @note If the API returns an empty response or JSON parsing fails, the function
#'   returns a data frame with NA values and logs a message.
annotate_one_module <- function(module_obj,
                                study_context,
                                provider,
                                llm_api_key,
                                model,
                                max_tokens  = 5000,
                                temperature = 0) {
  
  # 1) Build system prompt with global instructions
  system_prompt <- build_system_prompt()
  
  # 2) Build user prompt with module-specific data
  user_prompt <- build_module_prompt(module_obj, study_context)
  
  normalize_parsed <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.data.frame(x) && nrow(x) >= 1) return(as.list(x[1, , drop = FALSE]))
    if (is.list(x) && !is.null(names(x))) return(x)
    if (is.list(x) && length(x) >= 1 && is.list(x[[1]]) && !is.null(names(x[[1]]))) {
      return(x[[1]])
    }
    NULL
  }

  parse_json_response <- function(txt) {
    direct <- tryCatch(
      jsonlite::fromJSON(txt, simplifyVector = FALSE),
      error = function(e) NULL
    )
    direct <- normalize_parsed(direct)
    if (!is.null(direct)) return(direct)

    json_obj_txt <- extract_json_object(txt)
    if (!is.null(json_obj_txt)) {
      obj <- tryCatch(
        jsonlite::fromJSON(json_obj_txt, simplifyVector = FALSE),
        error = function(e) NULL
      )
      obj <- normalize_parsed(obj)
      if (!is.null(obj)) return(obj)
    }

    json_fence <- regmatches(txt, regexpr("```\\s*[\\s\\S]*?\\s*```", txt, perl = TRUE))
    if (length(json_fence) > 0 && nzchar(json_fence[1])) {
      inner <- sub("^```[a-zA-Z]*\\s*", "", json_fence[1], perl = TRUE)
      inner <- sub("\\s*```$", "", inner, perl = TRUE)
      fobj <- tryCatch(
        jsonlite::fromJSON(inner, simplifyVector = FALSE),
        error = function(e) NULL
      )
      fobj <- normalize_parsed(fobj)
      if (!is.null(fobj)) return(fobj)
    }

    NULL
  }

  parse_labeled_fallback <- function(txt) {
    if (is.null(txt) || !nzchar(txt)) return(NULL)
    lines <- unlist(strsplit(txt, "\n", fixed = TRUE))
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
    if (length(lines) == 0) return(NULL)

    get_value <- function(pattern) {
      hit <- grep(pattern, lines, perl = TRUE, ignore.case = TRUE)
      if (length(hit) == 0) return(NA_character_)
      val <- sub(pattern, "", lines[hit[1]], perl = TRUE, ignore.case = TRUE)
      val <- trimws(val)
      if (!nzchar(val)) NA_character_ else val
    }

    title <- get_value("^\\s*(title|module\\s*title)\\s*[:\\uFF1A]\\s*")
    summary <- get_value("^\\s*summary\\s*[:\\uFF1A]\\s*")
    rationale <- get_value("^\\s*rationale\\s*[:\\uFF1A]\\s*")
    conf_line <- get_value("^\\s*confidence(_label)?\\s*[:\\uFF1A]\\s*")

    confidence_label <- NA_character_
    if (!is.na(conf_line)) {
      if (grepl("\\bhigh\\b", conf_line, ignore.case = TRUE)) confidence_label <- "high"
      if (grepl("\\bmedium\\b", conf_line, ignore.case = TRUE)) confidence_label <- "medium"
      if (grepl("\\blow\\b", conf_line, ignore.case = TRUE)) confidence_label <- "low"
    }

    found_n <- sum(!is.na(c(title, summary, confidence_label, rationale)))
    if (found_n == 0) return(NULL)

    list(
      title = title,
      summary = summary,
      confidence_label = confidence_label,
      rationale = rationale
    )
  }

  parse_partial_json_fallback <- function(txt) {
    if (is.null(txt) || !nzchar(txt)) return(NULL)

    extract_json_string_value <- function(key) {
      pattern <- paste0('"', key, '"\\s*:\\s*"')
      m <- regexpr(pattern, txt, perl = TRUE)
      if (m[1] == -1) return(NA_character_)

      start <- m[1] + attr(m, "match.length")
      if (start > nchar(txt)) return(NA_character_)

      chars <- strsplit(substr(txt, start, nchar(txt)), "", fixed = TRUE)[[1]]
      if (length(chars) == 0) return(NA_character_)

      out <- character(0)
      esc <- FALSE
      closed <- FALSE
      for (ch in chars) {
        if (!esc && ch == "\\") {
          esc <- TRUE
          out <- c(out, ch)
          next
        }
        if (!esc && ch == "\"") {
          closed <- TRUE
          break
        }
        out <- c(out, ch)
        esc <- FALSE
      }

      val <- paste(out, collapse = "")
      val <- gsub("\\\\n", " ", val, perl = TRUE)
      val <- gsub("\\\\\"", "\"", val, perl = TRUE)
      val <- trimws(val)
      if (!nzchar(val)) return(NA_character_)
      if (!closed && nchar(val) < 5) return(NA_character_)
      val
    }

    title <- extract_json_string_value("title")
    summary <- extract_json_string_value("summary")
    rationale <- extract_json_string_value("rationale")
    confidence_label <- extract_json_string_value("confidence_label")

    if (!is.na(confidence_label)) {
      if (grepl("\\bhigh\\b", confidence_label, ignore.case = TRUE)) confidence_label <- "high"
      if (grepl("\\bmedium\\b", confidence_label, ignore.case = TRUE)) confidence_label <- "medium"
      if (grepl("\\blow\\b", confidence_label, ignore.case = TRUE)) confidence_label <- "low"
      if (!confidence_label %in% c("high", "medium", "low")) confidence_label <- NA_character_
    }

    found_n <- sum(!is.na(c(title, summary, confidence_label, rationale)))
    if (found_n == 0) return(NULL)

    list(
      title = title,
      summary = summary,
      confidence_label = confidence_label,
      rationale = rationale
    )
  }

  messages <- list(
    list(role = "system", content = system_prompt),
    list(role = "user",   content = user_prompt)
  )
  gemini_schema <- NULL
  if (identical(provider, "gemini")) {
    gemini_schema <- list(
      type = "OBJECT",
      properties = list(
        title = list(type = "STRING"),
        summary = list(type = "STRING"),
        confidence_label = list(type = "STRING", enum = c("high", "medium", "low")),
        rationale = list(type = "STRING")
      ),
      required = c("title", "summary", "confidence_label", "rationale")
    )
  }

  parsed <- NULL
  raw_text <- NA_character_
  current_max_tokens <- if (identical(provider, "gemini")) {
    max(4000, max_tokens)
  } else {
    max_tokens
  }

  for (attempt in 1:2) {
    res_call <- safe_chat_call(
      messages = messages,
      provider = provider,
      api_key = llm_api_key,
      model = model,
      max_tokens = current_max_tokens,
      temperature = temperature,
      response_schema = gemini_schema
    )

    if (is.null(res_call) ||
        is.null(res_call$content) ||
        !is.character(res_call$content) ||
        length(res_call$content) != 1 ||
        !nzchar(res_call$content)) {
      message("Empty or NULL response for module ", module_obj$module_id, " (attempt ", attempt, ")")
    } else {
      raw_text <- res_call$content
      parsed <- parse_json_response(raw_text)
      if (!is.null(parsed)) break
      message("No JSON found for module ", module_obj$module_id, " (attempt ", attempt, ")")

      if (identical(provider, "gemini") &&
          !is.null(res_call$finish_reason) &&
          identical(toupper(as.character(res_call$finish_reason)), "MAX_TOKENS")) {
        current_max_tokens <- min(max(3000, current_max_tokens * 2), 8192)
        message("Gemini output was truncated (MAX_TOKENS). Retrying with max_tokens = ", current_max_tokens)
      }
    }

    if (attempt == 1) {
      messages <- c(
        messages,
        list(list(
          role = "user",
          content = "Return ONLY a valid JSON object with keys: title, summary, confidence_label, rationale."
        ))
      )
    }
  }
  
  # Handle JSON parsing failure
  if (is.null(parsed)) {
    parsed <- parse_partial_json_fallback(raw_text)
  }

  if (is.null(parsed)) {
    parsed <- parse_labeled_fallback(raw_text)
  }

  if (is.null(parsed)) {
    return(
      data.frame(
        module_id        = module_obj$module_id,
        title            = NA_character_,
        summary          = NA_character_,
        confidence_label = NA_character_,
        rationale        = NA_character_,
        raw_response     = raw_text,
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Extract fields from parsed JSON, with fallback to NA
  title  <- parsed[["title"]]            %||% NA_character_
  summ   <- parsed[["summary"]]          %||% NA_character_
  clabel <- parsed[["confidence_label"]] %||% NA_character_
  rat    <- parsed[["rationale"]]        %||% NA_character_
  
  # Return structured data frame
  data.frame(
    module_id        = module_obj$module_id,
    title            = title,
    summary          = summ,
    confidence_label = clabel,
    rationale        = rat,
    raw_response     = raw_text,
    stringsAsFactors = FALSE
  )
}



#' Annotate All Modules Using LLM
#'
#' Iterates through all prepared modules and generates biological interpretations
#' for each using an LLM API. Results are combined into a single data frame.
#'
#' @param all_modules_ready List of prepared module objects (typically the output
#'   from \code{prepare_module()} for each module). Each element is a list with:
#'   \itemize{
#'     \item module_id (for example \code{"M1"})
#'     \item module_stats (Importance, Robustness, Annotated_fraction)
#'     \item top10_metabolites (per-metabolite data with functional summaries)
#'     \item other_metabolites (additional metabolites aggregated at class level)
#'   }
#' @param study_context List containing study metadata:
#'   - `biospecimen`: biological specimen type.
#'   - `phenotype`: phenotype of interest.
#'   - `design`: study design description (optional).
#' @param provider Character. LLM provider passed to `chat_api_call()`
#'   (e.g. `"openai"`, `"siliconflow"`, `"gemini"`).
#' @param llm_api_key Character. API key for the LLM service
#' @param model Character. Model name (e.g., "gpt-4", "gpt-3.5-turbo")
#' @param max_tokens Integer. Maximum tokens per response (default: 800)
#' @param temperature Numeric. Sampling temperature, 0 for deterministic (default: 0)
#' @param sleep_sec Numeric. Seconds to sleep between API calls to avoid rate limits (default: 0)
#'
#' @return A data frame with one row per module and columns
#'   `module_id`, `title`, `summary`, `confidence_label`, `rationale`,
#'   and `raw_response`.
#'
#' @details This function:
#'   \enumerate{
#'     \item Loops through each module in all_modules_ready
#'     \item Calls \code{annotate_one_module()} for each module
#'     \item Logs progress messages for each module
#'     \item Optionally adds delays between API calls to respect rate limits
#'     \item Combines all results into a single data frame
#'   }
#'
#' @section Rate Limiting:
#'   Use the sleep_sec parameter to add delays between API calls. For example:
#'   \itemize{
#'     \item sleep_sec = 0: No delay (default, may hit rate limits)
#'     \item sleep_sec = 1: 1 second delay (recommended for free tier APIs)
#'     \item sleep_sec = 0.5: Half-second delay (for paid APIs with higher limits)
#'   }
#'
#' @section Error Handling:
#'   If an individual module fails to annotate (e.g., API error, JSON parse error),
#'   that module's row will contain NA values but processing continues for remaining modules.
annotate_all_modules <- function(all_modules_ready,
                                 study_context,
                                 provider,
                                 llm_api_key,
                                 model,
                                 max_tokens  = 2000,
                                 temperature = 0,
                                 sleep_sec   = 0) {
  
  n_mod <- length(all_modules_ready)
  results_list <- lapply(seq_len(n_mod), function(i) {
    
    mod <- all_modules_ready[[i]]
    
    out <- annotate_one_module(
      module_obj    = mod,
      study_context = study_context,
      provider      = provider,
      llm_api_key   = llm_api_key,
      model         = model,
      max_tokens    = max_tokens,
      temperature   = temperature
    )
    
    if (sleep_sec > 0) {
      Sys.sleep(sleep_sec)
    }
    
    out
  })
  
  dplyr::bind_rows(results_list)
}
