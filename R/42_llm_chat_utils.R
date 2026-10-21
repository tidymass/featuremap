# ---- Shared helpers -------------------------------------------------------

# Seconds to wait for one LLM HTTP request. Models that "think" before
# answering (Gemini 2.5/3, reasoning models on SiliconFlow) regularly need
# more than a minute; without any limit a stuck connection would hang the
# analysis forever. Override with FEATUREMAP_LLM_TIMEOUT.
featuremap_llm_timeout <- function(default = 300) {
  value <- suppressWarnings(as.numeric(Sys.getenv("FEATUREMAP_LLM_TIMEOUT", unset = "")))
  if (length(value) != 1 || is.na(value) || value <= 0) default else value
}

# A failed call. `retryable` says whether trying again can help: network
# errors, timeouts, rate limits (429) and server errors (5xx) can; a wrong
# model ID or an invalid key (other 4xx) cannot, and retrying those only made
# users wait minutes before anything was reported.
llm_failure <- function(retryable, status = NA_integer_, message = "") {
  structure(
    list(status = status, message = message),
    class = "featuremap_llm_failure",
    retryable = isTRUE(retryable)
  )
}

is_llm_failure <- function(x) inherits(x, "featuremap_llm_failure")

llm_status_retryable <- function(status) {
  is.na(status) || status %in% c(408L, 429L) || status >= 500L
}

chat_openai_silicon_call <- function(
    messages,
    api_key,
    model,
    base_url,
    max_tokens = 2000,
    temperature = 0,
    token_parameter = c("max_tokens", "max_completion_tokens")
) {
  token_parameter <- match.arg(token_parameter)
  base_url <- sub("/+$", "", base_url)
  url <- paste0(base_url, "/v1/chat/completions")

  headers <- httr::add_headers(
    Authorization = paste("Bearer", api_key),
    "Content-Type" = "application/json"
  )

  build_body <- function(token_param, include_temperature = TRUE) {
    body <- list(
      model    = model,
      messages = messages
    )
    body[[token_param]] <- max_tokens
    if (isTRUE(include_temperature)) body$temperature <- temperature
    body
  }

  post_body <- function(body) {
    tryCatch({
      httr::POST(
        url,
        headers,
        body = jsonlite::toJSON(body, auto_unbox = TRUE),
        httr::timeout(featuremap_llm_timeout())
      )
    }, error = function(e) {
      message("HTTP error: ", e$message)
      return(NULL)
    })
  }

  network_failure <- llm_failure(TRUE, message = "network error or timeout")

  res <- post_body(build_body(token_parameter))

  if (is.null(res)) return(network_failure)

  if (httr::status_code(res) != 200) {
    response_text <- httr::content(res, "text", encoding = "UTF-8")
    if (identical(token_parameter, "max_tokens") &&
        grepl("max_tokens", response_text, fixed = TRUE) &&
        grepl("max_completion_tokens", response_text, fixed = TRUE)) {
      res <- post_body(build_body("max_completion_tokens"))
      if (!is.null(res) && httr::status_code(res) == 200) {
        token_parameter <- "max_completion_tokens"
      } else if (is.null(res)) {
        return(network_failure)
      } else {
        response_text <- httr::content(res, "text", encoding = "UTF-8")
      }
    }

    if (httr::status_code(res) != 200 &&
        grepl("temperature", response_text, fixed = TRUE) &&
        grepl("Unsupported", response_text, fixed = TRUE)) {
      retry_token_parameter <- if (identical(token_parameter, "max_tokens") &&
                                   grepl("max_completion_tokens", response_text, fixed = TRUE)) {
        "max_completion_tokens"
      } else {
        token_parameter
      }
      res <- post_body(build_body(retry_token_parameter, include_temperature = FALSE))
      if (is.null(res)) return(network_failure)
      response_text <- httr::content(res, "text", encoding = "UTF-8")
    }

    if (httr::status_code(res) != 200) {
      status <- httr::status_code(res)
      message("Chat HTTP status: ", status)
      message("Response: ", response_text)
      return(llm_failure(llm_status_retryable(status), status, response_text))
    }
  }

  cont <- tryCatch(httr::content(res, encoding = "UTF-8"), error = function(e) NULL)
  if (is.null(cont)) return(llm_failure(TRUE, 200L, "unreadable response body"))
  content_txt <- tryCatch(cont$choices[[1]]$message$content, error = function(e) NULL)
  finish_reason <- tryCatch(cont$choices[[1]]$finish_reason, error = function(e) NA_character_)
  if (is.null(content_txt)) {
    message("Failed to extract assistant content.")
    return(llm_failure(TRUE, 200L, "no assistant content in response"))
  }

  list(content = content_txt, raw = cont, finish_reason = finish_reason)
}

# SiliconFlow runs two independent sites, api.siliconflow.com (global) and
# api.siliconflow.cn (mainland China), with separate accounts: a key issued by
# one is rejected by the other. When the region is not configured we probe both
# and remember which one worked.
#
# That memory must be per API key, not per R process. One R process serves
# every Shiny session, so a process-wide cache (the previous
# options(featuremap.siliconflow.active_base_url = ...)) let whoever ran first
# pin the endpoint for everybody else, and users holding a key for the other
# site could never get a successful call.
.siliconflow_endpoint_cache <- new.env(parent = emptyenv())

# Identify a key without keeping the key itself around. digest is available in
# practice (htmltools depends on it), but the fallback keeps this working if it
# is not: the fingerprint only has to separate keys inside one process.
siliconflow_key_fingerprint <- function(api_key) {
  key <- if (is.null(api_key) || length(api_key) == 0) {
    ""
  } else {
    trimws(as.character(api_key[1]))
  }
  if (!nzchar(key)) return("__nokey__")

  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(key, algo = "sha256"))
  }

  codes <- utf8ToInt(key)
  checksum <- sum(as.numeric(codes) * seq_along(codes)) %% 2147483647
  paste0("len", nchar(key), "_", format(checksum, scientific = FALSE))
}

siliconflow_cached_base_url <- function(api_key) {
  fp <- siliconflow_key_fingerprint(api_key)
  if (!exists(fp, envir = .siliconflow_endpoint_cache, inherits = FALSE)) {
    return("")
  }
  value <- get(fp, envir = .siliconflow_endpoint_cache, inherits = FALSE)
  if (is.null(value) || length(value) == 0) return("")
  trimws(as.character(value[1]))
}

siliconflow_cache_base_url <- function(api_key, base_url) {
  fp <- siliconflow_key_fingerprint(api_key)
  assign(fp, base_url, envir = .siliconflow_endpoint_cache)
  invisible(base_url)
}

siliconflow_base_urls <- function(api_key = NULL) {
  # 1. An explicitly configured endpoint always wins.
  explicit <- getOption(
    "featuremap.siliconflow.base_url",
    Sys.getenv("SILICONFLOW_BASE_URL", unset = "")
  )
  if (is.null(explicit) || length(explicit) == 0) explicit <- ""
  explicit <- trimws(as.character(explicit[1]))
  if (nzchar(explicit)) return(sub("/+$", "", explicit))

  # 2. An explicitly configured region wins over anything auto-detected.
  region <- getOption(
    "featuremap.siliconflow.region",
    Sys.getenv("SILICONFLOW_REGION", unset = "auto")
  )
  if (is.null(region) || length(region) == 0) region <- "auto"
  region <- tolower(trimws(as.character(region[1])))

  if (region %in% c("cn", "china", "domestic", "mainland")) {
    return("https://api.siliconflow.cn")
  }
  if (region %in% c("com", "global", "international", "intl", "overseas", "foreign")) {
    return("https://api.siliconflow.com")
  }

  defaults <- c("https://api.siliconflow.com", "https://api.siliconflow.cn")

  # 3. Auto: start with whatever last worked for this key, but keep the other
  #    site as a fallback so a stale entry can never lock a user out.
  cached <- siliconflow_cached_base_url(api_key)
  if (nzchar(cached)) return(unique(c(sub("/+$", "", cached), defaults)))

  defaults
}

chat_siliconflow_call <- function(
    messages,
    api_key,
    model,
    max_tokens = 2000,
    temperature = 0
) {
  base_urls <- siliconflow_base_urls(api_key)
  last_failure <- llm_failure(TRUE, message = "no SiliconFlow endpoint configured")
  any_retryable <- FALSE
  for (i in seq_along(base_urls)) {
    base_url <- base_urls[[i]]
    if (i > 1) {
      message("SiliconFlow endpoint ", base_urls[[i - 1]],
              " failed; trying ", base_url, ".")
    }
    out <- chat_openai_silicon_call(
      messages = messages,
      api_key = api_key,
      model = model,
      base_url = base_url,
      max_tokens = max_tokens,
      temperature = temperature
    )
    if (!is_llm_failure(out)) {
      out$base_url <- base_url
      siliconflow_cache_base_url(api_key, base_url)
      return(out)
    }
    # A 4xx from one site (e.g. a key or model that only exists on the other
    # site) is still worth trying on the other one; only the overall result
    # decides whether the call is retried.
    last_failure <- out
    any_retryable <- any_retryable || isTRUE(attr(out, "retryable"))
  }
  attr(last_failure, "retryable") <- any_retryable
  last_failure
}

# "models/gemini-2.5-pro" and "gemini-2.5-pro" name the same model; the URL
# needs the bare ID.
gemini_model_id <- function(model) {
  sub("^models/", "", trimws(as.character(model %||% "")))
}

# Google strongly recommends leaving temperature at its default (1.0) for
# Gemini 3: lower values can make the model loop or degrade.
gemini_is_v3 <- function(model) {
  grepl("^gemini-3", gemini_model_id(model))
}

# Largest maxOutputTokens the model accepts. Gemini 1.x/2.0 stop at 8192;
# 2.5 and later allow 65536.
gemini_output_cap <- function(model) {
  if (grepl("^gemini-(1|2\\.0)", gemini_model_id(model))) 8192L else 65536L
}

# Thinking models (Gemini 2.5, 3) count their thinking tokens against
# maxOutputTokens, so a limit sized for the JSON answer alone can be spent
# entirely on thinking. Start with room for both.
gemini_initial_max_tokens <- function(max_tokens, model) {
  as.integer(min(max(8192L, max_tokens), gemini_output_cap(model)))
}

# Limit to retry with after a MAX_TOKENS stop; NULL when already at the cap.
gemini_next_max_tokens <- function(cur_max, model) {
  cap <- gemini_output_cap(model)
  if (cur_max >= cap) return(NULL)
  as.integer(min(max(3000L, cur_max * 2L), cap))
}

gemini_hit_max_tokens <- function(res_call) {
  !is.null(res_call) && !is_llm_failure(res_call) &&
    !is.null(res_call$finish_reason) &&
    identical(toupper(as.character(res_call$finish_reason)), "MAX_TOKENS")
}

# Gemini's schema fields `required` and `enum` are lists of strings, but
# jsonlite::toJSON(auto_unbox = TRUE) turns a one-element vector into a bare
# string ("required": "x"), which the API rejects. Mark them as arrays.
gemini_schema_arrays <- function(x) {
  if (!is.list(x)) return(x)
  nms <- names(x)
  for (i in seq_along(x)) {
    nm <- if (is.null(nms)) "" else nms[[i]]
    if (nm %in% c("required", "enum") && is.atomic(x[[i]])) {
      x[[i]] <- I(as.character(x[[i]]))
    } else if (is.list(x[[i]])) {
      x[[i]] <- gemini_schema_arrays(x[[i]])
    }
  }
  x
}

chat_gemini_call <- function(
    messages,
    api_key,
    model = "gemini-2.5-flash",
    max_tokens = 2000,
    temperature = 0,
    response_schema = NULL
) {
  role_vec <- vapply(messages, function(m) {
    if (is.null(m$role)) return("user")
    as.character(m$role)
  }, character(1))
  content_vec <- vapply(messages, function(m) {
    if (is.null(m$content)) return("")
    as.character(m$content)
  }, character(1))

  sys_prompt <- paste(content_vec[role_vec == "system"], collapse = "\n\n")
  user_prompt <- paste(content_vec[role_vec != "system"], collapse = "\n\n")

  url <- sprintf(
    "https://generativelanguage.googleapis.com/v1beta/models/%s:generateContent?key=%s",
    gemini_model_id(model), api_key
  )

  generation_config <- list(
    maxOutputTokens = as.integer(min(max_tokens, gemini_output_cap(model))),
    responseMimeType = "application/json"
  )
  if (!gemini_is_v3(model)) {
    generation_config$temperature <- temperature
  }

  body <- list(
    contents = list(list(
      role = "user",
      parts = list(list(text = user_prompt))
    )),
    generationConfig = generation_config
  )

  if (!is.null(response_schema)) {
    body$generationConfig$responseSchema <- gemini_schema_arrays(response_schema)
  }

  if (nzchar(sys_prompt)) {
    body$systemInstruction <- list(
      parts = list(list(text = sys_prompt))
    )
  }

  res <- tryCatch({
    httr::POST(
      url,
      httr::add_headers("Content-Type" = "application/json"),
      body = jsonlite::toJSON(body, auto_unbox = TRUE),
      httr::timeout(featuremap_llm_timeout())
    )
  }, error = function(e) {
    message("Gemini HTTP error: ", e$message)
    return(NULL)
  })

  if (is.null(res)) return(llm_failure(TRUE, message = "network error or timeout"))

  if (httr::status_code(res) != 200) {
    status <- httr::status_code(res)
    response_text <- httr::content(res, "text", encoding = "UTF-8")
    message("Gemini HTTP status: ", status)
    message("Response: ", response_text)
    return(llm_failure(llm_status_retryable(status), status, response_text))
  }

  cont <- tryCatch(httr::content(res, encoding = "UTF-8"), error = function(e) NULL)
  if (is.null(cont)) return(llm_failure(TRUE, 200L, "unreadable response body"))

  # Keep only the answer text: skip thought-summary parts and parts that carry
  # no text (e.g. only a thought signature). When thinking used up the whole
  # output budget, `parts` is empty and the result is "" with finish reason
  # MAX_TOKENS, which the callers answer with a larger limit.
  parts <- tryCatch(cont$candidates[[1]]$content$parts, error = function(e) NULL)
  txt <- paste(
    vapply(parts %||% list(), function(p) {
      if (isTRUE(p$thought) || is.null(p$text)) "" else as.character(p$text)[1]
    }, character(1)),
    collapse = ""
  )
  finish_reason <- tryCatch(cont$candidates[[1]]$finishReason, error = function(e) NA_character_)

  list(content = txt, raw = cont, finish_reason = finish_reason)
}

chat_api_call <- function(
    messages,
    provider = c("siliconflow", "openai", "gemini", "local"),
    api_key,
    model,
    base_url = NULL,
    max_tokens = 2000,
    temperature = 0,
    response_schema = NULL
) {
  provider <- match.arg(provider)

  if (provider == "siliconflow") {
    return(chat_siliconflow_call(
      messages = messages,
      api_key = api_key,
      model = model,
      max_tokens = max_tokens,
      temperature = temperature
    ))
  }

  if (provider == "openai") {
    return(chat_openai_silicon_call(
      messages = messages,
      api_key = api_key,
      model = model,
      base_url = "https://api.openai.com",
      max_tokens = max_tokens,
      temperature = temperature,
      token_parameter = "max_completion_tokens"
    ))
  }

  if (provider == "gemini") {
    return(chat_gemini_call(
      messages = messages,
      api_key = api_key,
      model = model,
      max_tokens = max_tokens,
      temperature = temperature,
      response_schema = response_schema
    ))
  }

  if (provider == "local") {
    if (is.null(base_url) || !nzchar(trimws(base_url))) {
      stop("provider = \"local\" requires a non-empty base_url.")
    }
    return(chat_openai_silicon_call(
      messages = messages,
      api_key = if (is.null(api_key) || !nzchar(trimws(api_key))) "EMPTY" else api_key,
      model = model,
      base_url = base_url,
      max_tokens = max_tokens,
      temperature = temperature
    ))
  }

  stop("Unsupported provider: ", provider)
}

safe_chat_call <- function(..., base_url = NULL) {
  attempts <- 5
  for (i in seq_len(attempts)) {
    out <- chat_api_call(..., base_url = base_url)
    if (!is.null(out) && !is_llm_failure(out)) return(out)

    if (is_llm_failure(out) && !isTRUE(attr(out, "retryable"))) {
      message(
        "Chat call failed (HTTP ", out$status, "); not retrying, since the ",
        "same request would fail again. Check the model ID and API key."
      )
      return(NULL)
    }
    if (i < attempts) {
      message("Chat call failed. Waiting ", i * 3, " seconds...")
      Sys.sleep(i * 3)
    }
  }
  NULL
}

extract_json_object <- function(txt) {
  if (is.null(txt) || !nzchar(txt)) return(NULL)

  m <- regmatches(txt, regexpr("```json\\s*[\\s\\S]*?\\s*```", txt, perl = TRUE))
  if (length(m) && nzchar(m)) {
    inner <- sub("^```json\\s*", "", m)
    inner <- sub("\\s*```$", "", inner)
    inner <- trimws(inner)
    if (nzchar(inner)) return(inner)
  }

  start <- regexpr("\\{", txt, perl = TRUE)
  if (start[1] == -1) return(NULL)
  s <- substring(txt, start[1])

  depth <- 0
  in_str <- FALSE
  esc <- FALSE
  for (i in seq_len(nchar(s))) {
    ch <- substr(s, i, i)
    if (esc) {
      esc <- FALSE
      next
    }
    if (ch == "\\") {
      if (in_str) esc <- TRUE
      next
    }
    if (ch == "\"") {
      in_str <- !in_str
      next
    }
    if (!in_str && ch == "{") depth <- depth + 1
    if (!in_str && ch == "}") {
      depth <- depth - 1
      if (depth == 0) return(substr(s, 1, i))
    }
  }
  NULL
}
