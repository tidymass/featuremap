chat_token_limit_body <- function(body, max_tokens, token_param = "max_tokens") {
  token_param <- match.arg(token_param, c("max_tokens", "max_completion_tokens"))
  body[[token_param]] <- max_tokens
  body
}

chat_openai_silicon_call <- function(
    messages,
    api_key,
    model,
    base_url,
    max_tokens = 2000,
    temperature = 0,
    token_param = c("max_tokens", "max_completion_tokens")
) {
  token_param <- match.arg(token_param)
  base_url <- sub("/+$", "", base_url)
  url <- paste0(base_url, "/v1/chat/completions")

  headers <- httr::add_headers(
    Authorization = paste("Bearer", api_key),
    "Content-Type" = "application/json"
  )

  body <- list(
    model       = model,
    messages    = messages,
    temperature = temperature
  )
  body <- chat_token_limit_body(body, max_tokens, token_param = token_param)

  res <- tryCatch({
    httr::POST(url, headers, body = jsonlite::toJSON(body, auto_unbox = TRUE))
  }, error = function(e) {
    message("HTTP error: ", e$message)
    return(NULL)
  })

  if (is.null(res)) return(NULL)

  if (httr::status_code(res) != 200) {
    message("Chat HTTP status: ", httr::status_code(res))
    message("Response: ", httr::content(res, "text", encoding = "UTF-8"))
    return(NULL)
  }

  cont <- tryCatch(httr::content(res, encoding = "UTF-8"), error = function(e) NULL)
  if (is.null(cont)) return(NULL)
  content_txt <- tryCatch(cont$choices[[1]]$message$content, error = function(e) NULL)
  finish_reason <- tryCatch(cont$choices[[1]]$finish_reason, error = function(e) NA_character_)
  if (is.null(content_txt)) {
    message("Failed to extract assistant content.")
    return(NULL)
  }

  list(content = content_txt, raw = cont, finish_reason = finish_reason)
}

siliconflow_base_urls <- function() {
  explicit <- getOption(
    "featuremap.siliconflow.base_url",
    Sys.getenv("SILICONFLOW_BASE_URL", unset = "")
  )
  if (is.null(explicit) || length(explicit) == 0) explicit <- ""
  explicit <- trimws(as.character(explicit[1]))
  if (nzchar(explicit)) return(sub("/+$", "", explicit))

  cached <- getOption("featuremap.siliconflow.active_base_url", "")
  if (!is.null(cached) && length(cached) > 0) {
    cached <- trimws(as.character(cached[1]))
    if (nzchar(cached)) return(sub("/+$", "", cached))
  }

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

  c("https://api.siliconflow.com", "https://api.siliconflow.cn")
}

chat_siliconflow_call <- function(
    messages,
    api_key,
    model,
    max_tokens = 2000,
    temperature = 0
) {
  base_urls <- siliconflow_base_urls()
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
      temperature = temperature,
      token_param = "max_tokens"
    )
    if (!is.null(out)) {
      out$base_url <- base_url
      options(featuremap.siliconflow.active_base_url = base_url)
      return(out)
    }
  }
  NULL
}

chat_gemini_call <- function(
    messages,
    api_key,
    model = "gemini-2.0-pro",
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
    model, api_key
  )

  body <- list(
    contents = list(list(
      role = "user",
      parts = list(list(text = user_prompt))
    )),
    generationConfig = list(
      temperature = temperature,
      maxOutputTokens = max_tokens,
      responseMimeType = "application/json"
    )
  )

  if (!is.null(response_schema)) {
    body$generationConfig$responseSchema <- response_schema
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
      httr::timeout(60)
    )
  }, error = function(e) {
    message("Gemini HTTP error: ", e$message)
    return(NULL)
  })

  if (is.null(res)) return(NULL)

  if (httr::status_code(res) != 200) {
    message("Gemini HTTP status: ", httr::status_code(res))
    message("Response: ", httr::content(res, "text", encoding = "UTF-8"))
    return(NULL)
  }

  cont <- tryCatch(httr::content(res, encoding = "UTF-8"), error = function(e) NULL)
  if (is.null(cont)) return(NULL)

  txt <- tryCatch({
    paste(
      vapply(
        cont$candidates[[1]]$content$parts,
        function(p) p$text,
        character(1)
      ),
      collapse = ""
    )
  }, error = function(e) NULL)
  finish_reason <- tryCatch(cont$candidates[[1]]$finishReason, error = function(e) NA_character_)

  if (is.null(txt)) {
    message("Failed to extract Gemini text.")
    return(NULL)
  }

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
      token_param = "max_completion_tokens"
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
      temperature = temperature,
      token_param = "max_tokens"
    ))
  }

  stop("Unsupported provider: ", provider)
}

safe_chat_call <- function(..., base_url = NULL) {
  for (i in 1:5) {
    out <- chat_api_call(..., base_url = base_url)
    if (!is.null(out)) return(out)

    message("Chat call failed. Waiting ", i * 3, " seconds...")
    Sys.sleep(i * 3)
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
