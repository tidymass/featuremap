


chat_openai_silicon_call <- function(
    messages,
    api_key,
    model,
    base_url,
    max_tokens = 2000,
    temperature = 0
) {
  url <- paste0(base_url, "/v1/chat/completions")
  
  headers <- httr::add_headers(
    Authorization = paste("Bearer", api_key),
    "Content-Type" = "application/json"
  )
  
  body <- list(
    model       = model,
    messages    = messages,
    max_tokens  = max_tokens,
    temperature = temperature
  )
  
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
  if (is.null(content_txt)) {
    message("Failed to extract assistant content.")
    return(NULL)
  }
  
  list(content = content_txt, raw = cont)
}

chat_gemini_call <- function(
    messages,
    api_key,
    model = "gemini-2.0-pro",
    max_tokens = 2000,
    temperature = 0
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
      temperature     = temperature,
      maxOutputTokens = max_tokens,
      responseMimeType = "application/json"
    )
  )

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
  
  # chat_gemini_call() -- FIXED VERSION
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
  
  
  if (is.null(txt)) {
    message("Failed to extract Gemini text.")
    return(NULL)
  }
  
  list(content = txt, raw = cont)
}


chat_api_call <- function(
    messages,
    provider = c("siliconflow", "openai", "gemini"),
    api_key,
    model,
    max_tokens = 2000,
    temperature = 0
) {
  
  provider <- match.arg(provider)
  
  if (provider == "siliconflow") {
    return(chat_openai_silicon_call(
      messages = messages,
      api_key = api_key,
      model = model,
      base_url = "https://api.siliconflow.cn",
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
      temperature = temperature
    ))
  }
  
  if (provider == "gemini") {
    return(chat_gemini_call(
      messages = messages,
      api_key = api_key,
      model = model,
      max_tokens = max_tokens,
      temperature = temperature
    ))
  }
  
  stop("Unsupported provider: ", provider)
}


safe_chat_call <- function(...) {
  for (i in 1:5) {
    out <- chat_api_call(...)
    if (!is.null(out)) return(out)
    
    message("Chat call failed. Waiting ", i * 3, " seconds...")
    Sys.sleep(i * 3)
  }
  NULL
}




#' Refine top-k PubMed articles for a single metabolite using SiliconFlow LLM
#'
#' @param res_single List. Output of embed_pubmed_for_metabolite().
#' @param provider Character. Chat provider (e.g. `"openai"`, `"siliconflow"`, `"gemini"`).
#' @param llm_api_key Character. SiliconFlow API key (for chat).
#' @param llm_model Character. Chat model name on SiliconFlow.
#' @param top_paper Integer. Number of final articles to keep (default 3).
#' @param max_chars_abstract Integer. Truncate abstracts in prompt (default 1500).
#' @param sleep_seconds Numeric. Delay in seconds before each LLM call.
#'
#' @return Same structure as res_single, but only top_paper articles kept.
#'         On failure, falls back to similarity top_paper.
refine_pubmed_for_metabolite_llm <- function(
    res_single,
    provider,
    llm_api_key,
    llm_model,
    top_paper = 3,
    max_chars_abstract = 1500,
    sleep_seconds = 1
) {

  if (is.null(res_single)) return(NULL)
  
  k <- length(res_single$pmids)
  if (k == 0) {
    return(res_single)
  }
  
  df <- data.frame(
    idx      = seq_len(k),
    pmid     = res_single$pmids,
    title    = res_single$titles,
    abstract = substr(res_single$abstracts, 1, max_chars_abstract),
    score    = res_single$scores,
    stringsAsFactors = FALSE
  )
  
  articles_txt <- paste(
    apply(df, 1, function(r) {
      paste0(
        "Article ", r[["idx"]], ":\n",
        "PMID: ", r[["pmid"]], "\n",
        "Title: ", r[["title"]], "\n",
        "Similarity_score: ", sprintf("%.4f", as.numeric(r[["score"]])), "\n",
        "Abstract: ", r[["abstract"]], "\n"
      )
    }),
    collapse = "\n\n"
  )
  
  user_prompt <- paste0(
    "You are an expert in metabolomics and biochemistry.\n\n",
    "We are studying the metabolite: ", res_single$met_name, ".\n",
    "We care especially about these aspects:\n",
    "1) Biochemical function of the metabolite;\n",
    "2) Endogenous metabolism and physiological role;\n",
    "3) Metabolic reactions involving this metabolite.\n\n",
    "I will give you ", k, " candidate PubMed articles that were pre-selected by embedding similarity.\n",
    "Please choose up to ", top_paper, " articles that are overall most relevant to these aspects.\n",
    "If none or only a few articles are relevant, return fewer (or even zero). \n",
    "Do NOT force selection if relevance is weak.\n",
    "Return ONLY a valid JSON object with the following schema (no extra text):\n",
    "{\n",
    '  \"selected_indices\": [i, j],\n',
    '  \"reason\": \"brief explanation\"\n',
    "}\n",
    "where indices are integers between 1 and ", k, ".\n\n",
    "Here are the candidate articles:\n\n",
    articles_txt
  )
  
  messages <- list(
    list(role = "system",
         content = "You rank scientific articles for metabolomics annotation. Always answer with valid JSON only."),
    list(role = "user", content = user_prompt)
  )
  
  
  parsed <- NULL
  llm_raw <- tryCatch({
    Sys.sleep(sleep_seconds)
    safe_chat_call(
      messages    = messages,
      provider    = provider,
      api_key     = llm_api_key,
      model       = llm_model,
      temperature = 0,
      max_tokens  = 2000
    )
    
  }, error = function(e) {
    message("LLM call failed: ", e$message)
    NULL
  })
  
  
  final_indices <- integer(0)
  if (!is.null(llm_raw) && !is.null(llm_raw$content)) {
    
    json_txt <- extract_json_object(llm_raw$content)
    parsed <- tryCatch({
      if (is.null(json_txt)) stop("No JSON object found")
      jsonlite::fromJSON(json_txt)
    }, error = function(e) NULL)
    
    
    if (!is.null(parsed) && !is.null(parsed$selected_indices)) {
      
      sel_idx <- as.integer(parsed$selected_indices)
      sel_idx <- sel_idx[sel_idx %in% seq_len(k)]
      sel_idx <- unique(sel_idx)
      
      if (length(sel_idx) > top_paper) {
        sel_idx <- sel_idx[1:top_paper]
      }
      
      final_indices <- sel_idx
    }
    
  }
  
  fallback_needed <- is.null(llm_raw) ||
    is.null(llm_raw$content) ||
    is.null(parsed) ||
    length(final_indices) == 0
  
  if (fallback_needed) {
    final_indices <- order(res_single$scores, decreasing = TRUE)[
      seq_len(min(top_paper, length(res_single$scores)))
    ]
  }
  

  res_single$pmids     <- res_single$pmids[final_indices]
  res_single$titles    <- res_single$titles[final_indices]
  res_single$abstracts <- res_single$abstracts[final_indices]
  res_single$scores    <- res_single$scores[final_indices]
  res_single$contexts    <- res_single$contexts[final_indices]
  
  res_single
}






#' Refine PubMed ranking for all modules & metabolites using SiliconFlow LLM
#'
#' @param module_pubmed_embed List. Output of your current embed pipeline:
#'   outer: modules; inner: metabolites; inner element is embed_pubmed_for_metabolite() result.
#' @param provider Character. Chat provider (e.g. `"openai"`, `"siliconflow"`, `"gemini"`).
#' @param llm_api_key Character. SiliconFlow API key for chat.
#' @param llm_model Character. Chat model name on SiliconFlow.
#' @param top_paper Integer. Number of articles to keep per metabolite (default 3).
#' @param max_chars_abstract Integer. Truncate abstracts for prompts.
#' @param sleep_seconds Numeric. Delay in seconds before each LLM call.
#'
#' @return Same structure as module_pubmed_embed, but each metabolite only keeps top_paper articles.
refine_module_pubmed_with_llm <- function(
    module_pubmed_embed,
    provider,
    llm_api_key,
    llm_model,
    top_paper = 3,
    max_chars_abstract = 1500,
    sleep_seconds = 1
) {
  
  out <- lapply(seq_along(module_pubmed_embed), function(mi) {
    
    message("--Reranking for ", names(module_pubmed_embed)[mi], " ...")
    mod_res <- module_pubmed_embed[[mi]]
    
    lapply(mod_res, function(res_single) {
      refined <- refine_pubmed_for_metabolite_llm(
        res_single        = res_single,
        provider          = provider,
        llm_api_key       = llm_api_key,
        llm_model         = llm_model,
        top_paper         = top_paper,
        max_chars_abstract = max_chars_abstract
      )
      refined
    })
  })

  names(out) <- names(module_pubmed_embed)
  
  return(out)
}




#' Summarize metabolite function from top PubMed articles using SiliconFlow LLM
#'
#' @param res_single List. Output of embed_pubmed_for_metabolite() *after* refinement:
#'   list(met_name, pmids, titles, abstracts, scores, ...)
#' @param provider Character. Chat provider (e.g. `"openai"`, `"siliconflow"`, `"gemini"`).
#' @param llm_api_key Character. SiliconFlow API key (for chat).
#' @param llm_model Character. Chat model name on SiliconFlow.
#' @param max_chars_abstract Integer. Truncate each abstract in the prompt (default 1200).
#'
#' @return Same list as res_single, with an extra element:
#'   $functional_summary (character string).
summarize_metabolite_pubmed_llm <- function(
    res_single,
    provider,
    llm_api_key,
    llm_model,
    max_chars_abstract = 1500
) {
  
  if (is.null(res_single)) return(NULL)
  
  if (length(res_single$pmids) == 0) {
    res_single$functional_summary <- NA_character_
    return(res_single)
  }
  
  k <- length(res_single$pmids)
  
  df <- data.frame(
    idx      = seq_len(k),
    pmid     = res_single$pmids,
    title    = res_single$titles,
    abstract = substr(res_single$abstracts, 1, max_chars_abstract),
    stringsAsFactors = FALSE
  )
  
  articles_txt <- paste(
    apply(df, 1, function(r) {
      paste0(
        "Article ", r[["idx"]], ":\n",
        "PMID: ", r[["pmid"]], "\n",
        "Title: ", r[["title"]], "\n",
        "Abstract (truncated): ", r[["abstract"]], "\n"
      )
    }),
    collapse = "\n\n"
  )
  
  user_prompt <- paste0(
    "You are an expert in metabolomics and biochemistry.\n\n",
    "The metabolite of interest is: ", res_single$met_name, ".\n",
    "You are given several PubMed abstracts that focus on this metabolite.\n",
    "Please synthesize the *overall* knowledge about this metabolite from these studies.\n\n",
    "Your tasks:\n",
    "1) Describe the main biochemical functions of this metabolite.\n",
    "2) Describe key metabolic pathways or reactions it participates in.\n",
    "3) Describe its physiological or pathophysiological roles (including disease relevance) if mentioned.\n",
    "4) Do NOT list individual studies separately; instead, merge them into a single coherent description.\n",
    "5) Avoid mentioning authors, years, sample sizes, or specific experimental details.\n\n",
    "Output ONLY a valid JSON object with the following format:\n",
    "{\n",
    '  \"metabolite_name\": \"', res_single$met_name, '\",\n',
    '  \"functional_summary\": \"<3-6 sentences summarizing functions, pathways, roles>\"\n',
    "}\n\n",
    "Here are the abstracts:\n\n",
    articles_txt
  )
  
  messages <- list(
    list(role = "system",
         content = "You summarize biomedical literature into concise functional descriptions. Always answer with valid JSON only."),
    list(role = "user", content = user_prompt)
  )
  
  llm_raw <- tryCatch({
    safe_chat_call(
      messages    = messages,
      provider    = provider,
      api_key     = llm_api_key,
      model       = llm_model,
      temperature = 0,
      max_tokens  = 2000
    )
  }, error = function(e) {
    message("LLM call failed in summarize_metabolite_pubmed_llm: ", e$message)
    NULL
  })
  
  res_single$functional_summary <- NA_character_
  
  if (!is.null(llm_raw) && !is.null(llm_raw$content)) {
    
    raw_txt <- llm_raw$content
    json_txt <- extract_json_object(raw_txt)
    
    parsed <- tryCatch({
      if (is.null(json_txt)) {
        stop("No JSON object could be extracted from LLM output")
      }
      jsonlite::fromJSON(json_txt)
    }, error = function(e) {
      message("Failed to parse JSON in summarize_metabolite_pubmed_llm: ", e$message)
      NULL
    })
    
    if (!is.null(parsed) && !is.null(parsed$functional_summary)) {
      res_single$functional_summary <- parsed$functional_summary
    } else {
      message("LLM JSON missing 'functional_summary'; keeping NA.")
    }
    
  } else {
    message("Empty LLM response in summarize_metabolite_pubmed_llm; keeping NA.")
  }
  
  res_single
}


#' Add functional summaries (from PubMed) to all metabolites in all modules
#'
#' @param module_pubmed_embed_refined List.
#'   Output of refine_module_pubmed_with_llm():
#'   outer: modules; inner: metabolites; inner element: res_single.
#' @param provider Character. Chat provider (e.g. `"openai"`, `"siliconflow"`, `"gemini"`).
#' @param llm_api_key Character.API key for chat.
#' @param llm_model Character. Chat model name.
#' @param max_chars_abstract Integer. Truncate abstracts for prompts.
#'
#' @return Same structure as module_pubmed_embed_refined, but each res_single
#'   has an extra element: $functional_summary.
add_metabolite_summaries_to_modules <- function(
    module_pubmed_embed_refined,
    provider,
    llm_api_key,
    llm_model, 
    max_chars_abstract = 1500
) {
  out <- lapply(seq_along(module_pubmed_embed_refined), function(mi) {
    
    message("--Summarizing for ", names(module_pubmed_embed_refined)[mi], " ...")
    
    mod_res <- module_pubmed_embed_refined[[mi]]
    
    lapply(mod_res, function(res_single) {
      summarize_metabolite_pubmed_llm(
        res_single        = res_single,
        provider          = provider,
        llm_api_key       = llm_api_key,
        llm_model         = llm_model,
        max_chars_abstract = max_chars_abstract
      )
    })
  })

  names(out) <- names(module_pubmed_embed_refined)
  out
}


extract_json_object <- function(txt) {
  if (is.null(txt) || !nzchar(txt)) return(NULL)
  
  # 1) Prefer fenced ```json ... ```
  m <- regmatches(txt, regexpr("```json\\s*[\\s\\S]*?\\s*```", txt, perl=TRUE))
  if (length(m) && nzchar(m)) {
    inner <- sub("^```json\\s*", "", m)
    inner <- sub("\\s*```$", "", inner)
    inner <- trimws(inner)
    if (nzchar(inner)) return(inner)
  }
  
  # 2) Bracket counting for first complete {...}
  start <- regexpr("\\{", txt, perl=TRUE)
  if (start[1] == -1) return(NULL)
  s <- substring(txt, start[1])
  
  depth <- 0
  in_str <- FALSE
  esc <- FALSE
  for (i in seq_len(nchar(s))) {
    ch <- substr(s, i, i)
    if (esc) { esc <- FALSE; next }
    if (ch == "\\") { if (in_str) esc <- TRUE; next }
    if (ch == "\"") { in_str <- !in_str; next }
    if (!in_str && ch == "{") depth <- depth + 1
    if (!in_str && ch == "}") {
      depth <- depth - 1
      if (depth == 0) return(substr(s, 1, i))
    }
  }
  NULL
}
