# ============================================================
# 1. 全局限速状态
# ============================================================

.emb_call_state <- new.env(parent = emptyenv())
.emb_call_state$last_time <- NULL


# ============================================================
# 2. 工具函数
# ============================================================

#' Calculate Cosine Similarity Between Two Vectors
#'
#' Computes the cosine similarity score between two numeric vectors.
#' Used internally to measure semantic closeness between query and
#' abstract embeddings.
#'
#' @param a Numeric vector. First embedding vector (e.g. query embedding).
#' @param b Numeric vector. Second embedding vector (e.g. abstract embedding).
#'   Must have the same length as \code{a}.
#'
#' @return Numeric scalar in \code{[-1, 1]}. Returns \code{0} when either
#'   vector has zero norm or when \code{length(a) != length(b)}.
cosine_sim <- function(a, b) {
  if (length(a) != length(b)) {
    warning("Vector dimensions do not match")
    return(0)
  }
  denom <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (denom == 0) return(0)
  sum(a * b) / denom
}


#' Build Metabolite Semantic Query Texts
#'
#' Generates a fixed set of three natural-language query sentences for a
#' given metabolite name. These sentences are embedded and used as the
#' query side of the cosine similarity ranking against PubMed abstracts.
#'
#' @param met_name Character string. Primary name of the metabolite
#'   (e.g. \code{"lactate"}).
#'
#' @return Named list of three character strings:
#'   \itemize{
#'     \item \code{biochem}  — biochemical mechanism query.
#'     \item \code{physio}   — endogenous metabolism and physiological role query.
#'     \item \code{reaction} — metabolic reaction query.
#'   }
build_metabolite_queries <- function(met_name) {
  list(
    biochem  = sprintf("%s is involved in biochemical mechanisms", met_name),
    physio   = sprintf("Endogenous metabolism and physiological role of %s", met_name),
    reaction = sprintf("Metabolic reactions involving %s", met_name)
  )
}


# ============================================================
# 3. Embedding 维度检测
# ============================================================

#' Look Up or Auto-Detect Embedding Dimension for a Model
#'
#' Returns the embedding output dimension for a given model name. For
#' well-known models the dimension is returned from a hard-coded lookup
#' table; for unknown models an API call is made to detect it
#' automatically.
#'
#' @param embed_model Character string. Model identifier
#'   (e.g. \code{"text-embedding-3-large"}, \code{"BAAI/bge-large-en-v1.5"}).
#' @param llm_api_key Character string or \code{NULL}. API key used only
#'   when the model is not in the built-in lookup table and auto-detection
#'   via \code{detect_embedding_dim} is required. Pass \code{NULL} to skip
#'   auto-detection (a warning is then issued).
#' @param provider Character string. Embedding API provider used during
#'   auto-detection. One of \code{"openai"}, \code{"siliconflow"}, or
#'   \code{"gemini"} (default: \code{"siliconflow"}).
#'
#' @return Integer scalar giving the embedding dimension, or \code{NULL}
#'   if the dimension cannot be determined (unknown model and no API key,
#'   or the test call failed).
#'
#' @details
#' Built-in lookup covers: \code{BAAI/bge-large-en-v1.5} (1024),
#' \code{BAAI/bge-base-en-v1.5} (768), \code{BAAI/bge-small-en-v1.5} (384),
#' \code{text-embedding-ada-002} (1536), \code{text-embedding-3-small} (1536),
#' \code{text-embedding-3-large} (3072), \code{embedding-001} (768).
get_embedding_dim <- function(embed_model,
                              llm_api_key = NULL,
                              provider    = "siliconflow") {
  
  model_dims <- list(
    "BAAI/bge-large-en-v1.5" = 1024,
    "BAAI/bge-base-en-v1.5"  = 768,
    "BAAI/bge-small-en-v1.5" = 384,
    "text-embedding-ada-002"  = 1536,
    "text-embedding-3-small"  = 1536,
    "text-embedding-3-large"  = 3072,
    "embedding-001"           = 768
  )
  
  if (embed_model %in% names(model_dims)) {
    return(model_dims[[embed_model]])
  }
  
  if (!embed_model %in% names(model_dims) && !is.null(llm_api_key)) {
    dim_detected <- detect_embedding_dim(
      llm_api_key = llm_api_key,
      embed_model = embed_model,
      provider    = provider
    )
    if (!is.null(dim_detected)) return(dim_detected)
  }
  
  warning("Could not determine embedding dimension for model: ", embed_model)
  warning("Provide llm_api_key to auto-detect, or add model to known_models list")
  return(NULL)
}


#' Auto-Detect Embedding Dimension via a Test API Call
#'
#' Sends a minimal test string (\code{"test"}) to the embedding API and
#' reads the number of columns in the returned matrix to infer the model's
#' output dimension. Called automatically by \code{get_embedding_dim} when
#' the model is not in the built-in lookup table.
#'
#' @param llm_api_key Character string. API key for the target provider.
#' @param embed_model Character string. Model identifier to probe
#'   (default: \code{"BAAI/bge-large-en-v1.5"}).
#' @param provider Character string. Embedding provider. One of
#'   \code{"openai"}, \code{"siliconflow"}, or \code{"gemini"}
#'   (default: \code{"siliconflow"}).
#'
#' @return Integer scalar giving the detected embedding dimension, or
#'   \code{NULL} if the test API call fails.
detect_embedding_dim <- function(llm_api_key,
                                 embed_model = "BAAI/bge-large-en-v1.5",
                                 provider    = "siliconflow") {

  test_emb <- get_embedding(
    text          = "test",
    llm_api_key   = llm_api_key,
    embed_model   = embed_model,
    provider      = provider,
    embedding_dim = NULL
  )
  
  if (is.null(test_emb)) {
    warning("Failed to detect embedding dimension - API call failed")
    return(NULL)
  }
  
  dim_detected <- ncol(test_emb)
  return(dim_detected)
}


# ============================================================
# 4. 核心 Embedding 调用
# ============================================================

#' Get Text Embeddings from an Embedding API
#'
#' Sends one or more text strings to a remote embedding endpoint and
#' returns a numeric matrix of embedding vectors. Supports OpenAI,
#' SiliconFlow, and Google Gemini providers. Includes global rate-limit
#' enforcement via a shared environment state.
#'
#' @param text Character vector. One or more non-empty strings to embed.
#'   \code{NA} and empty-string elements are silently removed before the
#'   API call. Returns \code{NULL} immediately if no valid text remains.
#' @param llm_api_key Character string. Bearer token or API key for the
#'   chosen provider. For Gemini this is appended as a URL query parameter;
#'   for OpenAI and SiliconFlow it is sent in the \code{Authorization}
#'   header.
#' @param embed_model Character string. Model identifier to use for
#'   embedding (default: \code{"text-embedding-3-large"}). Must be a model
#'   supported by the selected \code{provider}.
#' @param provider Character string. Embedding API provider. One of
#'   \code{"openai"} (default), \code{"siliconflow"}, or \code{"gemini"}.
#'   Controls endpoint URL, request body format, and response parsing.
#' @param embedding_dim Integer scalar or \code{NULL}. Expected number of
#'   dimensions in the returned embeddings. When provided, a warning is
#'   issued if the actual dimension differs. Pass \code{NULL} to skip the
#'   dimension check.
#'
#' @return Numeric matrix of shape \code{n_texts × embedding_dim}, where
#'   each row is the embedding vector for the corresponding input string.
#'   Returns \code{NULL} on any unrecoverable error (HTTP failure, empty
#'   response, parse error).
#'
#' @details
#' \strong{Rate limiting:} A provider-specific minimum inter-call interval
#' is enforced using \code{.emb_call_state$last_time}:
#' \code{gemini} = 2 s, \code{openai} / \code{siliconflow} = 0.3 s.
#'
#' \strong{Gemini restriction:} Gemini's \code{embedContent} endpoint
#' accepts only a single text per call; only \code{text[1]} is embedded
#' when \code{provider = "gemini"}.
#'
#' \strong{Error handling:} HTTP errors, JSON parse failures, and
#' API-level \code{error} / \code{error_msg} fields are caught and
#' reported via \code{message()}, returning \code{NULL} rather than
#' throwing.
get_embedding <- function(
    text,
    llm_api_key,
    embed_model   = "text-embedding-3-large",
    provider      = "openai",
    embedding_dim = NULL
) {
  
  ## ---- Global rate limit ----
  now <- Sys.time()
  if (!is.null(.emb_call_state$last_time)) {
    dt <- as.numeric(difftime(now, .emb_call_state$last_time, units = "secs"))
    min_interval <- switch(provider,
                           "gemini"      = 2,
                           "openai"      = 0.3,
                           "siliconflow" = 0.3,
                           1)
    if (dt < min_interval) Sys.sleep(min_interval - dt)
  }
  .emb_call_state$last_time <- Sys.time()
  
  ## ---- Input validation ----
  if (is.null(text) || length(text) == 0 || all(text == "")) {
    message("Empty text input")
    return(NULL)
  }
  text <- text[!is.na(text) & text != ""]
  if (length(text) == 0) return(NULL)
  
  ## ---- Configure provider ----
  if (provider == "siliconflow") {
    url  <- "https://api.siliconflow.cn/v1/embeddings"
    body <- list(model = embed_model,
                 input = if (length(text) == 1) text else as.list(text),
                 encoding_format = "float")
    headers <- httr::add_headers(Authorization = paste("Bearer", llm_api_key),
                                 "Content-Type" = "application/json")
    
  } else if (provider == "openai") {
    url  <- "https://api.openai.com/v1/embeddings"
    body <- list(model = embed_model,
                 input = if (length(text) == 1) text else as.list(text))
    headers <- httr::add_headers(Authorization = paste("Bearer", llm_api_key),
                                 "Content-Type" = "application/json")
    
  } else if (provider == "gemini") {
    url <- sprintf(
      "https://generativelanguage.googleapis.com/v1beta/models/%s:embedContent?key=%s",
      embed_model, llm_api_key)
    body    <- list(content = list(parts = list(list(text = text[1]))))
    headers <- httr::add_headers("Content-Type" = "application/json")
  }
  
  ## ---- Send request ----
  res_call <- tryCatch({
    httr::POST(url, headers,
               body = jsonlite::toJSON(body, auto_unbox = TRUE),
               httr::timeout(60))
  }, error = function(e) {
    message("HTTP error in get_embedding: ", e$message)
    return(NULL)
  })
  
  if (is.null(res_call)) return(NULL)
  
  if (httr::status_code(res_call) != 200) {
    message("HTTP status code: ", httr::status_code(res_call))
    message("Response: ", httr::content(res_call, "text", encoding = "UTF-8"))
    return(NULL)
  }
  
  ## ---- Parse response ----
  cont <- tryCatch(httr::content(res_call, encoding = "UTF-8"),
                   error = function(e) {
                     message("Failed to parse response: ", e$message)
                     return(NULL)
                   })
  
  if (is.null(cont)) return(NULL)
  
  if (!is.null(cont$error_msg)) {
    message("API error_msg: ", cont$error_msg)
    return(NULL)
  }
  if (!is.null(cont$error)) {
    message("API error: ",
            if (is.character(cont$error)) cont$error
            else jsonlite::toJSON(cont$error))
    return(NULL)
  }
  if (provider != "gemini" && (is.null(cont$data) || length(cont$data) == 0)) {
    message("No data field in response")
    return(NULL)
  }
  
  ## ---- Extract embeddings ----
  embs <- tryCatch({
    if (provider == "gemini") {
      list(unlist(cont$embedding$values))
    } else {
      lapply(cont$data, function(x) unlist(x$embedding))
    }
  }, error = function(e) {
    message("Failed to extract embeddings: ", e$message)
    return(NULL)
  })
  
  if (is.null(embs)) return(NULL)
  
  emb_matrix <- tryCatch(do.call(rbind, embs),
                         error = function(e) {
                           message("Failed to combine embeddings: ", e$message)
                           return(NULL)
                         })
  
  ## ---- Validate dimensions ----
  if (!is.null(emb_matrix)) {
    if (nrow(emb_matrix) != length(text))
      warning(sprintf("Dimension mismatch: %d texts but %d embeddings",
                      length(text), nrow(emb_matrix)))
    if (!is.null(embedding_dim) && ncol(emb_matrix) != embedding_dim)
      warning(sprintf("Unexpected embedding dimension: expected %d but got %d",
                      embedding_dim, ncol(emb_matrix)))
  }
  
  return(emb_matrix)
}


# ============================================================
# 5. 模块级摘要预计算
# ============================================================

#' Precompute PubMed Abstract Embeddings for an Entire Module
#'
#' Iterates over every metabolite in a module's PubMed result list and
#' converts each article context string into an embedding vector, storing
#' the results in a cache list for reuse by \code{embed_module_pubmed}.
#' Long contexts are split into fixed-size chunks and their chunk
#' embeddings are averaged into a single document vector.
#'
#' @param module_pub Named list. Each element is named by a metabolite ID
#'   and must contain:
#'   \itemize{
#'     \item \code{pmids}     — character vector of PubMed IDs.
#'     \item \code{titles}    — character vector of article titles.
#'     \item \code{abstracts} — character vector of abstract texts.
#'     \item \code{contexts}  — character vector of context strings to embed
#'       (typically title + abstract concatenated).
#'   }
#' @param llm_api_key Character string. API key forwarded to
#'   \code{get_embedding}.
#' @param embed_model Character string. Embedding model identifier
#'   (e.g. \code{"BAAI/bge-large-en-v1.5"}).
#' @param provider Character string. Embedding provider: \code{"openai"},
#'   \code{"siliconflow"}, or \code{"gemini"}.
#' @param max_chars_per_chunk Integer scalar. Maximum number of characters
#'   per embedding chunk when a context string exceeds this length
#'   (default: \code{1500}). Longer texts are split at this boundary and
#'   their chunk embeddings are averaged.
#' @param embedding_dim Integer scalar. Expected embedding dimension.
#'   \strong{Must not be \code{NULL}}; the function stops with an error if
#'   it is, because zero-padding of failed embeddings requires a known size.
#'
#' @return Named list (a subset of \code{names(module_pub)}) where each
#'   element contains:
#'   \itemize{
#'     \item \code{pmids}     — filtered PMIDs (non-empty contexts only).
#'     \item \code{titles}    — corresponding titles.
#'     \item \code{abstracts} — corresponding abstracts.
#'     \item \code{contexts}  — corresponding context strings.
#'     \item \code{emb}       — numeric matrix of shape
#'       \code{n_docs × embedding_dim}; rows with failed API calls are
#'       filled with zeros.
#'   }
#'   Metabolites with no valid (non-empty, non-NA) contexts are silently
#'   skipped and absent from the returned list.
#'
#' @details
#' For Gemini, each text is embedded in a separate API call (the endpoint
#' accepts only one input per request). For OpenAI and SiliconFlow,
#' multiple short texts are batched into a single vectorized call.
precompute_pubmed_embeddings_for_module <- function(
    module_pub,
    llm_api_key,
    embed_model,
    provider,
    max_chars_per_chunk = 1500,
    embedding_dim       = NULL
) {
  
  if (is.null(embedding_dim)) {
    stop("Cannot determine embedding dimension for model: ", embed_model,
         ". Please check your API connection.")
  }
  
  mets         <- names(module_pub)
  module_cache <- list()
  
  for (m in mets) {
    pubmed_obj <- module_pub[[m]]
    titles     <- pubmed_obj$titles
    abstracts  <- pubmed_obj$abstracts
    contexts   <- pubmed_obj$contexts
    pmids      <- pubmed_obj$pmids
    
    valid <- which(!is.na(contexts) & contexts != "")
    if (length(valid) == 0) next
    
    contexts  <- contexts[valid]
    abstracts <- abstracts[valid]
    titles    <- titles[valid]
    pmids     <- pmids[valid]
    
    n_docs       <- length(contexts)
    doc_emb_list <- vector("list", n_docs)
    
    n_chars   <- nchar(contexts, type = "chars")
    idx_short <- which(n_chars <= max_chars_per_chunk)
    idx_long  <- which(n_chars >  max_chars_per_chunk)
    
    ## -- Short contexts --
    if (length(idx_short) > 0) {
      if (provider == "gemini") {
        for (i in idx_short) {
          emb_i <- get_embedding(text = contexts[i], llm_api_key = llm_api_key,
                                 embed_model = embed_model, provider = provider,
                                 embedding_dim = embedding_dim)
          doc_emb_list[[i]] <- if (is.null(emb_i)) rep(0, embedding_dim) else as.numeric(emb_i[1, ])
        }
      } else {
        emb_short <- get_embedding(text = contexts[idx_short], llm_api_key = llm_api_key,
                                   embed_model = embed_model, provider = provider,
                                   embedding_dim = embedding_dim)
        if (is.null(emb_short)) {
          for (i in idx_short) doc_emb_list[[i]] <- rep(0, embedding_dim)
        } else {
          for (k in seq_along(idx_short)) {
            doc_emb_list[[idx_short[k]]] <- as.numeric(emb_short[k, ])
          }
        }
      }
    }
    
    ## -- Long contexts (chunked) --
    if (length(idx_long) > 0) {
      for (i in idx_long) {
        txt    <- contexts[i]
        n      <- nchar(txt, type = "chars")
        starts <- seq(1, n, by = max_chars_per_chunk)
        ends   <- pmin(starts + max_chars_per_chunk - 1, n)
        chunks <- substring(txt, starts, ends)
        
        if (provider == "gemini") {
          emb_list <- lapply(chunks, function(ch) {
            get_embedding(text = ch, llm_api_key = llm_api_key,
                          embed_model = embed_model, provider = provider,
                          embedding_dim = embedding_dim)
          })
          emb_list <- emb_list[!sapply(emb_list, is.null)]
          doc_emb_list[[i]] <- if (length(emb_list) == 0) rep(0, embedding_dim) else
            colMeans(do.call(rbind, emb_list))
        } else {
          chunk_emb <- get_embedding(text = chunks, llm_api_key = llm_api_key,
                                     embed_model = embed_model, provider = provider,
                                     embedding_dim = embedding_dim)
          doc_emb_list[[i]] <- if (is.null(chunk_emb)) rep(0, embedding_dim) else
            colMeans(chunk_emb)
        }
      }
    }
    
    emb_mat <- do.call(rbind, doc_emb_list)
    
    if (!is.null(emb_mat) && ncol(emb_mat) != embedding_dim)
      warning("Unexpected embedding dimension for ", m,
              ": expected ", embedding_dim, " but got ", ncol(emb_mat))
    
    module_cache[[m]] <- list(
      pmids     = pmids,
      titles    = titles,
      abstracts = abstracts,
      contexts  = contexts,
      emb       = emb_mat
    )
  }
  
  return(module_cache)
}


# ============================================================
# 6. 单个代谢物文章排名
# ============================================================

#' Embed and Rank PubMed Articles for a Single Metabolite
#'
#' Computes cosine similarity between a set of metabolite-specific query
#' embeddings and the embedding of each candidate PubMed abstract, then
#' returns the top-scoring articles that exceed a minimum similarity
#' threshold.
#'
#' @param met_name Character string. Metabolite database ID (e.g. HMDB or
#'   KEGG ID). Used only for log messages.
#' @param pubmed_obj Named list. PubMed data for the metabolite, containing:
#'   \itemize{
#'     \item \code{pmids}     — character vector of PubMed IDs.
#'     \item \code{titles}    — character vector of article titles.
#'     \item \code{abstracts} — character vector of abstract texts.
#'     \item \code{contexts}  — character vector of text strings to embed.
#'   }
#' @param met_real_name Character string. Human-readable metabolite name
#'   used to construct query sentences via \code{build_metabolite_queries}.
#' @param similarity Numeric scalar in \code{[0, 1]}. Minimum mean cosine
#'   similarity required for an article to be retained. Articles scoring
#'   below this threshold are discarded.
#' @param top_score Integer scalar. Maximum number of articles to return
#'   after filtering and ranking by similarity score (default: \code{5}).
#' @param llm_api_key Character string. API key forwarded to
#'   \code{get_embedding}.
#' @param embed_model Character string. Embedding model identifier
#'   (default: \code{"BAAI/bge-large-en-v1.5"}).
#' @param provider Character string. Embedding provider: \code{"openai"},
#'   \code{"siliconflow"}, or \code{"gemini"}
#'   (default: \code{"siliconflow"}).
#' @param abs_embeds Numeric matrix or \code{NULL}. Precomputed abstract
#'   embeddings (rows = articles). When supplied, the function skips the
#'   abstract embedding API call and uses this matrix directly. Row count
#'   must match \code{length(pubmed_obj$contexts)}; a mismatch triggers
#'   row-count truncation to the shorter of the two.
#' @param embedding_dim Integer scalar or \code{NULL}. Expected embedding
#'   dimension; used for a consistency check when \code{abs_embeds} is
#'   provided. Pass \code{NULL} to skip the check.
#'
#' @return Named list on success:
#'   \itemize{
#'     \item \code{met_name}  — value of \code{met_real_name}.
#'     \item \code{pmids}     — PMIDs of retained top articles.
#'     \item \code{titles}    — corresponding titles.
#'     \item \code{abstracts} — corresponding abstracts.
#'     \item \code{contexts}  — corresponding context strings.
#'     \item \code{scores}    — mean cosine similarity scores for each
#'       retained article, in descending order.
#'   }
#'   Returns \code{NULL} when no valid contexts exist, query embedding
#'   fails, abstract embedding fails, or no article meets the similarity
#'   threshold.
#'
#' @details
#' Three query sentences are generated by \code{build_metabolite_queries}
#' and each is independently embedded. The similarity score for an article
#' is the mean cosine similarity across all three query vectors. Articles
#' are ranked in descending score order and the top \code{top_score} are
#' returned.
embed_pubmed_for_metabolite <- function(
    met_name,
    pubmed_obj,
    met_real_name,
    similarity,
    top_score     = 5,
    llm_api_key,
    embed_model   = "BAAI/bge-large-en-v1.5",
    provider      = "siliconflow",
    abs_embeds    = NULL,
    embedding_dim = NULL
) {
  
  titles    <- pubmed_obj$titles
  contexts  <- pubmed_obj$contexts
  abstracts <- pubmed_obj$abstracts
  pmids     <- pubmed_obj$pmids
  
  valid <- which(!is.na(contexts) & contexts != "")
  if (length(valid) == 0) return(NULL)
  
  titles    <- titles[valid]
  abstracts <- abstracts[valid]
  pmids     <- pmids[valid]
  contexts  <- contexts[valid]
  
  ## 1) Query embeddings
  queries <- unname(unlist(build_metabolite_queries(met_real_name)))
  
  if (provider == "gemini") {
    q_emb_list <- lapply(queries, function(q)
      get_embedding(text = q, llm_api_key = llm_api_key,
                    embed_model = embed_model, provider = provider,
                    embedding_dim = embedding_dim))
    q_emb_list  <- q_emb_list[!sapply(q_emb_list, is.null)]
    if (length(q_emb_list) == 0) {
      message("Failed to get query embeddings for metabolite: ", met_name)
      return(NULL)
    }
    query_embeds <- do.call(rbind, q_emb_list)
  } else {
    query_embeds <- get_embedding(text = queries, llm_api_key = llm_api_key,
                                  embed_model = embed_model, provider = provider,
                                  embedding_dim = embedding_dim)
  }
  
  if (is.null(query_embeds)) {
    message("Failed to get query embeddings for metabolite: ", met_name)
    return(NULL)
  }
  
  ## 2) Abstract embeddings
  if (!is.null(abs_embeds)) {
    abs_embeds_use <- abs_embeds
    if (nrow(abs_embeds_use) != length(contexts)) {
      n_min          <- min(nrow(abs_embeds_use), length(contexts))
      abs_embeds_use <- abs_embeds_use[seq_len(n_min), , drop = FALSE]
      contexts       <- contexts[seq_len(n_min)]
      abstracts      <- abstracts[seq_len(n_min)]
      titles         <- titles[seq_len(n_min)]
      pmids          <- pmids[seq_len(n_min)]
    }
    if (!is.null(embedding_dim) && ncol(abs_embeds_use) != ncol(query_embeds))
      stop("Embedding dimension mismatch between abstracts and queries")
  } else {
    abs_embeds_use <- get_embedding(text = contexts, llm_api_key = llm_api_key,
                                    embed_model = embed_model, provider = provider,
                                    embedding_dim = embedding_dim)
    if (is.null(abs_embeds_use)) {
      message("Failed to get abstract embeddings for metabolite: ", met_name)
      return(NULL)
    }
  }
  
  ## 3) Similarity scoring
  sim_scores <- sapply(seq_len(nrow(abs_embeds_use)), function(i) {
    mean(apply(query_embeds, 1, function(qe) cosine_sim(abs_embeds_use[i, ], qe)))
  })
  sim_scores[is.na(sim_scores) | is.infinite(sim_scores)] <- 0
  
  keep_idx <- which(sim_scores >= similarity)
  if (length(keep_idx) == 0) {
    return(NULL)
  }
  
  ## 4) Top-k selection
  ordered_keep <- keep_idx[order(sim_scores[keep_idx], decreasing = TRUE)]
  top_id       <- ordered_keep[seq_len(min(top_score, length(ordered_keep)))]
  
  list(
    met_name  = met_real_name,
    pmids     = pmids[top_id],
    titles    = titles[top_id],
    abstracts = abstracts[top_id],
    contexts  = contexts[top_id],
    scores    = sim_scores[top_id]
  )
}


# ============================================================
# 7. 模块级批量排名（入口函数）
# ============================================================

#' Embed and Rank PubMed Articles for All Metabolites in a Module
#'
#' Iterates over every metabolite in a module PubMed result list, calls
#' \code{embed_pubmed_for_metabolite} for each one, and collects the
#' ranked results into a single named list. Optionally uses precomputed
#' abstract embeddings to avoid redundant API calls.
#'
#' @param module_pub Named list. Per-metabolite PubMed data as returned by
#'   the PubMed annotation pipeline. Each element must contain \code{pmids},
#'   \code{titles}, \code{abstracts}, and \code{contexts}.
#' @param llm_api_key Character string. API key forwarded to
#'   \code{get_embedding} for query and (if needed) abstract embedding.
#' @param refdata Character string. Reference database used to resolve the
#'   human-readable metabolite name from \code{module_df}. One of
#'   \code{"hmdb"} or \code{"kegg"}.
#' @param module_df Data frame or \code{NULL}. Source annotation table with
#'   at least a database ID column (\code{HMDB_ID} or \code{KEGG_ID}) and a
#'   \code{name} column. Used to look up the display name for each metabolite.
#'   When \code{NULL} or the ID is not found, the database ID itself is used
#'   as the metabolite name.
#' @param embed_model Character string. Embedding model identifier
#'   (default: \code{"BAAI/bge-large-en-v1.5"}).
#' @param provider Character string. Embedding provider: \code{"openai"},
#'   \code{"siliconflow"}, or \code{"gemini"}
#'   (default: \code{"siliconflow"}).
#' @param emb_cache Named list or \code{NULL}. Precomputed embedding cache
#'   as returned by \code{precompute_pubmed_embeddings_for_module}. When an
#'   entry for a metabolite is present in the cache, its stored embeddings
#'   and context data are used directly, bypassing a second abstract
#'   embedding API call. Pass \code{NULL} to recompute all embeddings
#'   on-the-fly.
#' @param top_score Integer scalar. Maximum number of top-ranked articles
#'   to retain per metabolite (default: \code{5}).
#' @param min_similarity Numeric scalar in \code{[0, 1]} or \code{NULL}.
#'   Minimum cosine similarity threshold passed to
#'   \code{embed_pubmed_for_metabolite}. Articles scoring below this value
#'   are discarded.
#' @param embedding_dim Integer scalar or \code{NULL}. Expected embedding
#'   dimension; forwarded to \code{embed_pubmed_for_metabolite} for
#'   consistency checking.
#'
#' @return Named list where each element is named by the metabolite database
#'   ID and contains the result of \code{embed_pubmed_for_metabolite}
#'   (\code{met_name}, \code{pmids}, \code{titles}, \code{abstracts},
#'   \code{contexts}, \code{scores}). Metabolites for which ranking returned
#'   \code{NULL} (no valid articles or errors) are silently excluded from
#'   the output.
#'
#' @details
#' Errors thrown inside \code{embed_pubmed_for_metabolite} are caught with
#' \code{tryCatch}; a message is printed and the metabolite is set to
#' \code{NULL} (then dropped from the result). This ensures that a single
#' failed metabolite does not abort the entire module.
embed_module_pubmed <- function(
    module_pub,
    llm_api_key,
    refdata        = c("hmdb", "kegg", "auto"),
    module_df      = NULL,
    embed_model    = "BAAI/bge-large-en-v1.5",
    provider       = "siliconflow",
    emb_cache      = NULL,
    top_score      = 5,
    min_similarity = NULL,
    embedding_dim  = NULL
) {
  refdata <- match.arg(refdata)
  
  mets <- names(module_pub)
  rak  <- list()
  
  for (m in mets) {
    
    ## Resolve human-readable name
    met_real_name <- m
    if (!is.null(module_df) && "Compound_ID" %in% colnames(module_df)) {
      hit <- module_df[module_df$Compound_ID == m, "name"]
      if (length(hit) > 0) met_real_name <- hit[1]
    } else if (refdata == "hmdb") {
      if (!is.null(module_df) && "HMDB_ID" %in% colnames(module_df)) {
        hit <- module_df[module_df$HMDB_ID == m, "name"]
        if (length(hit) > 0) met_real_name <- hit[1]
      }
    } else if (refdata == "kegg") {
      if (!is.null(module_df) && "KEGG_ID" %in% colnames(module_df)) {
        hit <- module_df[module_df$KEGG_ID == m, "name"]
        if (length(hit) > 0) met_real_name <- hit[1]
      }
    }
    
    ## Use cache if available, otherwise fall back to raw pubmed_obj
    if (!is.null(emb_cache) && !is.null(emb_cache[[m]])) {
      cache_m <- emb_cache[[m]]
      pubmed_obj_use <- list(pmids     = cache_m$pmids,
                             titles    = cache_m$titles,
                             abstracts = cache_m$abstracts,
                             contexts  = cache_m$contexts)
      abs_embeds_use <- cache_m$emb
    } else {
      pubmed_obj_use <- module_pub[[m]]
      abs_embeds_use <- NULL
    }
    
    rak[[m]] <- tryCatch({
      embed_pubmed_for_metabolite(
        met_name      = m,
        met_real_name = met_real_name,
        pubmed_obj    = pubmed_obj_use,
        llm_api_key   = llm_api_key,
        embed_model   = embed_model,
        provider      = provider,
        abs_embeds    = abs_embeds_use,
        top_score     = top_score,
        embedding_dim = embedding_dim,
        similarity    = min_similarity
      )
    }, error = function(e) {
      message("Error processing ", m, ": ", e$message)
      return(NULL)
    })
  }
  
  rak <- rak[!sapply(rak, is.null)]
  return(rak)
}
