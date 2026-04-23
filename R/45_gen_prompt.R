

build_system_prompt <- function() {
  glue::glue(
    "You are an expert in metabolomics and biochemical pathway interpretation.
Your task is to interpret one dysregulated metabolic module from an omics study.

For each module, you will receive:
- Module-level labels: Signal_concentration, Structural_stability, Annotation_completeness (each low | medium | high)
- Top metabolites: rank, metabolite name, chemical class, Priority_level (low | medium | high), and optional literature-derived reference text
- Module classes: class name, member metabolite names, proportion in the module, and Priority_level

Definitions:
- Signal_concentration: whether phenotype-associated Importance is concentrated in a few metabolites (high) or distributed across many metabolites (low).
- Structural_stability: average stability of module member metabolites across network inference procedures.
- Annotation_completeness: fraction of module member metabolites that are chemically annotated.
- Priority_level: interpretation priority within the module; high-priority metabolites and classes should drive the primary function.
- Literature-derived texts are supporting context only; they may be incomplete or indirect.

Using ONLY the provided metabolite/class information together with the module-level statistics, you must:

1) Propose a concise TITLE (2–5 words) describing the PRIMARY biochemical function of this module.
   - Prefer the most specific established canonical pathway or process supported by early-ranked high-priority metabolites.
   - If several early-ranked metabolites form one coherent biochemical sequence, branch, cycle, or precursor-product relationship, use that pathway/process as the title.
   - The title must reflect ONE dominant biochemical process. Do NOT combine multiple functions with 'and', 'or', or commas.
   - Name the process itself, not a description of it. Avoid phrases such as 'centered on', 'involving', 'related to', 'associated with', 'characterized by', or 'with secondary links to'.
   - Prefer a specific pathway over a broad parent category or chemical class.
   - Secondary signals may appear in the summary, but not in the title.
   - Bad title style: 'Lipid and Energy Metabolism Centered on Fatty Acids and Ketone Bodies'.
   - Good title style: 'Bile acid conjugation' or 'Sphingolipid biosynthesis'.
2) Provide a short SUMMARY (4-6 sentences) of the likely biological function of this module.
   - Focus on the biochemical function itself: the pathway, process, metabolite transformations, and physiological role suggested by the metabolite pattern.
   - Ground the interpretation in early-ranked high-priority metabolites, then use high-priority or high-proportion classes to support or qualify the pathway interpretation.
   - Explain what the module is doing biologically, not how you evaluated the module.
   - If secondary biochemical signals exist, mention them briefly as ancillary rather than expanding the title.

3) Assign a CONFIDENCE rating (high/medium/low).
   - If species is provided, you may use it to prefer biologically plausible pathway interpretations and pathway naming for that organism.
   - If biospecimen or phenotype is provided, use them only to judge plausibility of the inferred interpretation, not to infer the function itself.
   - Judge confidence using biological plausibility in the species and biospecimen, together with support from the module-level statistics.
   - Use Signal_concentration, Structural_stability, and Annotation_completeness mainly here, not as content for the summary.

4) Provide a short RATIONALE (1-2 sentences) explaining the confidence decision.
   - This is the right place to mention biological plausibility and module-level statistics if needed.
   - Do NOT mention missing literature as the main reason.
   - Do NOT make the rationale depend mainly on module-level statistics.

Return STRICT JSON only:
{{
  \"title\": \"...\",
  \"summary\": \"...\",
  \"confidence_label\": \"high\" | \"medium\" | \"low\",
  \"rationale\": \"1-2 short sentences explaining the confidence decision.\"
}}

Do NOT invent unsupported pathways or mechanisms."
  )
}




build_module_prompt <- function(module_obj, study_context = NULL) {
  
  scalar_text <- function(x, na = "NA") {
    if (is.null(x) || length(x) == 0) return(na)
    x <- as.character(x[1])
    x <- trimws(x)
    if (!nzchar(x) || is.na(x)) na else x
  }
  
  # ---- Module ID ----
  mid <- scalar_text(module_obj$module_id, na = "NA")
  
  # ---- Module-level stats (LEVELS only) ----
  stats <- module_obj$module_stats %||% list()
  
  signal_level   <- scalar_text(stats$Signal_concentration)
  stability_level <- scalar_text(stats$Structural_stability)
  annotation_level <- scalar_text(stats$Annotation_completeness)
  
  # ---- Top metabolites ----
  top_df <- module_obj$topN_metabolites
  
  if (is.null(top_df) || !is.data.frame(top_df) || nrow(top_df) == 0) {
    top_text <- "NA (no top metabolites available)"
  } else {
    
    if (!"functional_summary" %in% colnames(top_df)) {
      top_df$functional_summary <- NA_character_
    }
    
    top_blocks <- apply(top_df, 1, function(row) {
      row <- as.list(row)
      
      glue::glue(
        "Metabolite (top):
- Top_rank: {scalar_text(row$top_rank)}
- Name: {scalar_text(row$Compound_name)}
- Class: {scalar_text(row$class)}
- Priority_level: {scalar_text(row$priority_level)}
- Literature_reference_text: {scalar_text(row$functional_summary)}",
        .na = "NA"
      )
    })
    
    top_text <- paste(top_blocks, collapse = "\n\n")
  }
  
  # ---- Class-level structure (ALL metabolites, no non-top concept) ----
  class_df <- module_obj$other_metabolites
  
  if (is.null(class_df) || !is.data.frame(class_df) || nrow(class_df) == 0) {
    class_text <- "No class-level structure available."
  } else {
    
    class_blocks <- apply(class_df, 1, function(row) {
      row <- as.list(row)
      
      met_vec <- tryCatch(
        clean_field(row$metabolites),
        error = function(e) row$metabolites
      )
      
      met_str <- if (is.null(met_vec) || length(met_vec) == 0) {
        "NA"
      } else {
        paste(as.character(met_vec), collapse = ", ")
      }
      
      glue::glue(
        "Metabolite class:
- Class_name: {scalar_text(row$class)}
- Proportion_of_module: {scalar_text(row$n_metabolites)}
- Metabolite_names: {met_str}
- Priority_level: {scalar_text(row$priority_level)}",
        .na = "NA"
      )
    })
    
    class_text <- paste(class_blocks, collapse = "\n\n")
  }
  
  # ---- Final prompt ----
  has_context <- !is.null(study_context) && any(vapply(
    study_context[c("species", "biospecimen", "phenotype")],
    function(x) {
      if (is.null(x) || length(x) == 0) return(FALSE)
      val <- trimws(as.character(x[1]))
      nzchar(val) && !identical(val, "NA")
    },
    logical(1)
  ))

  context_text <- ""
  if (has_context) {
    context_text <- glue::glue(
      "
Study context:
- Species: {scalar_text(study_context$species)}
- Biospecimen: {scalar_text(study_context$biospecimen)}
- Phenotype: {scalar_text(study_context$phenotype)}

Use species to keep the interpretation biologically plausible for the organism when relevant.
Use biospecimen and phenotype mainly for plausibility and confidence assessment, not to infer the core function."
    )
  }

  prompt <- glue::glue("
Module information:
- Module ID: {mid}
- Signal_concentration: {signal_level}
- Structural_stability: {stability_level}
- Annotation_completeness: {annotation_level}

Top metabolites (strongest statistical signals in the module):
Listed in descending priority order, where smaller Top_rank means stronger evidence.

{top_text}

Class-level structure of the module (all metabolites grouped by class):

{class_text}{context_text}
")
  
  as.character(prompt)
}
