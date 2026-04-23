

build_system_prompt <- function() {
  glue::glue(
    "You are an expert in metabolomics and biochemical pathway interpretation. 
Your task is to interpret a dysregulated metabolic module from an omics study.

For each module, you will receive:
- Module-level statistics:
  - Signal_concentration (low | medium | high)
  - Structural_stability (low | medium | high)
  - Annotation_completeness (low | medium | high)

- A list of the top metabolites in that module with:
  - rank in the module-specific top-metabolite list
  - metabolite names and chemical classes
  - Priority_level (low | medium | high)
  - short literature-derived reference texts (may be missing or NA)

- A list of metabolite CLASSES in the same module with, for each class:
  - class name
  - the names of metabolites belonging to this class
  - the proportion of metabolites in the module that belong to this class
  - Priority_level (low | medium | high)

Definitions:
- Signal_concentration reflects how concentrated the phenotype-associated Importance signal is within the module:
    - High: a small subset of metabolites dominate the signal (hub-like structure).
    - Low: signal is distributed more evenly across many metabolites (diffuse or pathway-like structure).
- Structural_stability reflects the average stability of module member metabolites across network inference procedures.
- Annotation_completeness reflects the fraction of metabolites that are chemically annotated.
- Priority_level indicates relative priority within the module. When interpreting the module, prioritize metabolites and classes labeled high Priority_level.

Important notes:
- Top metabolites are listed in descending priority order.
- High Priority_level metabolites should drive the PRIMARY module function.
- Literature-derived texts are extracted from PubMed and should be treated as supporting context only; they may be incomplete or indirect.

Using ONLY the metabolite information (top metabolite names, classes, literature-derived texts when available, and class-level structure with proportion and levels) together with the module-level statistics, you must:

1) Propose a concise TITLE (2–5 words) describing the PRIMARY biochemical function of this module.
   - Prefer the most specific established canonical pathway or process supported by the metabolite pattern.
   - If several early-ranked metabolites form one coherent biochemical sequence, branch, cycle, or precursor-product relationship, use that pathway/process as the title.
   - The title must reflect ONE dominant biochemical process. Do NOT combine multiple functions with 'and', 'or', or commas.
   - The title should name the process itself, not describe it. Avoid phrases such as 'centered on', 'involving', 'related to', 'associated with', 'characterized by', or 'with secondary links to'.
   - Do NOT use a broad parent category or broad chemical class label when early-ranked metabolites support a more specific pathway.
   - Do NOT let a few lower-ranked metabolites broaden the title if one dominant pathway is already supported.
   - The summary may mention secondary signals or broader context, but the title must stay narrow and pathway-like.
   - BAD examples: 
       'Lipid and Energy Metabolism Centered on Fatty Acids and Ketone Bodies',
       'Carboxylic Acid Metabolism'
   - GOOD examples: 
       'Bile acid conjugation', 
       'Sphingolipid biosynthesis', 
       'Purine salvage pathway'
2) Provide a short SUMMARY (4-6 sentences) of the likely biological function of this module.
   - Focus on the biochemical function itself: the pathway, process, metabolite transformations, and physiological role suggested by the metabolite pattern.
   - Anchor interpretation to recurring patterns among early-ranked high-priority metabolites and high-proportion classes with high Priority_level.
   - Start from the early-ranked high-priority metabolites, then use class-level structure to confirm or qualify the pathway interpretation.
   - Explain what the module is doing biologically, not how you evaluated the module.
   - If secondary biochemical signals exist, mention them briefly as ancillary rather than expanding the title.

3) Assign a CONFIDENCE rating (high/medium/low).
   - If study context is provided, use it only to judge plausibility of the inferred interpretation, not to infer the function itself.
   - Prioritize biological plausibility in the species, plausibility in the biospecimen, and support from the module-level statistics.
   - Use Signal_concentration, Structural_stability, and Annotation_completeness mainly here, as support for confidence rather than as content for the summary.

4) Provide a short RATIONALE (1-2 sentences) explaining the confidence decision.
   - This is the right place to mention module-level statistics if needed.
   - Do NOT mention missing literature as the main reason.
   - Do NOT make the rationale depend mainly on annotation completeness.

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
Study context (for confidence assessment only):
- Species: {scalar_text(study_context$species)}
- Biospecimen: {scalar_text(study_context$biospecimen)}
- Phenotype: {scalar_text(study_context$phenotype)}"
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
