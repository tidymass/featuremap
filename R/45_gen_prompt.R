

build_system_prompt <- function() {
  glue::glue(
    "You are an expert in metabolomics and biochemical pathway interpretation. 
Your task is to interpret dysregulated metabolic network modules from an omics study.

For each module, you will receive:
- Study context (biospecimen, phenotype of interest, and species)
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
  - Importance_level (low | medium | high)
  - Robustness_level (low | medium | high)

Definitions:
- Signal_concentration reflects how concentrated the phenotype-associated Importance signal is within the module:
    - High: a small subset of metabolites dominate the signal (hub-like structure).
    - Low: signal is distributed more evenly across many metabolites (diffuse or pathway-like structure).
- Structural_stability reflects the average stability of module member metabolites across network inference procedures.
- Annotation_completeness reflects the fraction of metabolites that are chemically annotated.
- Importance reflects the strength of statistical association with the study phenotype.
- Robustness reflects structural stability across network inference or resampling.
- Class-level Importance_level and Robustness_level summarize the average node-level behavior within that class.

Important notes:
- Top metabolites are listed in descending priority order, with rank 1 being the strongest signal in the module.
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
2) Provide a short SUMMARY (5-8 sentences) of the likely biological function of this module.
   - Anchor interpretation to recurring patterns among top metabolites and high-proportion classes with high Importance_level.
   - Start from the early-ranked high-priority metabolites, then use class-level structure to confirm or qualify the pathway interpretation.
   - Use Structural_stability to assess reliability of the module structure.
   - If Signal_concentration is high, consider whether a small number of metabolites are likely driving the function.
   - If Signal_concentration is low, consider whether the module represents a coordinated pathway-level shift.
   - If secondary biochemical signals exist, mention them briefly as ancillary rather than expanding the title.

VERY IMPORTANT:
- When inferring the biological function and TITLE of the module, you must NOT use the phenotype, biospecimen and species information.

Then:

3) Assign a CONFIDENCE rating (high/medium/low) with this priority:

PRIMARY (most important): Biological plausibility in the study context
- Is the inferred module function biologically plausible in the given species?
- Is it plausible that this function/pathway would be observable or reflected in the given biospecimen?
  (e.g., plasma/serum reflects systemic metabolism and circulating lipids/amino acids; tissue-specific processes may be weakly reflected.)
- Is it plausible that this function could relate to the phenotype (ONLY as plausibility, not for inferring function)?

Rationale rules:
- 1-2 sentences only.
- Do NOT mention \"lack of literature\" or \"annotation completeness is low\" as the main reason.

Return STRICT JSON only:
{{
  \"title\": \"...\",
  \"summary\": \"...\",
  \"confidence_label\": \"high\" | \"medium\" | \"low\",
  \"rationale\": \"1-3 short sentences explaining the confidence decision.\"
}}

If evidence is weak or incoherent, clearly state this and assign LOW confidence.
Do NOT invent unsupported pathways or mechanisms."
  )
}




build_module_prompt <- function(module_obj, study_context) {
  
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
- Importance_level: {scalar_text(row$importance_level)}
- Robustness_level: {scalar_text(row$robustness_level)}",
        .na = "NA"
      )
    })
    
    class_text <- paste(class_blocks, collapse = "\n\n")
  }
  
  # ---- Study context ----
  sp  <- scalar_text(study_context$species)
  bio <- scalar_text(study_context$biospecimen)
  phe <- scalar_text(study_context$phenotype)
  
  # ---- Final prompt ----
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

{class_text}

Study context (for confidence assessment only):
- Species: {sp}
- Biospecimen: {bio}
- Phenotype: {phe}
")
  
  as.character(prompt)
}
