

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
- Priority_level: interpretation priority within the module; early-ranked high-priority metabolites should anchor the main interpretation.
- Literature-derived texts: supporting context only.

Using only the provided information, return four fields:

1) TITLE
   - Name the pathway in KEGG style: concise, standardized, and anchored to a specific biochemical pathway or process.
   - Follow KEGG naming conventions:
     - Use \"X metabolism\" when the module spans multiple steps of a compound class (e.g., \"Tyrosine metabolism\", \"Sphingolipid metabolism\").
     - Use \"X biosynthesis\" or \"X degradation\" when the directionality is clear from the metabolites (e.g., \"Sphingolipid biosynthesis\", \"Fatty acid degradation\").
     - Use a specific reaction or enzyme-level name when the module is anchored by a single well-defined step (e.g., \"Cysteine sulfinic acid oxidation\", \"Bilirubin glucuronidation\").
     - Use a physiological or system-level name when appropriate (e.g., \"Insulin signaling pathway\", \"MAPK signaling pathway\").
   - Choose the most specific level of granularity supported by the early-ranked high-priority metabolites — do not over-generalize (e.g., prefer \"Tyrosine metabolism\" over \"Amino acid metabolism\") but do not force artificial precision when the module is genuinely broad.
   - If the module's top metabolites point to a specific sub-pathway or reaction within a broader class, reflect that specificity in the title rather than defaulting to the parent pathway name (e.g., if a tyrosine module is anchored by catecholamine O-methylation products, prefer \"Catecholamine metabolism\" or \"Catecholamine O-methylation\" over the generic \"Tyrosine metabolism\").
   - If no canonical pathway fits, use the most specific neutral description of the dominant chemical transformation (e.g., \"Aromatic compound glucuronidation\", \"Phenolic sulfation\").

2) SUMMARY
   - Write 4-6 sentences as a factual biochemical annotation of the module.
   - Describe what the metabolites do biochemically: enzymatic conversions, pathway steps, metabolite roles, and physiological functions — as if writing an entry in a biochemical database.
   - Base the explanation mainly on early-ranked high-priority metabolites; use class structure and literature text only to support or qualify.
   - STRICTLY FORBIDDEN in the summary: do NOT state or imply whether metabolites are endogenous or exogenous, dietary, microbial, plant-derived, environmental, or whether the pathway is biologically plausible. Those judgments belong exclusively in RATIONALE.
   - Do NOT reference module-level statistics (signal concentration, annotation completeness, stability) in the summary.
   - Mention secondary signals only briefly at the end.

3) CONFIDENCE
   - Evaluate confidence across three dimensions using the Study context (species, biospecimen, phenotype):
     (a) Species validity: Does the study species possess this metabolic function? Is it an established pathway in that organism, or is it specific to other taxa (bacteria, plants, fungi)?
     (b) Biospecimen appropriateness: Is this pathway known to be active or detectable in the study biospecimen? Some pathways are tissue-specific (e.g., liver-specific) and may not be reflected in plasma or urine.
     (c) Phenotype relevance: Is this pathway plausibly perturbed in the study phenotype, or is it consistent with the biological question being asked?
   - Score holistically across all three dimensions:
   - High: ALL of the following must hold — (a) the pathway is an established endogenous function in the study species, (b) it is active and detectable in the study biospecimen, (c) it is plausibly relevant to the study phenotype, AND the early-ranked high-priority metabolites are themselves recognized endogenous substrates or products of this pathway in that species. If the module's top metabolites are primarily exogenous compounds (drugs, plant-derived alkaloids, industrial chemicals, environmental contaminants), the confidence CANNOT be high even if the named reaction type exists in the study species.
   - Medium: the pathway exists in the study species and fits the biospecimen, but at least one of the following applies — phenotype relevance is indirect or speculative; the dominant metabolites are primarily exogenous substrates being processed by endogenous enzymes (xenobiotic metabolism); or pathway support is mixed between endogenous and exogenous signals.
   - Low: the pathway fails on species validity — it is not a known function in the study species (e.g., bacterial cell wall biosynthesis, plant phenylpropanoid pathway), OR the module is overwhelmingly dominated by synthetic, industrial, or environmental compounds with no plausible endogenous metabolic interpretation.
   - Module-level statistics (signal concentration, structural stability, annotation completeness) are secondary; they may refine confidence within a tier but cannot override a biologically implausible assignment.

4) RATIONALE
   - Write 2-3 sentences explaining the confidence decision, explicitly addressing all three dimensions: species validity, biospecimen appropriateness, and phenotype relevance.
   - Mention module-level statistics only if they provide meaningful additional support or concern.

Return STRICT JSON only:
{{
  \"title\": \"...\",
  \"summary\": \"...\",
  \"confidence_label\": \"high\" | \"medium\" | \"low\",
  \"rationale\": \"2-3 sentences covering species validity, biospecimen appropriateness, and phenotype relevance.\"
}}

Do not invent unsupported pathways or mechanisms."
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

Use this context mainly to judge whether the titled pathway is biologically realistic in the study setting."
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
