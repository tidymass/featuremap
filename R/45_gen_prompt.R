

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
   - Write a 2-5 word name for one primary biochemical pathway, enzymatic process, or biological transformation.
   - Choose the most specific coherent pathway/process supported by the early-ranked high-priority metabolites.
   - Be precise: name the specific enzymatic step, metabolic route, or biological process — not just the compound class.
     Good examples: \"COMT-Mediated Catecholamine Clearance\", \"Sphingolipid De Novo Biosynthesis\", \"Branched-Chain Amino Acid Catabolism\", \"Cysteine Sulfinic Acid Oxidation\".
     Bad examples (too generic): \"Tyrosine Metabolism\", \"Aromatic Acid Metabolism\", \"Exogenous Compound Metabolism\".
   - Do NOT end the title with just \"Metabolism\" unless the module covers a genuinely broad, multi-step metabolic network.
   - If no canonical pathway is clear, describe the most specific chemical transformation or biotransformation observed (e.g., \"O-Methylation of Catechol Amines\", \"Phenolic Glucuronide Conjugation\").

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
   - High: the pathway is an established function in the study species (a), is appropriate for the biospecimen type (b), and is plausibly relevant to the study phenotype (c). Core module metabolites are recognized substrates or products of this pathway.
   - Medium: the pathway exists in the study species and fits the biospecimen, but one dimension is uncertain or indirect — e.g., phenotype relevance is speculative, substrates are primarily exogenous, or the pathway is active but not specific to the biospecimen.
   - Low: the pathway fails on the species dimension — it is not a known function in the study species (e.g., bacterial or plant biosynthesis), OR the module is dominated by synthetic, industrial, or drug-derived compounds that reflect exposure rather than active endogenous metabolism.
   - Module-level statistics (signal concentration, structural stability, annotation completeness) are secondary; they may refine confidence within a tier but cannot override a biologically implausible assignment.

4) RATIONALE
   - Write 1-2 sentences explaining the confidence decision, explicitly addressing the three dimensions: species validity, biospecimen appropriateness, and phenotype relevance.
   - Mention module-level statistics only if they provide meaningful additional support or concern.

Return STRICT JSON only:
{{
  \"title\": \"...\",
  \"summary\": \"...\",
  \"confidence_label\": \"high\" | \"medium\" | \"low\",
  \"rationale\": \"1-2 short sentences explaining the confidence decision.\"
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
