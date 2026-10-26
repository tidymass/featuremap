# featuremap

`featuremap` is an R package for identifying and interpreting functional
metabolic modules from untargeted metabolomics data. It connects metabolite
annotation, metabolic network extraction, module identification, and optional
literature- and LLM-assisted functional interpretation in one workflow.

## Installation

Install `BiocManager` and `remotes` first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
```

Install the upstream packages used by `featuremap`:

```r
remotes::install_github("tidymass/masstools")
remotes::install_github("tidymass/massdataset")
remotes::install_github("tidymass/metid")
remotes::install_github("tidymass/pcsf2featuremap")
```

Then install `featuremap` from GitHub:

```r
remotes::install_github(
  "tidymass/featuremap",
  dependencies = TRUE,
  repos = BiocManager::repositories(),
  upgrade = "ask"
)

library(featuremap)
```

Check the installed version:

```r
packageVersion("featuremap")
```

## Workflow

The main analysis follows four functions:

```r
res <- annotate_feature(...)
res <- extract_subnetwork(res = res, ...)
res <- identify_module(res = res, ...)
res <- annotate_module(res = res, ...)
```

The final annotation step can optionally use PubMed and an LLM provider. API
keys are not required for the first three functions.

## Tutorial

For input formats, complete examples, output interpretation, and the Shiny
workflow, see the [FeatureMAP tutorial](https://www.tidymass.org/featuremap_tutorial).

## Getting help

Open the package help page in R with:

```r
?featuremap
```

Issues and feature requests can be submitted through the
[FeatureMAP GitHub repository](https://github.com/tidymass/featuremap/issues).
