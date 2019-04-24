---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# GSEArot

Gene Set Enrichment Analysis (GSEA) with rotation testing.
This GSEA is similar to the original GSEA proposed by Subramanian et al. (2005).
The difference lies in the fact that rotation testing is replacing the permutation testing.
As shown by Dorum et al. (2009) rotation testing has higher power than permutation testing
in case of small sample sizes.

# References

Dorum, G., Snipen, L., Solheim, M. and Sabo (2009). Rotation Testing in Gene Set Enrichment
Analysis for Small Direct Comparison Experiments. \emph{Statistical Applications in Genetics
and Molecular Biology}, \bold{8}(1), article 34.

Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S.
and Mesirov, J. P (2005) Gene set enrichment analysis: A knowledge-based
approach for interpreting genome-wide expression pro?les, \emph{PNAS}, \bold{102},
15545-15550.

# Installation

```r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")
devtools::install_github("gdorum/relMix")
```
