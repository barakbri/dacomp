# Package dacomp
This is the GitHub page for the `dacomp` package. The `dacomp` package implements non parametric tests for differential abundance in microbiome counts data. 

## Installation
For lastest development version run:

```r
install.packages("devtools")
devtools::install_github('barakbri/dacomp')
```

## Where do I start?
The package vignette breifly goes over the model, assumptions and background. There are also code snippets showing how to analyze data for several study designs: 2- sample groups, $K$-sample groups, continuous phenotypes, paired study designs and multiple phenotypes.

```r
vignette('dacomp_main_vignette')
```


## Where can I read more?

The full paper is found on arXiv:
*Brill, Barak, Amnon Amir, and Ruth Heller. 2019. “Testing for Differential Abundance in Compositional Counts Data, with Application to Microbiome Studies.” arXiv Preprint arXiv:1904.08937.*

The introduction in the package vignette cites other highly relevant works.

## Reproducing paper results
A github repository for reproducing paper results paper is found [here](https://github.com/barakbri/CompositionalAnalysis_CodeBase)

## A note on version 1.1 of the package:
Version 1.1 of the package adds an addtional variant of the DACOMP method called "normalization by ratio". This variant of the test procedure is disabled by default. See the above paper for additional details on this variant of the `DACOMP` method.