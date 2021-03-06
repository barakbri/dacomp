# Package dacomp
This is the GitHub page for the `dacomp` package. The `dacomp` package implements non parametric tests for differential abundance in microbiome counts data. 

## Installation
For lastest development version run:

```r
install.packages("devtools")
devtools::install_github("barakbri/dacomp", build_opts = c("--no-resave-data", "--no-manual"))
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


## A note on version 1.2 of the package:
Version 1.2 of the package adds multiplicity adjusted P-values for the DS-FDR procedure. Adjusted P-values are returned from the `dacomp.test(...)` function, under the field `p.values.test.adjusted`. Users using the *normalization by rarefaction* option, (which is not the default), will find the DS-FDR adjusted P-values under the field `p.values.test.adjusted.ratio.normalization`.


## A note on version 1.24 of the package:
Version 1.24 improves the computation time of the Kruskal Wallis test statistic (test statistics should be identical up to machine error), and allows to compute the DSFDR multiplicity adjustment for only a part of the taxa,i.e. specify which taxa are tested, see documentation for `dacomp.test(...)` for additional details.

## A note on version 1.25 of the package:
the BioConductor package 'phyloseq' is now a suggested package, rather than a required package. It is required only for running the functions that generate a mock datasets in the example scripts.


## A note on version 1.26 of the package:
- The function `dacomp.test(...)` now returns effect sizes.
- The function `dacomp.select_references(...)` now has a default value of `0` for the parameter       `median_SD_threshold` meaning the reference taxa will be selected be default to be the lowest number such that the sequencing depth is at least `minTA` reads.
- The function  `dacomp.select_references(...)` can now run in parallel mode.
- The function  `dacomp.select_references(...)` now returns for each value cutoff value for reference selection the sample id for which the minimal number of reads was obtain in reference taxa. This can be used to determine if a specific sample is limiting the rarefaction depth, due to a small number of reads under the reference taxa.
- The function `dacomp.check_reference_set_is_valid.k_groups` has been removed. See next item for the recommended procedure for assessing the validity of the reference set.
- Added `dacomp.validate_references` a function for validating the reference set of taxa: the function checks if the reference taxa contain signal. If a signal is identified, taxa are removed from the reference set based on the order of their reference selection scores. The procedure continues iteratively, until no signal is found in the reference, or  The function was added to the documentation for reference selection and testing functions, and also to the vignette.