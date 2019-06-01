# Package dacomp
This is the GitHub page for the `dacomp` package. The `dacomp` package implements non parametric differential abundance tests for microbiome compositional counts data. 

## Installation
For lastest development version run:

```r
install.packages("devtools")
devtools::install_github('barakbri/dacomp', build_vignettes = TRUE)
```

## Where do I start?
The package vignette breifly goes over the model, assumptions and background. There are also code snippets showing how to analyze data for several study designs: 2- sample groups, $K$-sample groups, continuous phenotypes, paired study designs and multiple phenotypes.

```r
vignette('dacomp_main_vignette')
```


## Where can I read more?

The full paper is found here:
*Brill, Barak, Amnon Amir, and Ruth Heller. 2019. “Testing for Differential Abundance in Compositional Counts Data, with Application to Microbiome Studies.” arXiv Preprint arXiv:1904.08937.*

In addition, the following papers are highly relevant:

*Gloor, Gregory B, Jean M Macklaim, Vera Pawlowsky-Glahn, and Juan J Egozcue. 2017. “Microbiome Datasets Are Compositional: And This Is Not Optional.” Frontiers in Microbiology 8. Frontiers: 2224.*

*Kumar, M Senthil, Eric V Slud, Kwame Okrah, Stephanie C Hicks, Sridhar Hannenhalli, and Hector Corrada Bravo. 2018. “Analysis and Correction of Compositional Bias in Sparse Sequencing Count Data.” BMC Genomics 19 (1). BioMed Central: 799.*

*Mandal, Siddhartha, Will Van Treuren, Richard A White, Merete Eggesbø, Rob Knight, and Shyamal D Peddada. 2015. “Analysis of Composition of Microbiomes: A Novel Method for Studying Microbial Composition.” Microbial Ecology in Health and Disease 26 (1). Taylor & Francis: 27663.*
