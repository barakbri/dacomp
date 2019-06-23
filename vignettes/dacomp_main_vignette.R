## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE-----------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("barakbri/dacomp")

## ----eval = T,cache=FALSE--------------------------------------------------
library(dacomp)

set.seed(1)

data = dacomp.generate_example_dataset.two_sample(m1 = 100,
        n_X = 50,
        n_Y = 50,
        signal_strength_as_change_in_microbial_load = 0.1)

## ----eval = T--------------------------------------------------------------
#select references: (may take a minute)
result.selected.references = dacomp.select_references(
                      X = data$counts,
                      median_SD_threshold = 0.6, #APPLICATION SPECIFIC
                      verbose = F)


## ----eval = T--------------------------------------------------------------
print(result.selected.references)

