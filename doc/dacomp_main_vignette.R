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

## ----eval = T--------------------------------------------------------------
dacomp.plot_reference_scores(result.selected.references)

## ---- eval = T-------------------------------------------------------------
#multiplicity correction levels for the BH and DS-FDR methods
q_BH = q_DSFDR = 0.1

#Perform testing:
result.test = dacomp.test(X = data$counts, #counts data
                     y = data$group_labels, #phenotype in y argument
                     # obtained from dacomp.select_references(...):
                     ind_reference_taxa = result.selected.references, 
                     test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                     verbose = F,q = q_DSFDR) #multiplicity adjustment level

#These are the indices of taxa discoverted as differentially abundant:
# by applying a BH multiplicity adjustment on the P-values:
rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH) 
#by applying a DS-FDR multiplicity adjustment on the P-values:
rejected_DSFDR = result.test$dsfdr_rejected 

## ----eval = T--------------------------------------------------------------
print(result.test)

## ---- eval = F-------------------------------------------------------------
#  result.test = dacomp.test(X = data$counts,
#                       y = data$group_labels,
#                       ind_reference_taxa = result.selected.references,
#                       test = DACOMP.TEST.NAME.LOG_FOLD_DIFFERENCE_IN_MEANS,
#                       verbose = T,q = q_DSFDR)
#  
#  result.test = dacomp.test(X = data$counts,
#                       y = data$group_labels,
#                       ind_reference_taxa = result.selected.references,
#                       test = DACOMP.TEST.NAME.DIFFERENCE_IN_MEANS,
#                       verbose = T,q = q_DSFDR)
#  
#  result.test = dacomp.test(X = data$counts,
#                       y = data$group_labels,
#                       ind_reference_taxa = result.selected.references,
#                       test = DACOMP.TEST.NAME.TWO_PART_WILCOXON,
#                       verbose = T,q = q_DSFDR)

## ---- eval = F-------------------------------------------------------------
#  set.seed(1)
#  # Sample data:
#  # 30 is the number of samples, so we will have 60 rows.
#  # By default, 30 OTUs are differentially abundant
#  data = dacomp.generate_example_dataset_paired(30)
#  
#  # data$counts is matrix of counts:
#  # first 30 rows correspond to samples 1:30 under condition 1
#  # rows 31:60 correspond to samples 1:30 under condition 2
#  
#  #select references:
#  result.selected.references = dacomp.select_references(
#                                  X = data$counts,
#                                  median_SD_threshold = 0.6, #APPLICATION SPECIFIC
#                                  verbose = T)
#  
#  
#  length(result.selected.references$selected_references)
#  
#  #plot the reference selection scores:
#  #(can also be used to better set the median SD threshold)
#  dacomp.plot_reference_scores(result.selected.references)
#  
#  
#  #multiplicity correction levels for the BH and DS-FDR methods
#  q_BH = q_DSFDR = 0.1
#  
#  #Perform testing:
#  result.test = dacomp.test(
#                      X = data$counts, #counts matrix formated as required
#                      y = NULL, #supply a null phenotype
#                      ind_reference_taxa = result.selected.references,
#                      test = DACOMP.TEST.NAME.WILCOXON_SIGNED_RANK_TEST,
#                      verbose = T,q = q_DSFDR)
#  
#  #discoveries:
#  rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)
#  rejected_DSFDR = result.test$dsfdr_rejected

## ----eval = F--------------------------------------------------------------
#  set.seed(1)
#  data = dacomp.generate_example_dataset_continuous(n = 100,m1 = 30,
#  signal_strength_as_change_in_microbial_load = 0.1)
#  
#  #data$counts - matrix of counts
#  #data$covariate - a vector of 100 phenotype measurements,
#  #corresponding to the rows of X.
#  
#  
#  result.selected.references = dacomp.select_references(
#                            X = data$counts,
#                            median_SD_threshold = 0.6, #APPLICATION SPECIFIC
#                            verbose = T)
#  
#  #number of selected references
#  length(result.selected.references$selected_references)
#  
#  #plot the reference selection scores (can also be used to better set the median SD threshold)
#  dacomp.plot_reference_scores(result.selected.references)
#  
#  #multiplicity correction levels for the BH and DS-FDR methods
#  q_BH = q_DSFDR = 0.1
#  
#  #Perform testing:
#  result.test = dacomp.test(X = data$counts,
#                        y = data$covariate,test = DACOMP.TEST.NAME.SPEARMAN,
#                        ind_reference_taxa = result.selected.references,
#                        verbose = T,q = q_DSFDR)
#  
#  rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)
#  rejected_DSFDR = result.test$dsfdr_rejected
#  

## ----eval = F--------------------------------------------------------------
#  set.seed(1)
#  #generate data, with a multivariate phenotype:
#  
#  data = dacomp.generate_example_dataset_multivariate_example(
#    n = 100,
#    m1 = 30,
#    signal_strength_as_change_in_microbial_load = 0.1)
#  
#  #phenotype of dimensionality two, for each subject
#  head(data$covariate)
#  #            u1         u2
#  #[1,] 0.4820801 0.57487220
#  #[2,] 0.5995658 0.07706438
#  #[3,] 0.4935413 0.03554058
#  #[4,] 0.1862176 0.64279549
#  #[5,] 0.8273733 0.92861520
#  #[6,] 0.6684667 0.59809242
#  
#  #select references:
#  result.selected.references = dacomp.select_references(
#                                X = data$counts,
#                                median_SD_threshold = 0.5,
#                                verbose = T)
#  
#  
#  
#  #multiplicity correction levels for the BH and DS-FDR methods
#  q_BH = q_DSFDR = 0.1
#  
#  # The number of permutations performed for each test.
#  # Note that this number is passed as an argument to the function,
#  # AND must be exactly the number of permutations performed
#  # and returned by the supplied test function
#  nr_perm_to_perform = 1000
#  
#  # We will use the dcov test from package energy to
#  # to test for differential abundance
#  
#  library(energy)
#  
#  #this is the custom test function supplied by the user
#  # Input: array of rarefied reads, of length n
#  # Output: Array of test statistics, with right tailed alternative
#  # of length nr.perm +1. The first entry is the test statistic for the original data.
#  custom_test_function = function(X){
#    # compute test and permutations. Note that the phenotype
#    # is available to the test function
#    res = dcov.test(X, data$covariate, R=nr_perm_to_perform)
#    return(
#            c(
#              # first entry is the test statistic to the data:
#              res$statistic ,
#              # a vector of length nr_perm_to_perform containing
#              # test statistics computed for data with permuted
#              # phenotypes
#              res$replicates
#              )
#      )
#  }
#  
#  #Perform testing:
#  result.test = dacomp.test(X = data$counts,
#                            #note that y is NULL, phenotype is available to the test function:
#                            y = NULL,
#  
#                            # set test to be user defined:
#                            test = DACOMP.TEST.NAME.USER_DEFINED,
#                            ind_reference_taxa = result.selected.references,
#                            verbose = T,q = q_DSFDR,
#  
#                            #nr_perm must be identical to the number of
#                            #permutation returned from test function:
#                            nr_perm = nr_perm_to_perform,
#  
#                            #pass as argument the user defined test function:
#                            user_defined_test_function = custom_test_function)
#  
#  rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)
#  rejected_DSFDR = result.test$dsfdr_rejected
#  
#  

