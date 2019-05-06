test_that("multiplication works", {
  cat(paste0('\n\r'))

  
  if(!exists('DO_REFERENCE_SELECTION'))
    skip('DO_REFERENCE_SELECTION not defined, skipping')
  if(!DO_REFERENCE_SELECTION)
    skip('DO_REFERENCE_SELECTION is false, skipping')
  
  ###************************************************
  #check inputs
  ###************************************************
  
  set.seed(1)
  
  data = wcomp.generate_example_dataset(m1 = 100,
                                        n_X = 50,
                                        n_Y = 50,
                                        signal_strength_as_change_in_microbial_load = 0.1)
  
  expect_error(wcomp.select_references(X = data$counts,
                                       median_SD_threshold = -1, 
                                       verbose = F),info = "negative median SD threshold")
  
  expect_warning(wcomp.select_references(X = data$counts,
                                       median_SD_threshold = 0.4, 
                                       verbose = F),info = "medianSD under 0.5 should produce warning")
  
  expect_warning(wcomp.select_references(X = data$counts,
                                         median_SD_threshold = 1.6, 
                                         verbose = F),info = "medianSD above 1.5 should produce warning")
  
  expect_error(wcomp.select_references(X = data$counts,
                                       median_SD_threshold = 1.0, minimal_TA = -1,
                                       verbose = F),info = "minimal_TA cannot be negative")
  
  expect_error(wcomp.select_references(X = data$counts,
                                       median_SD_threshold = 1.0, minimal_TA = 50,maximal_TA = 10,
                                       verbose = F),info = "minimal_TA must be smaller than maximal TA")
  
  expect_error(wcomp.select_references(X = data$counts,
                                       median_SD_threshold = 1.0, minimal_TA = 10,maximal_TA = 50,Pseudo_Count_used = 0,
                                       verbose = F),info = "Pseudo_Count_used must be strictly greater than zero")
  
  expect_error(wcomp.select_references(X = data$counts,
                                       median_SD_threshold = 1.0, minimal_TA = 10,maximal_TA = 50,Pseudo_Count_used = 1,select_from = (1:(ncol(data$counts)+1)),
                                       verbose = F),info = "select_from must be a subset of 1:ncol(X)")
  
  
  ###************************************************
  #check returned class
  ###************************************************
  set.seed(1)
  
  data = wcomp.generate_example_dataset(m1 = 100,
                                        n_X = 50,
                                        n_Y = 50,
                                        signal_strength_as_change_in_microbial_load = 0.1)
  
  result.selected.references = wcomp.select_references(X = data$counts,
                                                       median_SD_threshold = 0.6, 
                                                       verbose = F)
  
  
  ###************************************************
  #check returned fields
  ###************************************************
  
  
  expect_is(result.selected.references,class = wcomp:::CLASS.LABEL.REFERENCE_SELECTION_OBJECT)
  
  
  expect_equal(names(result.selected.references),c("selected_references","mean_prevalence_over_the_sorted", "min_abundance_over_the_sorted", "ratio_matrix",                 
                             "scores","selected_MinAbundance","median_SD_threshold","minimal_TA", "maximal_TA"))
  
  
  ###************************************************
  # regression check on results
  ###************************************************
  
  library(digest)
  hash_computation_result = digest(result.selected.references, algo="md5")
  cat(paste0('Current MD5 of sum results: ',hash_computation_result,'\n\r'))
  hash_gold_standard = "61ed2ab005ef9c1735b7c722fa2feb41"
  expect_equal(hash_computation_result,hash_gold_standard)
  
})
