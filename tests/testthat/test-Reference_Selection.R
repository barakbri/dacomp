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
  
  data = dacomp.generate_example_dataset.two_sample(m1 = 100,
                                        n_X = 50,
                                        n_Y = 50,
                                        signal_strength_as_change_in_microbial_load = 0.1)
  
  expect_error(dacomp.select_references(X = data$counts,
                                       median_SD_threshold = -1, 
                                       verbose = F),info = "negative median SD threshold")
  
  expect_warning(dacomp.select_references(X = data$counts,
                                       median_SD_threshold = 0.4, 
                                       verbose = F),info = "medianSD under 0.5 should produce warning")
  
  expect_warning(dacomp.select_references(X = data$counts,
                                         median_SD_threshold = 1.6, 
                                         verbose = F),info = "medianSD above 1.5 should produce warning")
  
  expect_error(dacomp.select_references(X = data$counts,
                                       median_SD_threshold = 1.0, minimal_TA = -1,
                                       verbose = F),info = "minimal_TA cannot be negative")
  
  expect_error(dacomp.select_references(X = data$counts,
                                       median_SD_threshold = 1.0, minimal_TA = 50,maximal_TA = 10,
                                       verbose = F),info = "minimal_TA must be smaller than maximal TA")
  
  expect_error(dacomp.select_references(X = data$counts,
                                       median_SD_threshold = 1.0, minimal_TA = 10,maximal_TA = 50,Pseudo_Count_used = 0,
                                       verbose = F),info = "Pseudo_Count_used must be strictly greater than zero")
  
  expect_error(dacomp.select_references(X = data$counts,
                                       median_SD_threshold = 1.0, minimal_TA = 10,maximal_TA = 50,Pseudo_Count_used = 1,select_from = (1:(ncol(data$counts)+1)),
                                       verbose = F),info = "select_from must be a subset of 1:ncol(X)")
  
  
  ###************************************************
  #check returned class
  ###************************************************
  set.seed(1)
  
  data = dacomp.generate_example_dataset.two_sample(m1 = 100,
                                        n_X = 50,
                                        n_Y = 50,
                                        signal_strength_as_change_in_microbial_load = 0.1)
  
  result.selected.references = dacomp.select_references(X = data$counts,
                                                       median_SD_threshold = 0.6, 
                                                       verbose = F)
  
  result.selected.references.different.threshold = dacomp.select_references(X = data$counts,
                                                       median_SD_threshold = 0.5, 
                                                       verbose = F,Previous_Reference_Selection_Object = result.selected.references)
  
  
  ###************************************************
  #check returned fields
  ###************************************************
  
  
  expect_is(result.selected.references,class = dacomp:::CLASS.LABEL.REFERENCE_SELECTION_OBJECT)
  
  
  expect_equal(names(result.selected.references),c("selected_references","mean_prevalence_over_the_sorted", "min_abundance_over_the_sorted", "ratio_matrix",                 
                             "scores","selected_MinAbundance","median_SD_threshold","minimal_TA", "maximal_TA"))
  
  expect(result.selected.references$selected_MinAbundance >= result.selected.references.different.threshold$selected_MinAbundance,failure_message = 'selected_MinAbundance must be monotone non decreasing with threshold.')
  
  expect(!any(!(result.selected.references.different.threshold$selected_references %in% result.selected.references$selected_references)) ,failure_message = 'as medianSD threshold increases, taxa should only join the refernece set, not leave it.')
  
  result.selected.references.same.threshold = dacomp.select_references(X = data$counts,
                                                                           median_SD_threshold = 0.6, 
                                                                           verbose = F,Previous_Reference_Selection_Object = result.selected.references)
  expect_equal(
    digest::digest(result.selected.references, algo="md5"),
    digest::digest(result.selected.references.same.threshold, algo="md5"), info = "Selected references object by using previous computation of rations is not identical to original result"
               )
  
  ###************************************************
  # regression check on results
  ###************************************************
  
  dacomp:::compare_to_gold_standard(check_name = "Reference_Selection_VAL_Selected_References",obj_to_hash = result.selected.references)
  
})
