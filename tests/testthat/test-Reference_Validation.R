test_that("Reference Validation", {
  cat(paste0('\n\r'))
  
  if(!exists('DO_REFERENCE_VALIDATION'))
    skip('DO_REFERENCE_VALIDATION not defined, skipping')
  if(!DO_REFERENCE_VALIDATION)
    skip('DO_REFERENCE_VALIDATION is false, skipping')
  
  cat(paste0('\n\r'))

  ###************************************************
  # Run the example code
  ###************************************************
  
  set.seed(1)
  library(dacomp)
  #generate data with two study groups
  data = dacomp.generate_example_dataset.two_sample(n_X = 30,n_Y = 30,m1 = 50,signal_strength_as_change_in_microbial_load = 0.1)
  
  # select references. We purposely select reference taxa so that differentially abundant taxa enter the reference set. In general, select using median_SD_threshold=0, minimal_TA=100, see paper for discussion of this selection strategy
  result.selected.references = dacomp.select_references(X = data$counts,
                                                        median_SD_threshold = 1.3,
                                                        maximal_TA = 1000,
                                                        verbose = T)
  
  # some differentially abundant taxa entered the reference set:                                                    
  sum(result.selected.references$selected_references %in% data$select_diff_abundant)
  
  #run the sensitivity analysis.
  cleaned_references = dacomp.validate_references(X =  data$counts,
                                                  Y =  data$group_labels,
                                                  ref_obj = result.selected.references,
                                                  test =DACOMP.TEST.NAME.WILCOXON,
                                                  Q_validation = 0.1,
                                                  Minimal_Counts_in_ref_threshold = 50,
                                                  Reduction_Factor = 0.5,
                                                  Verbose = T,
                                                  disable_DSFDR = T,
                                                  NR_perm = 1000)
  
  
  
  ###************************************************
  #check inputs
  ###************************************************
  #now the reduced reference has no differentially abundant taxa inside....
  expect_equal(sum(cleaned_references %in% data$select_diff_abundant),0) 
  
  
  ###************************************************
  #check inputs
  ###************************************************
  
  ###************************************************
  # regression check on results
  ###************************************************
  dacomp:::compare_to_gold_standard(check_name = "Reference_Validation_VAL_result_ref_validity",obj_to_hash = cleaned_references)
  
})
