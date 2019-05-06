test_that("Test wcomp test function", {
  cat(paste0('\n\r'))
  
  
  if(!exists('DO_MAIN_TESTING'))
    skip('DO_MAIN_TESTING not defined, skipping')
  if(!DO_MAIN_TESTING   )
    skip('DO_MAIN_TESTING is false, skipping')
  
  set.seed(1)
  
  ###************************************************
  #generate data:
  ###************************************************
  
  data = wcomp.generate_example_dataset(m1 = 100,
                                        n_X = 50,
                                        n_Y = 50,
                                        signal_strength_as_change_in_microbial_load = 0.1)
  
  result.selected.references = wcomp.select_references(X = data$counts,
                                                       median_SD_threshold = 0.6, 
                                                       verbose = F)
  
  
  q_BH = q_DSFDR = 0.1
  
  
  
  ###************************************************
  #check inputs
  ###************************************************
  expect_error(wcomp.test(X = data$counts+0.5,
                          y = data$group_labels,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = q_DSFDR),info = "counts are not integers")
  
  expect_error(wcomp.test(X = NA,
                          y = data$group_labels,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = q_DSFDR),info = "counts are not matrix")
  
  expect_error(wcomp.test(X = data$counts,
                          y = data$group_labels[-1],
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = q_DSFDR),info = "labels not same length as counts")
  
  expect_error(wcomp.test(X = data$counts,
                          y = data$group_labels,
                          ind_reference_taxa = -1,verbose = F,q = q_DSFDR),info = "reference taxa must be subset of 1:ncol(X)")
  
  expect_error(wcomp.test(X = data$counts,
                          y = data$group_labels,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = NA,q = q_DSFDR),info = "verbose must be logical")
  
  expect_error(wcomp.test(X = data$counts,
                          y = data$group_labels,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = -1),info = "q must be between 0 and 1")
  
  expect_warning(wcomp.test(X = data$counts,
                          y = data$group_labels,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = 0.5),info = "abnormal value of q not detected")
  
  expect_error(wcomp.test(X = data$counts,
                            y = data$group_labels,
                            ind_reference_taxa = result.selected.references$selected_references,verbose = F,nr_perm = 50,q = q_DSFDR),info = "low number of permutations - nr_perm")
  
  expect_error(wcomp.test(X = data$counts,
                          y = data$group_labels,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,nr_perms_reference_validation = 50,q = q_DSFDR),info = "low number of permutations - nr_perms_reference_validation")
  
  expect_error(wcomp.test(X = data$counts,
                          y = data$group_labels,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,test =  'FOO_TEST',q = q_DSFDR),info = "invalid test")
  
  #check test on pairs
  #check two group test has two groups
  #check disable dsfdr
  
  
  ###************************************************
  #check returned class
  ###************************************************
  set.seed(1)
  result.test.with.class = wcomp.test(X = data$counts,
                           y = data$group_labels,
                           ind_reference_taxa = result.selected.references,verbose = F,q = q_DSFDR)
  set.seed(1)
  result.test = wcomp.test(X = data$counts,
                           y = data$group_labels,
                           ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = q_DSFDR) # can also use for example , test = 'TwoPartWilcoxon', show example
  #check results identical
  
  ###************************************************
  #check returned fields
  ###************************************************
  
  #regression test
})
