test_that("Test Plot function", {
  cat(paste0('\n\r'))
  
  
  if(!exists('DO_PLOTTING'))
    skip('DO_PLOTTING not defined, skipping')
  if(!DO_PLOTTING)
    skip('DO_PLOTTING is false, skipping')
  
  
  set.seed(1)
  
  ###************************************************
  #generate data:
  ###************************************************
  
  data = dacomp.generate_example_dataset(m1 = 100,
                                        n_X = 20,
                                        n_Y = 20,
                                        signal_strength_as_change_in_microbial_load = 0.5)
  result.selected.references = dacomp.select_references(X = data$counts,
                                                       median_SD_threshold = 0.6, 
                                                       verbose = F)
  
  
  ###************************************************
  # check bad inputs
  ###************************************************
  
  bad_input = result.selected.references
  class(bad_input) = 'SOMETHING_ELSE'
  expect_error(dacomp.plot_reference_scores(bad_input),info = 'dacomp.plot_reference_scores must except object of type CLASS.LABEL.REFERENCE_SELECTION_OBJECT but works despite bad input')
  
  expect_error(dacomp.plot_reference_scores(result.selected.references,label = 5),info = 'label must be character')
  
  expect_error(dacomp.plot_reference_scores(result.selected.references,label = NA),info = 'label must be character')
  
  expect_error(dacomp.plot_reference_scores(result.selected.references,quantiles_to_plot = NA),info = 'quantiles must be a vector of numbers, between 0 and 1')
  
  expect_error(dacomp.plot_reference_scores(result.selected.references,quantiles_to_plot = c(-5,1)),info = 'quantiles must be a vector of numbers, between 0 and 1')
  
  expect_error(dacomp.plot_reference_scores(result.selected.references,quantiles_to_plot = c(0.5,5)),info = 'quantiles must be a vector of numbers, between 0 and 1')
  
  
})
