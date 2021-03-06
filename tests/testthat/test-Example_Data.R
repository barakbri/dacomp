test_that("Data generation", {
  cat(paste0('\n\r'))
  
  
  if(!exists('DO_EXAMPLE_DATA'))
    skip('DO_EXAMPLE_DATA not defined, skipping')
  if(!DO_EXAMPLE_DATA  )
    skip('DO_EXAMPLE_DATA is false, skipping')
  
  ###************************************************
  #check inputs
  ###************************************************
  
  #expect error if m1, n_X,n_Y,signal_strength illegal
  expect_error({dacomp.generate_example_dataset.two_sample(m1 = -1,
                                               n_X = 50,
                                               n_Y = 50,
                                               signal_strength_as_change_in_microbial_load = 0.1)},info = "m1 = -1")
  expect_error({dacomp.generate_example_dataset.two_sample(m1 = NA,
                                               n_X = 50,
                                               n_Y = 50,
                                               signal_strength_as_change_in_microbial_load = 0.1)},info = "m1 = NA")
  expect_error({dacomp.generate_example_dataset.two_sample(m1 = NULL,
                                               n_X = 50,
                                               n_Y = 50,
                                               signal_strength_as_change_in_microbial_load = 0.1)},info = "m1 = NULL")
  expect_error({dacomp.generate_example_dataset.two_sample(m1 = 0,
                                               n_X = 0,
                                               n_Y = 50,
                                               signal_strength_as_change_in_microbial_load = 0.1)},info = "m1 = 0")
  
  expect_error({dacomp.generate_example_dataset.two_sample(m1 = 100,
                                               n_X = 0,
                                               n_Y = 50,
                                               signal_strength_as_change_in_microbial_load = 0.1)},info = "n_X = 0")
  
  expect_error({dacomp.generate_example_dataset.two_sample(m1 = 100,
                                               n_X = 4,
                                               n_Y = 50,
                                               signal_strength_as_change_in_microbial_load = 0.1)},info = "n_X = 4,small n_X")
  
  expect_error({dacomp.generate_example_dataset.two_sample(m1 = 100,
                                               n_X = 50,
                                               n_Y = 0,
                                               signal_strength_as_change_in_microbial_load = 0.1)},info = "n_Y = 0")
  
  expect_error({dacomp.generate_example_dataset.two_sample(m1 = 100,
                                               n_X = 50,
                                               n_Y = 4,
                                               signal_strength_as_change_in_microbial_load = 0.1)},info = "n_Y = 4,small n_X")
  
  expect_error({dacomp.generate_example_dataset.two_sample(m1 = 100,
                                               n_X = 50,
                                               n_Y = 50,
                                               signal_strength_as_change_in_microbial_load = -0.01)},info = "signal_strength_as_change_in_microbial_load < 0")
  
  expect_error({dacomp.generate_example_dataset.two_sample(m1 = 100,
                                               n_X = 50,
                                               n_Y = 50,
                                               signal_strength_as_change_in_microbial_load = 0.751)},info = "signal_strength_as_change_in_microbial_load > 0.75")
  
  # for continous:
  #expect error if m1, n_X,n_Y,signal_strength illegal
  expect_error({dacomp.generate_example_dataset_continuous(m1 = -1,
                                                n = 100,
                                                signal_strength_as_change_in_microbial_load = 0.1)},info = "m1 = -1")
  expect_error({dacomp.generate_example_dataset_continuous(m1 = NA,
                                                           n = 100,
                                                signal_strength_as_change_in_microbial_load = 0.1)},info = "m1 = NA")
  expect_error({dacomp.generate_example_dataset_continuous(m1 = NULL,
                                                           n = 100,
                                                signal_strength_as_change_in_microbial_load = 0.1)},info = "m1 = NULL")
  expect_error({dacomp.generate_example_dataset_continuous(m1 = 0,
                                                           n = 100,
                                                signal_strength_as_change_in_microbial_load = 0.1)},info = "m1 = 0")
  
  expect_error({dacomp.generate_example_dataset_continuous(m1 = 100,
                                                           n = 0,
                                                signal_strength_as_change_in_microbial_load = 0.1)},info = "n = 0")
  
  expect_error({dacomp.generate_example_dataset_continuous(m1 = 100,
                                                           n = 3,
                                                signal_strength_as_change_in_microbial_load = 0.1)},info = "small n")

  expect_error({dacomp.generate_example_dataset_continuous(m1 = 100,
                                                           n = 100,
                                                signal_strength_as_change_in_microbial_load = -0.01)},info = "signal_strength_as_change_in_microbial_load < 0")
  
  expect_error({dacomp.generate_example_dataset_continuous(m1 = 100,
                                                           n = 100,
                                                signal_strength_as_change_in_microbial_load = 0.51)},info = "signal_strength_as_change_in_microbial_load > 0.5")
  
  
  ###************************************************    
  #check returned class
  ###************************************************
  
  data = dacomp.generate_example_dataset.two_sample(m1 = 100,
                                        n_X = 50,
                                        n_Y = 50,
                                        signal_strength_as_change_in_microbial_load = 0.1);
  expect_is(data,class = class(list()))
  
  data_2 = dacomp.generate_example_dataset_continuous(m1 = 100,
                                         n = 100,
                                         signal_strength_as_change_in_microbial_load = 0.1);
  expect_is(data_2,class = class(list()))
  
  ###************************************************
  #check returned fields
  ###************************************************
  
  expect_equal(names(data),c("counts","group_labels","select_diff_abundant", "taxonomy"))
  expect_equal(names(data_2),c("counts" ,"covariate" ,"select_diff_abundant", "taxonomy"))
  
  #check dimensions of returned object
  expect_equal(dim(data$counts),c(100,1384))
  expect_equal(dim(data_2$counts),c(100,1384))
  
  expect_equivalent(table(data$group_labels),table(c(rep(0,50),rep(1,50))))
  expect(length(data_2$covariate)== 100,failure_message = "length of continuous covariate is not sample size")
  #check number of m1
  expect(length(data$select_diff_abundant),100)
  expect(length(data_2$select_diff_abundant),100)
  
})

