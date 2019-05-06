test_that("Data generation", {
  cat(paste0('\n\r'))
  
  
  if(!exists('DO_EXAMPLE_DATA  '))
    skip('DO_EXAMPLE_DATA not defined, skipping')
  if(!DO_EXAMPLE_DATA  )
    skip('DO_EXAMPLE_DATA is false, skipping')
  
  data = wcomp.generate_example_dataset(m1 = 100,
                                 n_X = 50,
                                 n_Y = 50,
                                 signal_strength_as_change_in_microbial_load = 0.1);
  expect_is(data,class = class(list()))
  expect_equal(names(data),c("counts","group_labels","select_diff_abundant"))
  
  #expect_equal(dim(data$counts),())
})

