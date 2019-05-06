test_that("Reference Validation", {
  cat(paste0('\n\r'))
  
  if(!exists('DO_REFERENCE_VALIDATION'))
    skip('DO_REFERENCE_VALIDATION not defined, skipping')
  if(!DO_REFERENCE_VALIDATION)
    skip('DO_REFERENCE_VALIDATION is false, skipping')
  
  cat(paste0('\n\r'))
  
  if(!exists('DO_SIMPLE_USE_CASE'))
    skip('DO_SIMPLE_USE_CASE not defined, skipping')
  if(!DO_SIMPLE_USE_CASE)
    skip('DO_SIMPLE_USE_CASE is false, skipping')
  
  
  library(wcomp)
  
  ###************************************************
  #Prepare data
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
  #check inputs
  ###************************************************
  expect_error(wcomp.check_reference_set_is_valid.k_groups(X_ref = 1,
                                                           Y = data$group_labels,
                                                           nr.perm = 10000,
                                                           verbose = T),info = "X_ref is not matrix")
  
  expect_error(wcomp.check_reference_set_is_valid.k_groups(X_ref = data$counts[,result.selected.references$selected_references],
                                                           Y = data$group_labels[-1],
                                                           nr.perm = 10000,
                                                           verbose = T),info = "X_ref and Y do not match in size")
  
  expect_error(wcomp.check_reference_set_is_valid.k_groups(X_ref = data$counts[,result.selected.references$selected_references],
                                                           Y = data$group_labels,
                                                           nr.perm = 'A',
                                                           verbose = T),info = "nr.perm is illegal, must be number")
  
  expect_error(wcomp.check_reference_set_is_valid.k_groups(X_ref = data$counts[,result.selected.references$selected_references],
                                                           Y = data$group_labels,
                                                           nr.perm = 50,
                                                           verbose = T),info = "nr.perm is illegal,must be at least 1000")
  
  expect_warning(wcomp.check_reference_set_is_valid.k_groups(X_ref = data$counts[,1,drop=F],
                                                           Y = data$group_labels,
                                                           nr.perm = 1000,
                                                           verbose = T),info = "reference with only one taxon must produce warning")
  
  ###************************************************
  #check returned class
  ###************************************************
  result.ref.validity = wcomp.check_reference_set_is_valid.k_groups(X_ref = data$counts[,result.selected.references$selected_references],Y = data$group_labels,nr.perm = 10000,verbose = T)
  
  expect_is(result.ref.validity,wcomp:::CLASS.LABEL.REFERENCE_VALIDATION_RESULT_OBJECT)
  
  ###************************************************
  #check returned fields
  ###************************************************
  expect_equal(names(result.ref.validity),c("p.value.HHG.L2","p.value.HHG.L1","p.value.HHG.BC","p.value.energy.L2","p.value.energy.L1","p.value.energy.BC","p.value.permanova_L2","p.value.permanova_L1", "p.value.permanova_BC"))
  
  
  ###************************************************
  # regression check on results
  ###************************************************
  
  library(digest)
  hash_computation_result = digest(result.ref.validity, algo="md5")
  cat(paste0('Current MD5 of sum results: ',hash_computation_result,'\n\r'))
  hash_gold_standard = "6a42b8624f524e322ebcdfbaff1f60b2"
  expect_equal(hash_computation_result,hash_gold_standard)
  
})
