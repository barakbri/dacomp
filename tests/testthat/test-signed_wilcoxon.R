test_that("Test Signed Wilcoxon", {
  cat(paste0('\n\r'))
  
  if(!exists('DO_TEST_SIGNED'))
    skip('DO_TEST_SIGNED not defined, skipping')
  if(!DO_TEST_SIGNED)
    skip('DO_TEST_SIGNED is false, skipping')
  
  library(dacomp)
  
  set.seed(1)
  data = dacomp.generate_example_dataset_paired(30)
  
  
  result.selected.references = dacomp.select_references(X = data$counts,
                                                        median_SD_threshold = 0.6, #APPLICATION SPECIFIC
                                                        verbose = F)
  
  #multiplicity correction levels for the BH and DS-FDR methods
  q_BH = q_DSFDR = 0.1
  
  #Perform testing:
  result.test = dacomp.test(X = data$counts,
                            y = NULL,
                            ind_reference_taxa = result.selected.references, test = DACOMP.TEST.NAME.WILCOXON_SIGNED_RANK_TEST,
                            verbose = F,q = q_DSFDR)
  
  # rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)
  # rejected_DSFDR = result.test$dsfdr_rejected
  # length(rejected_BH)
  # length(rejected_DSFDR)
  # sum(rejected_BH %in% data$select_diff_abundant)
  # sum(rejected_DSFDR %in% data$select_diff_abundant)
  
  #Regression test
  library(digest)
  hash_computation_result = digest::digest(result.test, algo="md5")
  cat(paste0('Current MD5 of sum results: ',hash_computation_result,'\n\r'))
  hash_gold_standard = "dfe76edb4ef768e2ae62d47c12b1f422"
  expect_equal(hash_computation_result,hash_gold_standard,label = 'Regression test failed on Wilcoxon signed test')
})


