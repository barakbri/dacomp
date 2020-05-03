test_that("Test dacomp test function", {
  cat(paste0('\n\r'))
  
  
  if(!exists('DO_MAIN_TESTING'))
    skip('DO_MAIN_TESTING not defined, skipping')
  if(!DO_MAIN_TESTING   )
    skip('DO_MAIN_TESTING is false, skipping')
  
  set.seed(1)
  
  ###************************************************
  # generate data:
  ###************************************************
  
  data = dacomp.generate_example_dataset.two_sample(m1 = 100,
                                        n_X = 50,
                                        n_Y = 50,
                                        signal_strength_as_change_in_microbial_load = 0.1)
  
  result.selected.references = dacomp.select_references(X = data$counts,
                                                       median_SD_threshold = 0.6, 
                                                       verbose = F)
  
  
  q_BH = q_DSFDR = 0.1
  
  
  
  ###************************************************
  # check inputs
  ###************************************************
  expect_error(dacomp.test(X = data$counts+0.5,
                          y = data$group_labels,
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = q_DSFDR),info = "counts are not integers")
  
  expect_error(dacomp.test(X = NA,
                          y = data$group_labels,
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = q_DSFDR),info = "counts are not matrix")
  
  expect_error(dacomp.test(X = data$counts,
                          y = data$group_labels[-1],
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = q_DSFDR),info = "labels not same length as counts")
  
  expect_error(dacomp.test(X = data$counts,
                          y = data$group_labels,
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = -1,verbose = F,q = q_DSFDR),info = "reference taxa must be subset of 1:ncol(X)")
  
  expect_error(dacomp.test(X = data$counts,
                          y = data$group_labels,
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = NA,q = q_DSFDR),info = "verbose must be logical")
  
  expect_error(dacomp.test(X = data$counts,
                          y = data$group_labels,
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = -1),info = "q must be between 0 and 1")
  
  expect_warning(dacomp.test(X = data$counts,
                          y = data$group_labels,
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = 0.5),info = "abnormal value of q not detected")
  
  expect_error(dacomp.test(X = data$counts,
                            y = data$group_labels,
                           test = DACOMP.TEST.NAME.WILCOXON,
                            ind_reference_taxa = result.selected.references$selected_references,verbose = F,nr_perm = 50,q = q_DSFDR),info = "low number of permutations - nr_perm")
  
  expect_error(dacomp.test(X = data$counts,
                          y = data$group_labels,
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,nr_perms_reference_validation = 50,q = q_DSFDR),info = "low number of permutations - nr_perms_reference_validation")
  
  expect_error(dacomp.test(X = data$counts,
                          y = data$group_labels,
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,test =  'FOO_TEST',q = q_DSFDR),info = "invalid test")
  
  expect_error(dacomp.test(X = data$counts[-1,],
                          y = NULL,
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,test =  'SignedWilcoxon',q = q_DSFDR),info = "odd number of rows for paired test")
  
  expect_error(dacomp.test(X = data$counts,
                          y = data$group_labels,
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,test =  'SignedWilcoxon',q = q_DSFDR),info = "Y not null for paired test")
  
  expect_error(dacomp.test(X = data$counts,
                          y = c(data$group_labels[-1],3),
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = q_DSFDR),info = "More than two groups for two group test")
  
  expect_error(dacomp.test(X = data$counts,
                          y = c(data$group_labels[-1],NA),
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = q_DSFDR),info = "Missing labels in y (NA)")
  
  expect_error(dacomp.test(X = data$counts,
                          y = c(data$group_labels[-1],NaN),
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = q_DSFDR),info = "Missing labels in y (NaN)")
  
  expect_error(dacomp.test(X = data$counts,
                          y = data$group_labels,
                          test = DACOMP.TEST.NAME.WILCOXON,
                          ind_reference_taxa = result.selected.references$selected_references,verbose = F,disable_DSFDR = 1,q = q_DSFDR),info = "disable_DSFDR is not valid logical")
  
  
  ###************************************************
  # check returned class
  ###************************************************
  set.seed(1)
  result.test.with.class = dacomp.test(X = data$counts,
                           y = data$group_labels,
                           test = DACOMP.TEST.NAME.WILCOXON,
                           ind_reference_taxa = result.selected.references,verbose = F,q = q_DSFDR)
  set.seed(1)
  result.test = dacomp.test(X = data$counts,
                           y = data$group_labels,
                           test = DACOMP.TEST.NAME.WILCOXON,
                           ind_reference_taxa = result.selected.references$selected_references,verbose = F,q = q_DSFDR) # can also use for example , test = 'TwoPartWilcoxon', show example
  
  #check results identical
  expect_identical(result.test.with.class,result.test,info = "dacomp.test results with reference object and vector of indices for references not identical")
  ###************************************************
  # check returned fields
  ###************************************************
  
  expect_identical(names(result.test.with.class),c("lambda","stats_matrix","p.values.test","p.values.test.adjusted","dsfdr_rejected","dsfdr_threshold" ))
  
  expect_identical(sort(which(is.na(result.test.with.class$lambda))),sort(result.selected.references$selected_references),info = "check missing lambda are only the given references")
  
  expect_identical(sort(which(is.na(result.test.with.class$p.values.test))),sort(result.selected.references$selected_references),info = "check missing p.values are only the given references")
  
  ###************************************************
  # regression test
  ###************************************************
  library(digest)
  hash_computation_result = digest::digest(result.test, algo="md5")
  cat(paste0('Current MD5 of sum results: ',hash_computation_result,'\n\r'))
  hash_gold_standard = "09b5ac0dadd6dd77633705b84484e535"
  expect_equal(hash_computation_result,hash_gold_standard)
  
  ###************************************************
  # check other variants work, need to expand this section...
  ###************************************************
  cat('sanity check variants: \n\r')
  cat('log difference in means: \n\r')
  result.test = dacomp.test(X = data$counts,
                            y = data$group_labels,
                            ind_reference_taxa = result.selected.references,
                            test = DACOMP.TEST.NAME.LOG_FOLD_DIFFERENCE_IN_MEANS,
                            verbose = T,q = q_DSFDR)
  cat(paste0('TP: ',sum(result.test$dsfdr_rejected%in%data$select_diff_abundant),'/',length(result.test$dsfdr_rejected),'\n\r'))
  
  cat(' difference in means: \n\r')
  result.test = dacomp.test(X = data$counts,
                            y = data$group_labels,
                            ind_reference_taxa = result.selected.references,
                            test = DACOMP.TEST.NAME.DIFFERENCE_IN_MEANS,
                            verbose = T,q = q_DSFDR)
  cat(paste0('TP: ',sum(result.test$dsfdr_rejected%in%data$select_diff_abundant),'/',length(result.test$dsfdr_rejected),'\n\r'))
  
  cat(' two part wilcoxon: \n\r')
  result.test = dacomp.test(X = data$counts,
                            y = data$group_labels,
                            ind_reference_taxa = result.selected.references,
                            test = DACOMP.TEST.NAME.TWO_PART_WILCOXON,
                            verbose = T,q = q_DSFDR)
  cat(paste0('TP: ',sum(result.test$dsfdr_rejected%in%data$select_diff_abundant),'/',length(result.test$dsfdr_rejected),'\n\r'))
  
  #KW was slow before v 1.24. As of V1.24, KW has C level implementation!
  
  cat(' KW-two groups: \n\r')
  result.test = dacomp.test(X = data$counts,
                            y = data$group_labels,
                            ind_reference_taxa = result.selected.references,
                            test = DACOMP.TEST.NAME.KRUSKAL_WALLIS,
                            verbose = T,q = q_DSFDR)
  cat(paste0('TP: ',sum(result.test$dsfdr_rejected%in%data$select_diff_abundant),'/',length(result.test$dsfdr_rejected),'\n\r'))
  
  cat(' KW-three groups: \n\r')
  result.test = dacomp.test(X = data$counts,
                            y = c(rep('A',33),rep('B',34),rep('C',33)),
                            ind_reference_taxa = result.selected.references,
                            test = DACOMP.TEST.NAME.KRUSKAL_WALLIS,
                            verbose = T,q = q_DSFDR)
  cat(paste0('TP: ',sum(result.test$dsfdr_rejected%in%data$select_diff_abundant),'/',length(result.test$dsfdr_rejected),'\n\r'))  
  
  
  
  ###************************************************
  # test continuous convariate
  ###************************************************
  set.seed(1)
  
  data = dacomp.generate_example_dataset_continuous(n = 100,m1 = 30,signal_strength_as_change_in_microbial_load = 0.1)
  
  
  result.selected.references = dacomp.select_references(X = data$counts,
                                                        median_SD_threshold = 0.6, #APPLICATION SPECIFIC
                                                        verbose = F)
  
  #multiplicity correction levels for the BH and DS-FDR methods
  q_BH = q_DSFDR = 0.1
  
  #Perform testing:
  result.test = dacomp.test(X = data$counts,
                            y = data$covariate,test = DACOMP.TEST.NAME.SPEARMAN,
                            ind_reference_taxa = result.selected.references,
                            verbose = F,q = q_DSFDR)
  
  rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)
  rejected_DSFDR = result.test$dsfdr_rejected
  
  #sum(rejected_BH %in% data$select_diff_abundant)
  #length(rejected_BH)
  
  library(digest)
  hash_computation_result = digest::digest(result.test, algo="md5")
  cat(paste0('Current MD5 of sum results: ',hash_computation_result,'\n\r'))
  hash_gold_standard_continous = "916ced97f0ec353a16b20e341cc67a43"
  expect_equal(hash_computation_result,hash_gold_standard_continous)
  
  ###************************************************
  # test adjusted P-values for DS-FDR
  ###************************************************
  set.seed(1)
  
  result.test = dacomp.test(X = data$counts,
                            y = data$covariate,test = DACOMP.TEST.NAME.SPEARMAN,
                            ind_reference_taxa = result.selected.references,
                            verbose = F,q = q_DSFDR,compute_ratio_normalization = T)
  
  expect_equal((result.test$p.values.test<=result.test$dsfdr_threshold),
               (result.test$p.values.test.adjusted<=q_DSFDR),
               label = 'DS-FDR adjusted P-values not equal to indicator for rejection (for regular tests)')
  
  expect_equal((result.test$p.values.test.ratio.normalization <= result.test$dsfdr_threshold_ratio_normalization),
               (result.test$p.values.test.adjusted.ratio.normalization <= q_DSFDR),
               label = 'DS-FDR adjusted P-values not equal to indicator for rejection (for regular tests)')
  
})
