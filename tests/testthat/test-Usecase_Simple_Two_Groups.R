test_that("Simple Use Case", {
  cat(paste0('\n\r'))
  
  if(!exists('DO_SIMPLE_USE_CASE'))
    skip('DO_SIMPLE_USE_CASE not defined, skipping')
  if(!DO_SIMPLE_USE_CASE)
    skip('DO_SIMPLE_USE_CASE is false, skipping')
  library(HHG)
  library(dacomp)
  
  set.seed(1)
  
  ###************************************************
  #generate data:
  ###************************************************
  
  data = dacomp.generate_example_dataset.two_sample(m1 = 100,
                                        n_X = 50,
                                        n_Y = 50,
                                        signal_strength_as_change_in_microbial_load = 0.1)
  
  ###************************************************
  #select references: (may take a minute)
  ###************************************************
  result.selected.references = dacomp.select_references(X = data$counts,
                                                       minimal_TA = 50, 
                                                       verbose = T)
  
  length(result.selected.references$selected_references)
  
  ###************************************************
  #plot the reference selection scores (can be used to better set the threshold...)
  ###************************************************
  dacomp.plot_reference_scores(result.selected.references)
  
  DIFF_ABUNDANT_TAXA_IN_REFERENCE_SET = sum(result.selected.references$selected_references %in% data$select_diff_abundant)
  cat(paste0('DIFF_ABUNDANT_TAXA_IN_REFERENCE_SET: ',DIFF_ABUNDANT_TAXA_IN_REFERENCE_SET,'\n\r'))
  
  q_BH = q_DSFDR = 0.1
  
  ###************************************************
  # Run dacomp
  ###************************************************
  result.test = dacomp.test(X = data$counts,
                           y = data$group_labels,
                           test = DACOMP.TEST.NAME.WILCOXON,
                           ind_reference_taxa = result.selected.references$selected_references,verbose = T,q = q_DSFDR) # can also use for example , test = 'TwoPartWilcoxon', show example
  
  rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)
  rejected_DSFDR = result.test$dsfdr_rejected
  
  
  TP = sum((rejected_BH %in% data$select_diff_abundant))
  TP_DSFDR = sum((rejected_DSFDR %in% data$select_diff_abundant))
  
  FDR = ifelse(length(rejected_BH)>0,
               sum(!(rejected_BH %in% data$select_diff_abundant))/length(rejected_BH),
               0)
  FDR_DSFDR = ifelse(length(rejected_DSFDR)>0,
                     sum(!(rejected_DSFDR %in% data$select_diff_abundant))/length(rejected_DSFDR),
                     0)
  
  cat(paste0('True positives: ',TP,', FDR: ',round(FDR,2),'\n\r'))
  
  
  FILTER_FOR_TESTING = rep(F,ncol(data$counts))
  FILTER_FOR_TESTING[data$select_diff_abundant] = T
  result.test.with.taxa.filter = dacomp.test(X = data$counts,
                            y = data$group_labels,
                            test = DACOMP.TEST.NAME.WILCOXON,
                            ind_reference_taxa = result.selected.references$selected_references,verbose = T,q = q_DSFDR,DSFDR_Filter = FILTER_FOR_TESTING)
  
  cat(paste0('True positives with filtering over tested hypotheses: ',length(result.test.with.taxa.filter$dsfdr_rejected),'\n\r'))
})
