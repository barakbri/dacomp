test_that("Regression tests to paper, power and FDR", {
  
  cat(paste0('\n\r'))
  
  if(!exists('DO_REGRESSION_TESTING'))
    skip('DO_REGRESSION_TESTING not defined, skipping')
  if(!DO_REGRESSION_TESTING)
    skip('DO_REGRESSION_TESTING is false, skipping')
  
  library(dacomp)
  
  ###************************************************
  # Function for estimating FDR and TP rate for a single scenario
  ###************************************************
  
  SCENARIO_run_backwards_compatability = function(scenario_ID,nr.reps = 10,q_BH = 0.1 , q_DSFDR = 0.1,verbose = T,CompositionalAnalysis_CodeBase_Location = 'E:/MCB3/CompositionalAnalysis_CodeBase/'){
    
    library(dacomp)
    mean_TP_BH = 0
    mean_FDR_BH = 0
    mean_TP_DSFDR = 0
    mean_FDR_DSFDR = 0
    for(b in 1:nr.reps){
      if(verbose){
        cat(paste0('Running scenario ID ',scenario_ID,', rep: ',b,' / ',nr.reps,'\n\r'))
      }
      data = REFSIM_generate_setting_wrapper(REFSIM_SETTINGS_LIST[[scenario_ID]])
      result.selected.references = dacomp.select_references(X = data$X,
                                                           median_SD_threshold = 1.3, 
                                                           verbose = F)
      
      
      DIFF_ABUNDANT_TAXA_IN_REFERENCE_SET = sum(result.selected.references$selected_references %in% data$select_diff_abundant)
      cat(paste0('DIFF_ABUNDANT_TAXA_IN_REFERENCE_SET: ',DIFF_ABUNDANT_TAXA_IN_REFERENCE_SET,'\n\r'))
      
      
      
      result.test = dacomp.test(X = data$X,
                               y = data$Y,
                               test = DACOMP.TEST.NAME.WILCOXON,
                               ind_reference_taxa = result.selected.references$selected_references,
                               verbose = F,q = q_DSFDR,nr_perm = 10000)
      
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
      
      mean_TP_BH = mean_TP_BH + TP
      mean_FDR_BH = mean_FDR_BH + FDR
      mean_TP_DSFDR = mean_TP_DSFDR + TP_DSFDR
      mean_FDR_DSFDR = mean_FDR_DSFDR + FDR_DSFDR
    }
    return(list(mean_TP_BH = mean_TP_BH/nr.reps,
                mean_FDR_BH = mean_FDR_BH/nr.reps,
                mean_TP_DSFDR = mean_TP_DSFDR/nr.reps,
                mean_FDR_DSFDR = mean_FDR_DSFDR/nr.reps))
  }
  
  ###************************************************
  # Run scenarios, print results and compare to hash value
  ###************************************************
  
  start_scenario = 1
  end_scenario = 25
  nr.reps = 1

  set.seed(1)
  res_list = list()
  for(s in start_scenario:end_scenario){
    cat(paste0('Regression test scenario: ',s,'\n\r'))
    res_list[[s]] = SCENARIO_run_backwards_compatability(s,nr.reps = nr.reps,q_BH = 0.1,q_DSFDR = 0.1,CompositionalAnalysis_CodeBase_Location = CompositionalAnalysis_CodeBase_Location)
  }
  result_table = unlist(res_list[[1]])
  for(i in 2:length(res_list)){
    result_table = rbind(result_table,unlist(res_list[[i]]))
  }
  
  print(result_table)
  
  dacomp:::compare_to_gold_standard(check_name = "Regression_and_Power_VAL_result_table",obj_to_hash = result_table)
})
