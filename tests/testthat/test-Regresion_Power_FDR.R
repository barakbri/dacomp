test_that("Regression tests to paper, power and FDR", {
  cat(paste0('\n\r'))
  
  if(!exists('DO_REGRESSION_TESTING'))
    skip('DO_REGRESSION_TESTING not defined, skipping')
  if(!DO_REGRESSION_TESTING)
    skip('DO_REGRESSION_TESTING is false, skipping')
  
  
  ###************************************************
  # Function for estimating FDR and TP rate for a single scenario
  ###************************************************
  
  SCENARIO_run_backwards_compatability = function(scenario_ID,nr.reps = 10,q_BH = 0.1 , q_DSFDR = 0.1,verbose = T){
    
    library(wcomp)
    
    set.seed(1)
    
    original_wd = getwd()
    setwd(CompositionalAnalysis_CodeBase_Location)
    source(paste0('MCB_Simulation/REFSIM_GenerateSettings_Index.R'))
    setwd(original_wd)
    
    mean_TP_BH = 0
    mean_FDR_BH = 0
    mean_TP_DSFDR = 0
    mean_FDR_DSFDR = 0
    for(b in 1:nr.reps){
      if(verbose){
        cat(paste0('Running scenario ID ',scenario_ID,', rep: ',b,' / ',nr.reps,'\n\r'))
      }
      data = REFSIM_generate_setting_wrapper(REFSIM_SETTINGS_LIST[[scenario_ID]])
      result.selected.references = wcomp.select_references(X = data$X,
                                                           median_SD_threshold = 1.3, 
                                                           verbose = T)
      
      
      
      
      
      
      DIFF_ABUNDANT_TAXA_IN_REFERENCE_SET = sum(result.selected.references$selected_references %in% data$select_diff_abundant)
      cat(paste0('DIFF_ABUNDANT_TAXA_IN_REFERENCE_SET: ',DIFF_ABUNDANT_TAXA_IN_REFERENCE_SET))
      
      
      
      result.test = wcomp.test(X = data$X,
                               y = data$Y,
                               ind_reference_taxa = result.selected.references$selected_references,
                               verbose = T,q = q_DSFDR,nr_perm = 10000,
                               return_results_also_on_reference_validation_fail = T)
      
      rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)
      rejected_DSFDR = result.test$rejected
      
      
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
  
  library(doRNG)
  library(doParallel)
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  start_scenario = 1
  end_scenario = 25 #25 #DEBUG
  nr.reps = 1 #10 #DEBUG
  scenario_results = foreach(s=start_scenario:end_scenario, .options.RNG=1234) %dorng% {
    SCENARIO_run_backwards_compatability(s,nr.reps = nr.reps,q_BH = 0.1,q_DSFDR = 0.1)
  }
  stopCluster(cl)
  result_table = unlist(scenario_results[[1]])
  for(i in 2:length(scenario_results)){
    result_table = rbind(result_table,unlist(scenario_results[[i]]))
  }  
  print(result_table)
  library(digest)
  hash_computation_result = digest(result_table, algo="md5")
  cat(paste0('Current MD5 of sum results: ',hash_computation_result,'\n\r'))
  hash_gold_standard = "c8bd697a35f4b892d10e779c46e21779"
  expect_equal(hash_computation_result,hash_gold_standard)
})
