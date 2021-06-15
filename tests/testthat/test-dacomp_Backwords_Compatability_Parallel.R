
test_that("Test DACOMP paper scenarios with multiple cores", {
  
  if(!exists('DO_BACKWORDS_COMPATABILITY_MULTIPLE_CORES'))
    skip('DO_BACKWORDS_COMPATABILITY_MULTIPLE_CORES not defined, skipping')
  if(!DO_BACKWORDS_COMPATABILITY_MULTIPLE_CORES)
    skip('DO_BACKWORDS_COMPATABILITY_MULTIPLE_CORES is false, skipping')
  
  
  library(dacomp)
  
  set.seed(1)
  
  SCENARIO_run_backwards_compatability = function(scenario_ID,nr.reps = 20,q_BH = 0.1 , q_DSFDR = 0.1,verbose = T){
    CompositionalAnalysis_CodeBase_Location = 'C:/MCB3/CompositionalAnalysis_CodeBase/Scripts/'
    source(paste0(CompositionalAnalysis_CodeBase_Location,'/REFSIM_GenerateSettings_Index.R'),chdir = T)

    library(dacomp)
    mean_TP_BH = 0
    mean_FDR_BH = 0
    mean_TP_DSFDR = 0
    mean_FDR_DSFDR = 0
    
    mean_TP_BH_ratio = 0
    mean_FDR_BH_ratio = 0
    mean_TP_DSFDR_ratio = 0
    mean_FDR_DSFDR_ratio = 0
    
    for(b in 1:nr.reps){
      if(verbose){
        cat(paste0('Running scenario ID ',scenario_ID,', rep: ',b,' / ',nr.reps,'\n\r'))
      }
      data = REFSIM_generate_setting_wrapper(REFSIM_SETTINGS_LIST[[scenario_ID]])
      result.selected.references = dacomp.select_references(X = data$X,
                                                            median_SD_threshold = 1.3, 
                                                            verbose = F)
      
      length(result.selected.references$selected_references)
      
      #plot the reference selection scores (can be used to better set the threshold...)
      dacomp.plot_reference_scores(result.selected.references)
      
      DIFF_ABUNDANT_TAXA_IN_REFERENCE_SET = sum(result.selected.references$selected_references %in% data$select_diff_abundant)
      cat(paste0('DIFF_ABUNDANT_TAXA_IN_REFERENCE_SET: ',DIFF_ABUNDANT_TAXA_IN_REFERENCE_SET))
      
      
      
      result.test = dacomp.test(X = data$X,
                                y = data$Y,
                                ind_reference_taxa = result.selected.references$selected_references,
                                test = DACOMP.TEST.NAME.WILCOXON,
                                verbose = F,q = q_DSFDR,nr_perm = 10000,compute_ratio_normalization = T)
      
      rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)
      rejected_DSFDR = result.test$dsfdr_rejected
      rejected_BH_ratio = which(p.adjust(result.test$p.values.test.ratio.normalization,method = 'BH')<=q_BH)
      rejected_DSFDR_ratio = result.test$dsfdr_rejected_ratio_normalization
      
      
      TP = sum((rejected_BH %in% data$select_diff_abundant))
      TP_DSFDR = sum((rejected_DSFDR %in% data$select_diff_abundant))
      
      TP_ratio = sum((rejected_BH_ratio %in% data$select_diff_abundant))
      TP_DSFDR_ratio = sum((rejected_DSFDR_ratio %in% data$select_diff_abundant))
      
      FDR = ifelse(length(rejected_BH)>0,
                   sum(!(rejected_BH %in% data$select_diff_abundant))/length(rejected_BH),
                   0)
      FDR_DSFDR = ifelse(length(rejected_DSFDR)>0,
                         sum(!(rejected_DSFDR %in% data$select_diff_abundant))/length(rejected_DSFDR),
                         0)
      
      FDR_ratio = ifelse(length(rejected_BH_ratio)>0,
                         sum(!(rejected_BH_ratio %in% data$select_diff_abundant))/length(rejected_BH_ratio),
                         0)
      FDR_DSFDR_ratio = ifelse(length(rejected_DSFDR_ratio)>0,
                               sum(!(rejected_DSFDR_ratio %in% data$select_diff_abundant))/length(rejected_DSFDR_ratio),
                               0)
      
      mean_TP_BH = mean_TP_BH + TP
      mean_FDR_BH = mean_FDR_BH + FDR
      mean_TP_DSFDR = mean_TP_DSFDR + TP_DSFDR
      mean_FDR_DSFDR = mean_FDR_DSFDR + FDR_DSFDR
      
      mean_TP_BH_ratio = mean_TP_BH_ratio + TP_ratio
      mean_FDR_BH_ratio = mean_FDR_BH_ratio + FDR_ratio
      mean_TP_DSFDR_ratio = mean_TP_DSFDR_ratio + TP_DSFDR_ratio
      mean_FDR_DSFDR_ratio = mean_FDR_DSFDR_ratio + FDR_DSFDR_ratio
      
      
      #cat(paste0('True positives: ',TP,', FDR: ',round(FDR,2),'\n\r'))
    }
    return(list(mean_TP_BH = mean_TP_BH/nr.reps,
                mean_FDR_BH = mean_FDR_BH/nr.reps,
                mean_TP_DSFDR = mean_TP_DSFDR/nr.reps,
                mean_FDR_DSFDR = mean_FDR_DSFDR/nr.reps,
                mean_TP_BH_ratio = mean_TP_BH_ratio/nr.reps,
                mean_FDR_BH_ratio = mean_FDR_BH_ratio/nr.reps,
                mean_TP_DSFDR_ratio = mean_TP_DSFDR_ratio/nr.reps,
                mean_FDR_DSFDR_ratio = mean_FDR_DSFDR_ratio/nr.reps))
  }
  
  
  library(doRNG)
  library(doParallel)
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  start_scenario = 1
  end_scenario = 25
  nr.reps = 1
  scenario_results = foreach(s=start_scenario:end_scenario, .options.RNG=1234) %dorng% { SCENARIO_run_backwards_compatability(s,nr.reps = nr.reps,q_BH = 0.1,q_DSFDR = 0.1) }
  stopCluster(cl)
  result_table = unlist(scenario_results[[1]])
  for(i in 2:length(scenario_results)){
    result_table = rbind(result_table,unlist(scenario_results[[i]]))
  }  
  print(result_table)
  
})
