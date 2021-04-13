
test_that("Test paired T1E and Power using multiple cores", {
  
  
  if(!exists('DO_PAIRED_T1E_AND_POWER_MULTIPLE_CORES'))
    skip('DO_PAIRED_T1E_AND_POWER_MULTIPLE_CORES not defined, skipping')
  if(!DO_PAIRED_T1E_AND_POWER_MULTIPLE_CORES)
    skip('DO_PAIRED_T1E_AND_POWER_MULTIPLE_CORES is false, skipping')
  
  
  ##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #T1E & power simulation:
  ##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  
  paired_setting_power_and_FDR = function(signal_strength=0.1,REPS = 1,n=30,m1 = 30){
    library(dacomp)
    fdr.bh = NULL
    fdr.dsfdr = NULL
    tp.bh = NULL
    tp.dsfdr = NULL
    
    for(r in 1:REPS){
      print(paste0('Doing rep: ',r,'/',REPS))
      data = dacomp.generate_example_dataset_paired(n = n,m1 = m1,signal_strength_as_change_in_microbial_load = signal_strength)  
      
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
      
      rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH') <= q_BH)
      rejected_DSFDR = result.test$dsfdr_rejected
      
      tp.bh = c(tp.bh,sum(rejected_BH %in% data$select_diff_abundant))
      tp.dsfdr = c(tp.dsfdr,sum(rejected_DSFDR %in% data$select_diff_abundant))
      
      current_fdr_dsfdr = 0
      if(length(rejected_DSFDR) >0)
        current_fdr_dsfdr = 1- sum(rejected_DSFDR %in% data$select_diff_abundant) / length(rejected_DSFDR)
      
      current_fdr_bh = 0
      if(length(rejected_BH) >0)
        current_fdr_bh = 1- sum(rejected_BH %in% data$select_diff_abundant) / length(rejected_BH)
      
      fdr.bh = c(fdr.bh,current_fdr_bh)  
      fdr.dsfdr = c(fdr.dsfdr,current_fdr_dsfdr)
    }
    ret= list(tp.bh = mean(tp.bh),fdr.bh = mean(fdr.bh),tp.dsfdr = mean(tp.dsfdr),fdr.dsfdr = mean(fdr.dsfdr))
    return(ret)
  }
  
  set.seed(1)
  library(doRNG)
  library(doParallel)
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  signal_vec = c(0,0.05,0.1,0.15,0.2,0.25,0.3)
  
  scenario_results = foreach(s=1:length(signal_vec), .options.RNG=1234,.combine = rbind) %dorng% {
    paired_setting_power_and_FDR(signal_strength = signal_vec[s],REPS = 10)
  }
  #scenario_results
  print(head(scenario_results,n = 10))
  stopCluster(cl)
  
  
})










