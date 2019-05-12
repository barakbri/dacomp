
Run_Test_Modular = function(report_file = 'E:/temp/wcomp_test_results.txt',
                            package_location = 'E:/wcomp/',
                            CompositionalAnalysis_CodeBase_Location = 'E:/MCB3/CompositionalAnalysis_CodeBase/MCB_Simulation/'){
  sink(file = report_file)
  
  set.seed(1)
  original_wd = getwd()
  setwd(CompositionalAnalysis_CodeBase_Location)
  source(paste0('REFSIM_GenerateSettings_Index.R'),chdir = T)
  setwd(package_location)
  
  DO_SIMPLE_USE_CASE <<- T
  DO_EXAMPLE_DATA <<- T
  DO_REFERENCE_SELECTION <<- T
  DO_REFERENCE_VALIDATION <<- T
  DO_PLOTTING <<- T
  DO_MAIN_TESTING <<- T
  DO_REGRESSION_TESTING <<- F
  devtools::test()
  sink()  
}

#Run_Test_Modular()