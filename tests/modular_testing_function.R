
Run_Test_Modular = function(report_directory = 'C:/Test_DACOMP/',
                            package_location = 'C:/dacomp/',
                            CompositionalAnalysis_CodeBase_Location = 'C:/MCB3/CompositionalAnalysis_CodeBase/Scripts/'){
  report_file = paste(report_directory,'dacomp_test_results.txt')
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
  DO_REGRESSION_TESTING <<- T
  DO_TEST_SIGNED <<-T
  DO_KW_REGRESSION_SWITCH <<-T
  DO_PAIRED_T1E_AND_POWER_MULTIPLE_CORES <<-T
  DO_BACKWORDS_COMPATABILITY_MULTIPLE_CORES<<-T
  
  devtools::test()
  sink()  
}

#Run_Test_Modular()