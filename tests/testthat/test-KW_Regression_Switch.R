test_that("Regression_KW_Switch", {
  cat(paste0('\n\r'))
  
  
  if(!exists('DO_KW_REGRESSION_SWITCH'))
    skip('DO_KW_REGRESSION_SWITCH not defined, skipping')
  if(!DO_KW_REGRESSION_SWITCH  )
    skip('DO_KW_REGRESSION_SWITCH is false, skipping')
  
  #generate data, and check that the ordering of our C level KW is the same as the R KW function, since we switched between stat functions as of version 1.24
  set.seed(1)
  B = 100
  stat_self = rep(NA,B)
  stat_R = rep(NA,B)
  for(i in 1:B){
    N=50
    X = rpois(N,5)
    X = rank(X)
    Y = sample(c('A','B','C'),size = N,replace = T)
    Y_pass = matrix(as.numeric(as.factor(Y))-1,ncol = 1)
    nr_groups = length(unique(Y_pass))
    stat_self[i] = dacomp:::rcpp_KW_PermTest_Given_Permutations(X,Y_pass,nr_groups)[[1]]
    stat_R[i] = (kruskal.test(X, as.factor(Y))$statistic)
    
  }
  
  expect(cor(stat_self,stat_R),1)
  
})



