
PermSumCount.test = function(X,Y,B = as.integer(c(40000)),DoWald = as.integer(c(1)),return_perms = F,disable.ties.correction = F){
  ranked_X = rank(X,ties.method = 'average')
  perms_val = as.integer(c(0))
  if(return_perms)
    perms_val = as.integer(c(1))
  res = subzero::rcpp_Wilcoxon_PermTest(ranked_X,Y,B,DoWald,perms_val)
  
  NR_Y0 = sum(Y == 0)
  N = length(Y)
  NR_Y1 = N - NR_Y0
  
  E_H0 = sum(ranked_X) * NR_Y1 / N
  
  unique_value_counts = table(ranked_X)
  which_unique_value_counts_are_ties = which(unique_value_counts>1)
  unique_value_counts_only_ties_correction = 0
  if(length(which_unique_value_counts_are_ties)>0 & !disable.ties.correction){
    t_r = unique_value_counts[which_unique_value_counts_are_ties]
    unique_value_counts_only_ties_correction  = sum(t_r * (t_r^2 - 1))
  }
  
  V_H0 = NR_Y0 * NR_Y1 * (N+1) / 12 - NR_Y0 * NR_Y1 * unique_value_counts_only_ties_correction / (12 * N * (N - 1))
  Z_score = (res[[1]] - E_H0)/sqrt(V_H0)
  pval.by.normal = 2*(1-pnorm(abs(Z_score)))
  
  ret = list(stat = res[[1]],p.value = res[[2]],b = res[[3]],pval.by.normal=pval.by.normal)
  if(return_perms)
    ret$perms = res[[7]]
  return(ret)
}



PermSumCount.multipleX.test = function(List_of_X,Y,B = as.integer(c(1000)),is.Asymptotic = F,disable.ties.correction = F){
  p = length(List_of_X)
  pval.vec = rep(NA,p)
  if(is.Asymptotic)
    B = 0
  
  for(i in 1:p){
    res = subzero::PermSumCount.test(List_of_X[[i]],Y,B,DoWald = as.integer(c(0)),return_perms = F,disable.ties.correction = disable.ties.correction)
    if(is.Asymptotic)
      pval.vec[i] = res$pval.by.normal
    else
      pval.vec[i] = res$p.value
  }
  
  p.value = mean(pval.vec)
  return(list(p.value=p.value))
}



TwoPartTest.test = function(X,Y,B = as.integer(c(40000)),DoWald = as.integer(c(1)),return_perms = F,disable.ties.correction = F){
  perms_val = as.integer(c(0))
  if(return_perms)
    perms_val = as.integer(c(1))
  
  X_ranked = X
  X_ranked[X>0] = rank(X_ranked[X>0],ties.method = 'average')
  
  
  
  res = subzero::rcpp_TwoPartTest(X_ranked,Y,B,DoWald,perms_val)
  additional_stats = res[[5]]
  
  
  Z_proportion = additional_stats[1];
  Z_wilcoxon = additional_stats[2] ;
  total_counts_1 = additional_stats[3];
  counts_zeros_0 = additional_stats[4];
  counts_zeros_1 = additional_stats[5];
  nr_samples_0_non_zero = additional_stats[6] ;
  nr_samples_1_non_zero = additional_stats[7];
  
  n_wilcoxon = nr_samples_0_non_zero + nr_samples_1_non_zero
  
  E_H0 = sum(X_ranked[X>0]) * nr_samples_1_non_zero / n_wilcoxon
  
  unique_value_counts = table(X_ranked[X>0])
  which_unique_value_counts_are_ties = which(unique_value_counts>1)
  unique_value_counts_only_ties_correction = 0
  if(length(which_unique_value_counts_are_ties)>0 & !disable.ties.correction){
    t_r = unique_value_counts[which_unique_value_counts_are_ties]
    unique_value_counts_only_ties_correction  = sum(t_r * (t_r^2 - 1))
  }
  
  V_H0 = nr_samples_0_non_zero * nr_samples_1_non_zero * (n_wilcoxon+1) / 12 - nr_samples_0_non_zero * nr_samples_1_non_zero * unique_value_counts_only_ties_correction / (12 * n_wilcoxon * (n_wilcoxon - 1))
  Z_Wilcoxon_Asymptotic = (total_counts_1 - E_H0)/sqrt(V_H0)
  
  
  
  if(is.na(Z_Wilcoxon_Asymptotic)){
    Z_Wilcoxon_Asymptotic = 0
  }
  
  
  if(is.na(Z_proportion)){
    Z_proportion = 0
  }
  
  stat_ties_corrected = Z_Wilcoxon_Asymptotic^2 + Z_proportion^2
  
  ret = list(stat = res[[1]],p.value = res[[2]],b = res[[3]],asymp.pval = 1-pchisq(stat_ties_corrected,2))  
  if(return_perms)
    ret$perms = res[[4]]
  return(ret)
}

TwoPartTest.multipleX.test = function(List_of_X,Y,B = as.integer(c(1000)),is.Asymptotic = F,disable.ties.correction = F){
  p = length(List_of_X)
  pval.vec = rep(NA,p)
  if(is.Asymptotic)
    B = 0
  
  for(i in 1:p){
    res = subzero::TwoPartTest.test(List_of_X[[i]],Y,B,DoWald = as.integer(c(0)),return_perms = F,disable.ties.correction = disable.ties.correction)
    if(is.Asymptotic)
      pval.vec[i] = res$asymp.pval
    else
      pval.vec[i] = res$p.value
  }
  p.value = mean(pval.vec)
  return(list(p.value=p.value))
  
}


SignedRankWilcoxon.statistic = function(X1,X2){
  differences = X1-X2
  second_is_bigger = 1*(differences > 0)
  differences[differences==0] = NA
  rank_differences = rank(abs(differences), ties.method = 'average',na.last = "keep")
  
  res = rcpp_Compute_Wilcoxon_Signed_Rank_Stat(rank_differences,second_is_bigger)
  #statistic:
  stat_mean = sum(rank_differences)/2
  stat_var = sum(rank_differences^2)/4
  stat = (res - stat_mean)/sqrt(stat_var)
  return(stat)
}
