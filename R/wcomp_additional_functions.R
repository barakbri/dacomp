
Compute.resample.test = function(X_matrix,Y_matrix,lambda,statistic = 'Wilcoxon'){
  disable.ties.correction = F
  nr_subsamples = ncol(X_matrix)
  nr_bootstraps  = ncol(Y_matrix)
  stats = rep(0,nr_bootstraps)
  current_stats = rep(0,nr_bootstraps)
  n = nrow(X_matrix)
  for(s in 1:nr_subsamples){
    ranked_X = rank(X_matrix[,s], ties.method = 'average')
    if(statistic == 'Wilcoxon'){
      current_stats = (subzero::rcpp_Wilcoxon_PermTest_Given_Permutations(ranked_X,Y_matrix)[[1]])
      current_stats = current_stats - sum(Y_matrix[,1])/nrow(Y_matrix) * sum(ranked_X)
      
      NR_Y0 = sum(Y_matrix[,1] == 0)
      N = nrow(Y_matrix)
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
      current_stats = current_stats/sqrt(V_H0) #convert to Z score
    }
    if(statistic == 'Avg.Diff'){
      current_stats = rep(0,nr_bootstraps)
      for(b in 1:ncol(Y_matrix)){
        current_stats[b] = mean(X_matrix[Y_matrix[,b]==0,s]) - mean(X_matrix[Y_matrix[,b]==1,s])
      }
    } 
    if(statistic == 'Log.Avg.Diff'){
      current_stats = rep(0,nr_bootstraps)
      for(b in 1:ncol(Y_matrix)){
        current_stats[b] = mean(log10(X_matrix[Y_matrix[,b]==0,s]+1)) - mean(log10(X_matrix[Y_matrix[,b]==1,s]+1))
      }
    }   
    if(statistic == 'TwoPartWilcoxon')
      current_stats = (subzero::rcpp_TwoPartTest_Given_Permutations(ranked_X,Y_matrix)[[1]])
    if(statistic == 'SignedWilcoxon'){
      
      for(j in 1:nr_bootstraps){
        g1 = X_matrix[Y_matrix[1:(n/2),j],s]
        g2 = X_matrix[Y_matrix[(n/2 +1):(n),j],s]
        add_stat = SignedRankWilcoxon.statistic(g1,g2)  
        current_stats[j] = current_stats[j] + add_stat
      }
      
    }
    if(statistic == 'KW'){
      for(j in 1:nr_bootstraps){
        add_stat = (kruskal.test(X_matrix[,s], as.factor(Y_matrix[,j]))$statistic)
        if(is.nan(add_stat))
          add_stat = 0
        current_stats[j] = current_stats[j] + add_stat
      }
    }
    stats = stats + (current_stats)
  }
  stat = stats/nr_subsamples
  if(statistic %in% c('Wilcoxon','Avg.Diff','Log.Avg.Diff','SignedWilcoxon')){
    stats = stats^2
  }
  return(stats)
}



dfdr_find_thresholds = function(stats,q=0.05,verbose = F){
  
  if(verbose){
    cat(paste0('computing rejection threshold for DS-FDR\n\r'))
  }
  
  #convert to pvalues
  for(j in 1:ncol(stats)){
    stats[,j] = (nrow(stats) +1+1 - rank(stats[,j],ties.method = 'min'))/(nrow(stats)+1) # one +1 is for permutations, one +1 is because rank is 1 based
  }
  
  C_possible_values = sort(unique(stats[1,]))
  C_possible_values = C_possible_values[C_possible_values<=q]
  FDR_hat_vec = rep(NA,length(C_possible_values))
  for(c in 1:length(C_possible_values)){
    C = C_possible_values[c]
    Vhat = sum(stats<=C)/nrow(stats)
    Rhat = sum(stats[1,]<=C)
    FDR_hat = Vhat/Rhat
    FDR_hat_vec[c] = FDR_hat
  }
  
  selected_c_ind = (which(FDR_hat_vec<= q))
  selected_c = -1 #no rejections are allowed, this is in P-value "units"
  if(length(selected_c_ind)>0){
    selected_c = max(C_possible_values[selected_c_ind])
  }
  return(selected_c)
}