
#internal function for computing test statistics for original data and data with permuted labels
# the function supports multiple tests (with the same permutations across tests)
# X_matrix is a matrix with a single column, with rarefied counts, across subjects
# coloumns of Y_matrix are different permutations of the group labels
# rows are samples
# the returned statistics are a sum of individual test statistics across columns of X (test the joint null hypothesis)
Compute.resample.test = function(X_matrix,Y_matrix,statistic = 'Wilcoxon'){
  
  if(ncol(X_matrix) != 1)
    stop('Error in Compute.resample.test: Function requires X_matrix with a single column')
  
  disable.ties.correction = F
  nr_subsamples = ncol(X_matrix)
  nr_bootstraps  = ncol(Y_matrix)
  stats = rep(0,nr_bootstraps)
  current_stats = rep(0,nr_bootstraps)
  n = nrow(X_matrix)
  
  # iterate over subsamples of the data
  
  ranked_X = rank(X_matrix[,1], ties.method = 'average') #break ties by average rank
  
  #the wilcoxon rank sum test
  if(statistic == 'Wilcoxon'){
    
    current_stats = (subzero::rcpp_Wilcoxon_PermTest_Given_Permutations(ranked_X,Y_matrix)[[1]]) #compute statistic and a sample of test statistic given from the null hypothesis
    current_stats = current_stats - sum(Y_matrix[,1])/nrow(Y_matrix) * sum(ranked_X)
    
    #compute test statistic with correction for ties:
    
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
  if(statistic == 'Avg.Diff'){ #test statistic is the difference in mean counts across two sample groups
    current_stats = rep(0,nr_bootstraps)
    for(b in 1:ncol(Y_matrix)){
      current_stats[b] = mean(X_matrix[Y_matrix[,b]==0,s]) - mean(X_matrix[Y_matrix[,b]==1,s])
    }
  } 
  if(statistic == 'Log.Avg.Diff'){ #test statistic is the log fold change across two sample groups
    current_stats = rep(0,nr_bootstraps)
    for(b in 1:ncol(Y_matrix)){
      current_stats[b] = mean(log10(X_matrix[Y_matrix[,b]==0,s]+1)) - mean(log10(X_matrix[Y_matrix[,b]==1,s]+1))
    }
  }   
  if(statistic == 'TwoPartWilcoxon'){#chi square score of a wilcoxon test and a two-sample test for equality of proportions (for zeroes in the data)
    X_ranked_without_zeroes = X_matrix[,1]
    X_ranked_without_zeroes[X_ranked_without_zeroes>0] = rank(X_ranked_without_zeroes[X_ranked_without_zeroes>0],ties.method = 'average')
    current_stats = (subzero::rcpp_TwoPartTest_Given_Permutations(X_ranked_without_zeroes,Y_matrix)[[1]])
  } 
    
  
  if(statistic == 'SignedWilcoxon'){ # signed wilcoxon
    
    for(j in 1:nr_bootstraps){
      g1 = X_matrix[Y_matrix[1:(n/2),j],s]
      g2 = X_matrix[Y_matrix[(n/2 +1):(n),j],s]
      current_stats[j] = SignedRankWilcoxon.statistic(g1,g2)  
    }
    
  }
  if(statistic == 'KW'){ #Kruskal Wallis
    for(j in 1:nr_bootstraps){
      add_stat = (kruskal.test(X_matrix[,s], as.factor(Y_matrix[,j]))$statistic)
      if(is.nan(add_stat))
        add_stat = 0
      current_stats[j] = add_stat
    }
  }
  stats = (current_stats)
  
  # tests with non negative test scores, and a two-sided hypothesis (with zero being mode of dist under )
  if(statistic %in% c('Wilcoxon','Avg.Diff','Log.Avg.Diff','SignedWilcoxon')){
    stats = stats^2
  }
  return(stats)
}

#wrapper for the signed rank test:
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


# A function for computing the rejection threshold for the DS-FDR procedure as described in Jiang et al. (2017)
# stats is matrix of test statistics, with a right sided alternative hypothesis.
# Each column is a different hypothesis, each row is a different resample of the group labels.
# The first row contains the test statistics for the original data
# q sets the requested FDR level for the DS-FDR procedure.
# The function returns a single number, the rejection threshold for the DS-FDR procedue in P-Value units (!!!).
# Hypotheses with P-values below or equal to computed threshold should be rejected.

dsfdr_find_thresholds = function(stats,q=0.05,verbose = F){
  
  if(verbose){
    cat(paste0('computing rejection threshold for DS-FDR\n\r'))
  }
  
  #convert to pvalues, right sided hypothesis. Note the treament of ties by ties.method
  for(j in 1:ncol(stats)){
    stats[,j] = (nrow(stats) +1+1 - rank(stats[,j],ties.method = 'min'))/(nrow(stats)+1) # one +1 is for permutations, one +1 is because rank is 1 based
  }
  #find the possible cutoff values
  C_possible_values = sort(unique(stats[1,]))
  C_possible_values = C_possible_values[C_possible_values<=q]
  
  # For each cutoff value, compute the estimated FDR
  FDR_hat_vec = rep(NA,length(C_possible_values))
  for(c in 1:length(C_possible_values)){
    C = C_possible_values[c]
    Vhat = sum(stats<=C)/nrow(stats)
    Rhat = sum(stats[1,]<=C)
    FDR_hat = Vhat/Rhat
    FDR_hat_vec[c] = FDR_hat
  }
  # select the highest threshold where FDR is maintained
  selected_c_ind = (which(FDR_hat_vec<= q))
  selected_c = -1 #no rejections are allowed, this is in P-value "units"
  if(length(selected_c_ind)>0){
    selected_c = max(C_possible_values[selected_c_ind])
  }
  return(selected_c)
}