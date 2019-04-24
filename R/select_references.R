#' Title
#'
#' @param X 
#' @param median_SD_threshold 
#' @param minimal_TA 
#' @param maximal_TA 
#' @param Psuedo_Count_used 
#' @param verbose 
#' @param select_from 
#'
#' @return
#' @export
#'
#' @examples
wcomp.select_references = function(X, median_SD_threshold = 1.3, 
                                                           minimal_TA = 10,
                                                           maximal_TA = 200,
                                                           Psuedo_Count_used = 1,
                                                           verbose = F,
                                                           select_from = NULL){
  factor_by_Median_Score = F
  adjustment = 0
  MIN_PREV = 0.0
  MAX_PREV = 1.0
  PSUEDOCOUNT_IS_PREVALENCE = F
  PSUEDOCOUNT_IS_ABUNDANCE_COMPLETION = F
  m = dim(X)[2]
  X_r  = X
  
  if(is.null(select_from)){
    select_from = 1:m
  }
  
  ratio_matrix = matrix(NA, ncol = m, nrow = m)
  PS_VEC = as.numeric(apply(X_r>0,2,mean))
  if(!PSUEDOCOUNT_IS_PREVALENCE){
    PS_VEC = rep(Psuedo_Count_used,ncol(X_r))
  }
  TOTAL_COUNTS_PER_SUBJECT = as.numeric(apply(X_r,1,sum))
  X_r_PS = X_r
  for(i in 1:nrow(X_r_PS)){
    X_r_PS[i,] = X_r_PS[i,]/TOTAL_COUNTS_PER_SUBJECT[i]
  }
  p_by_taxa = as.numeric(apply(X_r_PS,2,mean))
  
  for( i in 1:(m-1) ){
    if(verbose)
      print(paste0('Computing ratio for taxa ',i))
    
    for( j in (i+1):m ){
      
      X_i = X_r[,i]
      X_j = X_r[,j]
      r = log(((X_i+PS_VEC[i]) / (X_j+PS_VEC[j])))
      if(PSUEDOCOUNT_IS_ABUNDANCE_COMPLETION)
        r = log(((X_i+pmin(1,TOTAL_COUNTS_PER_SUBJECT*p_by_taxa[i])) / (X_j+pmin(1,TOTAL_COUNTS_PER_SUBJECT*p_by_taxa[j]))))
      ratio_matrix[i,j] = sd(r)
      ratio_matrix[j,i] = ratio_matrix[i,j]
    }
  }

  scores         = apply(ratio_matrix,2,function(x){median(x,na.rm = T)}) #mean
  prevalence_mat = 1*(X > 0)
  mean_prevalance = apply(prevalence_mat,2,mean)
  filter_cols = which(mean_prevalance>=MIN_PREV & mean_prevalance <= MAX_PREV)
  filter_cols = filter_cols[which(filter_cols %in% select_from)]
  scores = scores[filter_cols]
  sorted_columns = order(scores)
  original_ind = (filter_cols)[sorted_columns]
  sorted_scores = scores[sorted_columns]
    
  sorted_prevalence_mat = prevalence_mat[ , original_ind]
  sorted_X_mat = X[,original_ind]
  
  cummulative_sorted_prevalence_mat = sorted_prevalence_mat
  cummulative_sorted_X_mat = sorted_X_mat
  for(i in 1:nrow(cummulative_sorted_prevalence_mat)){
    cummulative_sorted_prevalence_mat[i,] = cummax(cummulative_sorted_prevalence_mat[i,])
    cummulative_sorted_X_mat[i,] = cumsum(cummulative_sorted_X_mat[i,])
  }
    
  
  mean_prevalence_over_the_sorted = as.numeric(apply(cummulative_sorted_prevalence_mat,2,mean))
  min_abundance_over_the_sorted = as.numeric(apply(cummulative_sorted_X_mat,2,min))
  
  possible_cut_points = which(min_abundance_over_the_sorted>=minimal_TA & min_abundance_over_the_sorted<=maximal_TA)
  if(length(possible_cut_points) == 0)
    possible_cut_points = min(which(min_abundance_over_the_sorted>=minimal_TA))
  if(length(possible_cut_points) == 1 & is.infinite(possible_cut_points[1]))
    possible_cut_points = min(which(min_abundance_over_the_sorted>=1))
  
  scores_possible_cut_points = sorted_scores[possible_cut_points]
  
  threshold_to_use = median_SD_threshold
  if(factor_by_Median_Score)
    threshold_to_use = median(sorted_scores)*median_SD_threshold
  
  possible_cut_points_above_threshold = which(scores_possible_cut_points >= threshold_to_use)
  if(length(possible_cut_points_above_threshold)>0)
    cut_point_selected = possible_cut_points[min(possible_cut_points_above_threshold)]
  else
    cut_point_selected = max(possible_cut_points)
  
  selected_references = original_ind[1:cut_point_selected]
  selected_MinAbundance = min_abundance_over_the_sorted[cut_point_selected]
  
  ret = list()
  ret$selected_references = selected_references
  ret$mean_prevalence_over_the_sorted = mean_prevalence_over_the_sorted
  ret$min_abundance_over_the_sorted = min_abundance_over_the_sorted
  ret$ratio_matrix = ratio_matrix
  ret$scores = scores
  ret$selected_MinAbundance = selected_MinAbundance
  ret$median_SD_threshold = median_SD_threshold
  ret$minimal_TA = minimal_TA
  ret$maximal_TA = maximal_TA
  return(ret)
}

