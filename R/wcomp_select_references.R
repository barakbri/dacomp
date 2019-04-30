CLASS.LABEL.REFERENCE_SELECTION_OBJECT = "wcomp.reference.selection.object"

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
wcomp.select_references = function(X, median_SD_threshold, 
                                                           minimal_TA = 10,
                                                           maximal_TA = 200,
                                                           Psuedo_Count_used = 1,
                                                           verbose = F,
                                                           select_from = NULL){
  input_check_result = check.input.wcomp.select_references(X, median_SD_threshold, minimal_TA,maximal_TA, Psuedo_Count_used, verbose, select_from)
  if(!input_check_result)
    stop('Input check failed on wcomp.select_references')
  
  adjustment = 0
  PSUEDOCOUNT_IS_PREVALENCE = F
  PSUEDOCOUNT_IS_ABUNDANCE_COMPLETION = F
  m = dim(X)[2]
  
  
  if(is.null(select_from)){
    select_from = 1:m
  }
  
  ratio_matrix = matrix(NA, ncol = m, nrow = m)
  
  for( i in 1:(m-1) ){
    if(verbose && i%%(floor(m/100))==1)
      cat(paste0('Computing pairwise ratios for taxon ',i,'/',m,'\n\r'))
    
    for( j in (i+1):m ){
      
      X_i = X[,i]
      X_j = X[,j]
      r = log(((X_i+Psuedo_Count_used) / (X_j+Psuedo_Count_used)))
      ratio_matrix[i,j] = sd(r)
      ratio_matrix[j,i] = ratio_matrix[i,j]
    }
  }

  scores         = apply(ratio_matrix,2,function(x){median(x,na.rm = T)}) 
  prevalence_mat = 1*(X > 0)
  mean_prevalance = apply(prevalence_mat,2,mean)
  filter_cols = 1:m
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
  class(ret) = CLASS.LABEL.REFERENCE_SELECTION_OBJECT
  return(ret)
}

check.input.wcomp.select_references = function(X, median_SD_threshold, minimal_TA,maximal_TA, Psuedo_Count_used, verbose, select_from){
  return(T)  
}

