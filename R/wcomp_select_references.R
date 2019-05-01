#object label for reference selection results
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
  #check inputs
  input_check_result = check.input.wcomp.select_references(X, median_SD_threshold, minimal_TA,maximal_TA, Psuedo_Count_used, verbose, select_from)
  if(!input_check_result)
    stop('Input check failed on wcomp.select_references')
  
  m = dim(X)[2]
  
  # for the default, null, all taxa can serve as reference
  if(is.null(select_from)){
    select_from = 1:m
  }
  
  #compute the SD_{j,k}, see paper for definition
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

  #compute the different median SD scores
  scores         = apply(ratio_matrix,2,function(x){median(x,na.rm = T)}) 
  #compute prevalance of taxa
  prevalence_mat = 1*(X > 0)
  mean_prevalance = apply(prevalence_mat,2,mean)
  
  #filter out columns that the user asked to not be considered
  filter_cols = 1:m
  filter_cols = filter_cols[which(filter_cols %in% select_from)]
  scores = scores[filter_cols]
  sorted_columns = order(scores)
  original_ind = (filter_cols)[sorted_columns]
  sorted_scores = scores[sorted_columns]
  
  #sort the prevalence matrix, by the order of medianSD scores.  
  sorted_prevalence_mat = prevalence_mat[ , original_ind]
  sorted_X_mat = X[,original_ind]
  
  # compute the cummulative prevalence and abundance of samples, for taxa from the lowest medianSD score up to any column on the matrix.
  cummulative_sorted_prevalence_mat = sorted_prevalence_mat
  cummulative_sorted_X_mat = sorted_X_mat
  for(i in 1:nrow(cummulative_sorted_prevalence_mat)){
    cummulative_sorted_prevalence_mat[i,] = cummax(cummulative_sorted_prevalence_mat[i,])
    cummulative_sorted_X_mat[i,] = cumsum(cummulative_sorted_X_mat[i,])
  }
    
  
  # minimum abundance across taxa, for each possible set of references, by including one additional taxon at a time in the reference set:
  mean_prevalence_over_the_sorted = as.numeric(apply(cummulative_sorted_prevalence_mat,2,mean))
  min_abundance_over_the_sorted = as.numeric(apply(cummulative_sorted_X_mat,2,min))
  
  #find the possible cut point by required abundance (these will be the fall back, if there are no point to find by requested threshold:
  
  possible_cut_points = which(min_abundance_over_the_sorted>=minimal_TA & min_abundance_over_the_sorted<=maximal_TA)
  # if no point are found, go the the firt place (effectivly increase the threshold) the the minimal number of counts is found in each subject
  if(length(possible_cut_points) == 0)
    possible_cut_points = min(which(min_abundance_over_the_sorted>=minimal_TA))
  
  # if not possible, find a cut point which is possible, with a lower number of counts
  if(length(possible_cut_points) == 1 & is.infinite(possible_cut_points[1]))
    possible_cut_points = min(which(min_abundance_over_the_sorted>=1))
  
  # out of all possible points that meet demands for abundance, we first consider those that also meet the demand for median SD:
  scores_possible_cut_points = sorted_scores[possible_cut_points] 
  
  threshold_to_use = median_SD_threshold
  
  possible_cut_points_above_threshold = which(scores_possible_cut_points >= threshold_to_use)
  if(length(possible_cut_points_above_threshold)>0) # if there are points which also meet the medianSD demand, we take the smallest one
    cut_point_selected = possible_cut_points[min(possible_cut_points_above_threshold)]
  else
    cut_point_selected = max(possible_cut_points) #we fall back to requiring only a minimal number of counts
  
  #the selected references
  selected_references = original_ind[1:cut_point_selected]
  
  #the minimum number of counts observed for reference taxa in a subject
  selected_MinAbundance = min_abundance_over_the_sorted[cut_point_selected] 
  
  #return results:
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
  
  # X - matrix of counts
  # median_SD_threshold - should be a valid, not positive number. Print out warnings
  # minimal_TA, - should be interger, reasonable range (warning), smaller then maximal_TA
  # maximal_TA,
  # Psuedo_Count_used,- should be valid integer
  # verbose, - should be logical
  # select_from - should be entirely contained in 1:m
  return(T)  
}

