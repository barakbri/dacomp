#object label for reference selection results
CLASS.LABEL.REFERENCE_SELECTION_OBJECT = "dacomp.reference.selection.object"

#' Select a set of reference taxa, used for testing differential abundance of other taxa.
#'
#' The function receives a table of microbiome counts, and selects a set of taxa used for normalization in  \code{\link{dacomp.test}}.
#' The counts matrix should be formated with taxa as columns, samples as rows. No rarefaction or perliminary normalization is required.
#' The first step of the computation consists of computing the standard deviation of the ratio of each pair of taxa, across subjects. Each taxon is involved in the computation of \eqn{m-1} pairwise standard deviations, with \eqn{m} being the number of taxa. The second step consists of finding the median pairwise standard deviation of each taxon, across all pairwise standard deviations computed.
#' This value, computed for each taxon, is named the `median SD` statistic. Taxa with values of the `median SD` statistic below \code{median_SD_threshold} are selected as the reference set of taxa.
#' See \code{vignette('dacomp_main_vignette')} for additional details and formulas.
#' The full algorithm and additional details are found in Subsection 3.1 in Brill et. al. (2019). 
#' 
#' @details 
#' Target Abundance (TA) limits: The function will attempt to use the argument \code{median_SD_threshold} as a threshold for classification. If needed, the function will increase the actual threshold used so that each sample has at least \code{minimal_TA} counts under taxa selected as reference. The function will lower the classification threshold if all samples have more than \code{maximal_Ta} reads under the selected set of reference taxa.
#' Computation may take up to a minute or two, for large datasets.
#' 
#' @param X Counts matrix, with rows representing samples, columns representing different taxa.
#' @param median_SD_threshold Critical value for the `median SD` statistic. Taxa with a `median SD` statistic will be taken.
#' @param minimal_TA The minimal number of counts required in each sample, in the taxa selected as reference. If the set of reference taxa has a sample with less than `minimal_TA` reads in the set of reference taxa selected, the function will increase the `median SD` value until all samples have at least \code{minimal_TA} reads in the selected set of reference taxa.
#' @param maximal_TA If all samples have more than \code{maximal_TA} reads available under the selected set of reference taxa, the `median SD` threshold value for classifying taxa as references will be lowered, until the condition is met. 
#' @param Pseudo_Count_used Pseudo count added to all count values, to avoid dividing by zero.
#' @param verbose If set to \code{TRUE}, messages will be displayed, as computation progresses.
#' @param select_from The default value of \code{NULL} indicates that all taxa are valid candidates for selection. The user may limit the set of possible candidates for the reference set, by supplying a list of candidates (by indices) using this argument.
#' @param Previous_Reference_Selection_Object If the user previously selected a set of reference taxa for the data with one threshold, and would like to select a new set of reference taxa with another threshold, the output of a previous run of \code{dacomp.select_references} can be supplied as an argument to speed up computations. See usage example in code snippet below.
#' @param run_in_parallel should computation be parallelized
#' @param Nr.Cores if computation is parallelised, how many parallel workers should be used
#' 
#' @return The function returns an object of type "dacomp.reference.selection.object", which is a list with the following fields:
#' \itemize{
#' \item{selected_references}{ - A vector with the indices of selected reference taxa.}
#' \item{mean_prevalence_over_the_sorted}{ - A vector, containing fraction of zero counts in the reference set of taxa, across samples, if: the lowest median SD are taken as reference, two lowest median SD are taken as reference, three lowest...}
#' \item{min_abundance_over_the_sorted}{ - A vector, containing the minimal number of counts observed in the reference set of taxa, across samples, if: the lowest median SD are taken as reference, two lowest median SD are taken as reference, three lowest...}
#' \item{ratio_matrix}{ - The matrix of SD_{j,k} as defined in the paper and the package vignette.}
#' \item{scores}{ - the median SD scores, S_j as defined in the package vignette and paper.}
#' \item{selected_MinAbundance}{ - The minimal number of counts, available under the reference set of taxa, across subjects.}
#' \item{median_SD_threshold}{ - The input supplied under the function argument with the same name, without modification by the Target Abundance feature (see 'details').}
#' \item{minimal_TA}{ - The input supplied under the function argument with the same name.}
#' \item{maximal_TA}{ - The input supplied under the function argument with the same name.}
#' }
#' 
#' @references 
#' Brill, Barak, Amnon Amir, and Ruth Heller. 2019. Testing for Differential Abundance in Compositional Counts Data, with Application to Microbiome Studies. arXiv Preprint arXiv:1904.08937.
#' 
#' @export
#'
#' @examples
#' #' \dontrun{
#' library(dacomp)
#' 
#' set.seed(1)
#' 
#' data = dacomp.generate_example_dataset.two_sample(m1 = 100,
#'        n_X = 50,
#'        n_Y = 50,
#'        signal_strength_as_change_in_microbial_load = 0.1)
#' 
#' # Select references: (may take a minute)
#' result.selected.references = dacomp.select_references(X = data$counts,
#'                                                      median_SD_threshold = 0.6, #APPLICATION SPECIFIC
#'                                                      verbose = T)
#' 
#' length(result.selected.references$selected_references)
#'
#' # Plot the reference selection scores (can also be used to better set the median SD threshold)
#' dacomp.plot_reference_scores(result.selected.references)
#' 
#' # Select a reference set with a different (lower) threshold. user can use the function argument
#' # 'Previous_Reference_Selection_Object' to provide a previous reference seleciton object for this data, to speed up computation.
#' 
#' result.selected.references.different.threshold = dacomp.select_references(X = data$counts,
#'           median_SD_threshold = 0.5, 
#'           verbose = F,
#'           Previous_Reference_Selection_Object = result.selected.references)
#'
#'
#' } 
dacomp.select_references = function(X, median_SD_threshold, 
                                                           minimal_TA = 10,
                                                           maximal_TA = 200,
                                                           Pseudo_Count_used = 1,
                                                           verbose = F,
                                                           select_from = NULL,
                                                           Previous_Reference_Selection_Object = NULL,
                                                           run_in_parallel = F,
                                                           Nr.Cores = 4
                                                           ){
  
  #we now call parallel_reference_select which supports additional options
  return(parallel_reference_select(X = X,
                                   median_SD_threshold = median_SD_threshold,
                                   minimal_TA = minimal_TA,
                                   maximal_TA = maximal_TA,
                                   Pseudo_Count_used = Pseudo_Count_used,
                                   verbose=verbose,
                                   select_from = select_from,
                                   Previous_Reference_Selection_Object = Previous_Reference_Selection_Object,run_in_parallel = F,Nr.Cores = 2))
}

check.input.select_references = function(X, median_SD_threshold, minimal_TA,maximal_TA, Pseudo_Count_used, verbose, select_from, Previous_Reference_Selection_Object){
  
  # X - matrix of counts
  MSG_X = 'X must be a valid counts matrix'
  if(!is.matrix(X))
    stop(MSG_X)
  if(any(X!=as.integer(X)))
    stop(MSG_X)
  if(any(X<0))
    stop(MSG_X)  
  
  # median_SD_threshold - should be a valid, not positive number. Print out warnings
  MSG_median_SD_threshold = "median_SD_threshold must be a valid numeric value, for most data types has values in the range [0.5,1.5]"
  if(!is.numeric(median_SD_threshold))
    stop(MSG_median_SD_threshold)
  if(median_SD_threshold <= 0)
    stop(MSG_median_SD_threshold)
  
  if(median_SD_threshold < 0.5 | median_SD_threshold > 1.5)
    warning('median_SD_threshold set to abnormal value, normally in range [0.5,1.5]')
  
  # minimal_TA,maximal_TA - should be interger, reasonable range (warning), smaller then maximal_TA
  MSG_minimal_TA = "minimal_TA, maximal_TA should be positive integers, valid range for the minimal number of counts in reference taxa, across subjects"
  if( (minimal_TA != as.integer(minimal_TA))| (maximal_TA != as.integer(maximal_TA)))
    stop(MSG_minimal_TA)
  if(minimal_TA>=maximal_TA | minimal_TA < 1)
    stop(MSG_minimal_TA)
  
  # Pseudo_Count_used,- should be valid integer
  MSG_Pseudo_Count_used = "Pseudo_Count_used must be a number, greater than zero."
  if(!is.numeric(Pseudo_Count_used))
    stop(MSG_Pseudo_Count_used)
  if(Pseudo_Count_used <= 0)
    stop(MSG_Pseudo_Count_used)
  
  # verbose, - should be logical
  if(!is.logical(verbose))
    stop('verbose should be a valid logical value')
  
  # select_from - should be entirely contained in 1:m
  if(!is.null(select_from))
    if(any(!(select_from %in% (1:ncol(X)))))
      stop('select_from must be a subset of 1:ncol(X)')
  
  if(!is.null(Previous_Reference_Selection_Object)){
    if(class(Previous_Reference_Selection_Object)!= CLASS.LABEL.REFERENCE_SELECTION_OBJECT)
      stop(paste0('Previous_Reference_Selection_Object is used to reselect a set of references, with a different Median SD threshold, given a previous object of type ',CLASS.LABEL.REFERENCE_SELECTION_OBJECT,'. This is used in order to save computation time. Otherwise, the argument of Previous_Reference_Selection_Object should be set to the default value of NULL'))
  }
  
  return(T)  
}


parallel_reference_select = function(X, median_SD_threshold, 
                                     minimal_TA = 10,
                                     maximal_TA = 200,
                                     Pseudo_Count_used = 1,
                                     verbose = F,
                                     select_from = NULL,
                                     Previous_Reference_Selection_Object = NULL,
                                     Nr.Cores = parallel::detectCores()-1,
                                     Precomputed_scores = NULL,
                                     run_in_parallel = T,
                                     skip_checks = F){
  
  if(!skip_checks & run_in_parallel){
    if(!require('doRNG') | !require('parallel') | !require('doParallel')){
      warning('doRNG, parallel and doParallel must be installed for parallelized computation in parallel_reference_select. At least one of the packages is not installed, so running on single core.')
      parallelize_computation = F
    }  
  }
  
  
  #check inputs
  if(!skip_checks){
    input_check_result = dacomp:::check.input.select_references(X, median_SD_threshold, minimal_TA,maximal_TA, Pseudo_Count_used, verbose, select_from, Previous_Reference_Selection_Object)
    if(!input_check_result)
      stop('Input check failed on dacomp.select_references')
  }
  
  if(is.null(Previous_Reference_Selection_Object) & is.null(Precomputed_scores)){
    m = dim(X)[2]
    
    #compute the SD_{j,k}, see paper for definition
    ratio_matrix = matrix(NA, ncol = m, nrow = m)
    
    
    SDs_for_fixed_i = function(i){
      ret = rep(NA,m)
      X_i = X[,i]
      for( j in (i+1):m ){
        X_j = X[,j]
        r = log(((X_i+Pseudo_Count_used) / (X_j+Pseudo_Count_used)))
        ret[j] = sd(r)
      }
      return(ret)
    }
    
    if(run_in_parallel){
      cl <- makeCluster(Nr.Cores)
      registerDoParallel(cl)
      
      ratio_matrix = foreach(i=1:(m-1), .options.RNG=1234,.combine = 'cbind') %dorng% {
        SDs_for_fixed_i(i)
      }
      stopCluster(cl)
      ratio_matrix = cbind(ratio_matrix,rep(NA,nrow(ratio_matrix)))
    }else{
      for(i in 1:(m-1)){
        ratio_matrix[,i] = SDs_for_fixed_i(i)
      }
    }
    
    
    for( i in 1:(m-1) ){
      for( j in (i+1):m ){
        #print(paste0("i ",i, "j ",j))
        ratio_matrix[i,j] = ratio_matrix[j,i]
      }
    }
    
  }else{
    if(is.null(Precomputed_scores)){
      m = ncol(Previous_Reference_Selection_Object$ratio_matrix)
      ratio_matrix = Previous_Reference_Selection_Object$ratio_matrix  
    }else{
      m = dim(X)[2]
    }
    
  }
  
  
  
  # for the default, null, all taxa can serve as reference
  if(is.null(select_from)){
    select_from = 1:m
  }
  
  
  #compute the different median SD scores
  if(is.null(Precomputed_scores))
    scores         = apply(ratio_matrix,2,function(x){median(x,na.rm = T)}) 
  else
    scores = Precomputed_scores         
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
  if(is.null(Precomputed_scores))
    ret$ratio_matrix = ratio_matrix
  ret$scores = scores
  ret$selected_MinAbundance = selected_MinAbundance
  ret$median_SD_threshold = median_SD_threshold
  ret$minimal_TA = minimal_TA
  ret$maximal_TA = maximal_TA
  class(ret) = dacomp:::CLASS.LABEL.REFERENCE_SELECTION_OBJECT
  return(ret)
}
