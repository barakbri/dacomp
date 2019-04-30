#class label for results of the wcomp test function
CLASS.LABEL.WCOMP_RESULT_OBJECT = "wcomp.reference.selection.object"

#' Title
#' @param X 
#' @param y 
#' @param ind_reference_taxa 
#' @param z 
#' @param test 
#' @param q 
#' @param nr_perm 
#' @param disable_DSFDR 
#' @param verbose 
#' 
#' @details 
#' 
#' @return
#'   
#' @references foo
#' @author foo
#'
#' @examples
wcomp.test = function(X,y,ind_reference_taxa,z = NULL,test = 'Wilcoxon', q=0.05,return_results_also_on_reference_validation_fail = F, nr_perm = 1/(q/(ncol(X)-length(ind_reference_taxa))),nr_perms_reference_validation = 10^4,T1E_reference_validation = 0.01, disable_DSFDR = F,verbose = F){
  
  #in case a user inserted a reference selection object, with take the indices from the object
  if(class(ind_reference_taxa) == CLASS.LABEL.REFERENCE_SELECTION_OBJECT){
    ind_reference_taxa = ind_reference_taxa$selected_references
  }
  
  #Check input validity
  input_check_result = check.input.wcomp.main(X,y,ind_reference_taxa,test, q,return_results_also_on_reference_validation_fail, nr_perm,nr_perms_reference_validation,T1E_reference_validation, disable_DSFDR,verbose)
  if(!input_check_result)
    stop('Input check failed on wcomp.test')
      
  #definitions and allocations:
  p = ncol(X)
  n = nrow(X)
  min_value_array = rep(NA,p)
  rarefied_data_length_array = rep(NA,p)
  pval_res = rep(NA,p)
  pval_res.CVM = rep(NA,p)
  taxa_nr_res = rep(NA,p)
  reference_values = rep(NA,n)
  lambda_selected = rep(NA,p)
  
  # Compute reference values
  if(length(ind_reference_taxa)>1){
    for(i in 1:n){
      reference_values[i] = sum(X[i,ind_reference_taxa])
    }
  }else{
    reference_values = X[,ind_reference_taxa] 
  }
  
  stats_matrix = matrix(NA, ncol = p, nrow = nr_perm+1)
  rarefaction_matrix = matrix(NA,nrow = n,ncol = 1)
  
  #We compute the permutation matrix. Note that this dependes on the type of test. If it is a signed rank test, permutations are only in each pair of subjects
  Y_matrix = matrix(NA, ncol = nr_perm+1, nrow = n)
  if(test == 'SignedWilcoxon'){
    Y_matrix[,1] = 1:n
    for( i in 2:ncol(Y_matrix)){
      #generate a permutation over signed subjects
      ind = 1:n
      for(t in 1:(n/2)){
        if(runif(1)<=0.5){
          ind[t] = t +(n/2)
          ind[t +(n/2)] = t
        }
      }
      Y_matrix[,i] = ind
    }
  }else{
    Y_matrix[,1] = y
    for( i in 2:ncol(Y_matrix)){
      Y_matrix[,i] = sample(y)
    }
  }

  #iterate over taxa and test
  for(i in 1:p){
    
    #no need to test reference taxa
    if(i %in% ind_reference_taxa){
      next
    }
    
    if(verbose)
      if(i%% ceiling(p/100) == 1)
        cat(paste0('Testing taxon : ',i,'/',p,' \n\r'))
    
    nom = X[,i]
    dnom = reference_values
    
    #choose rarefaction depth:
    total_reads_per_subject = nom+dnom
    
    min_value = min(total_reads_per_subject)
    min_value_array[i] = min_value
    to_keep = 1:length(y) # no filtration of samples allowed, see paper for details on how sample filtration at this point can induce bias
    
    z_keep = z
    
    #enforcing pairs to have two measurements for each subject, subjects with less or more measurements are excluded
    if(!is.null(z)){
      z_keep = z[to_keep]
      counts_for_subject = table(z_keep)
      z_to_keep = as.numeric(which(counts_for_subject == 2))
      z_values_to_keep = as.numeric(names(counts_for_subject)[z_to_keep])
      to_keep = which(z %in% z_values_to_keep) # we already filtered once by minimum value
      z_keep = z[to_keep]
    }
    
    nom_keep = nom[to_keep]
    dnom_keep = dnom[to_keep]
    y_keep = y[to_keep]
    rarefied_data_length_array[i] = length(nom_keep)
    
    # perform the actual subsample,
    nom_keep_original = nom_keep
    dnom_keep_original = dnom_keep
    
    #perform subsample and test
    
    temp_subsampled =  rhyper(n, nom_keep_original, dnom_keep_original,min_value)  
    rarefaction_matrix[,1] = temp_subsampled 
    stats_matrix[,i] = Compute.resample.test(rarefaction_matrix,Y_matrix,statistic = test)
    
  }
  
  #computes pvalues:
  stats = stats_matrix
  p.values = rep(NA,ncol(stats))
  for(i in 1:ncol(stats)){
    p.values[i] = mean(stats[,i]>=stats[1,i])
  }
  
  #compute DS-FDR:
  if(!disable_DSFDR){
    dsfdr_threshold = dsfdr_find_thresholds(stats[,-ind_reference_taxa,drop=F],q,verbose)  
  }
  
  p.values.test = p.values; p.values.test[ind_reference_taxa] = NA
  
  #test reference validity:
  if(verbose)
    cat(paste0('Running test to validate reference set\n\r'))
  
  test.reference.set.validity = wcomp.check_reference_set_is_valid(X_ref = X[,ind_reference_taxa],
                                                             Y = y, nr.perm = nr_perms_reference_validation,
                                                             verbose = verbose)
  # handle a case of signal detected in the reference set:
  possible_problem_in_reference_set_detected = F
  if(any(unlist(test.reference.set.validity)<=T1E_reference_validation)){
    possible_problem_in_reference_set_detected = T
    warning_msg = paste0('Warning: One are more tests for validating that no differentially abundant taxa have entered the reference set has rejected it\'s null hypothesis at level T1E_reference_validation = ',T1E_reference_validation,'. A different reference set of taxon or reference selection method may be considered. To return p.values and rejections toghether with this warning, set parameter \'return_results_also_on_reference_validation_fail\' to TRUE\n\r')
    warning(warning_msg)
  }
  
  #return results:
  
  ret = list()
  ret$test.reference.set.validity = test.reference.set.validity
  if(possible_problem_in_reference_set_detected){ #put warning if needed:
    ret$warning_msg = warning_msg
  }
  
  if(return_results_also_on_reference_validation_fail | !possible_problem_in_reference_set_detected){ #if reference validation procedure failed, dont put result. Only if the user specifically asked
    ret$lambda = min_value_array
    ret$stats_matrix = stats_matrix
    ret$p.values.test = p.values.test
    
    if(!disable_DSFDR){
      ret$rejected = which(p.values.test<=dsfdr_threshold)
      ret$dsfdr_threshold = dsfdr_threshold  
    }
  }
  class(ret) = CLASS.LABEL.WCOMP_RESULT_OBJECT
  return(ret)
}

#internal function for validating inputs on wcomp.test
check.input.wcomp.main = function(X, y, ind_reference_taxa, test, q, return_results_also_on_reference_validation_fail, nr_perm,nr_perms_reference_validation, T1E_reference_validation, disable_DSFDR, verbose){
  return(TRUE)
}

