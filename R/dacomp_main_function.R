#class label for results of the dacomp test function
CLASS.LABEL.DACOMP_RESULT_OBJECT = "dacomp.result.object"

#' Title
#'
#' @param X 
#' @param y 
#' @param ind_reference_taxa 
#' @param test 
#' @param q 
#' @param nr_perm 
#' @param disable_DSFDR 
#' @param verbose 
#'
#' @return
#' * lambda
#' * stats_matrix
#' * p.values.test
#' * dsfdr_rejected
#' * dsfdr_threshold
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' library(dacomp)
#' 
#' set.seed(1)
#' 
#' data = dacomp.generate_example_dataset(m1 = 100,
#'        n_X = 50,
#'        n_Y = 50,
#'        signal_strength_as_change_in_microbial_load = 0.1)
#' 
#' #select references: (may take a minute)
#' result.selected.references = dacomp.select_references(X = data$counts,
#'                                                      median_SD_threshold = 0.6, #APPLICATION SPECIFIC
#'                                                      verbose = T)
#' 
#' length(result.selected.references$selected_references)
#'
#' #plot the reference selection scores (can also be used to better set the median SD threshold)
#' dacomp.plot_reference_scores(result.selected.references)
#' 
#' #multiplicity correction levels for the BH and DS-FDR methods
#' q_BH = q_DSFDR = 0.1
#' 
#' #Perform testing:
#' result.test = dacomp.test(X = data$counts,
#'                       y = data$group_labels,
#'                       ind_reference_taxa = result.selected.references,
#'                       verbose = T,q = q_DSFDR)
#' 
#' rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)
#' rejected_DSFDR = result.test$dsfdr_rejected
#' 
#' # example with continous covariate
#' set.seed(1)
#'
#' data = dacomp.generate_example_dataset_continuous(n = 100,m1 = 30,
#' signal_strength_as_change_in_microbial_load = 0.1)
#'
#'
#' result.selected.references = dacomp.select_references(X = data$counts,
#'                                                      median_SD_threshold = 0.6, #APPLICATION SPECIFIC
#'                                                      verbose = T)
#' #number of selected references
#' length(result.selected.references$selected_references)
#'
#' #plot the reference selection scores (can also be used to better set the median SD threshold)
#' dacomp.plot_reference_scores(result.selected.references)
#'
#' #multiplicity correction levels for the BH and DS-FDR methods
#' q_BH = q_DSFDR = 0.1
#' 
#' #Perform testing:
#' result.test = dacomp.test(X = data$counts,
#'                          y = data$covariate,test = DACOMP.TEST.NAME.SPEARMAN,
#'                          ind_reference_taxa = result.selected.references,
#'                          verbose = T,q = q_DSFDR)
#'
#' rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)
#' rejected_DSFDR = result.test$dsfdr_rejected
#' }
dacomp.test = function(X,y,ind_reference_taxa,test = DACOMP.TEST.NAME.WILCOXON, q=0.05, nr_perm = 1/(q/(ncol(X)-length(ind_reference_taxa))), disable_DSFDR = F,user_defined_test_function = NULL, verbose = F ){
  
  #Preprocess inputs, before check:
  
  #in case a user inserted a reference selection object, with take the indices from the object
  if(class(ind_reference_taxa) == CLASS.LABEL.REFERENCE_SELECTION_OBJECT){
    ind_reference_taxa = ind_reference_taxa$selected_references
  }
  
  if(is.numeric(nr_perm))
    nr_perm = ceiling(nr_perm)
  
  #if test is in TEST.DEF.Y.IS.0.OR.1, convert Y to 0 and 1
  if(test %in% TEST.DEF.Y.IS.0.OR.1){
    y = as.numeric( as.factor(y) ) - 1
  }
     
  #Check input validity
  input_check_result = check.input.main(X,y,ind_reference_taxa,test, q, nr_perm, disable_DSFDR,user_defined_test_function,verbose)
  if(!input_check_result)
    stop('Input check failed on dacomp.test')
      
  #definitions and allocations:
  p = ncol(X)
  n = nrow(X)
  min_value_array = rep(NA,p)
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
  if(test %in% TEST.DEF.TESTS.ON.PAIRS){
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
  }else if (test %in% TEST.DEF.TESTS.ON.UNIVARIATE_CONTINOUS){
    Y_matrix[,1] = rank(y,ties.method = 'average')
    Y_matrix[,1] = Y_matrix[,1] - mean(Y_matrix[,1])
    for( i in 2:ncol(Y_matrix)){
      Y_matrix[,i] = sample(Y_matrix[,1])
    }
    
  }else if (test == DACOMP.TEST.NAME.USER_DEFINED){
     #user must plan and perform permutation independently.
  }else{
    Y_matrix[,1] = y
    for( i in 2:ncol(Y_matrix)){
      Y_matrix[,i] = sample(y)
    }
  }
  
  #iterate over taxa and test
  for(i in 1:p){
    
    if(verbose)
      if(i%% ceiling(p/10) == 1)
        cat(paste0('Testing taxon : ',i,'/',p,' \n\r'))
    
    #no need to test reference taxa
    if(i %in% ind_reference_taxa){
      next
    }
    
    nom = X[,i]
    dnom = reference_values
    
    #choose rarefaction depth:
    total_reads_per_subject = nom+dnom
    
    min_value = min(total_reads_per_subject)
    min_value_array[i] = min_value

    #perform subsample and test
    
    temp_subsampled =  rhyper(n, nom, dnom,min_value)  
    rarefaction_matrix[,1] = temp_subsampled 
    stats_matrix[,i] = Compute.resample.test(rarefaction_matrix,Y_matrix,statistic = test, user_defined_test_function = user_defined_test_function)
    
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
  
  #return results:
  
  ret = list()
  ret$lambda = min_value_array
  ret$stats_matrix = stats_matrix
  ret$p.values.test = p.values.test
  if(!disable_DSFDR){
    ret$dsfdr_rejected = which(p.values.test<=dsfdr_threshold)
    ret$dsfdr_threshold = dsfdr_threshold  
  }
  class(ret) = CLASS.LABEL.DACOMP_RESULT_OBJECT
  return(ret)
}

#internal function for validating inputs on dacomp.test
check.input.main = function(X, y, ind_reference_taxa, test, q, nr_perm, disable_DSFDR,user_defined_test_function, verbose){
  
  ##  dacomp.test: test X is numeric matrix of counts, 
  MSG_X = 'X must be a valid counts matrix'
  if(!is.matrix(X))
    stop(MSG_X)
  if(any(X!=as.integer(X)))
    stop(MSG_X)
  if(any(X<0))
    stop(MSG_X)
  
  if(!(test%in% dacomp::DACOMP.POSSIBLE.TEST.NAMES)){
    stop(paste0('Test ',test,' not part of list DACOMP.POSSIBLE.TEST.NAMES'))
  }
  
  ##test y is valid - 0 or 1 for the relevant tests, has two groups are more for KW, disregarded in signed rank test
  if(!(test %in% c(TEST.DEF.TESTS.ON.PAIRS,DACOMP.TEST.NAME.USER_DEFINED))){
    if(length(y)!= nrow(X))
      stop('length of y must be same as number of rows in X_ref')
    if(any(is.na(y))|any(is.nan(y)))
      stop('y has NA or NaNs - invalid observations')
    if(!(test %in% TEST.DEF.TESTS.ON.UNIVARIATE_CONTINOUS))
      if(min(table(y))<5)
        warning('Note: at least one sample group with less than 5 observations\n\r')
  }
  
  if(test %in% TEST.DEF.Y.IS.0.OR.1){
    MSG_two_groups_not_available = paste0('y has more than two levels - cannot use test = ',test)
    if(length(unique(y))!=2)
      stop(MSG_two_groups_not_available)
    if(!all.equal(sort(unique(y)),c(0,1))){
      stop(MSG_two_groups_not_available)
    }
  }
  if(test %in% TEST.DEF.TESTS.OVER.GROUPS){
    if(length(unique(y))<2)
      stop('y has only one level, cannot use test for groups')
  }
  
  #Y is null for signed rank test
  if(test %in% TEST.DEF.TESTS.ON.PAIRS){
    if(!is.null(y)){
      stop('For paired data, y has to be set to NULL, and the rows of X set such that the first N rows ')
    }
    if(nrow(X) %%2 ==1)
      stop('For paired data, X should have an even number of rows, but it has an odd number')
  }
  
  ##q is valid
  MSG_q = "q sets the required FDR level for the DS-FDR algorithm. Must be a number, >0"
  if(!is.numeric(q))
    stop(MSG_q)
  if(q <= 0)
    stop(MSG_q)
  
  if(q > 0.1)
    warning('Note: you have set q for DS-FDR to be greater than 0.1, irregular parameter setting.\n\r')
  
  ##nr_perm is valid,
  MSG_NR_PERM = 'nr.perm must be at integer, at least 1000'
  if(nr_perm != as.integer(nr_perm))
    stop(MSG_NR_PERM)
  if(nr_perm<1000)
    stop(MSG_NR_PERM)
  
  ##disable_DSFDR is valid 
  if(!is.logical(disable_DSFDR))
    stop('disable_DSFDR must be logical')
  
  ##verbose is valid
  if(!is.logical(verbose))
    stop('Verbose must be logical')
  
  if(any(!(ind_reference_taxa %in% 1:ncol(X))))
    stop('ind_reference_taxa must be subset of 1:ncol(X)')
  
  if(test == DACOMP.TEST.NAME.USER_DEFINED){
    if(typeof(user_defined_test_function) !=  "closure")
      stop('For user defined test, argument user_defined_test_function must be function returning P.value')
  }
  
  return(TRUE)
}

