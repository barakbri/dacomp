
leave_one_out_validation = function(X,Y,ref_obj_for_validation,test = dacomp::DACOMP.TEST.NAME.WILCOXON, Q = 0.1,NR_PERMS = 1000,Verbose = F,disable_DSFDR = T){
  if(class(ref_obj_for_validation) == dacomp:::CLASS.LABEL.REFERENCE_SELECTION_OBJECT){
    X_ref_for_validation = X[,ref_obj_for_validation$selected_references]
  }else{
    X_ref_for_validation = X[,ref_obj_for_validation]
  }
  
  
  temp_res = dacomp.test(X = X_ref_for_validation,y = Y,
                                 ind_reference_taxa = (1:ncol(X_ref_for_validation)),
                                 test = test,
                                 q = Q,
                                 nr_perm = NR_PERMS,
                                 disable_DSFDR = disable_DSFDR,
                                 compute_ratio_normalization = T,
                                 Test_All = T,
                                 verbose = F,
                                 return_rarefied_values = T)
  
  return(list(pval.vec = temp_res$p.values.test,pval.vec.ratio = temp_res$p.values.test.ratio.normalization,dacomp_res = temp_res))
}

select_references_by_n = function(ref_obj,maximal_TA,minimal_TA,select_from,X){
  ref_obj_to_return = dacomp.select_references(X = X, minimal_TA = minimal_TA,
                                               median_SD_threshold = 1E-6,maximal_TA = maximal_TA,
                                               Pseudo_Count_used = 1,verbose = F,
                                               select_from = select_from,
                                               Previous_Reference_Selection_Object = ref_obj)
  return(ref_obj_to_return$selected_references)
  #nr_ref_taxa_to_take = length(which(ref_obj$min_abundance_over_the_sorted<=n))+1
  #nr_ref_taxa_to_take = min(nr_ref_taxa_to_take,length(ref_obj$scores))
  #threshold_scrit = (sort(ref_obj$scores))[nr_ref_taxa_to_take]
  #return(which(ref_obj$scores<=threshold_scrit))
}


#' Diagnostic procedure for DACOMP reference set.
#' 
#' The function receives a data set (through \code{X} and \code{Y}), and a chosen reference set for the data (in \code{ref_obj}), evaluates if the reference set contains differentially abundant taxa, and if so, reselect the reference set.
#' 
#' @details At each iteration of the algorithm, all taxa in the reference set are tested for differential abundance in a leave one-out-manner: they are exluded from the reference set, and tested for differential abundance w.r.t to the remaining taxa in the reference set. The parameter \code{test} determines the test performed, and should be selected according to the values of \code{Y}, see additional details in \code{\link{dacomp.test}}. The \eqn{P}-values for \eqn{r} reference taxa are combined using the Simes test statistic (see 1986 paper): \eqn{p_{Simes} = min_{j=1, ..., r}\frac{r\times p_{(j)}}{j}}, where \eqn{p_{(j)}} are the ordered \eqn{P}-values. If \eqn{p_{Simes}} is smaller than \code{Q_validation}, we reselect the reference set set as follows.  If the current reference set has at least \eqn{\lambda_{min}} counts in all samples, the new reference set selected is the smallest set with at least \eqn{Reduction\_Factor\cdot \lambda_{min}} counts. Taxa are inserted into the reference set based on their reference scores (as given by the object \code{ref_obj}), until the target minimal abundance is reached.  The algorithm iterates over testing and reselection steps, until either \eqn{p_{Simes}} is larger than \code{Q_validation}, or the algorithm cannot select a reference set with at least \code{Minimal_Counts_in_ref_threshold} counts under the reference taxa for all samples.
#'
#' Setting the parameter \code{disable_DSFDR} to \code{F} uses an experimental feature, where combined of \eqn{P}-values is done using the DS-FDR multiple testing procedure, rather than the Simes test. 
#'
#' @param X matrix of 16S counts for data, rows are samples, columns are taxa. 
#' @param Y Vector of trait values. Entries in this vector should correspond to the rows of \code{X}
#' @param ref_obj A reference object, returned from \code{\link{dacomp.select_references}}.
#' @param test The type of the test to be used, passed to \code{\link{dacomp.test}} when testing reference taxa.
#' @param Q_validation The level of the Simes test (\eqn{\alpha}) used to determine if reference taxa show a signal.
#' @param Minimal_Counts_in_ref_threshold When shrinking the reference, the reference will not be shrinked past a point at least \code{Minimal_Counts_in_ref_threshold} counts are available in the reference set for each sample.
#' @param Reduction_Factor When a contamination is found in the reference set, a new reference set is selected with the following rule. If the tested reference set has at least \eqn{\lambda_{min}} counts in all samples, the new reference set selected is the smallest set with at least \eqn{Reduction_Factor \ctimes \lambda_{min}} counts. When constructing the new reference set, taxa are inserted into the reference set based on their reference scores (as given by the object \code{ref_obj}), until the target minimal abundance is reached. This reselection happens at each iteration of the algorithm.
#' @param NR_perm Number of permutations used for computing the \eqn{P}-values for the Simes intersection test. 
#' @param select_from Can be used to limit the set of taxa that the reference set is selected from. Receives a numeric vector detailing which taxa are valid. The default value of \code{NULL} means no taxa are excluded from the selection procedure.
#' @param Verbose Logical value indicating if messages should be printed. Messages include: the iteration of the algorithm; if a contamination has been found using the Simes test at each iteration; the number of taxa remaining in the reference set after reselection; and a message describing if the final reference included no signal, or iterations were halted due to the parameter \code{Minimal_Counts_in_ref_threshold}.
#' @param disable_DSFDR Logical value indicating if testing should be done using the Simes procedure (\code{True}, default value) or the DS-FDR procedure (\code{False}, experimental feature).
#' @return A numeric vector, containing the the indices of the taxa selected as the reference set by the function ( see description of algorithm under details). Indices are based on the columns of \code{X}. See example on how to use with testing procedure under \code{\link{dacomp.test}}.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' library(dacomp)
#' #generate data with two study groups
#' data = dacomp.generate_example_dataset.two_sample(n_X = 30,n_Y = 30,m1 = 50,signal_strength_as_change_in_microbial_load = 0.1)
#' 
#' # select references. We purposely select reference taxa so that differentially abundant taxa enter the reference set. In general, select using median_SD_threshold=0, minimal_TA=100, see paper for discussion of this selection strategy
#' result.selected.references = dacomp.select_references(X = data$counts,
#'                                                      median_SD_threshold = 1.3,
#'                                                      maximal_TA = 1000,
#'                                                      verbose = T)
#'                                                      
#' # some differentially abundant taxa entered the reference set:                                                    
#' sum(result.selected.references$selected_references %in% data$select_diff_abundant)
#' 
#' #run the sensitivity analysis.
#' cleaned_references = dacomp.validate_references(X =  data$counts,
#'                                                Y =  data$group_labels,
#'                                                ref_obj = result.selected.references,
#'                                                test =DACOMP.TEST.NAME.WILCOXON,
#'                                                Q_validation = 0.1,
#'                                                Minimal_Counts_in_ref_threshold = 50,
#'                                                Reduction_Factor = 0.5,
#'                                                Verbose = T,
#'                                                disable_DSFDR = T,
#'                                                NR_perm = 1000)
#'                                                
#' #now the reduced reference has no differentially abundant taxa inside....
#' sum(cleaned_references %in% data$select_diff_abundant)
#' 
#' # when testing, report results both for the original reference set selected, as well as the reduced reference.
#' }
#'  
dacomp.validate_references = function(X,Y,
                               ref_obj,
                               test,
                               Q_validation = 0.1,
                               Minimal_Counts_in_ref_threshold = 10,
                               Reduction_Factor = 0.9,
                               NR_perm = ceiling(1/(Q_validation/ncol(X))),
                               select_from = NULL,
                               Verbose = T,
                               disable_DSFDR = F
){
  if(!class(ref_obj) == dacomp:::CLASS.LABEL.REFERENCE_SELECTION_OBJECT){
    stop('ref_obj in dacomp.validate_references must be a valid reference selection object') 
  }
  if(length(ref_obj$selected_references)==1)
    return(ref_obj$selected_references)
  
  
  
  minimal_possible_reference_set = select_references_by_n(ref_obj = ref_obj,
                                                          minimal_TA = Minimal_Counts_in_ref_threshold,
                                                          maximal_TA  = Minimal_Counts_in_ref_threshold*10, #meaningless, we select by minimal_TA
                                                          select_from = select_from,
                                                          X = X)
  
  L1O_res = leave_one_out_validation(X =  X,
                                     Y =  Y,
                                     ref_obj_for_validation = ref_obj,
                                     test = test,
                                     Q = Q_validation,
                                     NR_PERMS = NR_perm,
                                     Verbose = F,
                                     disable_DSFDR = disable_DSFDR)
  nr_detected_as_contamination = which(L1O_res$dacomp_res$p.values.test.adjusted<=Q_validation)
  counts_in_ref_target = ref_obj$selected_MinAbundance
  Returned_references = ref_obj$selected_references
  
  #max_iter_debug = ceiling(log(Minimal_Counts_in_ref_threshold/counts_in_ref_target,
  #                             base = Reduction_Factor))
  
  check_identical = function(x,y){
    if(length(x) != length(y))
      return(F)
    return(all.equal(sort(as.numeric(x)),sort(as.numeric(y))))
  }
  
  while(length(nr_detected_as_contamination)>0 &
        counts_in_ref_target>=Minimal_Counts_in_ref_threshold){
    if(Verbose){
      print('###############')
      if(length(nr_detected_as_contamination) >0 ){
        print(paste0('Identified reference taxa with signal.'))    
      }else{
        print(paste0('No signal identified in reference taxa.'))
      }
    
    }
    new_references = Returned_references
    while(check_identical(new_references,Returned_references)){
      new_counts_in_ref_target = ceiling(counts_in_ref_target* Reduction_Factor)
      if(Verbose){
        print(paste0('Updating target rarefaction depth from ',counts_in_ref_target,' to ',new_counts_in_ref_target))
      }
      counts_in_ref_target = new_counts_in_ref_target
      
      new_references = select_references_by_n(ref_obj = ref_obj,
                                              minimal_TA = counts_in_ref_target,
                                              maximal_TA  = counts_in_ref_target*10, #meaningless, we select by minimal_TA
                                              select_from = select_from,
                                              X = X)
      if(Verbose){
        print(paste0('Length of reselected reference is ',length(new_references)))
        if(check_identical(new_references,Returned_references)){
          print(paste0('Same reference, continuing to lower target number of counts'))
        }
      }
      
      if(length(new_references)==0 |
         min(apply(X[,new_references,drop=F],1,sum)) < Minimal_Counts_in_ref_threshold |
         check_identical(new_references,minimal_possible_reference_set) |
         length(new_references) <= length(minimal_possible_reference_set)){
         S = 'Warning:Cannot decrease anymore. Returning reference that may be contaminated!'
         warning(S)
         return(Returned_references)
      }
    }
    
    if(Verbose){
      print(paste0('Retesting'))
    }
    L1O_res = leave_one_out_validation(X =  X,
                                       Y =  Y,
                                       ref_obj_for_validation = new_references,
                                       test = test,
                                       Q = Q_validation,
                                       NR_PERMS = NR_perm,
                                       Verbose = F,
                                       disable_DSFDR = disable_DSFDR)
    
    nr_detected_as_contamination = which(L1O_res$dacomp_res$p.values.test.adjusted<=Q_validation)
    if(Verbose){
      if(length(nr_detected_as_contamination) >0 ){
        print(paste0('Identified reference taxa with signal.'))    
      }else{
        print(paste0('No signal identified in reference taxa.'))
      }
    }
    Returned_references = new_references
  }
  return(Returned_references)
}


