#class label for results of the dacomp test function
CLASS.LABEL.DACOMP_RESULT_OBJECT = "dacomp.result.object"

#' Test taxa for differential abundance, using a given set of reference taxa used for normalization.
#' 
#' The function tests taxa for differential abundance, given a set of reference taxa used for normalization, as described in Brill et. al. (2019). The function supports several tests, by type of phenotype: two sample tests, K-sample tests, tests for association with a continuous phenotype and user defined tests. See `details` below on how the input argument \code{y} shold be formatted for different tests.
#'
#' @details
#' The function tests each taxon not in the reference set for differential abundance as follows. For each taxon, the following procedure is performed. First, an identical number of reads is taken from each sample, from the reads available under the reference set of taxa, and the taxon being tested. Next, a test of association is performed between the rarefied reads and the phenotype given by \code{y}. P-values are computed by permutations. Finally, after all P-values are computed, the DS-FDR threshold for rejection is computed. Hypotheses with P-values lower than the DS-FDR rejection threshold are rejected. Reference taxa are not tested for differential abundance. See \code{vignette('dacomp_main_vignette')} for additional details and examples.
#'
#' The function supports several tests. Different tests require different formats for the arguments \code{X} and \code{y}.
#' 
#' \itemize{
#' 
#' Two Sample Tests (treatment vs. control, healthy vs.sick): 
#' 
#' \item{DACOMP.TEST.NAME.WILCOXON}{ - (Preferred option for 2-sample testing) A Wilcoxon rank sum test between the rarefied reads from the two sample groups. \code{y} is any vector with two levels. For two sample testing this is the preferred test.}
#' \item{DACOMP.TEST.NAME.DIFFERENCE_IN_MEANS}{ - A two sample test for difference in means, based on permutation. \code{y} formated as for the Wilcoxon rank sum test.}
#' \item{DACOMP.TEST.NAME.LOG_FOLD_DIFFERENCE_IN_MEANS}{ - A two sample test for difference in means of the logatihm of counts, based on permutation. A pseudocount of 1 is added to the rarefied reads, so that the lograithm can be computed. \code{y} formated as for the Wilcoxon rank sum test.}
#' \item{DACOMP.TEST.NAME.TWO_PART_WILCOXON}{ - The two part test of Wagner et. al. (2011), for equality of distributions between two sample groups. \code{y} formated as for the Wilcoxon rank sum test.}
#' \item{DACOMP.TEST.NAME.WELCH}{ - Welch t.test over the rarefied counts. P-values are computed by permutations. \code{y} formated as for the Wilcoxon rank sum test.}
#' \item{DACOMP.TEST.NAME.WELCH_LOGSCALE}{ - Welch t.test over the logarithm of rarefied counts (after a pseudocount of 1 is added). P-values are computed by permutations. \code{y} formated as for the Wilcoxon rank sum test.}
#' 
#' Paired Design Tests:
#' 
#' \item{DACOMP.TEST.NAME.WILCOXON_SIGNED_RANK_TEST}{ - The Wilcoxon sign rank test for paired designs. For this test, \code{X} is formatted with the first \eqn{n/2} rows corresponding to the samples measured at condition 1, and the latter \eqn{n/2} rows measured at condition 2, i.e. rows \eqn{1} and \eqn{n/2 + 1} correspond to the same physical sample, under two conditions. For this test, \code{y} is set to \code{NULL}.}
#' 
#' Tests of Association for Ordinal/Continuous Phenotypes:
#' 
#' \item{DACOMP.TEST.NAME.SPEARMAN}{ - Test of association with a continuous, univariate phenotype. \code{y} is the measured phenotype across samples. The test performed is a permutation based test for the Spearman correlation coefficient with a two sided alternative.}
#' 
#' K-Sample Tests:
#' 
#' \item{DACOMP.TEST.NAME.KRUSKAL_WALLIS}{ - Test for equality of distributions between \eqn{K} sample groups. \code{y} should be contian the group labeling of different observations.}
#' 
#' User Defined Tests:
#' 
#' \item{DACOMP.TEST.NAME.USER_DEFINED}{ - Indicates that a custom test is supplied using the argument \code{user_defined_test_function}. The supplied function will receive a single argument, the vector of rarefied counts and will return an array of length \code{nr_perm +1}, containing the test statistic computed for the original data, along with test statistics computed for permuted phenotypes. Test statistics must have a right sided alternative. A complete example with code snippets is found in \code{vignette('dacomp_main_vignette')}.}
#' }
#' 
#' The parameter \code{compute_ratio_normalization} computes in addition to the test described above, a test of association between the phenotype \code{y} and the proportion of a tested taxon out of a subvector containing only the tested taxon and a set of reference taxa. This procedure, described as DACOMP-ratio in Brill et al. (2019), may enjoy higher power to detect differentially abundant taxa with a low number of counts, due to avoiding rarefaction. However, for scenarios with an extreme change in the microbial load of the measured ecology, e.g., a four-fold change in the microbial load between different study groups, this variant of DACOMP may have a slightly inflated FDR. See method description and simulation results in the paper for additional details.
#' 
#' We suggest using this function toghether with the function \code{\link{dacomp.validate_references}}. The function \code{\link{dacomp.validate_references}} is meant to detect if differentially abundant taxa have entered the reference set, and if so, reselect the reference set. We suggest that results be reported twice: once using the original reference selected, and once using the re-selected reference set. If results agree, than it is more likely that the initial reference set selected did not feature a contamination. If results show disagreement, this could be due to a signal entering the reference set. 
#' 
#' @param X Matrix of counts, with rows representing samples, columns representing taxa. See `details` for additional information on how to format this matrix for paired study designs.
#' @param y Vector of phenotype values by sample. See details on different formats used by different tests.
#' @param ind_reference_taxa One of two options: Object of type 'dacomp.reference.selection.object' returned from \code{\link{dacomp.select_references}}, or a vector of indices of taxa selected as a reference set for normalization. See package vignette for additional details.
#' @param test One of the values in the vector \code{DACOMP.POSSIBLE.TEST.NAMES}, or one of the constants available as \code{DACOMP.TEST.NAME.*} with asterisk representing the name of the test. See 'details' for additional details and inputs required, by test.
#' @param q The required FDR level for the DS-FDR algorithm, see Jiang et. al. (2017) for details.
#' @param nr_perm Number of permuations used for testing and computing P-values. The default number of permutations is set high enough to provide powerful inference after adjusting for multiplicity. Change the number of permutations only if you know how it effects power after correcting for multiplicity. 
#' @param disable_DSFDR Can be used to disable the DS-FDR computation (which may take up to a minute on large datasets). The default value is \code{FALSE}.
#' @param user_defined_test_function Argument for inputing a function for a custom test of association supplied by the user, see `details` below.
#' @param compute_ratio_normalization Argument determines if in addition to the test described above, the function should a run test based on normalization by division rather than rarefaction. See description under details.
#' @param verbose Should messages be printed to console, indicating computation progress. The default value is \code{FALSE}.
#' @param DSFDR_Filter A logical (boolean) vector specifying for which taxa should the DSFDR adjusted P-values be computed. Taxa whose corresponding entries in this vector are set to \code{} are excluded from testing, and will not have DSFDR adjusted Pvalues computed. If a taxon is set to be in the reference set, and its corresponding entry is set to \code{T} in this vector, it will still be excluded from testing. 
#' @param Test_All A logical value specifying if all test (including the reference taxa) should be tested for differential abundance. When testing a reference taxon for differential abundance, it is excluded from the reference set.
#' @param return_rarefied_values a logical value indicating if the rarefied draws used for testing should be returned as a matrix under the entry \code{rarefied_counts} in the returned object.
#' @return An object of type "dacomp.result.object", which is a list with the follow fields:
#' \itemize{
#' \item{lambda}{ - The subsampling depth for each tested taxon, as in step I under `details`.}
#' \item{stats_matrix}{ - A matrix of size \code{(nr_perm+1) X ncol(X)} containing the tests statistics for the data (first row) and test statistics computed for permuted values of \code{y} (other rows).}
#' \item{p.values.test}{ - A vector with P-values for the different tests of association, by taxa. P-values obtained by permutations. P-values for reference taxa will appear as \code{NA}.}
#' \item{p.values.test.adjusted}{ - A vector of P-values, adjusted for multiplicity, corresponding to \code{p.values.test}. By default, correction is done by the DS-FDR method. If \code{disable_DSFDR} is set to \code{TRUE}, the BH correction is performed. Entries lower than a value of \eqn{u} indicate a taxon declared differentially abundant when trying to control multiplicity at level \eqn{u}. Entries lower than the value defined for the parameter \code{q} will be indicated as discoveries under \code{dsfdr_rejected}}
#' \item{effect_size_estimates}{  - For correlation tests (Spearman), will give the spearman correlation coeffcient between rarefied counts and Y. For K-sample, 2-Sample and paired tests, will provide a string describing the ordering of mean ranks of rarefied counts, across levels of Y. }
#' \item{dsfdr_rejected}{ - A vector of taxa indices declared differentially abundant by the DS-FDR method for multiplicity adjustment. This field will not be available if \code{disable_DSFDR} is set to \code{TRUE}.}
#' \item{dsfdr_threshold}{ - The selected threshold, in terms of P-values, for declaring taxa as differentialy abundant. Taxa with P-values under this threshold will be declared diffentially abundant. This field will not be available if \code{disable_DSFDR} is set to \code{TRUE}.}
#' 
#' \item{p.values.test.ratio.normalization}{ - A vector with P-values for the different tests of association, by taxa, for the "normalization by ratio" variant of DACOMP.  This field will be available only if  \code{compute_ratio_normalization} is set to \code{TRUE}.}
#' \item{p.values.test.adjusted.ratio.normalization}{ - A vector of P-values, adjusted for multiplicity, for the tests performed using ratio normalization (given by \code{p.values.test.ratio.normalization}). By default, correction is done by the DS-FDR method. If \code{disable_DSFDR} is set to \code{TRUE}, the BH correction is performed. Entries lower than a value of \eqn{u} indicate a taxon declared differentially abundant when trying to control multiplicity at level \eqn{u}. Entries lower than the value defined for the parameter \code{q} will be indicated as discoveries under \code{dsfdr_rejected_ratio_normalization}}
#' \item{effect_size_estimates_ratio}{ - Similar to \code{effect_size_estimates}, but computed over the ratio between taxon counts and total number of reads available under both the reference taxa and the tested taxon.}
#' \item{dsfdr_rejected_ratio_normalization}{ -A vector of taxa indices declared differentially abundant by the DS-FDR method, similar to dsfdr_rejected, but using the P-values obtained for the "normalization by ratio" test. This field will be available only if  \code{compute_ratio_normalization} is set to \code{TRUE} and DS-FDR computation is not disabled.}
#' \item{dsfdr_threshold_ratio_normalization}{ - The selected threshold, in terms of P-values, for declaring taxa as differentialy abundant. Taxa with P-values, obtained with "normalization by ratio" type tests, under this threshold will be declared diffentially abundant. This field will be available only if  \code{compute_ratio_normalization} is set to \code{TRUE} and DS-FDR computation is not disabled.}
#' \item{rarefied_counts}{- The rarefied counts for each taxon-wise test, rows are samples, columns are taxa.}
#' }
#' 
#' @references 
#' Brill, Barak, Amnon Amir, and Ruth Heller. 2019. Testing for Differential Abundance in Compositional Counts Data, with Application to Microbiome Studies. arXiv Preprint arXiv:1904.08937.
#' 
#' Jiang, L, A Amir, JT Morton, R Heller, E Arias-Castro, and R Knight. 2017. Discrete False-Discovery Rate Improves Identification of Differentially Abundant Microbes. mSystems 2: E00092-17. Am Soc Microbiol.
#' 
#' Wagner, Brandie D, Charles E Robertson, and J Kirk Harris. 2011. Application of Two-Part Statistics for Comparison of Sequence Variant Counts. PloS One 6 (5). Public Library of Science: e20296.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' 
#' library(dacomp)
#' 
#' set.seed(1)
#' 
#' # example for a study with two groups:
#' 
#' data = dacomp.generate_example_dataset.two_sample(m1 = 100,
#'        n_X = 50,
#'        n_Y = 50,
#'        signal_strength_as_change_in_microbial_load = 0.1)
#' 
#' #select references: (may take a minute)
#' result.selected.references = dacomp.select_references(X = data$counts,
#'                                                      minimal_TA = 50, #Choosing the minimal number of reference taxa so that at least 50 reads are available under the reference for all samples
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
#'                       test = DACOMP.TEST.NAME.WILCOXON,nr_perm = 1000,
#'                       verbose = T,q = q_DSFDR)
#' 
#' rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)
#' rejected_DSFDR = result.test$dsfdr_rejected
#' 
#' #We note that our reference set may be contaminated, i.e. include some differentially abundant taxa. Therefore, we run a diagnostic check, trying to remove possible signals from the reference, and retesting:
#' Cleaned_references = dacomp.validate_references(X =  data$counts,
#'                                                Y =  data$group_labels,
#'                                                ref_obj = result.selected.references,
#'                                                test =DACOMP.TEST.NAME.WILCOXON,
#'                                                Q_validation = 0.1,
#'                                                Minimal_Counts_in_ref_threshold = 10,
#'                                                Reduction_Factor = 0.9,
#'                                                Verbose = T,
#'                                                disable_DSFDR = T,
#'                                                NR_perm = 10000)
#'                                                 
#' result.test.with.reduced.reference = dacomp.test(X = data$counts,
#'                        y = data$group_labels,
#'                       ind_reference_taxa = result.selected.references,
#'                       test = DACOMP.TEST.NAME.WILCOXON,nr_perm = 1000,
#'                       verbose = T,q = q_DSFDR)
#' 
#' rejected_BH = which(p.adjust(result.test.with.reduced.reference$p.values.test,method = 'BH')<=q_BH)
#' rejected_DSFDR = result.test.with.reduced.reference$dsfdr_rejected
#' 
#' # example with continous outcome
#' set.seed(1)
#'
#' data = dacomp.generate_example_dataset_continuous(n = 100,m1 = 30,
#' signal_strength_as_change_in_microbial_load = 0.2)
#'
#'
#' result.selected.references = dacomp.select_references(X = data$counts,
#'                                                      minimal_TA = 50, #Choosing the minimal number of reference taxa so that at least 50 reads are available under the reference for all samples
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
#'                          ind_reference_taxa = result.selected.references,nr_perm = 1000,
#'                          verbose = T,q = q_DSFDR)
#'
#' rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)
#' rejected_DSFDR = result.test$dsfdr_rejected
#' 
#' #We note that our reference set may be contaminated, i.e. include some differentially abundant taxa. Therefore, we run a diagnostic check, trying to remove possible signals from the reference, and retesting:
#' Cleaned_references = dacomp.validate_references(X =  data$counts,
#'                                                Y =  data$covariate,
#'                                                ref_obj = result.selected.references,
#'                                                test = DACOMP.TEST.NAME.SPEARMAN,
#'                                                Q_validation = 0.1,
#'                                                Minimal_Counts_in_ref_threshold = 10,
#'                                                Reduction_Factor = 0.9,
#'                                                Verbose = T,
#'                                                disable_DSFDR = T,
#'                                                NR_perm = 1000)
#' result.test.with.reduced.reference = dacomp.test(X = data$counts,
#'                                                 y = data$covariate,
#'                                                 ind_reference_taxa = Cleaned_references,
#'                                                 test = DACOMP.TEST.NAME.SPEARMAN,nr_perm = 1000,
#'                                                 verbose = T,q = q_DSFDR)
#' 
#' rejected_BH = which(p.adjust(result.test.with.reduced.reference$p.values.test,method = 'BH')<=q_BH)
#' rejected_DSFDR = result.test.with.reduced.reference$dsfdr_rejected
#' 
#' 
#' }
dacomp.test = function(X,y,ind_reference_taxa,test, q=0.05, nr_perm = max(1/(q/(ncol(X))),1000), disable_DSFDR = F,user_defined_test_function = NULL, compute_ratio_normalization = F, verbose = F,DSFDR_Filter = rep(T,ncol(X)),Test_All = F,return_rarefied_values = F ){
  
  #Preprocess inputs, before check:
  y_original = y #we keep a copy for k-sample effect size estimates
  #in case a user inserted a reference selection object, with take the indices from the object
  if(class(ind_reference_taxa) == CLASS.LABEL.REFERENCE_SELECTION_OBJECT){
    ind_reference_taxa = ind_reference_taxa$selected_references
  }
  
  if(is.numeric(nr_perm))
    nr_perm = ceiling(nr_perm)
  
  #if test is in TEST.DEF.Y.IS.0.OR.1, convert Y to 0 and 1
  if(test %in% c(TEST.DEF.TESTS.OVER.GROUPS,TEST.DEF.Y.IS.0.OR.1)){
    y = as.numeric( as.factor(y) ) - 1
  }
     
  #Check input validity
  input_check_result = check.input.main(X,y,ind_reference_taxa,test, q, nr_perm, disable_DSFDR,user_defined_test_function,verbose,DSFDR_Filter)
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
  effect_size_estimates  = rep(NA,p)
  effect_size_estimates_ratio  = rep(NA,p)
  
  if(return_rarefied_values){
    rarefied_counts =matrix(NA,nrow = n,ncol = p)
  }
  
  # Compute reference values
  if(length(ind_reference_taxa)>1){
    for(i in 1:n){
      reference_values[i] = sum(X[i,ind_reference_taxa])
    }
  }else{
    reference_values = X[,ind_reference_taxa] 
  }
  
  stats_matrix = matrix(NA, ncol = p, nrow = nr_perm+1)
  stats_matrix_ratio_normalization = matrix(NA, ncol = p, nrow = nr_perm+1)
  rarefaction_matrix = matrix(NA,nrow = n,ncol = 1)
  ratio_matrix = matrix(NA,nrow = n,ncol = 1)
  
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
    y_original = c(rep("First Group",n/2),rep("Second Group",n/2))
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
    if((!Test_All & i %in% ind_reference_taxa) | !DSFDR_Filter[i]){
      next
    }
    
    nom = X[,i]
    dnom = reference_values
    
    if(i %in% ind_reference_taxa){
      dnom  = dnom - nom
    }
    
    #choose rarefaction depth:
    total_reads_per_subject = nom+dnom
    
    min_value = min(total_reads_per_subject)
    min_value_array[i] = min_value

    #perform subsample and test
    
    temp_subsampled =  rhyper(n, nom, dnom,min_value)
    if(return_rarefied_values){
      rarefied_counts[,i] = temp_subsampled
    }
    rarefaction_matrix[,1] = temp_subsampled 
    ratio_matrix[,1] = nom/(dnom+nom)
    
    stats_matrix[,i] = Compute.resample.test(rarefaction_matrix,Y_matrix,statistic = test, user_defined_test_function = user_defined_test_function)
    
    if(compute_ratio_normalization){
      stats_matrix_ratio_normalization[,i] = Compute.resample.test(ratio_matrix,Y_matrix,statistic = test, user_defined_test_function = user_defined_test_function)  
    }else{
      stats_matrix_ratio_normalization[,i] = 0
    }
    
    #compute effect size estimates
    if(test %in% TEST.DEF.TESTS.ON.UNIVARIATE_CONTINOUS){
      effect_size_estimates[i] = suppressWarnings(as.character(cor(rank(rarefaction_matrix[,1]),rank(Y_matrix[,1])))) 
      if(compute_ratio_normalization){
        effect_size_estimates_ratio[i] = suppressWarnings(as.character(cor(rank(ratio_matrix[,1]),rank(Y_matrix[,1]))))
      }
    }else if(test != DACOMP.TEST.NAME.USER_DEFINED){
      effect_size_estimates[i] = description_for_KS(data.frame(X_rank = rank(rarefaction_matrix[,1]),
                                                               Y_perm = y_original))
      if(compute_ratio_normalization){
        effect_size_estimates_ratio[i] = description_for_KS(data.frame(X_rank = rank(ratio_matrix[,1]),
                                                                       Y_perm = y_original))
      }
    }
  }#end of iteration over taxa 
  stats = stats_matrix
  p.values = rep(NA,ncol(stats))
  p.values.ratio.normalization = rep(NA,ncol(stats))
  for(i in 1:ncol(stats)){
    p.values[i] = mean(stats[,i]>=stats[1,i])
    if(compute_ratio_normalization){
      p.values.ratio.normalization[i] = mean(stats_matrix_ratio_normalization[,i]>=stats_matrix_ratio_normalization[1,i])
    }
  }
  
  #compute DS-FDR:
  Taxa_for_DSFDR = DSFDR_Filter
  if(!Test_All)
    Taxa_for_DSFDR[ind_reference_taxa] = F
  if(!disable_DSFDR){
    dsfdr_obj = dsfdr_find_thresholds(stats[,Taxa_for_DSFDR,drop=F],q,verbose)  
    dsfdr_threshold = dsfdr_obj$selected_c
    Adj.P.value.DSFDR = dsfdr_obj$Adj.P.Value
    if(compute_ratio_normalization){
      dsfdr_obj_ratio_normalization = dsfdr_find_thresholds(stats_matrix_ratio_normalization[,Taxa_for_DSFDR,drop=F],q,F)  
      dsfdr_threshold_ratio_normalization = dsfdr_obj_ratio_normalization$selected_c
      Adj.P.value.DSFDR_ratio_normalization = dsfdr_obj_ratio_normalization$Adj.P.Value
    }
  }
  
  p.values.test = p.values
  if(!Test_All){
    p.values.test[ind_reference_taxa] = NA
  }
  p.values.test.adjusted = p.adjust(p.values.test,method = 'BH')
  if(!disable_DSFDR){
    p.values.test.adjusted[Taxa_for_DSFDR] = Adj.P.value.DSFDR
  }
  if(compute_ratio_normalization){
    p.values.test.ratio.normalization = p.values.ratio.normalization
    if(!Test_All){
      p.values.test.ratio.normalization[ind_reference_taxa] = NA
    }
    p.values.test.adjusted.ratio.normalization = p.adjust(p.values.test.ratio.normalization,method = 'BH')
    if(!disable_DSFDR){
      p.values.test.adjusted.ratio.normalization[Taxa_for_DSFDR] = Adj.P.value.DSFDR_ratio_normalization
    }
  }
  
  
  
  #return results:
  
  ret = list()
  ret$lambda = min_value_array
  ret$stats_matrix = stats_matrix
  ret$p.values.test = p.values.test
  ret$p.values.test.adjusted = p.values.test.adjusted
  ret$effect_size_estimates = effect_size_estimates
  if(!disable_DSFDR){
    ret$dsfdr_rejected = which(p.values.test<=dsfdr_threshold)
    ret$dsfdr_threshold = dsfdr_threshold  
  }
  
  if(compute_ratio_normalization){
    ret$p.values.test.ratio.normalization = p.values.test.ratio.normalization
    ret$p.values.test.adjusted.ratio.normalization = p.values.test.adjusted.ratio.normalization
    ret$effect_size_estimates_ratio = effect_size_estimates_ratio
    if(!disable_DSFDR){
      ret$dsfdr_rejected_ratio_normalization = which(p.values.test.ratio.normalization<=dsfdr_threshold_ratio_normalization)
      ret$dsfdr_threshold_ratio_normalization = dsfdr_threshold_ratio_normalization  
    }
  }
  if(return_rarefied_values){
    ret$rarefied_counts = rarefied_counts
  }
  class(ret) = CLASS.LABEL.DACOMP_RESULT_OBJECT
  return(ret)
}

#internal function for validating inputs on dacomp.test
check.input.main = function(X, y, ind_reference_taxa, test, q, nr_perm, disable_DSFDR,user_defined_test_function, verbose,DSFDR_Filter){
  
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
  
  if(ncol(X) != length(DSFDR_Filter)){
    stop('length of DSFDR_Filter not identical to ncol(X)')
  }
  return(TRUE)
}

