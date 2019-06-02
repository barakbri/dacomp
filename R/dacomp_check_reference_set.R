#Class label for object with results of a reference validation procedure
CLASS.LABEL.REFERENCE_VALIDATION_RESULT_OBJECT = "dacomp.reference.validation.result.object"

library(vegan)
library(HHG)
library(energy)


#' Test reference set of taxa for validity, by multivariate K-Sample testing
#' 
#' The function carries out the Reference Validation Procedure (RVP), described in Brill et. al. (2019), appendix A. The test verifies that no differentially abundant taxa have entered the reference set. The input to the function includes a counts matrix of taxa selected as a reference set (rows represent samples), and a vector \code{Y} of groups labels. The function first rarefies all samples to constant depth, then tests for multivariate equality of distributions across the differet study groups.
#' 
#' @details 
#' Rarefying samples to equal depth is mandatory for testing, see Brill et. al. (2019) Appendix A for discussion of the test and assumptions. The multivariate K-sample tests performed include the \code{hhg.test.k.sample} (Heller et. al. 2019) from package \code{HHG}, The \code{disco} test (Székely et. al. 2007 ) from package \code{energy} and the Permanova test (Anderson et, al. 2001), available through the function \code{adonis} in the \code{vegan} package.
#' All three tests require a distance (or dissimilarity) metric between samples. All three tests are evaluated with the L1,L2 distance metrics and Bray-Curtis dissimilarity metric.
#' This function supports categorical phenotypes. For continuous and multivariate phenotypes, a different procedure is required, making use of multivariate tests of independence such as \code{hhg.test} (from packge \code{HHG}), and \code{dcov.test} (from package \code{energy}).
#' 
#' @param X_ref Matrix of counts, containing the reference taxa alone. Rows represent samples, columns represent taxa.
#' @param Y Categorical variable of the measure phenotype, the group labeling.
#' @param nr.perm Number of permutations for computing P-values. 
#' @param verbose Should messages be printed to console, indicating which test is being computed.
#' 
#' @return An object of type "dacomp.reference.validation.result.object", is a list with the follows fields:
#' \itemize{
#' \item{p.value.HHG.L2}{ - P-value for the HHG K sample test, with the L2 distance metric.}
#' \item{p.value.HHG.L1}{ - P-value for the HHG K sample test, with the L1 distance metric.}
#' \item{p.value.HHG.BC}{ - P-value for the HHG K sample test, with the Bray-Curtis dissimilarity metric.}
#' \item{p.value.energy.L2}{ - P-value for the disco test, with the L2 distance metric.}
#' \item{p.value.energy.L1}{ - P-value for the disco test, with the L1 distance metric.}
#' \item{p.value.energy.BC}{ - P-value for the disco test, with the Bray-Curtis dissimilarity metric.}
#' \item{p.value.permanova_L2}{ - P-value for the permanova test, with the L2 distance metric.}
#' \item{p.value.permanova_L1}{ - P-value for the permanova test, with the L1 distance metric.}
#' \item{p.value.permanova_BC}{ - P-value for the permanova test, with the Bray-Curtis dissimilarity metric.}
#' }
#' 
#' @references
#' Brill, Barak, Amnon Amir, and Ruth Heller. 2019. Testing for Differential Abundance in Compositional Counts Data, with Application to Microbiome Studies. arXiv Preprint arXiv:1904.08937.
#' 
#' Szekely, Gabor J, Maria L Rizzo, Nail K Bakirov, and others. 2007. Measuring and Testing Dependence by Correlation of Distances. The Annals of Statistics 35 (6). Institute of Mathematical Statistics: 276994.
#' 
#' Ruth Heller, Yair Heller, and Malka Gorfine. A consistent multivariate test of association based on ranks of distances. Biometrika, 100(2):503 510, 2012.
#' 
#' Anderson, M.J. 2001. A new method for non-parametric multivariate analysis of variance. Austral Ecology, 26: 3246.
#' 
#' @examples
#' \dontrun{
#'
#' library(dacomp)
#' 
#' set.seed(1)
#' 
#' data = dacomp.generate_example_dataset.two_sample(m1 = 100, 
#'        n_X = 50, n_Y = 50,
#'        signal_strength_as_change_in_microbial_load = 0.1)
#' 
#' #select references: (may take a minute)
#' 
#' result.selected.references = dacomp.select_references(X = data$counts,
#'                                  median_SD_threshold = 0.6, #APPLICATION SPECIFIC
#'                                  verbose = T)
#' 
#' length(result.selected.references$selected_references)
#' 
#' #plot the reference selection scores (can be used to better set the median SD threshold...)
#' dacomp.plot_reference_scores(result.selected.references)
#' 
#' result.ref.validity = dacomp.check_reference_set_is_valid.k_groups(
#'                         X_ref = data$counts[,result.selected.references$selected_references],
#'                         Y = data$group_labels,
#'                         nr.perm = 10000,
#'                         verbose = T)
#' result.ref.validity
#'
#' } 
dacomp.check_reference_set_is_valid.k_groups = function(X_ref,Y,nr.perm=10^4,verbose = F){
  # check inputs
  input_check_result = check.input.check_reference_set_is_valid(X_ref,Y,nr.perm,verbose)
  if(!input_check_result)
    stop('Input check failed on dacomp.check_reference_set_is_valid')
  
  #the reference validation procedure works on reference sets with at least two taxa.
  if(dim(X_ref)[2] == 1){
    warning(' reference validation procedure cannot be run - a single reference taxon has been selected')
    ret = list()
    
    ret$p.value.HHG.L2 = 1
    ret$p.value.HHG.L1 = 1
    ret$p.value.HHG.BC = 1
    
    ret$p.value.energy.L2 = 1
    ret$p.value.energy.L1 = 1
    ret$p.value.energy.BC = 1
    
    ret$p.value.permanova_L2 = 1
    ret$p.value.permanova_L1 = 1
    ret$p.value.permanova_BC = 1
    return(ret)
  }
  
  #find the rarefaction depth, and rarefy all reference samples to this depth
  rarefy_depth = min(as.numeric(apply(X_ref,1,sum)))
  X_ref = vegan::rrarefy(X_ref,rarefy_depth)
  
  #compute distance matrices
  if(verbose)
    cat(paste0('Computing distance matrices\n\r'))
  
  dist_obj_L2 =dist(X_ref,method = "euclidean", diag = TRUE, upper = TRUE)
  dist_obj_L1 =dist(X_ref,method = "manhattan", diag = TRUE, upper = TRUE)
  dist_obj_BC = vegan::vegdist(X_ref,diag = TRUE, upper = TRUE)
  
  Dx_ref_L2 = as.matrix(dist_obj_L2)
  Dx_ref_L1 = as.matrix(dist_obj_L1)
  Dx_ref_BC = as.matrix(dist_obj_BC)
  # run the different tests: HHG, Permanova and DISCO
  # with L1,L2 and Bray-Curtis distance (BC actually not a distance...)
  
  if(verbose)
    cat(paste0('Testing reference set validity by HHG\n\r'))
  
  hhg.res.L2 = HHG::hhg.test.k.sample(Dx_ref_L2,Y,nr.threads = 1,nr.perm = nr.perm,perm.stats.wanted = T)
  hhg.res.L1 = HHG::hhg.test.k.sample(Dx_ref_L1,Y,nr.threads = 1,nr.perm = nr.perm,perm.stats.wanted = T)
  hhg.res.BC = HHG::hhg.test.k.sample(Dx_ref_BC,Y,nr.threads = 1,nr.perm = nr.perm,perm.stats.wanted = T)
  
  if(verbose)
    cat(paste0('Testing reference set validity by DISCO \n\r'))
  
  energy_L2 = energy::disco(Dx_ref_L2,factors = Y,distance = T,R = nr.perm)
  energy_L1 = energy::disco(Dx_ref_L1,factors = Y,distance = T,R = nr.perm)
  energy_BC = energy::disco(Dx_ref_BC,factors = Y,distance = T,R = nr.perm)
  
  if(verbose)
    cat(paste0('Testing reference set validity by Permanova \n\r'))
  
  permanova_L2 = vegan::adonis(dist_obj_L2~Y,permutations = nr.perm)
  permanova_L1 = vegan::adonis(dist_obj_L1~Y,permutations = nr.perm)
  permanova_BC = vegan::adonis(dist_obj_BC~Y,permutations = nr.perm)
  
  #return results
  ret = list()
  
  ret$p.value.HHG.L2 = hhg.res.L2$perm.pval.hhg.sc
  ret$p.value.HHG.L1 = hhg.res.L1$perm.pval.hhg.sc
  ret$p.value.HHG.BC = hhg.res.BC$perm.pval.hhg.sc
  
  ret$p.value.energy.L2 = energy_L2$p.value
  ret$p.value.energy.L1 = energy_L1$p.value
  ret$p.value.energy.BC = energy_BC$p.value
  
  ret$p.value.permanova_L2 = permanova_L2$aov.tab[1,6]
  ret$p.value.permanova_L1 = permanova_L1$aov.tab[1,6]
  ret$p.value.permanova_BC = permanova_BC$aov.tab[1,6]
  
  class(ret) = CLASS.LABEL.REFERENCE_VALIDATION_RESULT_OBJECT
  return(ret)
}

check.input.check_reference_set_is_valid = function(X_ref,Y,nr.perm,verbose){
  
  # X_ref, - check is numeric matrix
  MSG_X_REF = 'X_ref must be a valid counts matrix'
  
  if(!is.matrix(X_ref))
    stop(MSG_X_REF)
  if(any(X_ref!=as.integer(X_ref)))
    stop(MSG_X_REF)
  if(any(X_ref<0))
    stop(MSG_X_REF)
  
  # Y, - check is a group of values, same number of measurements as X_ref
  
  if(length(Y)!= nrow(X_ref))
    stop('length of Y must be same as number of rows in X_ref')
  if(any(is.na(Y))|any(is.nan(Y)))
    stop('Y has NA or NaNs - invalid observations')
  
  # nr.perm, number of permuatations should be valid
  MSG_NR_PERM = 'nr.perm must be at integer, at least 1000'
  if(nr.perm != as.integer(nr.perm))
    stop(MSG_NR_PERM)
  if(nr.perm<1000)
    stop(MSG_NR_PERM)
  
  # verbose - should be logical
  if(!is.logical(verbose))
    stop('Verbose must be logical')
  
  return(TRUE)
}
