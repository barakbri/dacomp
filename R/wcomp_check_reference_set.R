#Class label for object with results of a reference validation procedure
CLASS.LABEL.REFERENCE_VALIDATION_RESULT_OBJECT = "wcomp.reference.validation.result.object"

library(vegan)
library(HHG)
library(energy)


#' Title
#'
#' @param X_ref 
#' @param Y 
#' @param nr.perm 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
wcomp.check_reference_set_is_valid.k_groups = function(X_ref,Y,nr.perm=10^4,verbose = F){
  # check inputs
  input_check_result = check.input.wcomp.check_reference_set_is_valid(X_ref,Y,nr.perm,verbose)
  if(!input_check_result)
    stop('Input check failed on wcomp.check_reference_set_is_valid')
  
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

check.input.wcomp.check_reference_set_is_valid = function(X_ref,Y,nr.perm,verbose){
  
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
