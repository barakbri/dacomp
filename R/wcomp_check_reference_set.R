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
wcomp.check_reference_set_is_valid = function(X_ref,Y,nr.perm=10^4,verbose = F){
  
  input_check_result = check.input.wcomp.check_reference_set_is_valid(X_ref,Y,nr.perm,verbose)
  if(!input_check_result)
    stop('Input check failed on wcomp.check_reference_set_is_valid')
  
  if(dim(X_ref)[2] == 1){
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
  
  rarefy_depth = min(as.numeric(apply(X_ref,1,sum)))
  X_ref = vegan::rrarefy(X_ref,rarefy_depth)
  if(verbose)
    cat(paste0('Computing distance matrices:\n\r'))
  
  dist_obj_L2 =dist(X_ref,method = "euclidean", diag = TRUE, upper = TRUE)
  dist_obj_L1 =dist(X_ref,method = "manhattan", diag = TRUE, upper = TRUE)
  dist_obj_BC = vegan::vegdist(X_ref,diag = TRUE, upper = TRUE)
  
  Dx_ref_L2 = as.matrix(dist_obj_L2)
  Dx_ref_L1 = as.matrix(dist_obj_L1)
  Dx_ref_BC = as.matrix(dist_obj_BC)
  
  if(verbose)
    cat(paste0('HHG\n\r'))
  
  hhg.res.L2 = HHG::hhg.test.k.sample(Dx_ref_L2,Y,nr.threads = 1,nr.perm = nr.perm,perm.stats.wanted = T)
  hhg.res.L1 = HHG::hhg.test.k.sample(Dx_ref_L1,Y,nr.threads = 1,nr.perm = nr.perm,perm.stats.wanted = T)
  hhg.res.BC = HHG::hhg.test.k.sample(Dx_ref_BC,Y,nr.threads = 1,nr.perm = nr.perm,perm.stats.wanted = T)
  
  
  if(verbose)
    cat(paste0('DISCO \n\r'))
  
  energy_L2 = energy::disco(Dx_ref_L2,factors = Y,distance = T,R = nr.perm)
  energy_L1 = energy::disco(Dx_ref_L1,factors = Y,distance = T,R = nr.perm)
  energy_BC = energy::disco(Dx_ref_BC,factors = Y,distance = T,R = nr.perm)
  
  if(verbose)
    cat(paste0('Permanova \n\r'))
  
  permanova_L2 = vegan::adonis(dist_obj_L2~Y,permutations = nr.perm)
  permanova_L1 = vegan::adonis(dist_obj_L1~Y,permutations = nr.perm)
  permanova_BC = vegan::adonis(dist_obj_BC~Y,permutations = nr.perm)
  
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
  return(TRUE)
}
