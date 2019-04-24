
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
wcomp.test = function(X,y,ind_reference_taxa,z = NULL,test = 'Wilcoxon', q=0.05,nr_perm = 1/(q/(ncol(X)-length(ind_reference_taxa))), disable_DSFDR = F,verbose = F){
  
  nr_rarefactions_multiple_X = 1
  if(nr_rarefactions_multiple_X != 1){
    stop("Multiple rarefactions currently not supported")
  }
  p = ncol(X)
  n = nrow(X)
  min_value_array = rep(NA,p)
  rarefied_data_length_array = rep(NA,p)
  pval_res = rep(NA,p)
  pval_res.CVM = rep(NA,p)
  taxa_nr_res = rep(NA,p)
  reference_values = rep(NA,n)
  lambda_selected = rep(NA,p)
  remaining_Y = list()
  
  # Compute reference values
  if(length(ind_reference_taxa)>1){
    for(i in 1:n){
      reference_values[i] = sum(X[i,ind_reference_taxa])
    }
  }else{
    reference_values = X[,ind_reference_taxa] 
  }
  
  stats_matrix = matrix(NA, ncol = p, nrow = nr_perm+1)
  rarefaction_matrix = matrix(NA,nrow = n,ncol = nr_rarefactions_multiple_X)
  
  #We compute the permutation matrix
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

  for(i in 1:p){
    
    if(verbose)
      if(i%% ceiling(p/100) == 1)
        cat(paste0('Testing taxon : ',i,'/',m,' \n\r'))
    
    nom = X[,i]
    dnom = reference_values
    
    #
    if(i %in% ind_reference_taxa){
      nom = X[,i]
      dnom = reference_values - nom 
    }
      
    
    
    #choose rarefaction depth
    total_reads_per_subject = nom+dnom
    
    min_value = min(total_reads_per_subject)
    min_value_array[i] = min_value
    to_keep = which(total_reads_per_subject >= min_value) # no filtration of samples.
    
    z_keep = z
    
    #enforcing pairs to have both items
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
    remaining_Y[[i]] = y_keep
    rarefied_data_length_array[i] = length(nom_keep)
    
    # perform the actual subsample,
    nom_keep_original = nom_keep
    dnom_keep_original = dnom_keep
    
    #perform subsample and test
    for(s in 1:nr_rarefactions_multiple_X){
      temp_subsampled =  rhyper(n, nom_keep_original, dnom_keep_original,min_value)  
      rarefaction_matrix[,s] = temp_subsampled 
    }
    stats_matrix[,i] = Compute.resample.test(rarefaction_matrix,Y_matrix,min_value,statistic = test)
  }
  
  stats = stats_matrix
  p.values = rep(NA,ncol(stats))
  for(i in 1:ncol(stats)){
    p.values[i] = mean(stats[,i]>=stats[1,i])
  }
  
  if(!disable_DSFDR){
    C_test = dfdr_find_thresholds(stats[,-ind_reference_taxa,drop=F],q)  
  }
  
  p.values.test = p.values; p.values.test[ind_reference_taxa] = NA
  p.values.ref = p.values; p.values.ref[-ind_reference_taxa] = NA
  
  ret = list()
  ret$lambda = min_value_array
  ret$stats_matrix = stats_matrix
  ret$p.values.test = p.values.test
  ret$p.values.ref = p.values.ref
  
  if(!disable_DSFDR){
    ret$rejected = which(p.values.test<=C_test)
    ret$C_test = C_test  
  }
  return(ret)
}


