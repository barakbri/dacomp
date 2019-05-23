

#' Title
#'
#' @param n
#' @param m1 
#' @param signal_strength_as_change_in_microbial_load 
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' library(dacomp)
#'
#' set.seed(1)
#' data = dacomp.generate_example_dataset_continuous(n = 100,
#' m1 = 30,signal_strength_as_change_in_microbial_load = 0.1)
#'
#'
#'
#' } 
dacomp.generate_example_dataset_continuous = function(n, m1 = 30, signal_strength_as_change_in_microbial_load = 0.1){
  
  #check inputs
  input_check_result = check.input.dacomp.generate_example_dataset_continuous(n = n,m1 = m1, signal_strength_as_change_in_microbial_load)
  if(!input_check_result)
    stop('Input check failed on dacomp.generate_example_dataset')
  #load phyloseq and prepare data
  
  library(phyloseq)
  filepath = system.file("extdata", "study_1457_split_library_seqs_and_mapping.zip", package="phyloseq")
  sink(tempfile())
  kostic = (suppressWarnings(microbio_me_qiime(filepath)))
  sink()
  kostic = subset_samples(kostic, DIAGNOSIS == "Healthy") 
  kostic = prune_samples(sample_sums(kostic) > 500, kostic) #prune samples with a low number of counts
  kostic_counts = as.matrix(t(otu_table(kostic)))
  kostic_sample_data = sample_data(kostic)
  kostic_counts_data_healthy = kostic_counts
  #keep taxa that appear in  at least 2 subjects:
  kostic_counts_data_healthy= kostic_counts_data_healthy[ , apply(kostic_counts_data_healthy>0,2,sum) >= 2 ]
  kostic_N_reads = median(apply(kostic_counts_data_healthy,1,sum))
  
  
  p= ncol(kostic_counts_data_healthy)
  select_diff_abundant = sample(1:p,size = m1,replace = F)
  n_samples_to_sample_from = nrow(kostic_counts_data_healthy)
  N_reads = kostic_N_reads
  X = matrix(NA,nrow = n,ncol = p)
  u = runif(n)
  
  
  
  #sample 'sick':
  for(i in 1:n){
    sample_abundances_X = kostic_counts_data_healthy[sample(1:n_samples_to_sample_from,1),,drop = T]
    sample_abundances_X[1,select_diff_abundant] = sample_abundances_X[1,select_diff_abundant] + u[i] * rpois(m1,sum(sample_abundances_X) * signal_strength_as_change_in_microbial_load/m1)
    X[i,] = rmultinom(1,size = rpois(1,N_reads),prob = sample_abundances_X)  
  }
  #matrix of counts:
  counts = X
  
  ret = list(counts = counts,
             covariate = u,
             select_diff_abundant = select_diff_abundant)
  return(ret)
}

#internal function for checking inputs:
check.input.dacomp.generate_example_dataset_continuous = function(n,m1, signal_strength_as_change_in_microbial_load){
   
  # n - sample sizes
  if(n!= as.integer(n))
    stop('n must be integer')
  if(n<10)
    stop('Minimum number of 10 samples is required')
  
  # m1 - number of diff abundant taxa - should be between 0 and 100
  if(m1!=as.integer(m1))
    stop('m1 must be integer')
  if(m1<1 | m1>100)
    stop('m1 must be between 1 and 100')
  
  # signal_strength_as_change_in_microbial_load - should be between 0 and 0.5
  if(!is.numeric(signal_strength_as_change_in_microbial_load))
    stop('signal_strength_as_change_in_microbial_load must be between 0 and 0.5')
  if(signal_strength_as_change_in_microbial_load<0 | signal_strength_as_change_in_microbial_load>0.5)
    stop('signal_strength_as_change_in_microbial_load must be between 0 and 0.5')
  
  return(TRUE)
}


dacomp.generate_example_dataset_multivariate_example = function(n, m1 = 30, signal_strength_as_change_in_microbial_load = 0.1){
  
  #check inputs
  input_check_result = check.input.dacomp.generate_example_dataset_continuous(n = n,m1 = m1, signal_strength_as_change_in_microbial_load)
  if(!input_check_result)
    stop('Input check failed on dacomp.generate_example_dataset')
  #load phyloseq and prepare data
  
  library(phyloseq)
  filepath = system.file("extdata", "study_1457_split_library_seqs_and_mapping.zip", package="phyloseq")
  sink(tempfile())
  kostic = (suppressWarnings(microbio_me_qiime(filepath)))
  sink()
  kostic = subset_samples(kostic, DIAGNOSIS == "Healthy") 
  kostic = prune_samples(sample_sums(kostic) > 500, kostic) #prune samples with a low number of counts
  kostic_counts = as.matrix(t(otu_table(kostic)))
  kostic_sample_data = sample_data(kostic)
  kostic_counts_data_healthy = kostic_counts
  #keep taxa that appear in  at least 2 subjects:
  kostic_counts_data_healthy= kostic_counts_data_healthy[ , apply(kostic_counts_data_healthy>0,2,sum) >= 2 ]
  kostic_N_reads = median(apply(kostic_counts_data_healthy,1,sum))
  
  
  p= ncol(kostic_counts_data_healthy)
  select_diff_abundant = sample(1:p,size = m1,replace = F)
  n_samples_to_sample_from = nrow(kostic_counts_data_healthy)
  N_reads = kostic_N_reads
  X = matrix(NA,nrow = n,ncol = p)
  u1 = runif(n)
  u2 = runif(n)
  
  
  #sample 'sick':
  for(i in 1:n){
    sample_abundances_X = kostic_counts_data_healthy[sample(1:n_samples_to_sample_from,1),,drop = T]
    sample_abundances_X[1,select_diff_abundant] = sample_abundances_X[1,select_diff_abundant] + (u1[i] +u2[i]) * rpois(m1,sum(sample_abundances_X) * signal_strength_as_change_in_microbial_load/m1)
    X[i,] = rmultinom(1,size = rpois(1,N_reads),prob = sample_abundances_X)  
  }
  #matrix of counts:
  counts = X
  
  ret = list(counts = counts,
             covariate =cbind(u1,u2),
             select_diff_abundant = select_diff_abundant)
  return(ret)
}
