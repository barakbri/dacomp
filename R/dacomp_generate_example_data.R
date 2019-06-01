

#' Generate a simulated two sample dataset, based on data from the Phyloseq package
#'
#' This function generates a two-sample dataset, based on the \code{kostic} dataset (Kostic et. al. 2012) from the \code{phyloseq} package (McMurdie et. al. 2012). Simulated data is generated in a procedure similar to the one presented in Brill et. al. 2019, Subsection 4.1. See additionals details below.
#' 
#' @details 
#' Data is generated as follows.
#' In the first step, we generate a list of vectors of relative frequencies to sample from: only healthy subjects from the kostic colorectal dataset are selected. Samples with less than 500 reads are dropped. Only OTUs that appear in 2 or more subjects are retained.
#' In the second step, samples for group X are generated. For each sample, a vector of frequencies is chosen at random from the list generated in the first step. The observed sampled are multinomial random variables with a probability vector matching the selected frequencies, and a total number of reads realized from a Poisson distribution with a mean number of reads equal to the median number of reads across the samples listed in the first step.
#' In the third step, samples for group Y are generated. For each sample, a vector of frequencies is chosen at random, similar to group X. The frequencies of differentially abundant taxa is increased, with the increase realized from a poisson random variable, such that the total increase in microbial load across all differentially abundant taxa is equivlant to the signal strength specified by the user. Observed counts are sampled based on the updated frequncies.
#' @param n_X Number of samples from the first group
#' @param n_Y Number of samples from the second group
#' @param m1 Number of differentially abundant taxa
#' @param signal_strength_as_change_in_microbial_load A number in the range 0-0.75, indicating the fraction of the microbial load of group Y that is added due to the simulated condition. The complement of this fraction, is the fraction of the microbial load of group Y that is distribued across taxa as in group X.
#'
#' @return a list with the followig entries
#' \itemize{
#' \item{counts}{A counts matrix with \code{(n_X + n_Y)} rows, and 1384 columns, rows represent samples,columns represent taxa.}
#' \item{group_labels}{A vector of group labelings, with values 0 and 1}
#' \item{select_diff_abundant}{A vector containing the indices of taxa that are differentially abundant.}
#' }
#' 
#' @export
#'
#' @references 
#' Brill, Barak, Amnon Amir, and Ruth Heller. 2019. “Testing for Differential Abundance in Compositional Counts Data, with Application to Microbiome Studies.” arXiv Preprint arXiv:1904.08937.
#' Kostic, Aleksandar D, Dirk Gevers, Chandra Sekhar Pedamallu, Monia Michaud, Fujiko Duke, Ashlee M Earl, Akinyemi I Ojesina, et al. 2012. “Genomic Analysis Identifies Association of Fusobacterium with Colorectal Carcinoma.” Genome Research 22 (2). Cold Spring Harbor Lab: 292–98.
#' McMurdie, Paul J, and Susan Holmes. 2013. “Phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data.” PloS One 8 (4). Public Library of Science: e61217.
#' 
#' @examples
#' \dontrun{
#' library(dacomp)
#'
#' set.seed(1)
#' data = dacomp.generate_example_dataset.two_sample(m1 = 100,
#'        n_X = 50,
#'        n_Y = 50,
#'        signal_strength_as_change_in_microbial_load = 0.1)
#'
#'
#'
#' } 
dacomp.generate_example_dataset.two_sample = function(n_X = 30,n_Y = 30,m1 = 30, signal_strength_as_change_in_microbial_load = 0.1){
  
  #check inputs
  input_check_result = check.input.generate_example_dataset(n_X,n_Y,m1, signal_strength_as_change_in_microbial_load)
  if(!input_check_result)
    stop('Input check failed on dacomp.generate_example_dataset')
  #load phyloseq and prepare data
  
  kostic = get_kostic_data()
  kostic_counts_data_healthy = kostic$kostic_counts_data_healthy
  kostic_N_reads = kostic$kostic_N_reads  
  
  p= ncol(kostic_counts_data_healthy)
  select_diff_abundant = sample(1:p,size = m1,replace = F)
  n_samples_to_sample_from = nrow(kostic_counts_data_healthy)
  N_reads = kostic_N_reads
  X = matrix(NA,nrow = n_X,ncol = p)
  Y = matrix(NA,nrow = n_Y,ncol = p)
  
  #sample 'healthy'
  for(i in 1:n_X){
    sample_abundances_X = kostic_counts_data_healthy[sample(1:n_samples_to_sample_from,1),]
    X[i,] = rmultinom(1,size = rpois(1,N_reads),prob = sample_abundances_X)  
  }
  
  #sample 'sick':
  for(i in 1:n_Y){
    sample_abundances_Y = kostic_counts_data_healthy[sample(1:n_samples_to_sample_from,1),,drop = T]
    sample_abundances_Y[1,select_diff_abundant] = sample_abundances_Y[1,select_diff_abundant] + rpois(m1,sum(sample_abundances_Y) * signal_strength_as_change_in_microbial_load/m1)
    Y[i,] = rmultinom(1,size = rpois(1,N_reads),prob = sample_abundances_Y)  
  }
  #matrix of counts:
  counts = rbind(X,Y)
  #group labels:
  group_labels = c(rep(0,n_X),rep(1,n_Y))
  
  ret = list(counts = counts,
             group_labels = group_labels,
             select_diff_abundant = select_diff_abundant)
  return(ret)
}

#internal function for checking inputs:
check.input.generate_example_dataset = function(n_X,n_Y,m1, signal_strength_as_change_in_microbial_load){
  # n_X,n_Y - sample sizes
  if(n_X!= as.integer(n_X) | n_Y!= as.integer(n_Y))
    stop('n_X and n_Y must be integers')
  if(n_X<5 | n_Y<5 )
    stop('Minimum number of 5 samples in each study group is required')
  
  # m1 - number of diff abundant taxa - should be between 0 and 100
  if(m1!=as.integer(m1))
    stop('m1 must be integer')
  if(m1<1 | m1>100)
    stop('m1 must be between 1 and 100')
  
  # signal_strength_as_change_in_microbial_load - should be between 0 and 0.5
  if(!is.numeric(signal_strength_as_change_in_microbial_load))
    stop('signal_strength_as_change_in_microbial_load must be between 0 and 0.5')
  if(signal_strength_as_change_in_microbial_load<0 | signal_strength_as_change_in_microbial_load>0.75)
    stop('signal_strength_as_change_in_microbial_load must be between 0 and 0.75')
  
  return(TRUE)
}



#' Generate a simulated dataset with a continuous phenotype, based on data from the Phyloseq package
#' This function generates a dataset with a continuous Phenotype, based on the \code{kostic} dataset (Kostic et. al. 2012) from the \code{phyloseq} package (McMurdie et. al. 2012). Simulated data is generated in a procedure similar to the one presented in Brill et. al. 2019, Subsection 4.1. See additionals details below.
#' @details 
#' Data is generated as follows.
#' In the first step, we generate a list of frequency vectors to sample from: only healthy subjects from the kostic colorectal dataset are selected. Samples with less than 500 reads are dropped. Only OTUs that appear in 2 or more subjects are retained.
#' In the seccond step, a random phenotype is sampled for each sample, from a uniform(0,1) distribution.
#' In the third step, samples are generated. For each sample, a vector of frequencies is chosen at random, The differentially abundant taxa are increased, with the additions realized from a poisson random variable. The signal inserted is such that a phenotype with a value of 1 is equivlant to an increase in the microbial load, \code{signal_strength_as_change_in_microbial_load} in fraction of the original microbial load.
#' @param n Number of samples
#' @param m1 Number of differentially abundant taxa
#' @param signal_strength_as_change_in_microbial_load A number in the range 0-0.75, indicating the fraction of the microbial load that is added to the measured ecosystem, if the phenotype for the sample is equal to 1. For phenotypes with lower values, the change in the microbial load is proportional to the value of the measured phenotype.
#'
#' @return a list
#' \itemize{
#' \item{counts}{A counts matrix with \code{n} rows, and 1384 columns, rows represent samples,columns represent taxa.}
#' \item{covariate}{The measured phenotype}
#' \item{select_diff_abundant}{A vector containing the indices of taxa that are differentially abundant.}
#' }
#' 
#' @export
#' 
#' @references 
#' Brill, Barak, Amnon Amir, and Ruth Heller. 2019. “Testing for Differential Abundance in Compositional Counts Data, with Application to Microbiome Studies.” arXiv Preprint arXiv:1904.08937.
#' Kostic, Aleksandar D, Dirk Gevers, Chandra Sekhar Pedamallu, Monia Michaud, Fujiko Duke, Ashlee M Earl, Akinyemi I Ojesina, et al. 2012. “Genomic Analysis Identifies Association of Fusobacterium with Colorectal Carcinoma.” Genome Research 22 (2). Cold Spring Harbor Lab: 292–98.
#' McMurdie, Paul J, and Susan Holmes. 2013. “Phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data.” PloS One 8 (4). Public Library of Science: e61217.
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
  input_check_result = check.input.generate_example_dataset_continuous(n = n,m1 = m1, signal_strength_as_change_in_microbial_load)
  if(!input_check_result)
    stop('Input check failed on dacomp.generate_example_dataset')
  
  #load phyloseq and prepare data
  kostic = get_kostic_data()
  kostic_counts_data_healthy = kostic$kostic_counts_data_healthy
  kostic_N_reads = kostic$kostic_N_reads
  
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
check.input.generate_example_dataset_continuous = function(n,m1, signal_strength_as_change_in_microbial_load){
  
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


#' Generate an example dataset with a multivariate phenotype
#' Generate a simulated dataset, similar to \code{\link{dacomp.generate_example_dataset_continuous}} with the following difference: the generated dataset contains two phenotypes, instead of one.
#' The change observed in a sample, is monotone increasing with the values of each measured covariate.
#' 
#' @param n Number of samples.
#' @param m1 Number of differentially abundant taxa
#' @param signal_strength_as_change_in_microbial_load 
#'
#' @return a list
#' \itemize{
#' \item{counts}{A counts matrix with \code{n} rows, and 1384 columns, rows represent samples,columns represent taxa.}
#' \item{covariate}{The measured phenotype, a matrix of size \code{n X 2}, rows in this matrix correspond to the rows of \code{counts}}
#' \item{select_diff_abundant}{A vector containing the indices of taxa that are differentially abundant.}
#' }
#' 
#' @references 
#' Brill, Barak, Amnon Amir, and Ruth Heller. 2019. “Testing for Differential Abundance in Compositional Counts Data, with Application to Microbiome Studies.” arXiv Preprint arXiv:1904.08937.
#' Kostic, Aleksandar D, Dirk Gevers, Chandra Sekhar Pedamallu, Monia Michaud, Fujiko Duke, Ashlee M Earl, Akinyemi I Ojesina, et al. 2012. “Genomic Analysis Identifies Association of Fusobacterium with Colorectal Carcinoma.” Genome Research 22 (2). Cold Spring Harbor Lab: 292–98.
#' McMurdie, Paul J, and Susan Holmes. 2013. “Phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data.” PloS One 8 (4). Public Library of Science: e61217.
#' 
#' @export
#'
#' @examples
#' data = dacomp.generate_example_dataset_multivariate_example(30)
#' 
dacomp.generate_example_dataset_multivariate_example = function(n, m1 = 30, signal_strength_as_change_in_microbial_load = 0.1){
  
  #check inputs
  input_check_result = check.input.generate_example_dataset_continuous(n = n,m1 = m1, signal_strength_as_change_in_microbial_load)
  if(!input_check_result)
    stop('Input check failed on dacomp.generate_example_dataset')
  
  #load phyloseq and prepare data
  kostic = get_kostic_data()
  kostic_counts_data_healthy = kostic$kostic_counts_data_healthy
  kostic_N_reads = kostic$kostic_N_reads
  
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

# internal function for accessing the kostic dataset from the phyloseq package
get_kostic_data = function(){
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
  kostic_counts_data_healthy = kostic_counts_data_healthy[ , apply(kostic_counts_data_healthy>0,2,sum) >= 2 ]
  kostic_N_reads = median(apply(kostic_counts_data_healthy,1,sum))
  ret = list()
  ret$kostic_counts_data_healthy = kostic_counts_data_healthy
  ret$kostic_N_reads = kostic_N_reads
  return(ret)
}



#' Generate a simulated dataset for a paired study design.
#' 
#' The function generates a simulated data set for a paired study design. The first \code{n} rows correspond to ecosystems sampled, measured under condition 1. The next \code{n} rows correspond to the same ecosystems, measured under condition 2.
#' 
#' @details 
#' Data generation is similar to \code{\link{dacomp.generate_example_dataset.two_sample}} with the following difference: In the first step, samples for condition 1 are generated as in the two-sample function.
#' In the second step, the frequencies of taxa (used for sampling in the previous step) are increased to simulate the changing condition. Samples under condition 2 are multinomial samples, generated using the modified frequency vectors.
#'
#' @param n Number of ecosystems to sample. Actual number of samples will be twice this number.
#' @param m1 Number of differentially abundant taxa
#' @param signal_strength_as_change_in_microbial_load A number in the range 0-0.75, indicating the fractional increase in the microbial load of a sample, when changing from condition 1 to condition 2.
#'
#' @return a list
#' \itemize{
#' \item{counts}{A counts matrix with \code{(2 * n)} rows, and 1384 columns, rows represent samples,columns represent taxa.}
#' \item{select_diff_abundant}{A vector containing the indices of taxa that are differentially abundant.}
#' }
#' @export
#'
#' @references 
#' Brill, Barak, Amnon Amir, and Ruth Heller. 2019. “Testing for Differential Abundance in Compositional Counts Data, with Application to Microbiome Studies.” arXiv Preprint arXiv:1904.08937.
#' Kostic, Aleksandar D, Dirk Gevers, Chandra Sekhar Pedamallu, Monia Michaud, Fujiko Duke, Ashlee M Earl, Akinyemi I Ojesina, et al. 2012. “Genomic Analysis Identifies Association of Fusobacterium with Colorectal Carcinoma.” Genome Research 22 (2). Cold Spring Harbor Lab: 292–98.
#' McMurdie, Paul J, and Susan Holmes. 2013. “Phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data.” PloS One 8 (4). Public Library of Science: e61217.
#' 
#' @examples
#' data = dacomp.generate_example_dataset_paired(30)
#' 
dacomp.generate_example_dataset_paired = function(n, m1 = 30, signal_strength_as_change_in_microbial_load = 0.1){
  
  #check inputs
  input_check_result = check.input.generate_example_dataset_continuous(n = n,m1 = m1, signal_strength_as_change_in_microbial_load)
  if(!input_check_result)
    stop('Input check failed on dacomp.generate_example_dataset')
  #load phyloseq and prepare data
  
  kostic = get_kostic_data()
  kostic_counts_data_healthy = kostic$kostic_counts_data_healthy
  kostic_N_reads = kostic$kostic_N_reads
  
  p= ncol(kostic_counts_data_healthy)
  select_diff_abundant = sample(1:p,size = m1,replace = F)
  n_samples_to_sample_from = nrow(kostic_counts_data_healthy)
  N_reads = kostic_N_reads
  X = matrix(NA,nrow = n,ncol = p)
  Y = matrix(NA,nrow = n,ncol = p)
  
  
  for(i in 1:n){
    sample_abundances_X = kostic_counts_data_healthy[sample(1:n_samples_to_sample_from,1),,drop = T]
    sample_abundances_Y = sample_abundances_X
    sample_abundances_Y[1,select_diff_abundant] = sample_abundances_Y[1,select_diff_abundant] + rpois(m1,sum(sample_abundances_X) * signal_strength_as_change_in_microbial_load/m1)
    X[i,] = rmultinom(1,size = rpois(1,N_reads),prob = sample_abundances_X)  
    Y[i,] = rmultinom(1,size = rpois(1,N_reads),prob = sample_abundances_Y)  
  }
  #matrix of counts:
  counts = rbind(X,Y)
  
  ret = list(counts = counts,
             select_diff_abundant = select_diff_abundant)
  return(ret)
}
