


#' Title
#'
#' @param ref_select 
#' @param label 
#' @param quantiles_to_plot 
#' @param breaks_param 
#'
#' @return
#' @export
#'
#' @examples
wcomp.plot_reference_scores = function(ref_select,label="Histogram of medianSD statstic for reference selection",quantiles_to_plot = c(0.5,0.7,0.9), breaks_param = 30){
  
  input_check_result = check.input.wcomp.plot_reference_scores(ref_select,label,quantiles_to_plot,breaks_param)
  if(!input_check_result)
    stop('Input check failed on wcomp.plot_reference_scores')
  
  Target_MinAbundance_values = ref_select$target_abundance
  hist(ref_select$scores,breaks = breaks_param,main = label,xlab = "medianSD statistic")
  sorted_scores = sort(ref_select$scores)
  threshold_ind = 1
  for(threshold_ind in 1:length(Target_MinAbundance_values)){
    abline(v = sorted_scores[ref_select$arg_obj[threshold_ind]],col = threshold_ind+1,lwd = 3)
  }  
  qunatiles_in_scores = quantile(ref_select$scores,probs = quantiles_to_plot)
  for(i in 1:length(quantiles_to_plot)){
    abline(v = qunatiles_in_scores,col = 1,lty = 2,lwd = 3)
  }
}

check.input.wcomp.plot_reference_scores = function(ref_select,label,quantiles_to_plot,breaks_param){
  return(TRUE)  
}



