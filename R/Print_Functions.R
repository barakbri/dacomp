#' Print the results contained in an object from type dacomp.result.object
#'
#' @param obj Object returned from dacomp.test(...)
#'
#' @return
#' @export
#'
#' @examples
print.dacomp.result.object = function(obj){
  cat('DACOMP results object \n\r')
  cat(paste0('Nr. taxa tested: ',sum(!is.na(obj$lambda)),'\n\r'))
  cat(paste0('Nr. reference taxa : ',sum(is.na(obj$lambda)),'\n\r'))
  cat(paste0('Mean lambda : ',round(mean(obj$lambda, na.rm = T),2),'\n\r'))
  cat(paste0('Median lambda : ',median((obj$lambda),na.rm = T),'\n\r'))
  cat(paste0('Nr taxa identified by DS-FDR as differentially abundant: ',length(obj$dsfdr_rejected),'\n\r'))
}

#' Print information on the selected reference set, given an object from type reference.selection.object
#'
#' @param obj Object returned by dacomp.select.references(...)
#'
#' @return
#' @export
#'
#' @examples
print.dacomp.reference.selection.object = function(obj){
  cat('DACOMP reference selection object \n\r')
  cat(paste0('Thereshold for selecting reference taxa: ',obj$median_SD_threshold,'\n\r'))
  cat(paste0('For dacomp.select_references(), threshold given in units of medianSD score \n\r'))
  cat(paste0('For dacomp.select_references.by.split(), threshold given in units of -log(PV) \n\r'))
  cat(paste0('Nr. Selected references: ',length(obj$selected_references),'\n\r'))
  cat(paste0('Minimal number of counts observed in reference taxa (across subjects): ',(obj$selected_MinAbundance),'\n\r'))
}

#' Print the results of an object from type reference.validation.result.object
#'
#' @param obj Object returned by dacomp.check_reference_set_is_valid.k_groups
#'
#' @return
#' @export
#'
#' @examples
print.dacomp.reference.validation.result.object = function(obj){
  cat('DACOMP reference validation object \n\r')
  cat('Printing P-values by reference validation test: \n\r')
  cat('############################################### \n\r')
  for(i in 1:length(names(obj))){
    cat(paste0(names(obj)[i],' = ',round(obj[[i]],5),'\n\r'))
  }
}


