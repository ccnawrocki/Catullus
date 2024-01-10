#' Get differential expression data from an experiment.
#'
#' This function returns the data from specified differential expression testing
#' in a given \code{tiledbsoma} object. 
#' 
#' @param exp_object A \code{tiledbsoma} experiment object. 
#' @param de_name The string name of differential expression testing of interest 
#' within the varm slot of the given \code{tiledbsoma} object.
#' @export
GetDEData <- function(exp_object, 
                      de_name) {
  
  # Query the DE data from the SOMA. 
  de_df <- exp_object$ms$get("RNA")$varm$get(de_name)$read()$concat()$to_data_frame()
  
  # Return the data as a data frame.
  return(de_df)
}
