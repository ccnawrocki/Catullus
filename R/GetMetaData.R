#' Get the metadata from an experiment. 
#'
#' This function returns a data frame that contains the metadata from a given 
#' experiment. The user can optionally provide specific variables and/or cell 
#' IDs to make more specific queries from the metadata.
#'
#' @param exp_object A \code{tiledbsoma} experiment object. 
#' @param variables The string name of a singular variable of interest or a 
#' vector containing the string names of many variables of interest. \code{NULL} 
#' by default. 
#' @param cells The string name of a singular cell of interest or a vector 
#' containing the string names of many cells of interest. \code{NULL} by 
#' default. 
#' @param id_rownames A boolean indicating whether to make the row names of the 
#' resulting data frame the cell IDs. \code{TRUE} by default. 
#' @export
GetMetaData <- function(exp_object, 
                        variables = NULL, 
                        cells = NULL, 
                        id_rownames = T) {
  
  # Read the variables that the user wants. 
  if (is.null(variables) != T) {
    df <- exp_object$obs$read(column_names = c(variables,"obs_id"))$concat()$to_data_frame()
  }
  else {
    df <- exp_object$obs$read()$concat()$to_data_frame()
  }
  
  # Subset for the cells that the user wants.
  if (is.null(cells) != T) {
    df <- subset(df, obs_id %in% cells)
  }
  
  # Make the cell IDs the row names, if necessary. 
  if (id_rownames == T) {
    df <- tibble::column_to_rownames(df, var="obs_id")
  }
  
  # Return the data. 
  return(df)
}
