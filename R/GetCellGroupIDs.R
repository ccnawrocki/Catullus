#' Get the cell IDs of a group of interest from an experiment. 
#'
#' This function returns a vector that contains the cell IDs of a group of cells
#' that was specified by a query on an experiment's metadata. If no query is 
#' made, then all cell IDs will be returned.
#'
#' @param exp_object A \code{tiledbsoma} experiment object. 
#' @param cell_query_string A string representing a logical expression that can 
#' be used to filter for cells of a specific group. \code{NULL} by default. 
#' @export 
GetCellGroupIDs <- function(exp_object, 
                            cell_query_string = NULL) {
  
  # Do the querying. 
  v <- exp_object$obs$read(value_filter = cell_query_string)$concat()$obs_id |> as.vector()
  
  # Return the vector of cell IDs. 
  return(v)
}

