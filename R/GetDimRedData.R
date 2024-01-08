#' Get the dimensional reduction data from an experiment. 
#'
#' This function returns a data frame that contains the dimensional reduction 
#' data from a given experiment. The user can optionally provide specific 
#' reduction names and/or cell IDs to make more specific queries from the 
#' data.
#'
#' @param exp_object A \code{tiledbsoma} experiment object. 
#' @param red_name The string name of the dimensional reduction of interest. Can
#' be \code{"umap"} or \code{"pca"}. \code{"umap"} by default. 
#' @param cells The string name of a singular cell of interest or a vector 
#' containing the string names of many cells of interest. \code{NULL} by 
#' default. 
#' @export
GetDimRedData <- function(exp_object, 
                          red_name = "umap", 
                          cells = NULL) {
  
  # If there are cells specified, then check they are in the dataset, then 
  # establish a query string for them.
  if (is.null(cells) != T) {
    not_found <- cells[!(cells %in% exp_object$obs$read(column_names = "obs_id")$concat()$to_data_frame()$obs_id)]
    Catullus_cells <<- cells[!(cells %in% not_found)]
    if (length(not_found) > 0) {
      cat("Not in dataset: ", (str_flatten(not_found, collapse = ", ")), "\n")
    }
    cell_query_string <- "obs_id %in% Catullus_cells"
  }
  else {
    cell_query_string <- cells
  }
  
  # Form a query from the dataset.
  query <- exp_object$axis_query(
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(
      value_filter = cell_query_string
    )
  )
  
  # Execute the query on the dataset. 
  mat <- query$to_sparse_matrix(
    collection = "obsm",
    layer_name = paste("X_", red_name, sep = ""),
    obs_index = "obs_id"
  )
  
  # Convert the data to a data frame. 
  df <- data.frame(mat[,1], mat[,2])
  colnames(df) <- paste(c(red_name, red_name), c(1, 2), sep = "_")
  
  # Returning the data. 
  return(df)
}
