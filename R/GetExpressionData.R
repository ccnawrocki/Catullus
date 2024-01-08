#' Get the expression data from an experiment. 
#'
#' This function returns a sparse matrix in dgCMatrix format that contains the 
#' expression data from a given experiment. The user can optionally provide 
#' specific gene IDs and/or cell IDs to make more specific queries from the 
#' expression data. If this is done, the function has the side effect of 
#' creating global variables called \code{Catullus_genes}, which stores a 
#' vector of genes, and \code{Catullus_cells}, which stores a vector of cells.
#'
#' @param exp_object A \code{tiledbsoma} experiment object. 
#' @param genes The string name of a singular gene of interest or a vector 
#' containing the string names of many genes of interest. \code{NULL} by 
#' default. 
#' @param cells The string name of a singular cell of interest or a vector 
#' containing the string names of many cells of interest. \code{NULL} by 
#' default. 
#' @param X_slot The string name of the X layer from which the user wants to 
#' query expression data. Can be \code{"counts"}, \code{"data"}, or 
#' \code{"scale_data"}. \code{"counts"} by default.
#' @export
GetExpressionData <- function(exp_object,
                              genes = NULL,
                              cells = NULL,
                              X_slot = "counts") {
  
  # If there are genes specified, then check they are in the dataset, then 
  # establish a query string for them.
  if (is.null(genes) != T) {
    not_found <- genes[!(genes %in% exp_object$ms$get("RNA")$var$read(column_names="var_id")$concat()$to_data_frame()$var_id)]
    Catullus_genes <<- genes[!(genes %in% not_found)]
    if (length(not_found) > 0) {
      cat("Not in dataset: ", (str_flatten(not_found, collapse = ", ")), "\n")
    }
    gene_query_string <- "var_id %in% Catullus_genes"
  }
  else {
    gene_query_string <- genes
  }
  
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
    var_query = SOMAAxisQuery$new(
      value_filter = gene_query_string
    ), 
    obs_query = SOMAAxisQuery$new(
      value_filter = cell_query_string
    )
  )
  
  # Execute the query on the dataset. 
  mat <- query$to_sparse_matrix(
    collection = "X",
    layer_name = X_slot,
    obs_index = "obs_id",
    var_index = "var_id"
  ) |> as("CsparseMatrix")

  # Return the data. 
  return(mat)
}
