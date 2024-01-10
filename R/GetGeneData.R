#' Get the gene data from an experiment. 
#'
#' This function returns a data frame that contains the gene data from a given 
#' experiment. The user can optionally provide specific variables and/or gene 
#' names to make more specific queries from the gene data.
#'
#' @param exp_object A \code{tiledbsoma} experiment object. 
#' @param variables The string name of a singular variable of interest or a 
#' vector containing the string names of many variables of interest. \code{NULL}
#' by default. 
#' @param genes The string name of a singular gene of interest or a vector 
#' containing the string names of many genes of interest. \code{NULL} by 
#' default. 
#' @param id_rownames A boolean indicating whether to make the row names of the 
#' resulting data frame the gene names. \code{TRUE} by default. 
#' @export
GetGeneData <- function(exp_object, 
                        variables = NULL, 
                        genes = NULL, 
                        id_rownames = T) {
  
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
  
  # Query based on the genes and variables that the user wants.
  if (is.null(variables) != T) {
    df <- exp_object$ms$get("RNA")$var$read(value_filter = gene_query_string, 
                                            column_names = c(variables,"var_id"))$concat()$to_data_frame()
  }
  else {
    df <- exp_object$ms$get("RNA")$var$read(value_filter = gene_query_string, 
                                            column_names = c(variables))$concat()$to_data_frame()
  }
  
  # Make the cell IDs the row names, if necessary. 
  if (id_rownames == T) {
    df <- tibble::column_to_rownames(df, var="var_id")
  }
  
  # Return the gene data as a data frame. 
  return(df)
}
