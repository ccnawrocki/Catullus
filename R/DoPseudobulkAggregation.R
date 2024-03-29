#' Aggregate expression data in preparation for pseudobulk DE testing.
#' 
#' This function returns a sparse matrix in dgCMatrix format that contains the 
#' expression data from a given experiment, aggregated for replicates (usually  
#' samples or libraries) within DE testing conditions. The user can optionally 
#' provide specific gene IDs and/or cell IDs to make more specific queries from 
#' the expression data. If this is done, the function has the side effect of 
#' creating global variables called \code{Catullus_genes} that stores the vector 
#' of genes and \code{Catullus_cells} that stores the vector of cells.  
#' 
#' @param input_object A \code{tiledbsoma} experiment object or a list 
#' containing counts as the first entry and metadata as the second entry. 
#' @param condition_var The string name of a variable in the object's metadata 
#' that contains the condition labels that will be used for DE testing. 
#' @param replicate_var The string name of a variable in the object's metadata
#' that contains the replicate labels. These labels are usually the sample or 
#' library that the individual cells came from. 
#' @param separator A string value to separate the condition and the replicate 
#' labels in the column names of the resulting matrix. It is recommended to use 
#' a value that does not appear in any of the \code{condition_var} or 
#' \code{replicate_var} values. \code{":"} by default.
#' @param cells The string name of a singular cell of interest or a vector 
#' containing the string names of many cells of interest. \code{NULL} by 
#' default. Only applies if \code{input_object} is a \code{tiledbsoma} object. 
#' @param genes The string name of a singular gene of interest or a vector 
#' containing the string names of many genes of interest. \code{NULL} by 
#' default. Only applies if \code{input_object} is a \code{tiledbsoma} object. 
#' @export
DoPseudobulkAggregation <- function(input_object, 
                                    condition_var, 
                                    replicate_var,
                                    separator = ":",
                                    cells = NULL, 
                                    genes = NULL) {
  
  # Using Catullus functions to get the expression data and metadata, if needed.
  if (is(input_object, "SOMAExperiment")) {
    exp <- Catullus::GetExpressionData(exp_object=input_object, 
                                       cells=cells, 
                                       genes=genes)
    meta <- Catullus::GetMetaData(exp_object=input_object,
                                  variables = c(condition_var, replicate_var),
                                  cells=cells)
  }
  else if (is(input_object, "list")) {
    exp <- input_object[[1]]
    meta <- input_object[[2]]
  }
  else {
    cat("Not a valid input\n")
  }
  
  # Creating a design matrix based on the conditions and replicates. 
  des_mat <- Matrix::sparse.model.matrix(~ 0 + get(condition_var):get(replicate_var), data = meta) |> methods::as("CsparseMatrix")
  
  # Doing the following computation: t(exp) %*% des_mat. This is the same as the 
  # cross-product, and the crossprod function happens to be faster. This 
  # computation does the same thing as splitting the expression data into chunks
  # and taking the necessary sums, which is more intuitive. However, doing that 
  # is super memory-intensive and requires looping. This is much more efficient 
  # and it is vectorized (I think). Reading the Libra package's code helped me 
  # figure this out.
  agg_mat <- Matrix::crossprod(exp, des_mat) |> methods::as("CsparseMatrix")
  
  # Correcting the column names of the result.
  colnames(agg_mat) <- gsub(pattern = "(get\\(condition_var\\))|(get\\(replicate_var\\))", 
                            replacement = "",
                            x = colnames(agg_mat))
  if (separator != ":") {
    colnames(agg_mat) <- gsub(pattern = ":", 
                              replacement = separator, 
                              x = colnames(agg_mat))
  }
  
  # Removing any columns with zero counts.
  idx <- which(Matrix::colSums(agg_mat) == 0)
  agg_mat <- agg_mat[,-idx]
  
  # Returning the result.
  return(agg_mat)
}
