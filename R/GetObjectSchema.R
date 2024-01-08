#' Get the schema of an experiment object. 
#'
#' This function will print or return the basic schema of a given experiment 
#' object. Its main use is to allow users to examine the contents of an object 
#' and to see the various names of measurements throughout it.  
#'
#' @param exp_object A \code{tiledbsoma} experiment object. 
#' @param method A string name of the method that the user wants to employ. Can 
#' be \code{"print"}, which will simply print the schema. Can also be 
#' \code{"return"}, which will return a list object that the user can then parse 
#' through. \code{"print"} by default.   
#' @export
GetObjectSchema <- function(exp_object, 
                            method = "print") {
  
  schema_list <- list()
  
  ## Metadata
  meta_info <- exp_object$obs$colnames()
  if (method == "return") {
    schema_list[["METADATA"]] <- as.list(meta_info)
  }
  meta_info <- exp_object$obs$schema()$ToString()
  
  ## Expression Data
  express_slots <- exp_object$ms$get("RNA")$X$names() 
  if (method == "return") {
    schema_list[["EXPRESSION DATA"]][["slots"]] <- as.list(express_slots)
  }
  express_slots <- str_flatten(express_slots, collapse = ", ")
  express_shape <- exp_object$ms$get("RNA")$X$get("counts")$shape()
  express_shape <- paste(express_shape, c("cells", "genes"), sep = " ")
  if (method == "return") {
    schema_list[["EXPRESSION DATA"]][["shape"]] <- as.list(express_shape)
  }
  express_shape <- str_flatten(express_shape, collapse = ", ")
  
  ## Dimensional Reduction Data 
  dimred_info <- exp_object$ms$get("RNA")$obsm$names()
  if (method == "return") {
    schema_list[["DIMENSIONAL REDUCTION DATA"]] <- as.list(dimred_info)
  }
  dimred_info <- str_flatten(dimred_info, collapse = ", ")
  
  ## Gene Data
  gene_info <- exp_object$ms$get("RNA")$var$colnames()
  if (method == "return") {
    schema_list[["GENE DATA"]] <- as.list(gene_info)
  }
  gene_info <- exp_object$ms$get("RNA")$var$schema()$ToString()
  
  ## Printing or returning the schema for the user. 
  if (method == "return") {
    return(schema_list)
  }
  else {
    cat(paste("METADATA: ", deparse(substitute(exp_object)), "$obs", sep = ""), "\n", meta_info, "\n\n", 
        paste("EXPRESSION DATA: ", deparse(substitute(exp_object)), "$ms$get(\"RNA\")$X", sep = ""), "\n", "slots: ", express_slots, "\n", "shape: ", express_shape, "\n\n",
        paste("DIMESNIONAL REDUCTION DATA: ", deparse(substitute(exp_object)), "$ms$get(\"RNA\")$obsm", sep = ""), "\n", dimred_info, "\n\n", 
        paste("GENE DATA: ", deparse(substitute(exp_object)), "$ms$get(\"RNA\")$var", sep = ""), "\n", gene_info,
        sep = "")
  }
}
