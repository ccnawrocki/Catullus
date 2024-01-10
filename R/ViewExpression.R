#' View the expression of a gene or a gene list.
#' 
#' This function returns a plot in which a gene's expression values appear for 
#' every cell on an experiment's UMAP. The function can handle singular genes or 
#' a vector of genes. the function has the side effect of creating global 
#' variables called \code{Catullus_genes}, which stores a vector of genes, and 
#' \code{Catullus_cells}, which stores a vector of cells.
#' 
#' @param exp_object A \code{tiledbsoma} experiment object. 
#' @param genes The string name of a singular gene of interest or a vector 
#' containing the string names of many genes of interest. 
#' @param cells The string name of a singular cell of interest or a vector 
#' containing the string names of many cells of interest. \code{NULL} by 
#' default. 
#' @param X_slot The string name of the X layer for which the user wants to see 
#' expression. Can be \code{"counts"}, \code{"data"}, or \code{"scale_data"}. 
#' \code{"data"} by default. 
#' @param split_var The string name of the metadata column that the user wants 
#' to split the plot by. \code{NULL} by default. 
#' @param label_var The string name of the metadata column that the user wants 
#' to use to label the data. \code{NULL} by default.  
#' @param order_expression A boolean indicating whether to order the cells by 
#' decreasing expression. \code{TRUE} by default. 
#' @param colors A vector of a low color and a high color, in that order, 
#' to match the expression values to. \code{"snow2"} and \code{"blue"} by 
#' default. 
#' @param pt_size A float indicating the size of each point. 0.2 by default. 
#' @param label_size A float indicating the size the labels. \code{3} by 
#' default. 
#' @param plot_style A string name of the plot style that the user wants. Can be 
#' \code{"bare"}, \code{"standard"}, or \code{"detailed"}. \code{"bare"} by 
#' default. 
#' @param method A string name of the method that the user wants to employ. Can 
#' be \code{"print"}, which will simply print the plot or plots. Can also be 
#' \code{"return"}, which will return a list of plots that the user can then 
#' interact with later. \code{"print"} by default.      
#' @export
ViewExpression <- function(exp_object, 
                           genes,
                           cells = NULL, 
                           X_slot = "data",
                           split_var = NULL, 
                           label_var = NULL, 
                           order_expression = T,
                           colors = c("snow2","blue"), 
                           pt_size = 0.2, 
                           label_size = 3, 
                           plot_style = "bare",
                           method = "print") {
  
  # Querying the dimensional reduction data. 
  umap_df <- Catullus::GetDimRedData(exp_object = exp_object,
                                     red_name = "umap", 
                                     cells = cells)
  # Querying the expression data.
  exp_mat <- Catullus::GetExpressionData(exp_object = exp_object, 
                                         genes = genes, 
                                         cells = cells, 
                                         X_slot = X_slot)
  # Querying the metadata, if needed, and combining all data into a data frame.
  meta_vars <- c(split_var, label_var) |> stats::na.omit()
  if (is.null(meta_vars) != T) {
    meta_df <- Catullus::GetMetaData(exp_object = exp_object, 
                                     variables = meta_vars, 
                                     cells = cells, 
                                     id_rownames = T)
    umap_df <- cbind(umap_df, as.matrix(exp_mat), meta_df)
  }
  else {
    umap_df <- cbind(umap_df, as.matrix(exp_mat))
  }
  
  # Getting the position data for the labels, if needed. 
  if (is.null(label_var) != T) {
    lab_coords <- data.frame("umap_1" = tapply(umap_df[["umap_1"]], umap_df[[label_var]], median),
                             "umap_2" = tapply(umap_df[["umap_2"]], umap_df[[label_var]], median))
    lab_coords[[label_var]] <- rownames(lab_coords)
  }
  
  # Making the plots.
  plots <- list()
  for (gene in Catullus_genes) {
    
    # Ordering the data so the cells with high expression are in front, if needed.
    if (order_expression == T) {
      umap_df <- umap_df[order(umap_df[[gene]]),]
    }
    
    if (is.null(label_var) != T) {
      p <- ggplot2::ggplot(mapping=ggplot2::aes(x=umap_1, y=umap_2, color=get(gene), label=get(label_var))) + 
        ggplot2::geom_point(data=umap_df, size=pt_size) + 
        ggplot2::scale_color_gradient(low=colors[1], high=colors[2]) + 
        ggplot2::labs(color=gene) +
        ggplot2::geom_text(data=lab_coords, color="black", size=label_size)
    }
    else {
      p <- ggplot2::ggplot(mapping=ggplot2::aes(x=umap_1, y=umap_2, color=get(gene))) + 
        ggplot2::geom_point(data=umap_df, size=pt_size) + 
        ggplot2::scale_color_gradient(low=colors[1], high=colors[2]) + 
        ggplot2::labs(color=gene)
    }
    if (is.null(split_var) != T) {
      p <- p + ggplot2::facet_grid(cols=vars((get(split_var))))
    }
    if (plot_style == "bare") {
      theme_info <- ggplot2::theme(legend.position = "right",
                                   panel.grid = ggplot2::element_line(linetype = "blank"),
                                   panel.background = ggplot2::element_rect(fill="white", color="black"),
                                   strip.background = ggplot2::element_rect(fill = "white",  color="black"))
    }
    if (plot_style == "detailed") {
      theme_info <- ggplot2::theme(legend.position = "right",
                                   panel.grid.major = ggplot2::element_line(linetype = "dotted", color="black"),
                                   panel.grid.minor = ggplot2::element_line(linetype = "blank"),
                                   panel.background = ggplot2::element_rect(fill="grey", color="black"),
                                   strip.background = ggplot2::element_rect(fill = "grey", color = "black"))
    }
    if (plot_style == "standard") {
      theme_info <- ggplot2::theme()
    }
    p <- p + theme_info + ggplot2::coord_fixed()
    
    # Printing or saving and returning the plots.
    if (method == "print") {
      print(p)
    }
    if (method == "return") {
      plots[[gene]] <- p
    }
  }
  if (method == "return") {
    return(plots)
  }
}
