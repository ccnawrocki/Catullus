#' View one metadata variable or many.
#'
#' This function returns a plot in which each cell's status across one or more 
#' metadata variables is shown on an experiment's UMAP. The function can color, 
#' shape, and split cells based on metadata variables from the experiment. 
#'
#' @param exp_object A \code{tiledbsoma} experiment object. 
#' @param color_var The string name of the metadata column that the user wants 
#' to color the plot by. 
#' @param cells The string name of a singular cell of interest or a vector 
#' containing the string names of many cells of interest. \code{NULL} by 
#' default. 
#' @param labeling A boolean indicating whether to label the cells by the levels
#' in the metadata column that the plot is colored by. \code{FALSE} by default.  
#' @param legend A boolean indicating whether to include a legend for the levels 
#' in the metadata column that the plot is colored by. \code{TRUE} by default. 
#' @param split_var The string name of the metadata column that the user wants 
#' to split the plot by. \code{NULL} by default. 
#' @param shape_var The string name of the metadata column that the user wants 
#' to shape the plot by. \code{NULL} by default. 
#' @param pt_size A float indicating the size of each point. \code{0.2} by 
#' default. 
#' @param label_size A float indicating the size the labels. \code{3} by 
#' default. 
#' @param plot_style A string name of the plot style that the user wants. Can be 
#' \code{"bare"}, \code{"standard"}, or \code{"detailed"}. \code{"bare"} by 
#' default. 
#' @export
ViewMeta <- function(exp_object, 
                     color_var, 
                     cells = NULL, 
                     labeling = F, 
                     legend = T, 
                     split_var = NULL, 
                     shape_var = NULL, 
                     pt_size = 0.2, 
                     label_size = 3, 
                     plot_style = "bare") {
  
  # Querying the dimensional reduction data. 
  umap_df <- Catullus::GetDimRedData(exp_object = exp_object,
                                     red_name = "umap", 
                                     cells = cells)
  # Querying the metadata, and combining all data into a data frame.
  meta_vars <- c(color_var, split_var, shape_var) |> stats::na.omit()
  meta_df <- Catullus::GetMetaData(exp_object = exp_object, 
                                   variables = meta_vars, 
                                   cells = cells, 
                                   id_rownames = T)
  umap_df <- cbind(umap_df, meta_df)
  
  # Getting the position data for the labels, if needed. 
  if (labeling == T) {
    lab_coords <- data.frame("umap_1" = tapply(umap_df[["umap_1"]], umap_df[[color_var]], median),
                             "umap_2" = tapply(umap_df[["umap_2"]], umap_df[[color_var]], median))
    lab_coords[[color_var]] <- rownames(lab_coords)
  }
  
  # Making the plot. 
  if (is.null(shape_var) != T) {
    p <- ggplot2::ggplot(mapping = aes(x=umap_1, y=umap_2, color=get(color_var))) + 
      geom_point(data=umap_df, size=pt_size, aes(shape=get(shape_var))) + 
      labs(color=color_var, shape=shape_var) 
  }
  else {
    p <- ggplot2::ggplot(mapping = aes(x=umap_1, y=umap_2, color=get(color_var))) + 
      geom_point(data=umap_df, size=pt_size) + 
      labs(color=color_var) 
  }
  if (labeling == T) {
    p <- p + geom_text(data=lab_coords, color="black", size=label_size, aes(label=get(color_var)))
  }
  if (legend == F) {
    p <- p + guides(color = "none")
  }
  if (is.null(split_var) != T) {
    p <- p + facet_grid(cols=vars((get(split_var))))
  }
  if (plot_style == "bare") {
    theme_info <- theme(legend.position = "right",
                        panel.grid = element_line(linetype = "blank"),
                        panel.background = element_rect(fill="white", color="black"),
                        strip.background = element_rect(fill = "white",  color="black"))
  }
  if (plot_style == "detailed") {
    theme_info <- theme(legend.position = "right",
                        panel.grid.major = element_line(linetype = "dotted", color="black"),
                        panel.grid.minor = element_line(linetype = "blank"),
                        panel.background = element_rect(fill="grey", color="black"),
                        strip.background = element_rect(fill = "grey", color = "black"))
  }
  if (plot_style == "standard") {
    theme_info <- theme()
  }
  p <- p + theme_info + coord_fixed()
  
  # Returning the plot.
  return(p)
}
