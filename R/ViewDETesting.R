#' View differential expression testing results. 
#'
#' This function returns a volcano plot, displaying the results of specified 
#' differential expression testing within a \code{tiledbsoma} object. 
#' 
#' @param input_object A \code{tiledbsoma} experiment object, a data frame, or a
#' list returned by \code{DoDETesting}.
#' @param de_name The string name of differential expression testing of interest 
#' within the varm slot of the given \code{tiledbsoma} object. \code{NULL} by 
#' default. 
#' @param title A title for the plot. \code{NULL} by default.
#' @param repel_overlaps An integer specifying how many genes can be labeled. 
#' \code{50} by default. 
#' @param lfc_cutoff The log2FC cutoff for a gene to be labeled on the plot. 
#' \code{1} by default. 
#' @param p_adj_cutoff The p_adj cutoff for a gene to be labeled on the plot. 
#' \code{0.05} by default.
#' @param cols A vector of three colors corresponding to the how the colors will
#' appear on the plot, from left to right. \code{c("blue", "grey", "red")} by 
#' default.
#' @export
ViewDETesting <- function(input_object, 
                          de_name = NULL, 
                          title = NULL, 
                          repel_overlaps = 50, 
                          lfc_cutoff = 1, 
                          p_adj_cutoff = 0.05, 
                          cols = c("blue", "grey", "red")) {
  
  # Querying the DE data, if necessary.
  if (is(input_object, "SOMAExperiment")) {
    de_df <- Catullus::GetDEData(exp_object = input_object, 
                                 de_name = de_name)
  }
  else if (is(input_object, "data.frame")) {
    de_df <- input_object
  }
  else if (is(input_object, "list")) {
    de_df <- input_object[[3]]
  }
  else {
    cat("Not a valid input\n")
  }
  
  # Getting the comparison levels in the correct order.
  de_df <- de_df[order(de_df$log2FC, decreasing = F),]
  conds <- unique(de_df$group)
  
  # Adding information to the DE data for plotting. 
  de_df$delabel <- 'Neither'
  de_df$delabel[de_df$p_adj<p_adj_cutoff & de_df$log2FC < (-1*lfc_cutoff)] <- conds[1]
  de_df$delabel[de_df$p_adj<p_adj_cutoff & de_df$log2FC > lfc_cutoff] <- conds[2]
  de_df$lbl <- NA
  de_df$lbl[de_df$delabel != "Neither"] <- de_df$gene[de_df$delabel != 'Neither']
  
  # Making the plot.
  p <- ggplot2::ggplot(data = de_df, 
                       mapping = ggplot2::aes(x=log2FC, 
                                              y=-log10(p_adj), 
                                              col=delabel, 
                                              label=lbl)) +
    ggplot2::geom_point() + 
    ggplot2::scale_color_manual(values = cols, 
                                labels = c(conds[1], "Neither", conds[2])) +
    ggplot2::geom_vline(xintercept=c(-1*lfc_cutoff, lfc_cutoff), col="black") + 
    ggplot2::geom_hline(yintercept=-log10(p_adj_cutoff), col="black") + 
    ggrepel::geom_text_repel(max.overlaps = repel_overlaps, show.legend = F)
  if (is.null(title) != T) {
    p <- p +
      ggplot2::labs(title=title)
  }
  else {
    p <- p + 
      ggplot2::labs(title=paste(conds[1], "vs", conds[2], sep = " "))
  }
  theme_info <- ggplot2::theme(legend.position = "right",
                               panel.grid = ggplot2::element_line(linetype = "blank"),
                               panel.background = ggplot2::element_rect(fill="white", color="black"),
                               strip.background = ggplot2::element_rect(fill = "white",  color="black"))
  p <- p + theme_info
  
  # Returning the plot.
  return(p)
}
