#' Perform differential expression testing.
#'
#' This function returns the results of testing for differentially-expressed 
#' genes across two conditions, using various methodologies/models. The function
#' has pseudobulking capabilities built in and supports the following tests: 
#' Wilcoxon Rank Sum (with and without \code{presto}), DESeq2, limma, edgeR, and 
#' glmGamPoi.
#'
#' @param exp_object A \code{tiledbsoma} experiment object. 
#' @param condition_var The string name of the metadata column that contains the
#' condition that is being tested.
#' @param comparison Either a singular value that appears in the 
#' \code{condition_var} column or a vector of two values that appear in the 
#' \code{condition_var} column. If a singular value is given, cells/samples 
#' matching that value will be compared to all other cells/samples. If a vector
#' of two values in given, cells/samples in the two groups will be compared.
#' @param test The string name of the test to use. Can be \code{"wilcox"}, 
#' \code{"wilcox_presto"}, \code{"DESeq2"}, \code{"limma"}, \code{"edgeR"}, or 
#' \code{"glmGamPoi"}. \code{"wilcox"} by default. Note: a very basic design is 
#' used for tests other than \code{"wilcox"} or \code{"wilcox_presto"}. To use 
#' more complex designs, revert back to the packages from where these tests 
#' come. 
#' @param replicate_var The string name of the metadata column that will be used
#' to make replicates for pseudobulking purposes. If \code{NULL}, then 
#' pseudobulking will not be performed. \code{NULL} by default.
#' @param method The string name of the method that the user wants to employ. 
#' Can be \code{"save"}, which will save the results of the testing inside of 
#' the given \code{exp_object}, or \code{"return"}, which will return the 
#' results of the testing as a data frame. \code{"return"} by default. 
#' @param separator A string value to separate the condition and the replicate 
#' labels in the column names of the resulting matrix. It is recommended to use 
#' a value that does not appear in any of the \code{condition_var} or 
#' \code{replicate_var} values. \code{":"} by default.
#' @param cells The string name of a singular cell of interest or a vector 
#' containing the string names of many cells of interest. \code{NULL} by 
#' default. 
#' @param genes The string name of a singular gene of interest or a vector 
#' containing the string names of many genes of interest. 
#' @export
DoDETesting <- function(exp_object, 
                        condition_var,
                        comparison, 
                        test = "wilcox",
                        replicate_var = NULL,
                        method = "return", 
                        separator = ":",
                        cells = NULL,
                        genes = NULL) {
  
  # Make the comparison levels easier to deal with.
  if (length(comparison) == 1) {
    comparison <- c(comparison, "Other")
  }

  # Pseudobulk, if necessary.
  if (is.null(replicate_var) != T) {
    cat("Pseudobulking...")
    cts <- Catullus::DoPseudobulkAggregation(exp_object = exp_object, 
                                             condition_var = condition_var, 
                                             replicate_var = replicate_var, 
                                             separator = separator, 
                                             cells = cells, 
                                             genes = genes)
    meta <- data.frame(row.names = colnames(cts))
    meta[[condition_var]] <- ifelse(test = (stringr::str_split(rownames(meta), separator, simplify = T)[,1] == comparison[1]), 
                                    yes = comparison[1], 
                                    no = comparison[2])
    cat(" DONE\n")
  }
  
  # If pseudobulking is not needed, just query the data. 
  else {
    cts <- Catullus::GetExpressionData(exp_object = exp_object, 
                                       genes = genes, 
                                       cells = cells, 
                                       X_slot = "counts") |> Matrix::t()
    meta <- GetMetaData(exp_object = exp_object, 
                        variables = c(condition_var), 
                        cells = cells)
    meta[[condition_var]] <- ifelse(test = (meta[[condition_var]] == comparison[1]), 
                                    yes = comparison[1], 
                                    no = comparison[2])
  }
  
  # Next, do the testing method indicated by the user. 
  if (test == "wilcox") {
    
    cat("Performing Wilcoxon Rank Sum Tests...")
    
    # Normalizing the data.
    cts_n <- Matrix::t(Matrix::t(cts) / Matrix::colSums(cts))
    cts_n <- log1p(cts_n*10000)
    
    # Getting the normalized counts for each group.
    g1_samps <- subset(meta, get(condition_var) == comparison[1]) |> rownames()
    g2_samps <- subset(meta, get(condition_var) == comparison[2]) |> rownames()
    g1_n <- cts_n[,g1_samps] 
    g2_n <- cts_n[,g2_samps]
    
    # Doing the test.
    pvals <- matrixTests::row_wilcoxon_twosample(x = as.matrix(g1_n), 
                                                 y = as.matrix(g2_n),
                                                 null = 0, exact = F)$pvalue
    pvals_adj <- stats::p.adjust(p = pvals, 
                                 method = "fdr")
    
    # Getting LFC. 
    sizefactor <- Matrix::colSums(cts)/median(Matrix::colSums(cts))
    g1_sizefactor <- sizefactor[g1_samps]
    g2_sizefactor <- sizefactor[g2_samps]
    g1_cts <- cts[,g1_samps]
    g2_cts <- cts[,g2_samps]
    g1_gene_cts <- Matrix::rowSums(g1_cts)
    g2_gene_cts <- Matrix::rowSums(g2_cts)
    l2fc <- log2((1+g2_gene_cts)/(1+sum(g2_sizefactor)))-log2((1+g1_gene_cts)/(1+sum(g1_sizefactor)))
    
    # Creating the output table. 
    de_df <- data.frame("log2FC"=l2fc, 
                        "p"=pvals,
                        "p_adj"=pvals_adj, 
                        "group"=ifelse(test = (l2fc < 0), 
                                       yes = comparison[1],
                                       no = comparison[2]), 
                        "gene"=names(l2fc))
    
    # Getting proportions, if necessary.
    if (is.null(replicate_var)) {
      de_df$prop_1 <- Matrix::rowSums(g1_cts > 0)/ncol(g1_cts)
      de_df$prop_2 <- Matrix::rowSums(g2_cts > 0)/ncol(g2_cts)
    }
    
    cat(" DONE\n")
    
  }
  else if (test == "wilcox_presto") {
    
    cat("Performing Wilcoxon Rank Sum Tests using presto...")
    
    # Normalizing the data.
    cts_n <- Matrix::t(Matrix::t(cts) / Matrix::colSums(cts))
    cts_n <- log1p(cts_n*10000)
    
    # Doing the test.
    de <- presto::wilcoxauc(cts_n, t(meta))
    
    # Creating the output table.
    de <- de[-(1:(0.5*nrow(de))),]
    de_df <- data.frame("log2FC"=de$logFC, 
                        "p"=de$pval, 
                        "p_adj"=de$padj, 
                        row.names=de$feature)
    de_df$group <- ifelse(test = (de_df$log2FC < 0), 
                          yes = comparison[1], 
                          no = comparison[2])
    de_df$gene <- rownames(de_df)
    
    # Getting proportions, if necessary.
    g1_samps <- subset(meta, get(condition_var) == comparison[1]) |> rownames()
    g2_samps <- subset(meta, get(condition_var) == comparison[2]) |> rownames()
    g1_cts <- cts[,g1_samps]
    g2_cts <- cts[,g2_samps]
    if (is.null(replicate_var)) {
      de_df$prop_1 <- Matrix::rowSums(g1_cts > 0)/ncol(g1_cts)
      de_df$prop_2 <- Matrix::rowSums(g2_cts > 0)/ncol(g2_cts)
    }
    
    cat(" DONE\n")
    
  }
  else if (test == "DESeq2") {
    
    cat("Performing DESeq2...\n")
    cat("DESeq2 Messages:\n")
    
    # Setting up the object.
    colnames(meta) <- c("condition")
    meta$condition <- as.factor(meta$condition)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, 
                                          colData = meta, 
                                          design = ~condition)
    
    # Filtering and correcting the levels.
    idx <- Matrix::rowSums(DESeq2::counts(dds) >= 10) >= 3
    dds <- dds[idx,]
    dds$condition <- factor(dds$condition, levels = comparison)
    
    # Running the model.
    dds <- DESeq2::DESeq(dds)
    
    # Creating the output table. 
    de_df <- DESeq2::results(dds) |> as.data.frame() |> na.omit()
    colnames(de_df) <- c("baseMean", "log2FC", "lfcSE", "stat", "p", "p_adj")
    de_df$group <- ifelse(test = (de_df$log2FC < 0), 
                          yes = comparison[1], 
                          no = comparison[2])
    de_df$gene <- rownames(de_df)
    
    cat("DONE\n")
    
  }
  else if (test == "limma") {
    
    cat("Performing limma...")
    
    # Creating the object, normalizing, and filtering. 
    de0 <- edgeR::DGEList(counts=cts,group=factor(meta[[condition_var]])) |> 
      edgeR::calcNormFactors()
    cutoff <- 1
    drop <- which(apply(edgeR::cpm(de0), 1, max) < cutoff)
    de <- de0[-drop,]
    
    # Creating the necessary design matrix. 
    group <- as.factor(meta[[condition_var]])
    des_m <- model.matrix(~ 0 + group)

    # Running the limma model. 
    y <- limma::voom(de, des_m)
    fit <- limma::lmFit(y, des_m)
    
    # Doing the tests and calculating statistics.
    x <- paste("group", comparison[2], "-", "group", comparison[1], sep = "")
    contr <- limma::makeContrasts(contrasts = c(x), 
                                  levels = colnames(stats::coef(fit)))
    de <- limma::contrasts.fit(fit, contr) |> 
      limma::eBayes()
    
    # Creating the output table. 
    de_df <- limma::topTable(de, sort.by = "P", n = Inf)
    colnames(de_df) <- c("log2FC", "AveExpr", "t", "p", "p_adj", "B")
    de_df$group <- ifelse(test = (de_df$log2FC < 0), 
                          yes = comparison[1], 
                          no = comparison[2])
    de_df$gene <- rownames(de_df)
    
    cat(" DONE\n")
    
  }
  else if (test == "edgeR") {
    
    cat("Performing edgeR...")
    
    # Creating the object and normalizing. 
    de0 <- edgeR::DGEList(counts=cts,group=factor(meta[[condition_var]])) |> 
      edgeR::calcNormFactors()
    
    # Creating the necessary design matrix. 
    des_m <- model.matrix(~ 0 + de0$samples$group)
    colnames(des_m) <- levels(de0$samples$group)
    
    # Running the dispersion estimations and the exact test.
    de <- edgeR::estimateGLMCommonDisp(de0, des_m) |> 
      edgeR::estimateGLMTrendedDisp(des_m, method="power") |> 
      edgeR::estimateGLMTagwiseDisp(des_m) |>
      edgeR::exactTest(pair=c(comparison[1],comparison[2]))
    
    # Creating the output table. 
    de_df <- edgeR::topTags(de, n = Inf)$table
    colnames(de_df) <- c("log2FC", "logCPM", "p", "p_adj")
    de_df$group <- ifelse(test = (de_df$log2FC < 0), 
                          yes = comparison[1], 
                          no = comparison[2])
    de_df$gene <- rownames(de_df)
    
    cat(" DONE\n")
    
  }
  else if (test == "glmGamPoi") {
    
    cat("Performing glmGamPoi...")
    
    # Fitting the glmGamPoi model. 
    fit <- glmGamPoi::glm_gp(data = as.matrix(cts), 
                             col_data = meta, 
                             design = ~ get(condition_var))
    
    # Doing the tests. 
    de_df <- glmGamPoi::test_de(fit, contrast = get(colnames(fit$model_matrix)[2]))
    
    # Creating the output table.
    colnames(de_df) <- c("gene", "p", "p_adj", "F_stat", "df1", "df2", "log2FC")
    de_df$group <- ifelse(test = (de_df$log2FC < 0), 
                          yes = comparison[1], 
                          no = comparison[2])
    rownames(de_df) <- de_df$gene
    
    cat(" DONE\n")
    
  }
  else {
    cat("Not a valid test.\n")
  }
  
  # Returning or saving the results.
  if (method == "save") {
    
    # Naming the result.
    de_name <- paste(comparison[1], "vs", comparison[2], test, "de", sep = "_")
    
    # Writing and registering the result to the SOMA object on disk.
    tmp <- tiledbsoma::SOMAExperimentOpen(uri = exp_object$uri, mode = "WRITE")
    sdf <- tiledbsoma::write_soma(de_df, uri = de_name, soma_parent = tmp$ms$get("RNA")$varm)
    tmp$ms$get("RNA")$varm$set(sdf, name = de_name)
    cat("Saved ", de_name, " in varm.\n", sep = "")
    
    # Returning the updated SOMA object.
    uri_to_open <- exp_object$uri
    tmp2 <- tiledbsoma::SOMAExperimentOpen(uri_to_open)
    remove(tmp)
    return(tmp2)
  }
  else {
    
    # Simply returning the df. 
    return(de_df)
  }
}
