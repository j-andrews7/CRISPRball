# js code to enable/disable tabs.
.utils.js <- function() {
    out <- "
shinyjs.disableTab = function(name) {
  var tab = $('.nav.navbar-nav li a[data-value=\"' + name + '\"]');
  tab.bind('click.tab', function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass('disabled');
}

shinyjs.enableTab = function(name) {
  var tab = $('.nav.navbar-nav li a[data-value=\"' + name + '\"]');
  tab.unbind('click.tab');
  tab.removeClass('disabled');
}
"
    out
}


#' Parse gene summary data for easier plotting and display
#' @param df data.frame of gene summary data. Gene IDs should be in the first column.
#' @param sig.thresh Numeric scalar for significance threshold to consider a gene a hit.
#' @param lfc.thresh Numeric scalar for absolute log fold change threshold to consider a gene a hit.
#' @param positive.ctrl.genes Character vector of gene identifiers to label as positive controls.
#' @param essential.genes Character vector of gene identifiers to label as essential genes.
#' @param depmap.genes data.frame of DepMap gene summary data.
#'
#' @return A data.frame of gene summary with additional, easier to plot, columns added.
#' @author Jared Andrews
#' @export
#' @examples
#' library(CRISPRball)
#' d1.genes <- read.delim(system.file("extdata", "esc1.gene_summary.txt",
#'     package = "CRISPRball"
#' ), check.names = FALSE)
#' out.df <- gene_ingress(d1.genes, 0.05, 0.5)
gene_ingress <- function(df, sig.thresh, lfc.thresh, positive.ctrl.genes = NULL,
                         essential.genes = NULL, depmap.genes = NULL) {
    if (!is.null(essential.genes)) {
        df$essential <- df[[1]] %in% essential.genes
    }

    if (!is.null(positive.ctrl.genes)) {
        df$Positive_Control <- df[[1]] %in% positive.ctrl.genes
    }

    if (!is.null(depmap.genes)) {
        df$DepMap_CRISPR_Essential <- df[[1]] %in%
            depmap.genes$gene_name[depmap.genes$dataset %in%
                c("Chronos_Combined", "Chronos_Score", "Chronos_Achilles") &
                depmap.genes$common_essential == TRUE]

        df$DepMap_CRISPR_Selective <- df[[1]] %in%
            depmap.genes$gene_name[depmap.genes$dataset %in%
                c("Chronos_Combined", "Chronos_Score", "Chronos_Achilles") &
                depmap.genes$strongly_selective == TRUE]

        df$DepMap_RNAi_Essential <- df[[1]] %in% depmap.genes$gene_name[depmap.genes$dataset == "RNAi_merged" &
            depmap.genes$common_essential == TRUE]
        df$DepMap_RNAi_Selective <- df[[1]] %in% depmap.genes$gene_name[depmap.genes$dataset == "RNAi_merged" &
            depmap.genes$strongly_selective == TRUE]
    }

    if ("neg|score" %in% colnames(df)) {
        df <- .RRA_ingress(df, sig.thresh, lfc.thresh)
    }

    return(df)
}


#' Read and parse MAGeCK MLE output gene summary file
#' 
#' This function reads the gene summary file output by \code{mageck mle} and
#' parses it into a list of data.frames, one for each sample. The sample names
#' are extracted from the column names of the input file and used as the names
#' of the list elements.
#' 
#' @param filepath Path to the gene summary file output by \code{mageck mle}.
#' 
#' @return A named list of data.frames containing MAGeCK MLE output,
#'   one for each sample contained in the file.
#' @author Jared Andrews
#' @importFrom utils read.delim
#' @export
#' 
#' @examples
#' library(CRISPRball)
#' mle_gene_summary <- file.path(system.file("extdata", "beta_leukemia.gene_summary.txt", package = "CRISPRball"))
#' gene_data <- read_mle_gene_summary(mageck_mlmle_gene_summarye_file)
read_mle_gene_summary <- function(filepath) {
    # Read the table
    data <- read.delim(filepath, check.names = FALSE)

    # Get the sample names
    samples <- unique(gsub("\\|.*", "", names(data)[-(1:2)]))

    # Create a list to hold the data frames
    dataframes <- list()

    # For each sample, create a data frame and add it to the list
    for (sample in samples) {
        # Find the columns for this sample
        sample_cols <- grep(paste0("^", sample, "\\|"), names(data))

        # Add the 'Gene' and 'sgRNA' columns
        sample_cols <- c(1, 2, sample_cols)

        # Extract the data for this sample
        df <- data[, sample_cols]

        # Remove the sample name from the header
        names(df) <- gsub(paste0("^", sample, "\\|"), "", names(df))

        # Fix the column names to not use '-'
        names(df) <- gsub("-", ".", names(df))
        dataframes[[sample]] <- df
    }

    return(dataframes)
}


#' Rename columns and add additional columns to gene summary data
#' @param df data.frame of gene summary data
#'
#' @return data.frame of gene summary data with renamed columns.
#' @author Jared Andrews
#' @rdname INTERNAL_RRA_ingress
.RRA_ingress <- function(df, fdr.thresh, lfc.thresh) {
    df$LFC <- as.numeric(df$`neg|lfc`)

    df$RRAscore <- apply(df, 1, function(x) {
        x <- split(unname(x), names(x))
        as.numeric(ifelse(as.numeric(x$`neg|score`) < as.numeric(x$`pos|score`), x$`neg|score`, x$`pos|score`))
    })

    df$FDR <- apply(df, 1, function(x) {
        x <- split(unname(x), names(x))
        as.numeric(ifelse(as.numeric(x$`neg|fdr`) < as.numeric(x$`pos|fdr`), x$`neg|fdr`, x$`pos|fdr`))
    })

    df$hit_type <- apply(df, 1, function(x) {
        x <- split(unname(x), names(x))
        ifelse(as.numeric(x$`neg|fdr`) < fdr.thresh & as.numeric(x$LFC) < -as.numeric(lfc.thresh), "neg",
            ifelse(as.numeric(x$`pos|fdr`) < fdr.thresh & as.numeric(x$LFC) > as.numeric(lfc.thresh), "pos", NA)
        )
    })

    df$pval <- apply(df, 1, function(x) {
        x <- split(unname(x), names(x))
        as.numeric(ifelse(as.numeric(x$`neg|p-value`) < as.numeric(x$`pos|p-value`),
            x$`neg|p-value`, x$`pos|p-value`
        ))
    })

    df$goodsgrna <- apply(df, 1, function(x) {
        x <- split(unname(x), names(x))
        as.integer(ifelse(as.numeric(x$`neg|score`) < as.numeric(x$`pos|score`),
            x$`neg|goodsgrna`, x$`pos|goodsgrna`
        ))
    })

    df$Rank <- rank(df$LFC)
    df$RandomIndex <- sample(seq(nrow(df)), nrow(df))

    return(df)
}


#' Read user-uploaded gene summary files and assign sample names
#' @param fileList A list of gene summary files uploaded by the user.
#'
#' @importFrom utils read.delim
#'
#' @return A named list of data.frames read from gene summary files uploaded by the user.
#' @author Jared Andrews
#' @rdname INTERNAL_gene_summ_ingress
.gene_summ_ingress <- function(fileList) {
    out <- lapply(fileList$datapath, read.delim, check.names = FALSE)
    names(out) <- sapply(fileList$name, gsub,
        pattern = ".gene_summary.txt",
        replacement = "", fixed = TRUE
    )
    return(out)
}


#' Read user-uploaded sgrna summary files and assign sample names
#' @param fileList A list of sgrna summary files uploaded by the user.
#'
#' @return A named list of data.frames read from sgrna summary files uploaded by the user.
#'
#' @importFrom utils read.delim
#'
#' @author Jared Andrews
#' @rdname INTERNAL_sgrna_summ_ingress
.sgrna_summ_ingress <- function(fileList) {
    out <- lapply(fileList$datapath, read.delim, check.names = FALSE)
    names(out) <- sapply(fileList$name, gsub,
        pattern = ".sgrna_summary.txt",
        replacement = "", fixed = TRUE
    )
    return(out)
}
