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
#' @param es.thresh Numeric scalar for absolute log fold change threshold to consider a gene a hit.
#' @param es.col Character scalar for the column name of the effect size value.
#' @param sig.col Character scalar for the column name of the significance value.
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
#' out.df <- gene_ingress(d1.genes, 0.05, 0.5, es.col = "LFC", sig.col = "fdr")
gene_ingress <- function(df, sig.thresh, es.thresh, es.col, sig.col, positive.ctrl.genes = NULL,
                         essential.genes = NULL, depmap.genes = NULL) {

    if ("neg|score" %in% colnames(df)) {
        # These handle initial load where UI has not populated the sig.col and es.col values yet.
        if (is.null(es.col)) {
            es.col <- "LFC"
        }

        if (is.null(sig.col)) {
            sig.col <- "fdr"
        }
        df <- .rra_ingress(df, sig.thresh, es.thresh, es.col = es.col, sig.col = sig.col)
    } else if (any(grepl("beta", colnames(df)))) {
        # These handle initial load where UI has not populated the sig.col and es.col values yet.
        if (is.null(es.col)) {
            es.col <- "beta"
        }

        if (is.null(sig.col)) {
            sig.col <- "fdr"
        }
        df <- .mle_ingress(df, sig.thresh, es.thresh, es.col = es.col, sig.col = sig.col)
    } else {
        stop("Unknown gene summary format.")
    }

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
#' mle_gene_summary <- file.path(system.file("extdata", "beta_leukemia.gene_summary.txt",
#'     package = "CRISPRball"
#' ))
#' gene_data <- read_mle_gene_summary(mle_gene_summary)
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
#' @param sig.thresh Numeric scalar for significance threshold to consider a gene a hit.
#' @param es.thresh Numeric scalar for absolute effect size threshold to consider a gene a hit.
#' @param sig.col Character string for the column name of the significance value.
#'   Should be one of "fdr" or "pval".
#' @param es.col Character string for the column name of the effect size value.
#'
#' @return data.frame of gene summary data with renamed columns.
#' @author Jared Andrews
#' @rdname INTERNAL_rra_ingress
.rra_ingress <- function(df, sig.thresh, es.thresh, sig.col = "fdr", es.col = "LFC") {
    df$LFC <- as.numeric(df$`neg|lfc`)

    df$RRAscore <- apply(df, 1, function(x) {
        x <- split(unname(x), names(x))
        as.numeric(ifelse(as.numeric(x$`neg|score`) < as.numeric(x$`pos|score`), x$`neg|score`, x$`pos|score`))
    })

    # Ease of use for plotting RRAscore as rank and interpreting directionality
    df$signed_log10_RRAscore <- apply(df, 1, function(x) {
        x <- split(unname(x), names(x))
        as.numeric(ifelse(as.numeric(x$`neg|score`) < as.numeric(x$`pos|score`),
            -(-log10(as.numeric(x$`neg|score`))), -log10(as.numeric(x$`pos|score`))
        ))
    })

    df$fdr <- apply(df, 1, function(x) {
        x <- split(unname(x), names(x))
        as.numeric(ifelse(as.numeric(x$`neg|fdr`) < as.numeric(x$`pos|fdr`), x$`neg|fdr`, x$`pos|fdr`))
    })

    df$pval <- apply(df, 1, function(x) {
        x <- split(unname(x), names(x))
        as.numeric(ifelse(as.numeric(x$`neg|p-value`) < as.numeric(x$`pos|p-value`),
            x$`neg|p-value`, x$`pos|p-value`
        ))
    })

    if (!sig.col %in% colnames(df)) {
        message(paste0("Column '", sig.col, "' not found in data. Using 'fdr' instead."))
        sig.col <- "fdr"
    }

    if (!es.col %in% colnames(df)) {
        message(paste0("Column '", es.col, "' not found in data. Using 'LFC' instead."))
        es.col <- "LFC"
    }

    df$hit_type <- apply(df, 1, function(x) {
        x <- split(unname(x), names(x))
        ifelse(as.numeric(x[[sig.col]]) < sig.thresh & as.numeric(x[[es.col]]) < -as.numeric(es.thresh), "neg",
            ifelse(as.numeric(x[[sig.col]]) < sig.thresh & as.numeric(x[[es.col]]) > as.numeric(es.thresh), "pos", NA)
        )
    })

    df$goodsgrna <- apply(df, 1, function(x) {
        x <- split(unname(x), names(x))
        as.integer(ifelse(as.numeric(x$`neg|score`) < as.numeric(x$`pos|score`),
            x$`neg|goodsgrna`, x$`pos|goodsgrna`
        ))
    })

    df$Rank <- rank(df[[es.col]])
    df$RandomIndex <- sample(seq(nrow(df)), nrow(df))

    return(df)
}


#' Rename columns and add additional columns to gene summary data from MAGeCK mle
#' @param df data.frame of gene summary data.
#' @param sig.thresh Numeric scalar for significance threshold to consider a gene a hit.
#' @param es.thresh Numeric scalar for effect size threshold to consider a gene a hit.
#' @param sig.col Column name for significance values in \code{df}. Should be one of
#'   "fdr", "p.value", "wald.p.value", or "wald.fdr".
#' @param es.col Column name for effect size values in \code{df}. Should be one of
#'   "beta" or "z".
#'
#' @return data.frame of gene summary data from MAGeCK mle with renamed columns.
#' @author Jared Andrews
#' @rdname INTERNAL_mle_ingress
.mle_ingress <- function(df, sig.thresh, es.thresh, sig.col = "fdr", es.col = "beta") {
    if (!sig.col %in% colnames(df)) {
        message(paste0("Column '", sig.col, "' not found in data. Using 'fdr' instead."))
        sig.col <- "fdr"
    }

    if (!es.col %in% colnames(df)) {
        message(paste0("Column '", es.col, "' not found in data. Using 'beta' instead."))
        es.col <- "beta"
    }

    # To match RRA output.
    df$id <- df$Gene

    df$hit_type <- apply(df, 1, function(x) {
        x <- split(unname(x), names(x))
        ifelse(as.numeric(x[[sig.col]]) < sig.thresh & as.numeric(x[[es.col]]) < -as.numeric(es.thresh), "neg",
            ifelse(as.numeric(x[[sig.col]]) < sig.thresh & as.numeric(x[[es.col]]) > as.numeric(es.thresh), "pos", NA)
        )
    })

    df$Rank <- rank(df[[es.col]])
    df$RandomIndex <- sample(seq(nrow(df)), nrow(df))

    return(df)
}


#' Read user-uploaded gene summary files and assign sample names
#' @param fileList A list of gene summary files uploaded by the user.
#'   MAGeCK mle format will be autodetected.
#'
#' @importFrom utils read.delim
#'
#' @return A named list of data.frames read from gene summary files uploaded by the user.
#' @author Jared Andrews
#' @rdname INTERNAL_gene_summ_ingress
.gene_summ_ingress <- function(fileList) {
    checker <- read.delim(fileList$datapath[[1]], check.names = FALSE)

    if ("neg|score" %in% colnames(checker)) {
        out <- lapply(fileList$datapath, read.delim, check.names = FALSE)
        names(out) <- sapply(fileList$name, gsub,
            pattern = ".gene_summary.txt",
            replacement = "", fixed = TRUE
        )
    } else if (any(grepl("beta", colnames(checker)))) {
        out <- read_mle_gene_summary(fileList$datapath[[1]])
    } else {
        stop("Unknown gene summary format.")
    }

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
