#' Define observers for QC tabs
#'
#' Define a series of observers to run PCA based on new data upload or button click.
#'
#' @param input The Shiny input object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return
#' Observers are set up to monitor the reactives necessary for the PCA.
#' A \code{NULL} is invisibly returned.
#'
#' @author Jared Andrews
#'
#' @rdname INTERNAL_create_qc_observers
#' @importFrom shiny observeEvent
#' @importFrom matrixStats rowVars rowMaxs rowMins
#' @importFrom PCAtools pca
.create_qc_observers <- function(input, robjects) {
    # nocov start
    # This is used so that the matrix and metadata input for PCA are updated when new data is uploaded.
    observeEvent(req(colnames(robjects$norm.counts[, c(-1, -2)]) == gsub("-", ".", robjects$count.summary$Label)), {
        slmed <- robjects$norm.counts
        slmat <- as.matrix(slmed[, c(-1, -2)])
        mat <- log2(slmat + 1)
        rownames(mat) <- slmed$sgRNA

        meta <- robjects$count.summary

        robjects$pca.mat <- mat
        robjects$pca.meta <- meta
    })

    observeEvent(input$pca.update, {
        pca.meta <- robjects$pca.meta
        pca.mat <- robjects$pca.mat

        # Filter samples from QC table.
        if (!is.null(input$count.summary_rows_all) & input$meta.filt) {
            pca.meta <- pca.meta[input$count.summary_rows_all, ]
            pca.mat <- pca.mat[, input$count.summary_rows_all]
        }

        rownames(pca.meta) <- gsub("-", ".", pca.meta$Label)

        # Remove guides with no variance in counts, as they break the PCA.
        pca.mat <- pca.mat[(rowMaxs(pca.mat) - rowMins(pca.mat) > 0), ]

        # If input to use top N features instead rather than percent-based feature removal, account for that
        if (input$keep.top.n) {
            pca.mat <- pca.mat[order(rowVars(pca.mat), decreasing = TRUE), ]
            pca.mat <- pca.mat[1:input$var.n.keep, ]
            var.remove <- 0
        } else {
            var.remove <- input$var.remove
        }

        pca.meta <- pca.meta[colnames(pca.mat), ]

        if (ncol(pca.mat) > 1) {
            robjects$pc <- pca(pca.mat,
                metadata = pca.meta,
                removeVar = var.remove,
                scale = input$scale,
                center = input$center
            )
        }
    })
    # nocov end

    invisible(NULL)
}
