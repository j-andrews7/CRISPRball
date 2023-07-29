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
#' @importFrom shiny observeEvent
#' @importFrom matrixStats rowVars rowMaxs rowMins
#' @importFrom PCAtools pca
#'
#' @author Jared Andrews
#'
#' @rdname INTERNAL_create_qc_observers
.create_qc_observers <- function(input, robjects) {
    # This is used so that the matrix and metadata input for PCA are updated when new data is uploaded.
    # nocov start
    observeEvent(robjects$norm.counts, {
        slmed <- robjects$norm.counts
        slmat <- as.matrix(slmed[, c(-1, -2)])
        mat <- log2(slmat + 1)
        rownames(mat) <- slmed$sgRNA

        if (!is.null(robjects$count.summary)) {
            meta <- robjects$count.summary
        } else {
            meta <- NULL
        }

        robjects$pca.mat <- mat
        robjects$pca.meta <- meta
    })
    # nocov end

    # nocov start
    observeEvent(robjects$count.summary, {
        robjects$pca.meta <- robjects$count.summary
    })
    # nocov end

    # nocov start
    observeEvent(input$pca.update, {
        pca.meta <- robjects$pca.meta
        pca.mat <- robjects$pca.mat

        # Filter samples from summary table.
        if (!is.null(input$count.summary_rows_all) & input$meta.filt & !is.null(pca.meta)) {
            pca.meta <- pca.meta[input$count.summary_rows_all, ]
            pca.mat <- pca.mat[, input$count.summary_rows_all]
        }

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

        # Check that metadata rows and matrix column names match.
        # If not, set the metadata to NULL so PCA still runs.
        # This keeps the app from crashing if the user uploads a metadata file with a different set of samples.
        rownames(pca.meta) <- pca.meta$Label
        if (!is.null(pca.meta) & any(!colnames(pca.mat)[c(-1, -2)] %in% pca.meta$Label)) {
            pca.meta <- NULL
        }

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
