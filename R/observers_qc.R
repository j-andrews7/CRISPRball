#' Define observers for PCA
#'
#' Define a series of observers to run PCA based on new data upload or button click.
#'
#' @param input The Shiny input object from the server function.
#' @param output The Shiny output object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return
#' Observers are set up to monitor the reactives necessary for the PCA.
#' A \code{NULL} is invisibly returned.
#'
#' @author Jared Andrews
#'
#' @rdname INTERNAL_create_pca_observers
#' @importFrom shiny observeEvent updateSelectizeInput
#' @importFrom shinyWidgets updatePickerInput
#' @importFrom shinyjs js
.create_pca_observers <- function(input, output, robjects) {
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

    observeEvent(robjects$pc, {
        output$qc.pca <- renderPlotly({
            req(input$dim1, input$dim2, input$dim3)

            pc.res <- robjects$pc

            pl.cols <- NULL
            pl.shapes <- NULL
            pl.col <- "black"
            hov.text <- NULL

            # Get marker aesthetics mappings.
            # Drop unused factor levels if possible.
            if (isolate(input$bip.color) != "") {
                pl.cols <- pc.res$metadata[, isolate(input$bip.color), drop = TRUE]
                if (is.factor(pl.cols)) {
                    pl.cols <- droplevels(pl.cols)
                }
                pl.col <- dittoColors()[seq_along(unique(pc.res$metadata[, isolate(input$bip.color), drop = TRUE]))]
            }

            if (isolate(input$bip.shape) != "") {
                pl.shapes <- pc.res$metadata[, isolate(input$bip.shape), drop = TRUE]
                if (is.factor(pl.shapes)) {
                    pl.shapes <- droplevels(pl.shapes)
                }
            }

            # Just throw label on hover for now.
            hov.text <- paste0("</br><b>Label:</b> ", pc.res$metadata$Label)

            # Check if 2D is wanted.
            if (isolate(input$bip.twod)) {
                fig <- plot_ly(pc.res$rotated,
                    x = as.formula(paste0("~", isolate(input$dim1))),
                    y = as.formula(paste0("~", isolate(input$dim2))),
                    type = "scatter",
                    mode = "markers",
                    marker = list(size = 15),
                    color = pl.cols,
                    colors = pl.col,
                    symbol = pl.shapes,
                    symbols = c(
                        "circle", "square", "diamond", "cross",
                        "diamond-open", "circle-open", "square-open", "x"
                    ),
                    text = hov.text,
                    hoverinfo = "text"
                ) %>%
                    layout(
                        xaxis = list(
                            showgrid = FALSE, showline = TRUE, mirror = TRUE, zeroline = FALSE,
                            title = paste0(
                                isolate(input$dim1),
                                " (", format(round(pc.res$variance[isolate(input$dim1)], 2), nsmall = 2), "%)"
                            )
                        ),
                        yaxis = list(
                            showgrid = FALSE, showline = TRUE, mirror = TRUE, zeroline = FALSE,
                            title = paste0(
                                isolate(input$dim2),
                                " (", format(round(pc.res$variance[isolate(input$dim2)], 2), nsmall = 2), "%)"
                            )
                        )
                    )

                fig <- fig %>% toWebGL()

                # Plot loadings.
                if (isolate(input$bip.loadings)) {
                    lengthLoadingsArrowsFactor <- 1.5

                    # Get number of loadings to display.
                    xidx <- order(abs(pc.res$loadings[, isolate(input$dim1)]), decreasing = TRUE)
                    yidx <- order(abs(pc.res$loadings[, isolate(input$dim2)]), decreasing = TRUE)
                    vars <- unique(c(
                        rownames(pc.res$loadings)[xidx][seq_len(isolate(input$bip.n.loadings))],
                        rownames(pc.res$loadings)[yidx][seq_len(isolate(input$bip.n.loadings))]
                    ))

                    # get scaling parameter to match between variable loadings and rotated loadings
                    # This is cribbed almost verbatim from PCAtools code.
                    r <- min(
                        (max(pc.res$rotated[, isolate(input$dim1)]) - min(pc.res$rotated[, isolate(input$dim1)]) /
                            (max(pc.res$loadings[, isolate(input$dim1)]) - min(pc.res$loadings[, isolate(input$dim1)]))),
                        (max(pc.res$rotated[, isolate(input$dim2)]) - min(pc.res$rotated[, isolate(input$dim2)]) /
                            (max(pc.res$loadings[, isolate(input$dim2)]) - min(pc.res$loadings[, isolate(input$dim2)])))
                    )

                    fig <- fig %>%
                        add_segments(
                            x = 0, xend = pc.res$loadings[vars, isolate(input$dim1)] * r * lengthLoadingsArrowsFactor,
                            y = 0, yend = pc.res$loadings[vars, isolate(input$dim2)] * r * lengthLoadingsArrowsFactor,
                            line = list(color = "black"), inherit = FALSE, showlegend = FALSE, hoverinfo = "text"
                        ) %>%
                        add_annotations(
                            x = pc.res$loadings[vars, isolate(input$dim1)] * r * lengthLoadingsArrowsFactor,
                            y = pc.res$loadings[vars, isolate(input$dim2)] * r * lengthLoadingsArrowsFactor,
                            ax = 0, ay = 0, text = vars, xanchor = "center", yanchor = "bottom"
                        )
                }
            } else {
                # Generate plot.
                fig <- plot_ly(pc.res$rotated,
                    x = as.formula(paste0("~", isolate(input$dim1))),
                    y = as.formula(paste0("~", isolate(input$dim2))),
                    z = as.formula(paste0("~", isolate(input$dim3))),
                    type = "scatter3d",
                    mode = "markers",
                    color = pl.cols,
                    colors = pl.col,
                    symbol = pl.shapes,
                    symbols = c(
                        "circle", "square", "diamond", "cross", "diamond-open",
                        "circle-open", "square-open", "x"
                    ),
                    text = hov.text,
                    hoverinfo = "text"
                ) %>%
                    layout(scene = list(
                        xaxis = list(title = paste0(
                            isolate(input$dim1), " (",
                            format(round(pc.res$variance[isolate(input$dim1)], 2), nsmall = 2), "%)"
                        )),
                        yaxis = list(title = paste0(
                            isolate(input$dim2), " (",
                            format(round(pc.res$variance[isolate(input$dim2)], 2), nsmall = 2), "%)"
                        )),
                        zaxis = list(title = paste0(
                            isolate(input$dim3), " (",
                            format(round(pc.res$variance[isolate(input$dim3)], 2), nsmall = 2), "%)"
                        )),
                        camera = list(eye = list(x = 1.5, y = 1.8, z = 0.4))
                    ))
            }
            fig <- fig %>%
                config(
                    edits = list(
                        annotationPosition = TRUE,
                        annotationTail = FALSE
                    ),
                    toImageButtonOptions = list(format = "svg"),
                    displaylogo = FALSE,
                    plotGlPixelRatio = 7
                )

            fig
        })
    })
    
    # nocov end

    invisible(NULL)
}
