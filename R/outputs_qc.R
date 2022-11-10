#' Render general output for the QC tab
#'
#' Create rendering expressions for the QC tab outputs.
#'
#' @param input The Shiny input object from the server function.
#' @param output The Shiny output object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return A \linkS4class{NULL} is invisibly returned
#' and rendering expressions for QC tab features are added to \code{output}.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny renderUI renderPlot tagList column selectInput
#' @importFrom plotly renderPlotly ggplotly layout config plot_ly toWebGL add_segments add_annotations
#' @importFrom MAGeCKFlute BarView MapRatesView
#' @importFrom shinyWidgets updatePickerInput
#' @importFrom shinyjs js
#' @importFrom dittoSeq dittoColors
#' @importFrom grid grid.newpage grid.text
#' @importFrom DT renderDT datatable formatStyle
#' @importFrom graphics hist legend lines
#' @importFrom stats cor as.formula
#' @import ggplot2
#' @rdname INTERNAL_create_qc_output
.create_qc_output <- function(input, output, robjects) {
    # nocov start
    output$pca.comps <- renderUI({
      req(robjects$pc)
      pcs <- robjects$pc

      tagList(
        fluidRow(
          column(4, selectInput("dim1", "Dim1:", choices = pcs$components, selected = "PC1")),
          column(4, selectInput("dim2", "Dim2:", choices = pcs$components, selected = "PC2")),
          if (length(pcs$components) > 2) {
            column(4, selectInput("dim3", "Dim3:", choices = pcs$components, selected = "PC3"))
          }
        )
      )
    })
    # nocov end

    # nocov start
    output$qc.pca <- renderPlotly({
        req(input$dim1, input$dim2)

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
    # nocov end

    # nocov start
    output$qc.gini <- renderPlotly({
      gg <- BarView(robjects$count.summary,
        x = "Label",
        y = "GiniIndex",
        ylab = "Gini index",
        main = "sgRNA Read Distribution"
      )

      gg + theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)
      )

      ggplotly(gg, tooltip = c("y")) %>%
        layout(
          yaxis = list(range = list(0, max(robjects$count.summary$GiniIndex) + .05)),
          xaxis = list(tickangle = 315)
        ) %>%
        config(
          toImageButtonOptions = list(format = "svg"),
          displaylogo = FALSE,
          plotGlPixelRatio = 7
        )
    })
    # nocov end

    # nocov start
    output$qc.missed <- renderPlotly({
      gg <- BarView(robjects$count.summary,
        x = "Label", y = "Zerocounts", fill = "#394E80",
        ylab = "Zero Count sgRNAs", main = "Fully Depleted sgRNAs"
      )

      gg + theme_classic() + theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)
      ) +
        ylim(0, max(robjects$count.summary$Zerocounts) + 5)

      ggplotly(gg, tooltip = c("y")) %>%
        layout(
          yaxis = list(range = list(0, max(robjects$count.summary$Zerocounts) + 5)),
          xaxis = list(tickangle = 315)
        ) %>%
        config(
          toImageButtonOptions = list(format = "svg"),
          displaylogo = FALSE,
          plotGlPixelRatio = 7
        )
    })
    # nocov end

    # nocov start
    output$qc.map <- renderPlot({
      MapRatesView(robjects$count.summary) + scale_x_discrete(guide = guide_axis(angle = 45))
    })
    # nocov end

    # TODO: rewrite this.
    # nocov start
    output$qc.histplot <- renderPlot({
      colors <- dittoColors()

      slmed <- robjects$norm.counts
      tabsmat <- as.matrix(log2(slmed[, c(-1, -2)] + 1))
      colnames(tabsmat) <- colnames(slmed)[c(-1, -2)]
      samplecol <- colors[((1:ncol(tabsmat)) %% length(colors))]
      tgz <- hist(tabsmat, breaks = 40)

      if (ncol(tabsmat) >= 1) {
        histlist <- lapply(1:ncol(tabsmat), function(X) {
          return(hist(tabsmat[, X], plot = FALSE, breaks = tgz$breaks))
        })
        xrange <- range(unlist(lapply(histlist, function(X) {
          X$mids
        })))
        yrange <- range(unlist(lapply(histlist, function(X) {
          X$counts
        })))
        hst1 <- histlist[[1]]
        plot(hst1$mids, hst1$counts,
          type = "b", pch = 20, xlim = c(0, xrange[2] * 1.2),
          ylim = c(0, yrange[2] * 1.2), xlab = "log2(counts)", ylab = "Frequency",
          main = "Distribution of read counts", col = samplecol[1]
        )
      }

      if (ncol(tabsmat) >= 2) {
        for (i in 2:ncol(tabsmat)) {
          hstn <- histlist[[i]]
          lines(hstn$mids, hstn$counts, type = "b", pch = 20, col = samplecol[i])
        }
      }

      legend("topright", colnames(tabsmat), pch = 20, lwd = 1, col = samplecol)
    })
    # nocov end

    # TODO: rewrite this, add color min/max/mid selectors.
    # nocov start
    output$qc.corr <- renderPlot({
      slmed <- robjects$norm.counts
      slmat <- as.matrix(slmed[, c(-1, -2)])
      slmat.log <- log2(slmat + 1)

      if (ncol(slmat.log) > 1) {
        ComplexHeatmap::pheatmap(cor(slmat.log),
          heatmap_legend_param = list(title = "Pearson\nCorr."),
          main = "Correlation Matrix"
        )
      } else {
        grid.newpage()
        grid.text("Only one sample, no correlation possible.")
      }
    })
    # nocov end

    # nocov start
    output$count.summary <- renderDT(server = FALSE, {
      DT::datatable(robjects$count.summary,
        rownames = FALSE,
        filter = "top",
        extensions = c("Buttons", "Scroller"),
        options = list(
          search = list(regex = TRUE),
          lengthMenu = list(c(10, 25, 50, -1), c("10", "25", "50", "all")),
          dom = "Blfrtip",
          buttons = c("copy", "csv", "excel", "pdf", "print"),
          scrollX = TRUE,
          deferRender = TRUE,
          scrollY = 600,
          scroller = TRUE
        )
      ) %>% DT::formatStyle(0, target = "row", lineHeight = "80%")
    })
    # nocov end
    invisible(NULL)
}