#' Create interactive bar plot
#'
#' Create an interactive bar plot for specific summary information for each sample in the dataset.
#'
#' @param count.summary data.frame containing the count summary information for each sample.
#' @param x Character scalar for column of count.summary to plot along the x-axis.
#' @param y Character scalar for column of count.summary to plot along the y-axis.
#' @param title Character scalar for title of the plot.
#' @param xlab Character scalar for label of the x-axis.
#' @param ylab Character scalar for label of the y-axis.
#' @param fill Character scalar for bar fill color in hex.
#' @param yaxis.addition Numeric scalar for additional space to add to the y-axis.
#'
#' @return An interactive bar chart showing the specified summary information for each sample.
#'      The axis and plot title are editable.
#'
#' @author Jared Andrews
#' @importFrom plotly ggplotly layout config
#' @importFrom MAGeCKFlute BarView
#' @importFrom ggplot2 theme element_text
#' @export
#'
#' @seealso \code{\link[MAGeCKFlute]{BarView}}, for a static bar plot from the same count summary data.
#' @examples
#' library(CRISPRball)
#' count.summ <- read.delim(system.file("extdata", "escneg.countsummary.txt",
#'     package = "CRISPRball"
#' ), check.names = FALSE)
#' # Gini Index plot
#' plot_bar(count.summ)
#'
#' # Zero count sgRNAs plot
#' plot_bar(count.summ,
#'     x = "Label", y = "Zerocounts", title = "Fully Depleted sgRNAs",
#'     fill = "#394E80", ylab = "Zero Count sgRNAs", yaxis.addition = 10
#' )
plot_bar <- function(count.summary,
                     x = "Label",
                     y = "GiniIndex",
                     title = "sgRNA Read Distribution",
                     xlab = NULL,
                     ylab = "Gini Index",
                     fill = "#E69F00",
                     yaxis.addition = 0.05) {

    gg <- BarView(count.summary,
        x = x,
        y = y,
        ylab = ylab,
        xlab = xlab,
        main = title,
        fill = fill
    )

    gg + theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)
    )

    ggplotly(gg, tooltip = c("y")) %>%
        layout(
            yaxis = list(range = list(0, max(count.summary[[y]]) + yaxis.addition)),
            xaxis = list(tickangle = 315)
        ) %>%
        config(
            toImageButtonOptions = list(format = "svg"),
            displaylogo = FALSE,
            plotGlPixelRatio = 7,
            edits = list(
                axisTitleText = TRUE,
                titleText = TRUE
            )
        )
}


#' Create a plotly plot from a distribution of values
#'
#' This function creates a plotly plot with the distribution of values for each
#' column in the `mat` matrix, using different colors for each column.
#' The legend will show the column names from `mat`, and the plot will have
#' the title "Distribution of read counts".
#'
#' @param mat A matrix with the data to plot.
#' @param xlab Character scalar for label of the x-axis.
#' @param ylab Character scalar for label of the y-axis.
#' @param title A character scalar for the title of the plot.
#' @param show.grid A boolean for whether to show the grid lines.
#' @return A plotly plot with the distribution of read counts.
#'
#' @importFrom plotly plot_ly add_trace layout config ggplotly
#' @importFrom graphics hist
#' @export
#' @author Jared Andrews
#' @examples
#' library(CRISPRball)
#' cts <- read.delim(system.file("extdata", "escneg.count_normalized.txt",
#'     package = "CRISPRball"
#' ), check.names = FALSE)
#' cts.log <- as.matrix(log2(cts[, c(-1, -2)] + 1))
#' colnames(cts.log) <- colnames(cts)[c(-1, -2)]
#'
#' plot_hist(cts.log,
#'     title = "Distribution of read counts",
#'     xlab = "log2(counts + 1)", ylab = "Frequency"
#' )
plot_hist <- function(mat,
                      title = NULL,
                      xlab = "Values",
                      ylab = "Frequency",
                      show.grid = FALSE) {
    
    histo <- hist(mat, breaks = 40)

    if (ncol(mat) >= 1) {
        histlist <- lapply(1:ncol(mat), function(x) {
            return(hist(mat[, x], plot = FALSE, breaks = histo$breaks))
        })
        xrange <- range(unlist(lapply(histlist, function(x) {
            x$mids
        })))
        yrange <- range(unlist(lapply(histlist, function(x) {
            x$counts
        })))
        hst1 <- histlist[[1]]
        fig <- plot_ly(
            x = hst1$mids, y = hst1$counts,
            type = "scatter", mode = "lines+markers",
            name = colnames(mat)[1]
        )
    }

    if (ncol(mat) >= 2) {
        for (i in 2:ncol(mat)) {
            hstn <- histlist[[i]]
            fig <- fig %>% add_trace(
                x = hstn$mids, y = hstn$counts,
                type = "scatter", mode = "lines+markers",
                name = colnames(mat)[i]
            )
        }
    }

    ay <- list(
        showline = TRUE,
        mirror = TRUE,
        linecolor = toRGB("black"),
        linewidth = 0.5,
        showgrid = show.grid,
        layer = "below traces",
        zeroline = FALSE,
        ticks = "outside",
        zerolinewidth = 0.5,
        title = ylab,
        range = c(0, yrange[2] * 1.1)
    )

    ax <- list(
        showline = TRUE,
        mirror = TRUE,
        linecolor = toRGB("black"),
        linewidth = 0.5,
        zeroline = FALSE,
        showgrid = show.grid,
        layer = "below traces",
        ticks = "outside",
        zerolinewidth = 0.5,
        title = xlab,
        range = c(0, xrange[2] * 1.1)
    )

    fig <- fig %>%
        layout(
            title = title,
            xaxis = ax,
            yaxis = ay
        ) %>%
        config(
            toImageButtonOptions = list(format = "svg"),
            displaylogo = FALSE,
            plotGlPixelRatio = 7,
            edits = list(
                axisTitleText = TRUE,
                titleText = TRUE
            )
        )

    fig
}


#' Plot a Correlation Heatmap
#'
#' This function creates a heatmap using ComplexHeatmap to display the correlation values
#' in a matrix. The color of each cell in the heatmap is determined by the
#' corresponding correlation value, using a color ramp that ranges from the
#' minimum value color to a maximum value color.
#'
#' @param mat A matrix containing the correlation values.
#' @param min.color Character scalar for the hexadecimal color code for
#'   the minimum values in the heatmap.
#' @param max.color Character scalar for the hexadecimal color code for
#'   the maximum values in the heatmap.
#' @param plot.title Character scalar for title of the plot.
#' @param legend.title Character scalar for title of the legend.
#'
#' @return A Heatmap object.
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @examples
#' library(CRISPRball)
#' norm.counts <- read.delim(system.file("extdata", "escneg.count_normalized.txt",
#'     package = "CRISPRball"
#' ), check.names = FALSE)
#' norm.counts <- as.matrix(norm.counts[, c(-1, -2)])
#' norm.counts.log <- log2(norm.counts + 1)
#' cor.mat <- cor(norm.counts.log)
#' plot_correlation_heatmap(cor.mat)
#' @author Jared Andrews
#' @export
plot_correlation_heatmap <- function(mat,
                                     min.color = "#FF0000",
                                     max.color = "#0000FF",
                                     legend.title = "Pearson Corr.",
                                     plot.title = "Correlation Matrix") {

    Heatmap(mat,
        name = legend.title,
        column_title = plot.title,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        col = colorRamp2(c(min(mat), max(mat)), c(min.color, max.color))
    )
}


#' Plot a biplot from a PCAtools `pca` object
#'
#' This function plots a biplot from a PCAtools \code{\link[PCAtools]{pca}} object.
#'
#' @param pca.res A \code{\link[PCAtools]{pca}} generated by the PCAtools package.
#' @param dim.x Character scalar for the principal component to plot on the x-axis.
#' @param dim.y Character scalar for the principal component to plot on the y-axis.
#' @param dim.z Character scalar for the principal component to plot on the z-axis.
#' @param plot.title Character scalar for the title of the plot.
#' @param color.by Character scalar for the column name to use for coloring points.
#' @param shape.by Character scalar for the column name to use for shaping points.
#' @param hover.info Character vector of column name(s) to include in the hover info for each point.
#' @param show.loadings Boolean indicating whether to show component loadings on the plot.
#' @param n.loadings Integer for number of loadings to show.
#' @param pt.size Numeric size of the plotted points.
#' @return A plotly plot of the PCA biplot, or a text grob indicating no PCA was provided.
#' @export
#'
#' @author Jared Andrews
#'
#' @importFrom plotly plot_ly layout config add_segments add_annotations ggplotly
#' @importFrom stats as.formula
#' @importFrom dittoSeq dittoColors
#'
#' @seealso \code{\link[PCAtools]{pca}}
#'
#' @examples
#' library("PCAtools")
#' library("CRISPRball")
#' col <- 10
#' row <- 2000
#' mat <- matrix(
#'     rexp(col * row, rate = 0.1),
#'     ncol = col
#' )
#' rownames(mat) <- paste0("gene", 1:nrow(mat))
#' colnames(mat) <- paste0("sample", 1:ncol(mat))
#'
#' metadata <- data.frame(row.names = colnames(mat))
#' metadata$Group <- rep(NA, ncol(mat))
#' metadata$Group[seq(1, 10, 2)] <- "A"
#' metadata$Group[seq(2, 10, 2)] <- "B"
#'
#' p <- pca(mat, metadata = metadata, removeVar = 0.1)
#' plot_pca_biplot(p, color.by = "Group")
plot_pca_biplot <- function(pca.res,
                            dim.x = "PC1",
                            dim.y = "PC2",
                            dim.z = NULL,
                            plot.title = "PCA Biplot",
                            color.by = NULL,
                            shape.by = NULL,
                            hover.info = NULL,
                            show.loadings = FALSE,
                            n.loadings = 3,
                            pt.size = 12) {
    pl.cols <- NULL
    pl.shapes <- NULL
    pl.col <- "black"
    hov.text <- NULL

    # Warn if color.by, shape.by, and hover.info are provided with no metadata
    if (is.null(pca.res$metadata)) {
        if (!is.null(color.by) || !is.null(shape.by) || !is.null(hover.info)) {
            warning("No metadata was provided to the pca function, so color.by, shape.by, and hover.info will be ignored.")
        }
    }

    # Define axes layouts.
    az <- NULL
    zizzy <- NULL
    if (!is.null(dim.z)) {
        az <- list(
            showgrid = FALSE, showline = TRUE, mirror = TRUE, zeroline = FALSE,
            title = paste0(
                dim.z,
                " (", format(round(pca.res$variance[dim.z], 2), nsmall = 2), "%)"
            )
        )

        # Create formula for z axis, can't be done with ifelse statement
        zizzy <- as.formula(paste0("~", dim.z))
    }

    ax <- list(
        showgrid = ifelse(is.null(dim.z), FALSE, TRUE), showline = TRUE, mirror = TRUE, zeroline = FALSE,
        title = paste0(
            dim.x,
            " (", format(round(pca.res$variance[dim.x], 2), nsmall = 2), "%)"
        )
    )

    ay <- list(
        showgrid = ifelse(is.null(dim.z), FALSE, TRUE), showline = TRUE, mirror = TRUE, zeroline = FALSE,
        title = paste0(
            dim.y,
            " (", format(round(pca.res$variance[dim.y], 2), nsmall = 2), "%)"
        )
    )

    # Get marker aesthetics mappings.
    # Drop unused factor levels if possible.
    if (!is.null(color.by) && !is.null(pca.res$metadata)) {
        pl.cols <- pca.res$metadata[, color.by, drop = TRUE]
        if (is.factor(pl.cols)) {
            pl.cols <- droplevels(pl.cols)
        }
        pl.col <- dittoColors()[seq_along(unique(pca.res$metadata[, color.by, drop = TRUE]))]
    }

    if (!is.null(shape.by) & !is.null(pca.res$metadata)) {
        pl.shapes <- pca.res$metadata[, shape.by, drop = TRUE]
        if (is.factor(pl.shapes)) {
            pl.shapes <- droplevels(pl.shapes)
        }
    }

    # Get hover info.
    if (!is.null(hover.info) & !is.null(pca.res$metadata)) {
        for (n in hover.info) {
            hov.text <- paste0(hov.text, "</br><b>", n, ":</b>", pca.res$metadata[[n]])
        }
    }


    fig <- plot_ly(pca.res$rotated,
        x = as.formula(paste0("~", dim.x)),
        y = as.formula(paste0("~", dim.y)),
        z = zizzy,
        type = ifelse(is.null(dim.z), "scatter", "scatter3d"),
        mode = "markers",
        marker = list(size = pt.size),
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
            scene = list(
                xaxis = ax,
                yaxis = ay,
                zaxis = az
            ),
            xaxis = ax,
            yaxis = ay,
            title = plot.title
        )

    fig <- fig %>% toWebGL()

    # Plot loadings.
    # plotly doesn't support easily adding segments in 3D, so skip if 3D is requested.
    if (show.loadings & is.null(dim.z)) {
        lengthLoadingsArrowsFactor <- 1.5

        # Get number of loadings to display.
        xidx <- order(abs(pca.res$loadings[, dim.x]), decreasing = TRUE)
        yidx <- order(abs(pca.res$loadings[, dim.y]), decreasing = TRUE)
        vars <- unique(c(
            rownames(pca.res$loadings)[xidx][seq_len(n.loadings)],
            rownames(pca.res$loadings)[yidx][seq_len(n.loadings)]
        ))

        # get scaling parameter to match between variable loadings and rotated loadings
        # This is identical to PCAtools approach.
        r <- min(
            (max(pca.res$rotated[, dim.x]) - min(pca.res$rotated[, dim.x]) /
                (max(pca.res$loadings[, dim.x]) - min(pca.res$loadings[, dim.x]))),
            (max(pca.res$rotated[, dim.y]) - min(pca.res$rotated[, dim.y]) /
                (max(pca.res$loadings[, dim.y]) - min(pca.res$loadings[, dim.y])))
        )

        fig <- fig %>%
            add_segments(
                x = 0, xend = pca.res$loadings[vars, dim.x] * r * lengthLoadingsArrowsFactor,
                y = 0, yend = pca.res$loadings[vars, dim.y] * r * lengthLoadingsArrowsFactor,
                line = list(color = "black"), inherit = FALSE, showlegend = FALSE, hoverinfo = "text"
            ) %>%
            add_annotations(
                x = pca.res$loadings[vars, dim.x] * r * lengthLoadingsArrowsFactor,
                y = pca.res$loadings[vars, dim.y] * r * lengthLoadingsArrowsFactor,
                ax = 0, ay = 0, text = vars, xanchor = "center", yanchor = "bottom"
            )
    }

    fig <- fig %>%
        config(
            edits = list(
                annotationPosition = TRUE,
                annotationTail = FALSE,
                axisTitleText = TRUE,
                titleText = TRUE
            ),
            toImageButtonOptions = list(format = "svg"),
            displaylogo = FALSE,
            plotGlPixelRatio = 7
        )

    fig
}
