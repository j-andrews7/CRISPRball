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
plot_bar <- function(count.summary, x = "Label",
                     y = "GiniIndex", title = "sgRNA Read Distribution", xlab = NULL,
                     ylab = "Gini Index", fill = "#E69F00", yaxis.addition = 0.05) {
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
#' @importFrom plotly plot_ly add_trace layout config
#' @export
#' @author Jared Andrews
#' @examples
#' library(CRISPRball)
#' counts <- read.delim(system.file("extdata", "escneg.count_normalized.txt",
#'     package = "CRISPRball"
#' ), check.names = FALSE)
#' slmed <- robjects$norm.counts
#' tabsmat <- as.matrix(log2(slmed[, c(-1, -2)] + 1))
#' colnames(tabsmat) <- colnames(slmed)[c(-1, -2)]
#'
#' plot_hist(tabsmat,
#'     title = "Distribution of read counts",
#'     xlab = "log2(counts + 1)", ylab = "Frequency"
#' )
plot_hist <- function(mat, title = NULL, xlab = "Values",
                      ylab = "Frequency", show.grid = FALSE) {
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
#' @return A Heatmap object. If only one column is present in the `data` matrix,
#'   returns a text grob indicating only one sample was provided.
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
plot_correlation_heatmap <- function(mat, min.color = "#FF0000",
                                     max.color = "#0000FF", legend.title = "Pearson Corr.",
                                     plot.title = "Correlation Matrix") {
    Heatmap(mat,
        name = legend.title,
        column_title = plot.title,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        col = colorRamp2(c(min(mat), max(mat)), c(min.color, max.color))
    )
}
