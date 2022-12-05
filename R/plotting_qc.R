#' Create interactive bar plot
#'
#' Create an interactive bar plot for specific summary information for each sample in the dataset.
#'
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
