#' Make empty heatmap with text
#'
#' @param text Character scalar to show in plot area.
#' @return A heatmap with no data and text in the plot area.
#'
#' @author Jared Andrews
#' @rdname INTERNAL_empty_heatmap
#' @importFrom ComplexHeatmap Heatmap draw
.empty_heatmap <- function(text) {
    data <- matrix(0, nrow = 1, ncol = 1)

    # Create labels for the axes
    colnames(data) <- c("a")
    rownames(data) <- c("1")

    # Create the heatmap
    heatmap <- Heatmap(data,
        col = "white",
        name = "empty", show_row_names = FALSE,
        show_column_names = FALSE, show_heatmap_legend = FALSE,
        column_title = text
    )

    # Draw the heatmap
    ht <- draw(heatmap)
    return(ht)
}

#' Create an empty ggplot2 plot or plotly plot with input text
#'
#' This function creates an empty ggplot2 or plotly plot and places a user-provided text
#' string in the middle of the plot.
#'
#' @param text Character scalar to show in plot area.
#' @param plotly Boolean indicating whether to return a plotly object.
#' @return Either a ggplot object or a plotly object if \code{plotly = TRUE}.
#'
#' @author Jared Andrews
#'
#' @rdname INTERNAL_empty_plot
#' @seealso \code{\link[ggplot2]{geom_text}}, \code{\link[ggplot2]{theme_void}}
#' @importFrom ggplot2 theme_void geom_text theme margin ggplot
#' @importFrom plotly ggplotly layout
.empty_plot <- function(text = NULL, plotly = FALSE) {
    plot <- ggplot() +
        theme_void() +
        theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
        geom_text(aes(x = 0.5, y = 0.5, label = text),
            inherit.aes = FALSE, check_overlap = TRUE
        )

    if (plotly) {
        plot <- ggplotly(plot)
        plot <- plot %>% layout(
            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, showline = FALSE),
            yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, showline = FALSE),
            plot_bgcolor = "white",
            showlegend = FALSE,
            autosize = TRUE,
            margin = list(l = 0, r = 0, b = 0, t = 0)
        )
    }

    return(plot)
}
