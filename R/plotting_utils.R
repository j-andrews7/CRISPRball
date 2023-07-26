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

#' Plot text on empty plotly plot
#'
#' @param title Character scalar to show in plot area.
#' @return A plotly plot with no data and text in the plot area.
#' 
#' @author Jared Andrews
#' @rdname INTERNAL_empty_plot
#' @importFrom plotly plotly_empty config layout
.empty_plot <- function(title = NULL) {
    p <- plotly_empty(type = "scatter", mode = "markers") %>%
        config(
            displayModeBar = FALSE
        ) %>%
        layout(
            title = list(
                text = title,
                yref = "paper",
                y = 0.5
            )
        )
    return(p)
}
