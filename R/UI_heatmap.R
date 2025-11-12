#' Create a tabPanel for the gene heatmap tab
#'
#' Create a \code{\link{tabPanel}} with UI elements for the gene heatmap tab.
#'
#' @return
#' A \code{\link{tabPanel}} with UI elements for the gene heatmap tab.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny tabPanel sidebarLayout sidebarPanel mainPanel fluidRow column hr textAreaInput selectInput actionButton div span downloadButton icon
#' @importFrom shinyWidgets prettyCheckbox
#' @importFrom shinyBS popify
#' @importFrom shinycssloaders withSpinner
#' @importFrom DT DTOutput
#' @importFrom InteractiveComplexHeatmap InteractiveComplexHeatmapOutput
#' @rdname INTERNAL_create_tab_heatmap
.create_tab_heatmap <- function() {
    # nocov start
    tabPanel(
        title = "Heatmaps",
        id = "heatmap",
        sidebarLayout(
            sidebarPanel(
                width = 3,
                h4("Heatmap Controls"),
                hr(),
                textAreaInput(
                    "heatmap.genes",
                    label = "Genes",
                    value = "",
                    placeholder = "Enter space, comma, or newline delimited genes",
                    rows = 4
                ),
                selectInput(
                    "heatmap.value.term",
                    label = "Value to display",
                    choices = character(0)
                ),
                prettyCheckbox(
                    "heatmap.scale",
                    label = "Scale genes (rows)",
                    value = TRUE,
                    animation = "smooth",
                    status = "success",
                    bigger = TRUE,
                    icon = icon("check")
                ),
                prettyCheckbox(
                    "heatmap.cluster.columns",
                    label = "Cluster cell lines",
                    value = TRUE,
                    animation = "smooth",
                    status = "success",
                    bigger = TRUE,
                    icon = icon("check")
                ),
                div(actionButton("heatmap.update", "Update Heatmap"), align = "center")
            ),
            mainPanel(
                width = 9,
                fluidRow(
                    column(
                        width = 12,
                        span(
                            popify(
                                icon("circle-info", style = "font-size: 20px"),
                                title = "Gene Heatmap",
                                content = c(
                                    "Visualize uploaded gene summary data for a custom set of genes across the available datasets.",
                                    "Enter one or more genes, choose the value to display, and optionally scale rows or turn off column clustering.",
                                    "Click a heatmap cell to view the underlying value for the selected dataset."
                                ),
                                placement = "bottom",
                                trigger = "hover",
                                options = list(container = "body")
                            ),
                            div(downloadButton("dl_plot.heatmap", "Download Heatmap", class = "btn-dl"), style = "display:inline-block; float:right"),
                            withSpinner(InteractiveComplexHeatmapOutput(
                                heatmap_id = "gene_heatmap",
                                title1 = "Gene Heatmap",
                                layout = "1-2",
                                width1 = 950,
                                height1 = 550,
                                response = "click",
                                output_ui = div(DTOutput("heatmap.info"), style = "font-size:70%;")
                            ))
                        )
                    )
                )
            )
        )
    )
    # nocov end
}
