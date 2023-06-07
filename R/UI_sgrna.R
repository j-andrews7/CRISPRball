#' Create a tabPanel for the sgrna tab
#'
#' Create a \code{\link{tabPanel}} with UI elements for the sgrna tab.
#'
#' @param sgrna.choices A character vector containing dataset names.
#' @param sgrna.gene A character vector containing gene identifiers.
#'
#' @return
#' A \code{\link{tabPanel}} with UI elements for the sgrna tab.
#'
#' @author Jared Andrews
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom colourpicker colourInput
#' @importFrom shinyBS tipify popify bsCollapse bsCollapsePanel
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyjqui jqui_resizable
#' @importFrom plotly plotlyOutput
#'
#' @rdname INTERNAL_create_tab_sgrna
.create_tab_sgrna <- function(sgrna.choices, sgrna.gene) {
    # nocov start
    tabPanel(
        title = "sgRNA",
        id = "sgrna",
        sidebarLayout(
            sidebarPanel(
                width = 2,
                h4("Plot Controls"),
                hr(),
                div(
                    fluidRow(
                        column(
                            6,
                            tipify(selectizeInput("sgrna.sel1", "Dataset 1:", choices = sgrna.choices),
                                "Dataset shown in top row.", "right",
                                options = list(container = "body")
                            )
                        ),
                        column(
                            6,
                            tipify(selectizeInput("sgrna.sel2", "Dataset 2:", choices = sgrna.choices),
                                "Dataset shown in bottom row.", "right",
                                options = list(container = "body")
                            )
                        )
                    ),
                    fluidRow(
                        column(
                            12,
                            pickerInput("sgrna.gene", "Choose gene:",
                                choices = sgrna.gene,
                                multiple = FALSE, options = list(`live-search` = TRUE, `actions-box` = TRUE)
                            ),
                            prettyCheckbox("sgrna.rank.ascending",
                                label = "Ascending Rank", value = TRUE,
                                animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")
                            ),
                            uiOutput("sgrna.rank.options"),
                        )
                    ),
                    style = "background-color: #FFFFFF; padding: 3px; margin-bottom: 3px; border: 1px solid #bce8f1; "
                ),
            ),
            mainPanel(
                width = 10,
                fluidRow(
                    column(
                        width = 2,
                        span(
                            popify(icon("circle-info", style = "font-size: 20px"),
                                title = "Counts Plot",
                                c(
                                    "This rank plot shows the normalized counts for each sgRNA for the ",
                                    "specified gene across the samples that make up the dataset comparision. ",
                                    "Hover over a point for additional info."
                                ),
                                placement = "bottom", trigger = "hover", options = list(container = "body")
                            ),
                            div(downloadButton("dl_plot.sgrna1.counts", "Download Counts Plot", class = "btn-dl"), style = "display:inline-block; float:right"),
                            withSpinner(jqui_resizable(plotlyOutput("sgrna1.counts")))
                        )
                    ),
                    column(
                        width = 4,
                        span(
                            popify(icon("circle-info", style = "font-size: 20px"),
                                title = "Rank Plot",
                                c(
                                    "This rank plot shows the 'rank by' term on the y-axis and the sgRNA rank on the x-axis. ",
                                    "sgRNAs for the selected gene will be highlighted. ",
                                    "Click and drag to zoom in. Hover over a point for additional info."
                                ),
                                placement = "bottom", trigger = "hover", options = list(container = "body")
                            ),
                            div(downloadButton("dl_plot.sgrna1.rank", "Download Rank Plot", class = "btn-dl"), style = "display:inline-block; float:right"),
                            withSpinner(jqui_resizable(plotlyOutput("sgrna1.rank")))
                        )
                    ),
                    column(
                        width = 6,
                        jqui_resizable(div(DT::dataTableOutput("sgrna1.detail"), style = "font-size:80%;"))
                    )
                ),
                hr(),
                fluidRow(
                    column(
                        width = 2,
                        withSpinner(jqui_resizable(plotlyOutput("sgrna2.counts"))),
                        div(downloadButton("dl_plot.sgrna2.counts", "Download Counts Plot", class = "btn-dl"), style = "display:inline-block; float:right; height:18px"),
                    ),
                    column(
                        width = 4,
                        withSpinner(jqui_resizable(plotlyOutput("sgrna2.rank"))),
                        div(downloadButton("dl_plot.sgrna2.rank", "Download Rank Plot", class = "btn-dl"), style = "display:inline-block; float:right"),
                    ),
                    column(
                        width = 6,
                        jqui_resizable(div(DT::dataTableOutput("sgrna2.detail"), style = "font-size:80%;"))
                    )
                )
            )
        )
    )
    # nocov end
}


#' Create a tabPanel for the sgrna summary tab
#'
#' Create a \code{\link{tabPanel}} with UI elements for the sgrna summary tab.
#'
#' @return
#' A \code{\link{tabPanel}} with UI elements for the sgrna summary tab.
#'
#' @author Jared Andrews
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom shinycssloaders withSpinner
#'
#' @rdname INTERNAL_create_tab_sgrna_summary
.create_tab_sgrna_summary <- function() {
    # nocov start
    tabPanel(
        title = "sgRNA Summary Tables",
        id = "sgrna-tables",
        br(),
        div(withSpinner(DT::dataTableOutput("sgrna1.summary")), style = "font-size:80%;"),
        br(),
        div(withSpinner(DT::dataTableOutput("sgrna2.summary")), style = "font-size:80%;")
    )
    # nocov end
}
