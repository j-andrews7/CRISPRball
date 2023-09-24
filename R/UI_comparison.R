#' Create a tabPanel for the comparison tab
#'
#' Create a \code{\link{tabPanel}} with UI elements for the comparison tab.
#'
#' @param datasets A character vector containing dataset names.
#'
#' @return
#' A \code{\link{tabPanel}} with UI elements for the comparison tab.
#'
#' @author Jared Andrews
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom shinyBS tipify
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyjqui jqui_resizable
#' @importFrom InteractiveComplexHeatmap InteractiveComplexHeatmapOutput
#'
#' @rdname INTERNAL_create_tab_comparison
.create_tab_comparison <- function(datasets) {
    # nocov start
    tabPanel(
        title = "Comparisons",
        id = "comparisons",
        sidebarLayout(
            sidebarPanel(
                width = 2,
                h4("Plot Controls"),
                hr(),
                div(
                    fluidRow(
                        column(
                            12,
                            tipify(
                                div(
                                    id = "comp_select",
                                    selectizeInput("comp.sets", "Datasets:",
                                        # Upset plots only allow 30 sets max.
                                        if (length(datasets) > 30) {
                                            selected <- datasets[seq_len(30)]
                                        } else {
                                            selected <- datasets
                                        },
                                        choices = datasets,
                                        multiple = TRUE,
                                        options = list(maxItems = 30)
                                    )
                                ),
                                "Datasets to compare.", "right",
                                options = list(container = "body")
                            )
                        ),
                        uiOutput("comp.term.options"),
                        column(
                            6,
                            numericInput("comp.sig.th", "Significance threshold:",
                                min = 0, max = 1, step = 0.01, value = 0.05
                            )
                        ),
                        column(
                            6,
                            numericInput("comp.es.th", "Effect size threshold:",
                                min = 0, max = Inf, step = 0.05, value = 0.5
                            )
                        )
                    ),
                    fluidRow(
                        column(
                            12,
                            tipify(
                                prettyCheckbox("comp.rem.ess",
                                    label = "Remove essential genes", value = FALSE,
                                    animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")
                                ),
                                "Remove essential genes if any provided to function.", "right",
                                options = list(container = "body")
                            ),
                            tipify(
                                prettyCheckbox("comp.dep.crispr.ess",
                                    label = "Remove DepMap CRISPR essential genes", value = FALSE,
                                    animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")
                                ),
                                "Remove DepMap Chronos Combined, Score, and Achilles common essential genes from latest release.",
                                "right",
                                options = list(container = "body")
                            ),
                            tipify(
                                prettyCheckbox("comp.dep.rnai.ess",
                                    label = "Remove DepMap RNAi essential genes", value = FALSE,
                                    animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")
                                ),
                                "Remove RNAi common essential genes from latest DepMap release.", "right",
                                options = list(container = "body")
                            ),
                            tipify(
                                prettyCheckbox("comp.dep.crispr.sel",
                                    label = "Remove DepMap CRISPR selective genes", value = FALSE,
                                    animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")
                                ),
                                "Remove DepMap Chronos Combined, Score, and Achilles strongly selective genes from latest release.", "right",
                                options = list(container = "body")
                            ),
                            tipify(
                                prettyCheckbox("comp.dep.rnai.sel",
                                    label = "Remove DepMap RNAi selective genes", value = FALSE,
                                    animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")
                                ),
                                "Remove DepMap RNAi strongly selective genes from latest release.", "right",
                                options = list(container = "body")
                            ),
                            tipify(
                                prettyCheckbox("comp.rem.pos",
                                    label = "Remove positive control genes", value = FALSE,
                                    animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")
                                ),
                                "Remove positive control genes if provided.", "right",
                                options = list(container = "body")
                            )
                        )
                    ),
                    div(actionButton("comp.update", "Update Upset Plot"), align = "center"),
                    style = "background-color: #FFFFFF; padding: 3px; margin-bottom: 3px; border: 1px solid #bce8f1; "
                ),
            ),
            mainPanel(
                width = 10,
                fluidRow(
                    column(
                        width = 12,
                        span(
                            popify(icon("circle-info", style = "font-size: 20px"),
                                title = "Upset Plot",
                                c(
                                    "This upset plot shows shared and distinct positive hits between the chosen datasets. ",
                                    "Click a set for additional information in the output table."
                                ),
                                placement = "bottom", trigger = "hover", options = list(container = "body")
                            ),
                            withSpinner(InteractiveComplexHeatmapOutput(
                                heatmap_id = "overlap_pos",
                                title1 = "Positively Selected Hits", layout = "1-2-3",
                                title3 = "Shared Positively Selected Hits for Clicked Set",
                                width1 = 950, width3 = 600, height1 = 300,
                                response = "click", # Removes the sub-heatmap
                                output_ui = div(DT::dataTableOutput("comp.pos.info"), style = "font-size:70%;")
                            ))
                        ),
                        hr(),
                        span(
                            popify(icon("circle-info", style = "font-size: 20px"),
                                title = "Upset Plot",
                                c(
                                    "This upset plot shows shared and distinct negative hits between the chosen datasets. ",
                                    "Click a set for additional information in the output table."
                                ),
                                placement = "bottom", trigger = "hover", options = list(container = "body")
                            ),
                            withSpinner(InteractiveComplexHeatmapOutput(
                                heatmap_id = "overlap_neg",
                                title1 = "Negatively Selected Hits", layout = "1-2-3",
                                title3 = "Shared Negatively Selected Hits for Clicked Set",
                                width1 = 950, width3 = 600, height1 = 300,
                                response = "click", # Removes the sub-heatmap
                                output_ui = div(DT::dataTableOutput("comp.neg.info"), style = "font-size:70%;")
                            ))
                        )
                    )
                )
            )
        )
    )
    # nocov end
}
