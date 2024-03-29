#' Create a tabPanel for the QC tab
#'
#' Create a \code{\link{tabPanel}} with UI elements for the QC tab.
#'
#' @param meta.choices A character vector containing the metadata columns that 
#'   can be used to color and/or shape points in the PCA biplot.
#'
#' @return
#' A \code{\link{tabPanel}} with UI elements for the QC tab.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny sidebarLayout sidebarPanel mainPanel fluidRow column hr br div numericInput 
#'   selectizeInput h4 h3 uiOutput downloadButton tabPanel span actionButton icon conditionalPanel plotOutput
#' @importFrom shinyBS tipify popify bsCollapse bsCollapsePanel
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyjqui jqui_resizable
#' @importFrom plotly plotlyOutput
#'
#' @rdname INTERNAL_create_tab_qc
.create_tab_qc <- function(meta.choices) {
    # nocov start
    tabPanel(
        title = "QC",
        id = "qc",
        sidebarLayout(
            sidebarPanel(
                width = 2,
                h4("Plot Controls"),
                hr(),
                bsCollapse(
                    open = "pca.settings",
                    bsCollapsePanel(
                        title = span(icon("plus"), "PCA Settings"), value = "pca.settings", style = "info",
                        uiOutput("pca.comps"),
                        conditionalPanel(
                            condition = "input['keep.top.n'] == false",
                            numericInput("var.remove", "Remove this proportion of features ranked by variance:",
                                min = 0, max = 1, step = 0.01, value = 0
                            )
                        ),
                        conditionalPanel(
                            condition = "input['keep.top.n'] == true",
                            numericInput("var.n.keep", "Number of features to retain by variance:",
                                min = 2, max = Inf, step = 1, value = 500
                            )
                        ),
                        fluidRow(
                            column(
                                6,
                                tipify(
                                    prettyCheckbox("center", strong("Center data"), TRUE,
                                        bigger = FALSE,
                                        animation = "smooth", status = "success",
                                        icon = icon("check"), width = "100%"
                                    ),
                                    "Zero center the data before performing PCA.", "right",
                                    options = list(container = "body")
                                ),
                                tipify(
                                    prettyCheckbox("keep.top.n", strong("Limit by top N features"), FALSE,
                                        bigger = FALSE,
                                        animation = "smooth", status = "success",
                                        icon = icon("check"), width = "100%"
                                    ),
                                    "Limit PCA to top N features ranked by variance.", "right",
                                    options = list(container = "body")
                                ),
                                tipify(
                                    prettyCheckbox("meta.filt", strong("Filter via metadata table"), TRUE,
                                        bigger = FALSE,
                                        animation = "smooth", status = "success",
                                        icon = icon("check"), width = "100%"
                                    ),
                                    "Filter PCA samples to those in metadata table.", "right",
                                    options = list(container = "body")
                                )
                            ),
                            column(
                                6,
                                tipify(
                                    prettyCheckbox("scale", strong("Scale data"), TRUE,
                                        bigger = FALSE,
                                        animation = "smooth", status = "success",
                                        icon = icon("check"), width = "100%"
                                    ),
                                    "Scale the data to have unit variance before performing PCA.", "right",
                                    options = list(container = "body")
                                ),
                                tipify(
                                    numericInput("pca.pt.size", "Point size:",
                                        min = 1, max = Inf, step = 0.1, value = 12
                                    ),
                                    "PCA biplot point size.", "right",
                                    options = list(container = "body")
                                )
                            )
                        ),
                        fluidRow(
                            column(6, tipify(
                                selectizeInput("bip.color", "Color by:",
                                    # ifelse apparently can't return NULL.
                                    choices = c("", meta.choices), selected = ifelse("Label" %in% meta.choices, "Label", "")
                                ),
                                "Metadata variable by which samples are colored.", "right",
                                options = list(container = "body")
                            )),
                            column(6, tipify(
                                selectizeInput("bip.shape", "Shape by:",
                                    choices = c("", meta.choices)
                                ),
                                "Metadata variable by which samples are shaped.", "right",
                                options = list(container = "body")
                            ))
                        ),
                        fluidRow(
                            column(
                                6,
                                tipify(
                                    prettyCheckbox("bip.twod", strong("Limit to 2D"), TRUE,
                                        bigger = FALSE,
                                        animation = "smooth", status = "success",
                                        icon = icon("check"), width = "100%"
                                    ),
                                    "Limit PCA biplot to 2D.", "right",
                                    options = list(container = "body")
                                )
                            ),
                            column(
                                6,
                                tipify(
                                    prettyCheckbox("bip.loadings", strong("Plot Loadings"), FALSE,
                                        bigger = FALSE,
                                        animation = "smooth", status = "success",
                                        icon = icon("check"), width = "100%"
                                    ),
                                    "Plot top PCA loadings for each PC.", "right",
                                    options = list(container = "body")
                                )
                            ),
                            column(
                                12,
                                tipify(
                                    numericInput("bip.n.loadings", "Loadings:",
                                        min = 0, max = 100, step = 1, value = 3
                                    ),
                                    "Number of PCA loadings to plot (if checked).", "right",
                                    options = list(container = "body")
                                )
                            ),
                            div(actionButton("pca.update", "Update PCA"), align = "center")
                        )
                    ),
                    bsCollapsePanel(
                        title = span(icon("plus"), "Correlation Plot Settings"), value = "corr.settings", style = "info",
                        fluidRow(
                            column(
                                12,
                                tipify(colourInput("corr.max.col", "Maximum color:", value = "#FF0000"),
                                    "Maximum correlation color.", "right",
                                    options = list(container = "body")
                                ),
                                tipify(colourInput("corr.min.col", "Maximum color:", value = "#FFFFFF"),
                                    "Minimum correlation color.", "right",
                                    options = list(container = "body")
                                )
                            ),
                            div(actionButton("corr.update", "Update Correlation Plot"), align = "center")
                        )
                    )
                )
            ),
            mainPanel(
                width = 10,
                fluidRow(
                    column(
                        width = 4,
                        span(
                            popify(icon("circle-info", style = "font-size: 20px"), "Gini Index",
                                c(
                                    "This plot shows the gini index for each sample, ",
                                    "which is a measure of read count inequality between guides. ",
                                    "For initial timepoints and early passages, this should usually be <0.1, ",
                                    "but should increase as selection occurs."
                                ),
                                placement = "bottom", trigger = "hover", options = list(container = "body")
                            ),
                            div(downloadButton("dl_plot.qc.gini", "Download Gini Plot", class = "btn-dl"), style = "display:inline-block; float:right"),
                            withSpinner(jqui_resizable(plotlyOutput("qc.gini")))
                        ),
                        span(
                            popify(icon("circle-info", style = "font-size: 20px"), "Read Distributions",
                                c(
                                    "This plot shows the read distribution across guides for each sample. ",
                                    "The distribution should be normal and relatively tight for initial ",
                                    "timepoints and early passages. As selection occurs, the ",
                                    "distribution will begin to spread, and there may be buildup at the extremes."
                                ),
                                placement = "top", trigger = "hover", options = list(container = "body")
                            ),
                            withSpinner(jqui_resizable(plotlyOutput("qc.histplot"))),
                            div(downloadButton("dl_plot.qc.hist", "Download Read Distribution Plot", class = "btn-dl"), style = "display:inline-block; float:right"),
                        )
                    ),
                    column(
                        width = 4,
                        span(
                            popify(icon("circle-info", style = "font-size: 20px"), "Zero Count sgRNAs",
                                c(
                                    "This plot shows the number of guides with zero reads for each sample. ",
                                    "For early passages and initial timepoints, this count should ideally be zero ",
                                    "but should increase as selection occurs. Useful for assessing library quality."
                                ),
                                placement = "bottom", trigger = "hover", options = list(container = "body")
                            ),
                            div(downloadButton("dl_plot.qc.missed", "Download sgRNA Count Plot", class = "btn-dl"), style = "display:inline-block; float:right"),
                            withSpinner(jqui_resizable(plotlyOutput("qc.missed")))
                        ),
                        span(
                            popify(icon("circle-info", style = "font-size: 20px"), "Sample Correlations",
                                c(
                                    "This plot shows correlation between samples. Typically, initial timepoints ",
                                    "and early passages, even from different tissues or conditions, correlate well, ",
                                    "but diverge more and more as selection occurs."
                                ),
                                placement = "top", trigger = "hover", options = list(container = "body")
                            ),
                            jqui_resizable(plotOutput("qc.corr")),
                            div(downloadButton("dl_plot.qc.corr", "Download Correlation Matrix", class = "btn-dl"), style = "display:inline-block; float:right")
                        )
                    ),
                    column(
                        width = 4,
                        span(
                            popify(icon("circle-info", style = "font-size: 20px"), "Mapping Rates",
                                c("This plot shows read mapping rates for each sample. 50-75% mapped is typical."),
                                placement = "bottom", trigger = "hover", options = list(container = "body")
                            ),
                            div(downloadButton("dl_plot.qc.map", "Download Mapping Rates", class = "btn-dl"), style = "display:inline-block; float:right"),
                            jqui_resizable(plotOutput("qc.map"))
                        ),
                        span(
                            popify(icon("circle-info", style = "font-size: 20px"), "Principal Componenet Analysis",
                                c(
                                    "This biplot shows two (or three) PCs from a principal component analysis. ",
                                    "For initial timepoints, libraries should be nearly identical, even for different tissues or cell lines. ",
                                    "As selection occurs, samples should diverge."
                                ),
                                placement = "top", trigger = "hover", options = list(container = "body")
                            ),
                            withSpinner(jqui_resizable(plotlyOutput("qc.pca"))),
                            div(downloadButton("dl_plot.qc.pca", "Download PCA Plot", class = "btn-dl"), style = "display:inline-block; float:right")
                        ),
                    )
                )
            )
        )
    )
    # nocov end
}


#' Create a tabPanel for the qc summary tab
#'
#' Create a \code{\link{tabPanel}} with UI elements for the qc summary tab.
#'
#' @return
#' A \code{\link{tabPanel}} with UI elements for the qc summary tab.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny tabPanel br
#' @importFrom DT DTOutput
#'
#' @rdname INTERNAL_create_tab_qc_summary
.create_tab_qc_summary <- function() {
    # nocov start
    tabPanel(
        title = "QC Table",
        id = "qc-table",
        br(),
        DTOutput("count.summary")
    )
    # nocov end
}
