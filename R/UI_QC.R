.tab_qc <- function(summ.choices) {
    tabPanel(
        title = "QC",
        id = "qc",
        sidebarLayout(
            sidebarPanel(
                width = 2,
                h4("Plot Controls"),
                hr(),
                bsCollapse(open = "pca.settings",
                        bsCollapsePanel(
                            title = span(icon("plus"), "PCA Settings"), value = "pca.settings", style = "info",
                            uiOutput("pca.comps"),
                            conditionalPanel(
                            condition = "input['keep.top.n'] == false",
                            numericInput("var.remove", "Remove this proportion of features ranked by variance:",
                                            min = 0, max = 1, step = 0.01, value = 0)
                            ),
                            conditionalPanel(
                            condition = "input['keep.top.n'] == true",
                            numericInput("var.n.keep", "Number of features to retain by variance:",
                                            min = 2, max = Inf, step = 1, value = 500)
                            ),
                            fluidRow(
                            column(6,
                                    tipify(prettyCheckbox("center", strong("Center data"), TRUE, bigger = FALSE,
                                                            animation = "smooth", status = "success",
                                                            icon = icon("check"), width = "100%"),
                                            "Zero center the data before performing PCA.", "right", options = list(container = "body")),
                                    tipify(prettyCheckbox("keep.top.n", strong("Limit by top N features"), FALSE, bigger = FALSE,
                                                            animation = "smooth", status = "success",
                                                            icon = icon("check"), width = "100%"),
                                            "Limit PCA to top N features ranked by variance.", "right", options = list(container = "body"))
                            ),
                            column(6,
                                    tipify(prettyCheckbox("scale", strong("Scale data"), TRUE, bigger = FALSE,
                                                            animation = "smooth", status = "success",
                                                            icon = icon("check"), width = "100%"),
                                            "Scale the data to have unit variance before performing PCA.", "right", options = list(container = "body"))
                            )
                            ),
                            tipify(prettyCheckbox("meta.filt", strong("Filter via metadata table"), TRUE, bigger = FALSE,
                                                animation = "smooth", status = "success",
                                                icon = icon("check"), width = "100%"),
                                    "Filter PCA samples to those in metadata table.", "right", options = list(container = "body")),
                            fluidRow(
                            column(6, tipify(selectizeInput("bip.color", "Color by:",
                                                            choices = summ.choices),
                                                "Metadata variable by which samples are colored.", "right", options = list(container = "body"))),
                            column(6, tipify(selectizeInput("bip.shape", "Shape by:",
                                                            choices = summ.choices),
                                                "Metadata variable by which samples are shaped.", "right", options = list(container = "body")))
                            ),
                            fluidRow(
                            column(6,
                                    tipify(prettyCheckbox("bip.twod", strong("Limit to 2D"), TRUE, bigger = FALSE,
                                                            animation = "smooth", status = "success",
                                                            icon = icon("check"), width = "100%"),
                                            "Limit PCA biplot to 2D.", "right", options = list(container = "body"))
                            ),
                            column(6,
                                    tipify(prettyCheckbox("bip.loadings", strong("Plot Loadings"), FALSE, bigger = FALSE,
                                                            animation = "smooth", status = "success",
                                                            icon = icon("check"), width = "100%"),
                                            "Plot top PCA loadings for each PC.", "right", options = list(container = "body")
                                    )
                            ),
                            column(12,
                                tipify(numericInput("bip.n.loadings", "Loadings:",
                                                    min = 0, max = 100, step = 1, value = 5),
                                        "Number of PCA loadings to plot (if checked).", "right", options = list(container = "body")
                                )
                            ),
                            div(actionButton("pca.update", "Update PCA"), align = "center")
                            )
                        )
                )
            ),
            mainPanel(
                width = 10,
                fluidRow(
                    column(width = 4,
                            span(popify(icon("info-circle", style="font-size: 20px"), "Gini Index",
                                        c("This plot shows the gini index for each sample, ",
                                        "which is a measure of read count inequality between guides. ",
                                        "For initial timepoints and early passages, this should usually be <0.1, ",
                                        "but should increase as selection occurs."),
                                        placement = "bottom", trigger = "hover", options = list(container = "body")),
                                jqui_resizable(plotlyOutput("qc.gini"))),
                            span(popify(icon("info-circle", style="font-size: 20px"), "Read Distributions",
                                        c("This plot shows the read distribution across guides for each sample. ",
                                        "The distribution should be normal and relatively tight for initial ",
                                        "timepoints and early passages. As selection occurs, the ",
                                        "distribution will begin to spread, and there may be buildup at the extremes."),
                                        placement = "top", trigger = "hover", options = list(container = "body")),
                                jqui_resizable(plotOutput("qc.histplot")))
                    ),
                    column(width = 4,
                            span(popify(icon("info-circle", style="font-size: 20px"), "Zero Count gRNAs",
                                        c("This plot shows the number of guides with zero reads for each sample. ",
                                        "For early passages and initial timepoints, this count should ideally be zero ",
                                        "but should increase as selection occurs. Useful for assessing library quality."),
                                        placement = "bottom", trigger = "hover", options = list(container = "body")),
                                jqui_resizable(plotlyOutput("qc.missed"))),
                            span(popify(icon("info-circle", style="font-size: 20px"), "Sample Correlations",
                                        c("This plot shows correlation between samples. Typically, initial timepoints ",
                                        "and early passages, even from different tissues or conditions, correlate well, ",
                                        "but diverge more and more as selection occurs."),
                                        placement = "top", trigger = "hover", options = list(container = "body")),
                                jqui_resizable(plotOutput("qc.corr")))
                    ),
                    column(width = 4,
                            span(popify(icon("info-circle", style="font-size: 20px"), "Mapping Rates",
                                        c("This plot shows read mapping rates for each sample. 50-75% mapped is typical."),
                                        placement = "bottom", trigger = "hover", options = list(container = "body")),
                                jqui_resizable(plotOutput("qc.map"))),
                            span(popify(icon("info-circle", style="font-size: 20px"), "Principal Componenet Analysis",
                                        c("This biplot shows two (or three) PCs from a principal component analysis. ",
                                        "For initial timepoints, libraries should be nearly identical, even for different tissues or cell lines. ",
                                        "As selection occurs, samples should diverge."),
                                        placement = "top", trigger = "hover", options = list(container = "body")),
                                withSpinner(jqui_resizable(plotlyOutput("qc.pca")))),
                    )
                )
            )
        )
    )
}