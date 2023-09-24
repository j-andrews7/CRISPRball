#' Render general outputs for the Comparisons tab
#'
#' Create rendering expressions for the Comparisons tab outputs.
#'
#' @param input The Shiny input object from the server function.
#' @param output The Shiny output object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return A \linkS4class{NULL} is invisibly returned
#' and rendering expressions for Gene tabs features are added to \code{output}.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny isolate selectInput tagList renderUI
#' @importFrom plotly renderPlotly
#' @importFrom DT renderDT datatable formatStyle
#' @rdname INTERNAL_create_comparisons_outputs
.create_comparisons_outputs <- function(input, output, robjects) {
    # nocov start
    output$comp.term.options <- renderUI({
        req(robjects$set1.genes)
        df <- robjects$set1.genes
        # Get only numeric variables.
        choices <- names(df)[vapply(df, is.numeric, logical(1))]

        tagList(
            column(
                6,
                selectInput("comp.esterm", "Effect size term:",
                    choices = choices,
                    selected = ifelse("LFC" %in% choices, "LFC",
                        ifelse("beta" %in% choices, "beta", choices[1])
                    )
                )
            ),
            column(
                6,
                selectInput("comp.sigterm", "Significance term:",
                    choices = choices,
                    selected = ifelse("fdr" %in% choices, "fdr",
                        ifelse("pval" %in% choices, "pval", choices[1])
                    )
                )
            )
        )
    })
    # nocov end
}
