#' Create download handlers for all interactive plots.
#'
#' @param output The Shiny output object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return A \linkS4class{NULL} is invisibly returned
#' and handlers for plotly plot downloads are added to \code{output}.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny downloadHandler
#' @importFrom htmlwidgets saveWidget
#' @importFrom shinyjqui jqui_resizable
#' @rdname INTERNAL_create_dl_outputs
.create_dl_outputs <- function(output, robjects) {
    plotters <- list(
        "plot.qc.pca", "plot.qc.missed", "plot.qc.gini", "plot.qc.hist",
        "plot.gene1.vol", "plot.gene1.rank", "plot.gene1.lawn",
        "plot.gene2.vol", "plot.gene2.rank", "plot.gene2.lawn",
        "plot.sgrna1.counts", "plot.sgrna1.rank",
        "plot.sgrna2.counts", "plot.sgrna2.rank",
        "plot.depmap.essplot", "plot.depmap.expplot", "plot.depmap.cnplot",
        "plot.depmap.lineages", "plot.depmap.sublineage"
    )

    # nocov start
    lapply(plotters, function(x) {
        output[[paste0("dl_", x)]] <- downloadHandler(
            filename = function() {
                paste(x, "-", Sys.Date(), ".html", sep = "")
            },
            content = function(file) {
                # export plotly html widget as a temp file to download.
                saveWidget(jqui_resizable(robjects[[x]]),
                    file,
                    selfcontained = TRUE
                )
            }
        )
    })
    # nocov end

    invisible(NULL)
}
