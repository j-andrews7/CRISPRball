#' Define observers for UI element enabling/disabling/hiding
#'
#' Define a series of observers to adjust UI usability/visibility based on data availability.
#'
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return
#' Observers are set up to monitor the reactives necessary for the dynamic UI elements.
#' A \code{NULL} is invisibly returned.
#'
#' @importFrom shinyjs show hide disable enable
#' @importFrom shiny observeEvent
#'
#' @author Jared Andrews
#'
#' @rdname INTERNAL_create_ui_observers
.create_ui_observers <- function(robjects) {
    # nocov start
    observe({
        if (is.null(robjects$depmap.gene)) {
            shinyjs::hide("dep.crispr.ess")
            shinyjs::hide("dep.crispr.sel")
            shinyjs::hide("dep.rnai.ess")
            shinyjs::hide("dep.rnai.sel")

            shinyjs::hide("comp.dep.crispr.ess")
            shinyjs::hide("comp.dep.crispr.sel")
            shinyjs::hide("comp.dep.rnai.ess")
            shinyjs::hide("comp.dep.rnai.sel")

            # Hide depmap tab if database not provided.
            # Tried disable, still looks/feels selectable which may be confusing.
            shinyjs::hide(selector = '.navbar-nav a[data-value="DepMap"')
        }
    })
    # nocov end

    # nocov start
    observe({
        if (is.null(robjects$essential.genes)) {
            shinyjs::hide("rem.ess")
            shinyjs::hide("comp.rem.ess")
        }
    })
    # nocov end

    # nocov start
    observe({
        if (is.null(robjects$positive.ctrl.genes)) {
            shinyjs::hide("rem.pos")
            shinyjs::hide("comp.rem.pos")
        }
    })
    # nocov end

    # Disable certain inputs if only one dataset provided.
    # nocov start
    observeEvent(robjects$gene.data, {
        if (length(robjects$gene.data) == 1) {
            shinyjs::disable("gene.sel2")
            shinyjs::hide("highlight.common")
            shinyjs::hide("dl_plot.gene2.vol")
            shinyjs::hide("dl_plot.gene2.rank")
            shinyjs::hide("dl_plot.gene2.lawn")
        } else {
            shinyjs::enable("gene.sel2")
            shinyjs::show("highlight.common")
            shinyjs::show("dl_plot.gene2.vol")
            shinyjs::show("dl_plot.gene2.rank")
            shinyjs::show("dl_plot.gene2.lawn")
        }
    })
    # nocov end

    # nocov start
    observeEvent(robjects$sgrna.data, {
        if (length(robjects$sgrna.data) == 1) {
            shinyjs::disable("sgrna.sel2")
            shinyjs::hide("dl_plot.sgrna2.counts")
            shinyjs::hide("dl_plot.sgrna2.rank")
        } else {
            shinyjs::enable("sgrna.sel2")
            shinyjs::show("dl_plot.sgrna2.counts")
            shinyjs::show("dl_plot.sgrna2.rank")
        }
    })
    # nocov end

    invisible(NULL)
}
