#' Define observers for sgRNA tab
#'
#' Define a series of observers to parse sgRNA data based 
#'   on new data upload or dataset selection.
#'
#' @param input The Shiny input object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return
#' Observers are set up to monitor the reactives necessary for the sgRNA plots.
#' A \code{NULL} is invisibly returned.
#'
#' @author Jared Andrews
#'
#' @rdname INTERNAL_create_sgrna_observers
.create_sgrna_observers <- function(input, robjects) {
    # nocov start
    observe({
        df <- robjects$sgrna.data[[input$sgrna.sel1]]
        df$Rank <- rank(df$LFC)
        robjects$set1.sgrnas <- df

        if (length(robjects$sgrna.data) > 1) {
            df <- robjects$sgrna.data[[input$sgrna.sel2]]
            df$Rank <- rank(df$LFC)
            robjects$set2.sgrnas <- df
        }
    })
    # nocov end

    invisible(NULL)
}
