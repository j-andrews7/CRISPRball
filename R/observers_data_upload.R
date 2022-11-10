#' Define observers for data upload
#'
#' Define a series of observers to track additional new data being uploaded. These observers
#' also enable tabs as data is made available.
#'
#' @param input The Shiny input object from the server function.
#' @param session The Shiny session object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return
#' Observers are set up to monitor the data upload UI elements that can accept new data.
#' A \code{NULL} is invisibly returned.
#'
#' @author Jared Andrews
#'
#' @rdname INTERNAL_create_upload_observers
#' @importFrom shiny observeEvent updateSelectizeInput
#' @importFrom shinyWidgets updatePickerInput
#' @importFrom shinyjs js
.create_upload_observers <- function(input, session, robjects) {
    # nocov start
    observeEvent(input$geneSummaryFiles, {
        robjects$gene.data <- .gene_summ_ingress(input$geneSummaryFiles)
        if (!is.null(robjects$gene.data)) {
            js$enableTab("Gene (Overview)")
            js$enableTab("Gene Summary Tables")
            updateSelectizeInput(session, "gene.sel1", choices = names(robjects$gene.data), server = TRUE)
            updateSelectizeInput(session, "gene.sel2", choices = names(robjects$gene.data), server = TRUE)
        }
    })

    observeEvent(input$sgrnaSummaryFiles, {
        robjects$sgrna.data <- .sgrna_summ_ingress(input$sgrnaSummaryFiles)
        if (!is.null(robjects$sgrna.data)) {
            js$enableTab("sgRNA")
            js$enableTab("sgRNA Summary Tables")
            updateSelectizeInput(session, "sgrna.sel1", choices = names(robjects$sgrna.data), server = TRUE)
            updateSelectizeInput(session, "sgrna.sel2", choices = names(robjects$sgrna.data), server = TRUE)
            updatePickerInput(session, "sgrna.gene", choices = unique(c(robjects$sgrna.data[[1]]$Gene)))
        }
    })

    observeEvent(input$countSummary, {
        robjects$count.summary <- read.delim(input$countSummary$datapath)
        if (!is.null(robjects$count.summary)) {
            js$enableTab("QC")
            js$enableTab("QC Table")
            updateSelectizeInput(session, "bip.color", choices = c("", colnames(robjects$count.summary)), server = TRUE)
            updateSelectizeInput(session, "bip.shape", choices = c("", colnames(robjects$count.summary)), server = TRUE)
        }
    })

    observeEvent(input$countNormFile, {
        robjects$norm.counts <- read.delim(input$countNormFile$datapath)
        if (!is.null(robjects$norm.counts)) {
            js$enableTab("QC")
        }
    })
    # nocov end

    invisible(NULL)
}
