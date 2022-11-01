#' Define observers for data upload
#'
#' Define a series of observers to track additional new data being uploaded.
#'
#' @param input The Shiny input object from the server function.
#' @param session The Shiny session object from the server function.
#' @param robjects A reactive list of values generated in the \code{\link{iSEE}} app.
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
    observeEvent(input$geneSummaryFiles, {
        new.data <- .gene_summ_ingress(input$geneSummaryFiles)
        robjects$gene.data <- new.data
        if (!is.null(robjects$gene.data)) {
        js$enableTab('Gene (Overview)')
        js$enableTab('Gene Summary Tables')
        updateSelectizeInput(session, 'gene.sel1', choices = names(robjects$gene.data), server = TRUE)
        updateSelectizeInput(session, 'gene.sel2', choices = names(robjects$gene.data), server = TRUE)
        }
    })

    observeEvent(input$sgrnaSummaryFiles, {
        new.data <- .sgrna_summ_ingress(input$sgrnaSummaryFiles)
        robjects$sgrna.data <- new.data
        if (!is.null(robjects$sgrna.data)) {
        js$enableTab('sgRNA')
        js$enableTab('sgRNA Summary Tables')
        updateSelectizeInput(session, 'sgrna.sel1', choices = names(robjects$sgrna.data), server = TRUE)
        updateSelectizeInput(session, 'sgrna.sel2', choices = names(robjects$sgrna.data), server = TRUE)
        updatePickerInput(session, 'sgrna.gene', choices = unique(c(robjects$sgrna.data[[1]]$Gene)))
        }
    })

    observeEvent(input$countSummary, {
        new.data <- read.delim(input$countSummary$datapath)
        robjects$count.summary <- new.data
        if (!is.null(robjects$count.summary)) {
        js$enableTab('QC')
        js$enableTab('QC Table')
        updateSelectizeInput(session, 'bip.color', choices = c('', colnames(robjects$count.summary)), server = TRUE)
        updateSelectizeInput(session, 'bip.shape', choices = c('', colnames(robjects$count.summary)), server = TRUE)
        }
    })

    observeEvent(input$countNormFile, {
        new.data <- read.delim(input$countNormFile$datapath)
        robjects$norm.counts <- new.data
        if (!is.null(robjects$norm.counts)) {
        js$enableTab('QC')
        }
    })
}