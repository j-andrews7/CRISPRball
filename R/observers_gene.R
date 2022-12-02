#' Define observers for Gene tabs
#'
#' Define a series of observers to parse gene data based on new data upload or dataset selection.
#'
#' @param input The Shiny input object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return
#' Observers are set up to monitor the reactives necessary for the Gene plots.
#' A \code{NULL} is invisibly returned.
#'
#' @author Jared Andrews
#'
#' @rdname INTERNAL_create_gene_observers
.create_gene_observers <- function(input, robjects) {
    # nocov start
    observeEvent(c(input$gene.sel1, input$gene.update), {
        df <- robjects$gene.data[[input$gene.sel1]]
        robjects$set1.genes <- .gene_ingress(df,
            sig.thresh = input$gene.fdr.th, lfc.thresh = input$gene.lfc.th,
            positive.ctrl.genes = robjects$positive.ctrl.genes, essential.genes = robjects$essential.genes,
            depmap.genes = robjects$depmap.gene
        )
    })
    # nocov end

    # nocov start
    observeEvent(c(input$gene.sel2, input$gene.update), {
        if (length(robjects$gene.data) > 1) {
            df <- robjects$gene.data[[input$gene.sel2]]
            if (!is.null(df)) {
                robjects$set2.genes <- .gene_ingress(df,
                    sig.thresh = input$gene.fdr.th, lfc.thresh = input$gene.lfc.th,
                    positive.ctrl.genes = robjects$positive.ctrl.genes, essential.genes = robjects$essential.genes,
                    depmap.genes = robjects$depmap.gene
                )

                # Get overlapping hits between sets if needed.
                s1 <- robjects$set1.genes
                s2 <- robjects$set2.genes

                set1.hits <- s1$id[s1$hit_type %in% c("neg", "pos")]
                set2.hits <- s2$id[s2$hit_type %in% c("neg", "pos")]

                robjects$common.hits <- set1.hits[set1.hits %in% set2.hits]
            }
        }
    })
    # nocov end

    # On click, the key field of the event data contains the gene symbol.
    # Add that gene to the set of all "selected" genes. Double click will clear all labels.
    # nocov start
    plot.suf <- list("volc1", "volc2", "rank1", "rank2", "lawn1", "lawn2")
    lapply(plot.suf, function(x) {
        obj <- paste0("clicked.", x)
        observeEvent(event_data("plotly_click", source = paste0(robjects$h.id, "_", x)), {
            gene <- event_data("plotly_click", source = paste0(robjects$h.id, "_", x))
            gene_old_new <- rbind(robjects[[obj]], gene)
            keep <- gene_old_new[gene_old_new$customdata %in% names(which(table(gene_old_new$customdata) == 1)), ]

            if (nrow(keep) == 0) {
                robjects[[obj]] <- NULL
            } else {
                robjects[[obj]] <- keep
            }
        })

        observeEvent(event_data("plotly_doubleclick", source = paste0(robjects$h.id, "_", x)), {
            robjects[[obj]] <- NULL
        })
    })
    # nocov end

    invisible(NULL)
}
