#' Define observers for the Comparison Tab
#'
#' Define observers to generate UpSet plots and display clicked data in DataTable outputs.
#'
#' @param input The Shiny input object from the server function.
#' @param session The Shiny session object from the server function.
#' @param output The Shiny output object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return
#' Observers are set up to monitor the dataset comparison UI element and generate UpSet plots.
#' A \code{NULL} is invisibly returned.
#'
#' @author Jared Andrews
#'
#' @rdname INTERNAL_create_comparison_observers
#' @importFrom shiny observeEvent updateSelectizeInput
#' @importFrom shinyWidgets updatePickerInput
#' @importFrom shinyjs js
.create_comparisons_observers <- function(input, session, output, robjects) {
    # nocov start
    if (length(robjects$gene.data) > 1) {
        observeEvent(input$comp.update,
            {
                # Get hits for each dataset based on thresholds.
                robjects$comps <- lapply(input$comp.sets, function(x) {
                    df <- robjects$gene.data[[x]]
                    .gene_ingress(df,
                        sig.thresh = input$comp.fdr.th, lfc.thresh = input$comp.lfc.th,
                        positive.ctrl.genes = robjects$positive.ctrl.genes, essential.genes = robjects$essential.genes, depmap.genes = depmap.gene
                    )
                })

                names(robjects$comps) <- input$comp.sets

                robjects$comp.pos.genes <- lapply(robjects$comps, .remove_unwanted_hits,
                    robjects = robjects, input = input
                )

                names(robjects$comp.pos.genes) <- input$comp.sets

                robjects$comp.neg.genes <- lapply(robjects$comps, .remove_unwanted_hits,
                    robjects = robjects, input = input, pos = FALSE
                )

                names(robjects$comp.neg.genes) <- input$comp.sets

                # Get the combination matrices.
                robjects$pos.m <- make_comb_mat(robjects$comp.pos.genes)
                robjects$neg.m <- make_comb_mat(robjects$comp.neg.genes)

                # Make upset plots.
                ht.pos <- draw(UpSet(robjects$pos.m))
                ht.neg <- draw(UpSet(robjects$neg.m))

                #TODO: See if there's a way to move these out of this function.
                .pos_click_action <- function(df, output) {
                    if (!is.null(df)) {
                        # Get combination members and shared hits.
                        comb.name <- comb_name(robjects$pos.m[, df$column_index])
                        comb.genes <- extract_comb(robjects$pos.m, comb.name)
                        comb.members <- set_name(robjects$pos.m)[unlist(strsplit(comb.name, "")) == "1"]

                        # Make output df from genes using info from all combination members.
                        member.dfs <- lapply(comb.members, function(x) {
                            og.df <- robjects$comps[[x]]
                            og.df <- og.df[og.df$id %in% comb.genes, ]
                            rownames(og.df) <- og.df$id
                            og.df <- og.df[, c("FDR", "LFC")]
                            colnames(og.df) <- paste0(colnames(og.df), "_", x)
                            og.df
                        })

                        comp_pos_df <- Reduce(
                            function(dtf1, dtf2) {
                                out.d <- merge(dtf1, dtf2, by = "row.names", all.x = TRUE)
                                rownames(out.d) <- out.d$Row.names
                                out.d$Row.names <- NULL
                                out.d
                            },
                            member.dfs
                        )
                    } else {
                        comp_pos_df <- data.frame(
                            Gene = character()
                        )
                        comb.members <- NULL
                    }

                    output$comp.pos.info <- renderDT(server = FALSE, {
                        DT::datatable(comp_pos_df,
                            rownames = TRUE,
                            filter = "top",
                            extensions = c("Buttons"),
                            caption = paste("Positively Selected Hits in:", paste0(comb.members, collapse = ", ")),
                            options = list(
                                search = list(regex = TRUE),
                                pageLength = 10,
                                dom = "Blfrtip",
                                buttons = c("copy", "csv", "excel", "pdf", "print")
                            )
                        ) %>% DT::formatStyle(0, target = "row", lineHeight = "30%")
                    })
                }

                .neg_click_action <- function(df, output) {
                    if (!is.null(df)) {
                        # Get combination members and shared hits.
                        comb.name <- comb_name(robjects$neg.m[, df$column_index])
                        comb.genes <- extract_comb(robjects$neg.m, comb.name)
                        comb.members <- set_name(robjects$neg.m)[unlist(strsplit(comb.name, "")) == "1"]

                        # Make output df from genes using info from all combination members.
                        member.dfs <- lapply(comb.members, function(x) {
                            og.df <- robjects$comps[[x]]
                            og.df <- og.df[og.df$id %in% comb.genes, ]
                            rownames(og.df) <- og.df$id
                            og.df <- og.df[, c("FDR", "LFC")]
                            colnames(og.df) <- paste0(colnames(og.df), "_", x)
                            og.df
                        })

                        comp_neg_df <- Reduce(
                            function(dtf1, dtf2) {
                                out.d <- merge(dtf1, dtf2, by = "row.names", all.x = TRUE)
                                rownames(out.d) <- out.d$Row.names
                                out.d$Row.names <- NULL
                                out.d
                            },
                            member.dfs
                        )
                    } else {
                        comp_neg_df <- data.frame(
                            Gene = character()
                        )
                        comb.members <- NULL
                    }

                    output$comp.neg.info <- renderDT(server = FALSE, {
                        DT::datatable(comp_neg_df,
                            rownames = TRUE,
                            filter = "top",
                            extensions = c("Buttons"),
                            caption = paste("Negatively Selected Hits in:", paste0(comb.members, collapse = ", ")),
                            options = list(
                                search = list(regex = TRUE),
                                pageLength = 10,
                                dom = "Blfrtip",
                                buttons = c("copy", "csv", "excel", "pdf", "print")
                            )
                        ) %>% DT::formatStyle(0, target = "row", lineHeight = "30%")
                    })
                }

                # Make the output.
                makeInteractiveComplexHeatmap(input, output, session, ht.pos,
                    heatmap_id = "overlap_pos", click_action = .pos_click_action
                )

                makeInteractiveComplexHeatmap(input, output, session, ht.neg,
                    heatmap_id = "overlap_neg", click_action = .neg_click_action
                )
            },
            ignoreInit = TRUE
        )
    }
    # nocov end

    invisible(NULL)
}


#' Remove unwanted hits in Comparisons UpSet plots and tables
#'
#' @param df A data.frame containing gene hits.
#' @param robjects A reactive list of values generated in the server function.
#' @param input The Shiny input object from the server function.
#' @param pos A boolean indicating whether to get genes to remove for positive selection.
#'   FALSE indicates negative selection.
#'
#' @return
#' Character vector of gene identifiers to remove.
#'
#' @author Jared Andrews
#' @rdname INTERNAL_remove_unwanted_hits
.remove_unwanted_hits <- function(df, robjects, input, pos = TRUE) {
    to.remove <- c()

    if (input$comp.rem.ess) {
        to.remove <- c(to.remove, robjects$essential.genes)
    }

    if (input$comp.rem.pos) {
        to.remove <- c(to.remove, robjects$positive.ctrl.genes)
    }

    if (input$comp.dep.crispr.ess) {
        to.remove <- c(to.remove, df$id[df$DepMap_CRISPR_Essential == TRUE])
    }

    if (input$comp.dep.rnai.ess) {
        to.remove <- c(to.remove, df$id[df$DepMap_RNAi_Essential == TRUE])
    }

    if (input$comp.dep.crispr.sel) {
        to.remove <- c(to.remove, df$id[df$DepMap_CRISPR_Selective == TRUE])
    }

    if (input$comp.dep.crispr.sel) {
        to.remove <- c(to.remove, df$id[df$DepMap_RNAi_Selective == TRUE])
    }

    if (pos) {
        out <- df$id[!df$id %in% unique(to.remove) & df$FDR < input$comp.fdr.th & df$LFC > input$comp.lfc.th]
    } else {
        out <- df$id[!df$id %in% unique(to.remove) & df$FDR < input$comp.fdr.th & df$LFC < -input$comp.lfc.th]
    }
}
