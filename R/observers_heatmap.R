#' Define observers for the gene heatmap tab
#'
#' @param input The Shiny input object from the server function.
#' @param session The Shiny session object from the server function.
#' @param output The Shiny output object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return
#' Observers are set up to build the gene heatmap and respond to click events.
#' A \code{NULL} is invisibly returned.
#'
#' @author Jared Andrews
#'
#' @rdname INTERNAL_create_heatmap_observers
#' @importFrom shiny observeEvent isolate req
#' @importFrom DT renderDT datatable formatStyle
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom InteractiveComplexHeatmap makeInteractiveComplexHeatmap
#' @importFrom circlize colorRamp2
.create_heatmap_observers <- function(input, session, output, robjects) {
    # nocov start
    output$heatmap.info <- renderDT(server = FALSE, {
        datatable(
            data.frame(Info = "Upload gene summary data, enter genes, and update the heatmap."),
            rownames = FALSE,
            options = list(dom = "t", paging = FALSE, ordering = FALSE)
        ) %>% formatStyle(0, target = "row", lineHeight = "50%")
    })

    observeEvent(robjects$gene.data, {
        choices <- .heatmap_value_choices(robjects$gene.data)
        selected <- if (length(choices) == 0) {
            character(0)
        } else if ("LFC" %in% choices) {
            "LFC"
        } else {
            choices[1]
        }
        shiny::updateSelectInput(session,
            "heatmap.value.term",
            choices = choices,
            selected = selected
        )
    }, ignoreNULL = FALSE)

    observeEvent(input$heatmap.update, {
        req(!is.null(robjects$gene.data), length(robjects$gene.data) > 0)

        genes <- .parse_gene_input(isolate(input$heatmap.genes))
        value.term <- isolate(input$heatmap.value.term)

        if (length(genes) == 0) {
            ht <- .empty_heatmap("Enter one or more genes to draw the heatmap.")
            robjects$plot.heatmap <- ht
            makeInteractiveComplexHeatmap(
                input, output, session, ht,
                heatmap_id = "gene_heatmap",
                click_action = function(df, output) {
                    output$heatmap.info <- renderDT(server = FALSE, {
                        datatable(
                            data.frame(Info = "No genes entered."),
                            rownames = FALSE,
                            options = list(dom = "t", paging = FALSE, ordering = FALSE)
                        ) %>% formatStyle(0, target = "row", lineHeight = "50%")
                    })
                }
            )
            return(invisible(NULL))
        }

        if (is.null(value.term) || !nzchar(value.term)) {
            ht <- .empty_heatmap("Select a value to display and update the heatmap.")
            robjects$plot.heatmap <- ht
            makeInteractiveComplexHeatmap(
                input, output, session, ht,
                heatmap_id = "gene_heatmap",
                click_action = function(df, output) {
                    output$heatmap.info <- renderDT(server = FALSE, {
                        datatable(
                            data.frame(Info = "No value selected."),
                            rownames = FALSE,
                            options = list(dom = "t", paging = FALSE, ordering = FALSE)
                        ) %>% formatStyle(0, target = "row", lineHeight = "50%")
                    })
                }
            )
            return(invisible(NULL))
        }

        heatmap.data <- .get_user_heatmap_data(
            genes = genes,
            value.term = value.term,
            gene.data = robjects$gene.data
        )

        if (is.null(heatmap.data)) {
            ht <- .empty_heatmap("No data available for the requested genes.")
            robjects$plot.heatmap <- ht
            makeInteractiveComplexHeatmap(
                input, output, session, ht,
                heatmap_id = "gene_heatmap",
                click_action = function(df, output) {
                    output$heatmap.info <- renderDT(server = FALSE, {
                        datatable(
                            data.frame(Info = "No data available for the requested genes."),
                            rownames = FALSE,
                            options = list(dom = "t", paging = FALSE, ordering = FALSE)
                        ) %>% formatStyle(0, target = "row", lineHeight = "50%")
                    })
                }
            )
            return(invisible(NULL))
        }

        mat <- heatmap.data$matrix
        legend.title <- heatmap.data$legend_title
        missing.genes <- heatmap.data$missing_genes

        scaled.rows <- isolate(input$heatmap.scale)

        if (scaled.rows) {
            mat <- .scale_heatmap_rows(mat)
            legend.title <- paste0(legend.title, " (row Z-score)")
        }

        col.fun <- .heatmap_color_fun(mat, scaled.rows)
        column.title <- .heatmap_column_title(value.term, mat, missing.genes)

        ht <- Heatmap(
            mat,
            name = legend.title,
            col = col.fun,
            cluster_rows = TRUE,
            cluster_columns = isolate(input$heatmap.cluster.columns),
            column_title = column.title,
            row_names_side = "left",
            column_names_side = "bottom",
            na_col = "#f5f5f5"
        )

        ht <- draw(ht)
        robjects$plot.heatmap <- ht

        output$heatmap.info <- renderDT(server = FALSE, {
            datatable(
                data.frame(Info = "Click a heatmap cell to view gene and dataset details."),
                rownames = FALSE,
                options = list(dom = "t", paging = FALSE, ordering = FALSE)
            ) %>% formatStyle(0, target = "row", lineHeight = "50%")
        })

        .heatmap_click_action <- function(df, output) {
            info.df <- .build_heatmap_click_table(df, mat, value.term)

            output$heatmap.info <- renderDT(server = FALSE, {
                datatable(
                    info.df,
                    rownames = FALSE,
                    options = list(dom = "t", paging = FALSE, ordering = FALSE)
                ) %>% formatStyle(0, target = "row", lineHeight = "50%")
            })
        }

        makeInteractiveComplexHeatmap(
            input, output, session, ht,
            heatmap_id = "gene_heatmap",
            click_action = .heatmap_click_action
        )
    })
    # nocov end

    invisible(NULL)
}

.parse_gene_input <- function(x) {
    if (is.null(x) || identical(x, "")) {
        return(character())
    }

    genes <- strsplit(x, ",|\n|\r|\t| ")[[1]]
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    unique(genes)
}

.heatmap_value_choices <- function(gene.data) {
    if (is.null(gene.data) || length(gene.data) == 0) {
        return(character())
    }

    terms <- unique(unlist(lapply(gene.data, function(df) {
        if (!"Gene" %in% colnames(df)) {
            return(character())
        }
        names(df)[vapply(df, is.numeric, logical(1))]
    })))

    terms <- terms[terms != ""]
    sort(unique(terms))
}

.get_user_heatmap_data <- function(genes, value.term, gene.data) {
    if (length(gene.data) == 0) {
        return(NULL)
    }

    dataset.names <- names(gene.data)
    value.frames <- lapply(gene.data, function(df) {
        if (!"Gene" %in% colnames(df)) {
            return(NULL)
        }

        if (!value.term %in% colnames(df)) {
            data.frame(Gene = df$Gene, Value = NA_real_, stringsAsFactors = FALSE)
        } else {
            vals <- suppressWarnings(as.numeric(df[[value.term]]))
            data.frame(Gene = df$Gene, Value = vals, stringsAsFactors = FALSE)
        }
    })

    names(value.frames) <- dataset.names
    value.frames <- value.frames[!vapply(value.frames, is.null, logical(1))]

    if (length(value.frames) == 0) {
        return(NULL)
    }

    all.genes <- unique(unlist(lapply(value.frames, `[[`, "Gene")))
    genes <- unique(genes)
    requested.genes <- genes[genes %in% all.genes]

    if (length(requested.genes) == 0) {
        return(NULL)
    }

    gene.order <- requested.genes
    mat <- matrix(NA_real_, nrow = length(gene.order), ncol = length(value.frames),
        dimnames = list(gene.order, names(value.frames))
    )

    for (i in seq_along(value.frames)) {
        df <- value.frames[[i]]
        idx <- match(df$Gene, gene.order, nomatch = 0)
        keep <- idx > 0 & !is.na(df$Value)
        if (any(keep)) {
            mat[idx[keep], i] <- df$Value[keep]
        }
    }

    mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
    mat <- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]

    if (!nrow(mat) || !ncol(mat)) {
        return(NULL)
    }

    list(
        matrix = mat,
        legend_title = value.term,
        missing_genes = setdiff(genes, rownames(mat))
    )
}

.scale_heatmap_rows <- function(mat) {
    res <- mat
    for (i in seq_len(nrow(mat))) {
        x <- mat[i, ]
        keep <- !is.na(x)
        if (!any(keep)) {
            next
        }
        m <- mean(x[keep])
        s <- stats::sd(x[keep])
        if (is.na(s) || s == 0) {
            res[i, keep] <- 0
        } else {
            res[i, keep] <- (x[keep] - m) / s
        }
    }
    res
}

.heatmap_color_fun <- function(mat, scaled) {
    rng <- range(mat, na.rm = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2])) {
        return(colorRamp2(c(-1, 0, 1), c("#313695", "#FFFFBF", "#A50026")))
    }
    if (scaled) {
        limit <- max(abs(rng), na.rm = TRUE)
        if (!is.finite(limit) || limit == 0) {
            limit <- 1
        }
        rng <- c(-limit, 0, limit)
    } else if (rng[1] == rng[2]) {
        delta <- if (rng[1] == 0) 1 else abs(rng[1]) * 0.1
        rng <- c(rng[1] - delta, rng[1], rng[1] + delta)
    } else {
        rng <- c(rng[1], mean(rng), rng[2])
    }
    colorRamp2(rng, c("#313695", "#FFFFBF", "#A50026"))
}

.heatmap_column_title <- function(value.term, mat, missing.genes) {
    title <- sprintf("%s (%d genes x %d datasets)", value.term, nrow(mat), ncol(mat))

    if (length(missing.genes) > 0) {
        title <- paste0(title, "\nMissing: ", paste(missing.genes, collapse = ", "))
    }

    title
}

.build_heatmap_click_table <- function(df, mat, value.term) {
    if (is.null(df) || nrow(df) == 0) {
        return(data.frame(Info = "Click a heatmap cell to view gene and dataset details."))
    }

    row.idx <- df$row_index
    col.idx <- df$column_index
    gene <- rownames(mat)[row.idx]
    dataset <- colnames(mat)[col.idx]
    value <- mat[row.idx, col.idx]

    out <- data.frame(
        Gene = gene,
        Dataset = dataset,
        Value = round(value, 4),
        stringsAsFactors = FALSE
    )

    if (!is.null(value.term) && nzchar(value.term)) {
        names(out)[names(out) == "Value"] <- value.term
    }

    out
}
