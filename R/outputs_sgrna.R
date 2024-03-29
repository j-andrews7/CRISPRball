#' Render general outputs for the sgRNA tabs
#'
#' Create rendering expressions for the sgRNA tab outputs.
#'
#' @param input The Shiny input object from the server function.
#' @param output The Shiny output object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return A \linkS4class{NULL} is invisibly returned
#' and rendering expressions for sgRNA tabs features are added to \code{output}.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny isolate
#' @importFrom plotly renderPlotly
#' @importFrom DT renderDT datatable formatStyle %>%
#' @rdname INTERNAL_create_sgrna_outputs
.create_sgrna_outputs <- function(input, output, robjects) {
    # nocov start
    output$sgrna.rank.options <- renderUI({
        req(robjects$set1.sgrnas)
        df <- robjects$set1.sgrnas

        tagList(
            selectInput("sgrna.rankby", "Rank by:", choices = names(df), selected = ifelse("LFC" %in% names(df), "LFC", NULL))
        )
    })
    # nocov end

    # nocov start
    output$sgrna1.summary <- renderDT(server = FALSE, {
        req(robjects$set1.sgrnas)

        df <- robjects$set1.sgrnas

        datatable(df,
            rownames = FALSE,
            filter = "top",
            extensions = c("Buttons"),
            caption = paste0(input$sgrna.sel1, " sgRNA Summary"),
            options = list(
                search = list(regex = TRUE),
                pageLength = 10,
                dom = "Blfrtip",
                buttons = c("copy", "csv", "excel", "pdf", "print")
            )
        ) %>% formatStyle(0, target = "row", lineHeight = "50%")
    })
    # nocov end

    # nocov start
    output$sgrna1.counts <- renderPlotly({
        req(robjects$set1.sgrnas, input$sgrna.gene)

        df <- robjects$set1.sgrnas
        df <- df[df$Gene == input$sgrna.gene, ]

        fig <- .make_sgrna_pairplot(df)

        robjects$plot.sgrna1.counts <- fig
        fig
    })
    # nocov end

    # nocov start
    output$sgrna1.rank <- renderPlotly({
        req(robjects$set1.sgrnas, input$sgrna.rankby)
        input$rank.update

        df <- robjects$set1.sgrnas

        hov.info <- c("Gene")

        highlight <- NULL
        highlight <- df$sgrna[df$Gene == input$sgrna.gene]

        fig <- plot_rank(
            res = df,
            ylim = list(
                min(df[[isolate(input$sgrna.rankby)]]) - 0.5,
                max(df[[isolate(input$sgrna.rankby)]]) + 0.5
            ),
            y.thresh = 0,
            y.lines = FALSE,
            sig.thresh = 0,
            h.id = robjects$h.id,
            h.id.suffix = "_sgrank1",
            sig.term = "FDR",
            rank.term = isolate(input$sgrna.rankby),
            rank.ascending = isolate(input$sgrna.rank.ascending),
            feat.term = "sgrna",
            hover.info = hov.info,
            fs = NULL,
            up.color = "#A6A6A6",
            down.color = "#A6A6A6",
            insig.color = "#A6A6A6",
            sig.opacity = 1,
            insig.opacity = 1,
            sig.size = 5,
            insig.size = 5,
            label.size = 8,
            webgl = TRUE,
            webgl.ratio = 7,
            show.counts = FALSE,
            show.hl.counts = FALSE,
            counts.size = 8,
            highlight.featsets = NULL,
            highlight.feats = highlight,
            featsets = NULL,
            highlight.feats.color = "red",
            highlight.feats.size = 8,
            highlight.feats.opac = 1,
            highlight.feats.linecolor = "black",
            highlight.feats.linewidth = 0.5,
            highlight.featsets.color = "#A6A6A6",
            highlight.featsets.size = 7,
            highlight.featsets.opac = 1,
            highlight.featsets.linecolor = "black",
            highlight.featsets.linewidth = 0.5,
            highlight.feats.label = FALSE
        )

        robjects$plot.sgrna1.rank <- fig
        fig
    })
    # nocov end

    # nocov start
    output$sgrna1.detail <- renderDT({
        req(robjects$set1.sgrnas, input$sgrna.gene)

        df <- robjects$set1.sgrnas
        df <- df[df$Gene == input$sgrna.gene, ]

        target <- which(names(df) %in% c(
            "control_mean", "treat_mean",
            "control_var", "adj_var", "high_in_treatment", "p.low",
            "p.high", "p.twosided", "score"
        )) - 1

        datatable(df,
            rownames = FALSE,
            filter = "top",
            extensions = c("Buttons"),
            caption = paste0(input$sgrna.sel1, " ", input$sgrna.gene, " sgRNA Details"),
            options = list(
                pageLength = 10,
                dom = "Blfrtip",
                buttons = c("copy", "csv", "excel", "pdf", "print"),
                columnDefs = list(list(visible = FALSE, targets = target))
            )
        ) %>% formatStyle(0, target = "row", lineHeight = "50%")
    })
    # nocov end

    # nocov start
    output$sgrna2.summary <- renderDT(server = FALSE, {
        req(robjects$set2.sgrnas)

        df <- robjects$set2.sgrnas

        datatable(df,
            rownames = FALSE,
            filter = "top",
            extensions = c("Buttons"),
            caption = paste0(input$sgrna.sel2, " sgRNA Summary"),
            options = list(
                search = list(regex = TRUE),
                pageLength = 10,
                dom = "Blfrtip",
                buttons = c("copy", "csv", "excel", "pdf", "print")
            )
        ) %>% formatStyle(0, target = "row", lineHeight = "50%")
    })
    # nocov end

    # nocov start
    output$sgrna2.counts <- renderPlotly({
        req(robjects$set2.sgrnas, input$sgrna.gene)

        df <- robjects$set2.sgrnas
        df <- df[df$Gene == input$sgrna.gene, ]

        fig <- .make_sgrna_pairplot(df)

        robjects$plot.sgrna2.counts <- fig
        fig
    })
    # nocov end

    # nocov start
    output$sgrna2.rank <- renderPlotly({
        req(robjects$set2.sgrnas, input$sgrna.rankby)
        input$rank.update

        df <- robjects$set2.sgrnas

        hov.info <- c("Gene")

        highlight <- NULL
        highlight <- df$sgrna[df$Gene == input$sgrna.gene]

        fig <- plot_rank(
            res = df,
            ylim = list(
                min(df[[isolate(input$sgrna.rankby)]]) - 0.5,
                max(df[[isolate(input$sgrna.rankby)]]) + 0.5
            ),
            y.thresh = 0,
            y.lines = FALSE,
            sig.thresh = 0,
            h.id = robjects$h.id,
            h.id.suffix = "_sgrank1",
            sig.term = "FDR",
            rank.term = isolate(input$sgrna.rankby),
            rank.ascending = isolate(input$sgrna.rank.ascending),
            feat.term = "sgrna",
            hover.info = hov.info,
            fs = NULL,
            up.color = "#A6A6A6",
            down.color = "#A6A6A6",
            insig.color = "#A6A6A6",
            sig.opacity = 1,
            insig.opacity = 1,
            sig.size = 5,
            insig.size = 5,
            label.size = 8,
            webgl = TRUE,
            webgl.ratio = 7,
            show.counts = FALSE,
            show.hl.counts = FALSE,
            counts.size = 8,
            highlight.featsets = NULL,
            highlight.feats = highlight,
            featsets = NULL,
            highlight.feats.color = "red",
            highlight.feats.size = 8,
            highlight.feats.opac = 1,
            highlight.feats.linecolor = "black",
            highlight.feats.linewidth = 0.5,
            highlight.featsets.color = "#A6A6A6",
            highlight.featsets.size = 7,
            highlight.featsets.opac = 1,
            highlight.featsets.linecolor = "black",
            highlight.featsets.linewidth = 0.5,
            highlight.feats.label = FALSE
        )

        robjects$plot.sgrna2.rank <- fig
        fig
    })
    # nocov end

    # nocov start
    output$sgrna2.detail <- renderDT({
        req(robjects$set2.sgrnas, input$sgrna.gene)

        df <- robjects$set2.sgrnas
        df <- df[df$Gene == input$sgrna.gene, ]

        target <- which(names(df) %in% c(
            "control_mean", "treat_mean",
            "control_var", "adj_var", "high_in_treatment", "p.low",
            "p.high", "p.twosided", "score"
        )) - 1

        datatable(df,
            rownames = FALSE,
            filter = "top",
            extensions = c("Buttons"),
            caption = paste0(input$sgrna.sel2, " ", input$sgrna.gene, " sgRNA Details"),
            options = list(
                pageLength = 10,
                dom = "Blfrtip",
                buttons = c("copy", "csv", "excel", "pdf", "print"),
                columnDefs = list(list(visible = FALSE, targets = target))
            )
        ) %>% formatStyle(0, target = "row", lineHeight = "50%")
    })
    # nocov end
    invisible(NULL)
}
