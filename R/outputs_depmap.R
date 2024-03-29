#' Render general outputs for the DepMap tab
#'
#' Create rendering expressions for the DepMap tab outputs.
#'
#' @param input The Shiny input object from the server function.
#' @param output The Shiny output object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return A \linkS4class{NULL} is invisibly returned
#' and rendering expressions for DepMap tab features are added to \code{output}.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny isolate renderUI
#' @importFrom plotly renderPlotly
#' @importFrom DT renderDT datatable formatStyle
#' @rdname INTERNAL_create_depmap_outputs
.create_depmap_outputs <- function(input, output, robjects) {
    # nocov start
    output$depmap.deplines <- renderUI({
        req(input$depmap.gene, robjects$depmap.gene)
        input$dm.dep.update

        dep.info <- get_depmap_essentiality(input$depmap.gene, robjects$depmap.gene)
        .make_dependency_tag(
            dep.info = dep.info,
            dep.release = robjects$depmap.release,
            crispr.color = isolate(input$dep.crispr.color),
            rnai.color = isolate(input$dep.rnai.color)
        )
    })
    # nocov end

    # nocov start
    output$depmap.essplot <- renderPlotly({
        req(input$depmap.gene, robjects$depmap.meta, robjects$pool)
        input$dm.dep.update

        df <- get_depmap_plot_data(
            gene = input$depmap.gene, data.type = "dependency",
            depmap.meta = robjects$depmap.meta, depmap.pool = robjects$pool
        )

        fig <- plot_depmap_dependency(
            df = df,
            crispr.color = isolate(input$dep.crispr.color),
            rnai.color = isolate(input$dep.rnai.color),
            depline = isolate(input$dep.depline),
            plot.grid = isolate(input$dep.plot.grid)
        )

        robjects$plot.depmap.essplot <- fig
        fig
    })
    # nocov end

    # nocov start
    output$depmap.expplot <- renderPlotly({
        req(input$depmap.gene, robjects$depmap.meta, robjects$pool)
        input$dm.exp.update

        df <- get_depmap_plot_data(
            gene = input$depmap.gene, data.type = "ccle_tpm",
            depmap.meta = robjects$depmap.meta, depmap.pool = robjects$pool
        )

        fig <- plot_depmap_expression(
            df = df,
            color = isolate(input$exp.color),
            plot.grid = isolate(input$exp.plot.grid)
        )

        robjects$plot.depmap.expplot <- fig
        fig
    })
    # nocov end

    # nocov start
    output$depmap.cnplot <- renderPlotly({
        req(input$depmap.gene, robjects$depmap.meta, robjects$pool)
        input$dm.cn.update

        df <- get_depmap_plot_data(
            gene = input$depmap.gene, data.type = "cn",
            depmap.meta = robjects$depmap.meta, depmap.pool = robjects$pool
        )

        fig <- plot_depmap_cn(
            df = df,
            color = isolate(input$cn.color),
            plot.grid = isolate(input$cn.plot.grid)
        )

        robjects$plot.depmap.cnplot <- fig
        fig
    })
    # nocov end

    # nocov start
    output$depmap.lineages <- renderPlotly({
        req(input$depmap.gene, robjects$depmap.meta, robjects$pool)
        input$dm.lineage.update

        plot.val <- isolate(input$lin.data)

        switch(plot.val,
            dependency = {
                plot.val <- "dependency"
            },
            crispr = {
                plot.val <- "dependency"
            },
            rnai = {
                plot.val <- "dependency"
            },
            cn = {
                plot.val <- "log_copy_number"
            },
            ccle_tpm = {
                plot.val <- "rna_expression"
            }
        )

        df <- get_depmap_plot_data(
            gene = input$depmap.gene, data.type = isolate(input$lin.data),
            depmap.meta = robjects$depmap.meta, depmap.pool = robjects$pool
        )

        fig <- plot_depmap_lineages(
            df = df,
            plot.val = plot.val,
            group.by = isolate(input$lin.group),
            label.size = isolate(input$lin.label.size),
            pt.color = isolate(input$lin.pt.color),
            pt.size = isolate(input$lin.pt.size),
            boxplot.fill = isolate(input$lin.box.fill),
            boxplot.line.color = isolate(input$lin.box.color),
            depline = isolate(input$lin.depline)
        )

        robjects$plot.depmap.lineages <- fig
        fig
    })
    # nocov end

    # nocov start
    output$depmap.sublineage <- renderPlotly({
        req(input$depmap.gene, robjects$depmap.meta, robjects$pool)
        input$dm.sublineage.update

        plot.val <- isolate(input$lin.data)

        switch(plot.val,
            dependency = {
                plot.val <- "dependency"
            },
            crispr = {
                plot.val <- "dependency"
            },
            rnai = {
                plot.val <- "dependency"
            },
            cn = {
                plot.val <- "log_copy_number"
            },
            ccle_tpm = {
                plot.val <- "rna_expression"
            }
        )

        df <- get_depmap_plot_data(
            gene = input$depmap.gene, data.type = isolate(input$lin.data),
            depmap.meta = robjects$depmap.meta, depmap.pool = robjects$pool
        )

        fig <- plot_depmap_lineages(
            df = df,
            plot.val = plot.val,
            group.by = "lineage_subtype",
            lineage = isolate(input$sub.lineage),
            label.size = isolate(input$sub.label.size),
            pt.color = isolate(input$sub.pt.color),
            pt.size = isolate(input$sub.pt.size),
            boxplot.fill = isolate(input$sub.box.fill),
            boxplot.line.color = isolate(input$sub.box.color),
            depline = isolate(input$sub.depline)
        )

        robjects$plot.depmap.sublineage <- fig
        fig
    })
    # nocov end

    # nocov start
    output$depmap.geneinfo <- renderUI({
        req(input$depmap.gene, robjects$depmap.gene)

        .make_gene_tag(input$depmap.gene)
    })
    # nocov end
    invisible(NULL)
}
