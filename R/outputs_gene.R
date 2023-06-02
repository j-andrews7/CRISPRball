#' Render general outputs for the Gene tabs
#'
#' Create rendering expressions for the gene tab outputs.
#'
#' @param input The Shiny input object from the server function.
#' @param output The Shiny output object from the server function.
#' @param robjects A reactive list of values generated in the server function.
#'
#' @return A \linkS4class{NULL} is invisibly returned
#' and rendering expressions for Gene tabs features are added to \code{output}.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny isolate selectInput tagList renderUI
#' @importFrom plotly renderPlotly
#' @importFrom DT renderDT datatable formatStyle
#' @rdname INTERNAL_create_gene_outputs
.create_gene_outputs <- function(input, output, robjects) {
    # nocov start
    output$gene.rank.options <- renderUI({
        req(robjects$set1.genes)
        df <- robjects$set1.genes

        tagList(
            selectInput("gene.rankby", "Rank by:", choices = names(df), selected = ifelse("LFC" %in% names(df), "LFC", NULL))
        )
    })
    # nocov end

    # nocov start
    output$gene1.summary <- renderDT(server = FALSE, {
        req(robjects$set1.genes)
        # Remove columns that are redundant or confusing.
        target <- which(names(robjects$set1.genes) %in% c(
            "neg|score", "neg|p-value", "neg|rank",
            "neg|lfc", "pos|score", "pos|p-value", "pos|rank",
            "pos|lfc", "RandomIndex", "Rank", "goodsgrna"
        )) - 1

        df <- robjects$set1.genes

        if (!is.null(robjects$common.hits)) {
            df$Overlap <- df$id %in% robjects$common.hits
        }

        DT::datatable(df,
            rownames = FALSE,
            filter = "top",
            extensions = c("Buttons"),
            caption = paste0(input$gene.sel1, " Gene Summary"),
            options = list(
                search = list(regex = TRUE),
                pageLength = 10,
                dom = "Blfrtip",
                buttons = c("copy", "csv", "excel", "pdf", "print"),
                columnDefs = list(list(visible = FALSE, targets = target))
            )
        ) %>% DT::formatStyle(0, target = "row", lineHeight = "50%")
    })
    # nocov end

    # nocov start
    output$gene1.vol <- renderPlotly({
        req(robjects$set1.genes)
        input$vol.update

        df <- robjects$set1.genes

        hov.info <- c("hit_type", "num", "goodsgrna")

        # Remove common essential genes if needed.
        if (isolate(input$rem.ess) & !is.null(df$essential)) {
            df <- df[!df$essential, ]
        }

        # Remove positive control genes if needed.
        if (isolate(input$rem.pos) & !is.null(df$Positive_Control)) {
            df <- df[!df$Positive_Control, ]
        }

        # Remove DepMap stuff if requested.
        if (!is.null(robjects$depmap.gene)) {
            if (isolate(input$dep.crispr.ess)) {
                df <- df[!df$DepMap_CRISPR_Essential, ]
            }

            if (isolate(input$dep.crispr.sel)) {
                df <- df[!df$DepMap_CRISPR_Selective, ]
            }

            if (isolate(input$dep.rnai.ess)) {
                df <- df[!df$DepMap_RNAi_Essential, ]
            }

            if (isolate(input$dep.rnai.sel)) {
                df <- df[!df$DepMap_RNAi_Selective, ]
            }
        }

        highlight <- NULL
        if (!is.null(isolate(input$hl.genes)) & isolate(input$hl.genes) != "") {
            highlight.feats <- strsplit(isolate(input$hl.genes), ",|\\s|,\\s")[[1]]
            highlight <- highlight.feats[highlight.feats != ""]
        }

        # Add common hits to highlight.
        if (isolate(input$highlight.common)) {
            highlight <- unique(c(robjects$common.hits, highlight))
        }


        fig <- plot_volcano(
            res = df,
            xlim = isolate(input$vol.x),
            ylim = isolate(input$vol.y),
            fc.thresh = isolate(input$gene.lfc.th),
            fc.lines = isolate(input$vol.fcline),
            sig.thresh = isolate(input$gene.fdr.th),
            sig.line = isolate(input$vol.sigline),
            h.id = robjects$h.id,
            h.id.suffix = "_volc1",
            sig.term = "FDR",
            lfc.term = "LFC",
            feat.term = "id",
            hover.info = hov.info,
            fs = robjects$clicked.volc1,
            up.color = isolate(input$up.color),
            down.color = isolate(input$down.color),
            insig.color = isolate(input$insig.color),
            sig.opacity = isolate(input$sig.opa),
            insig.opacity = isolate(input$insig.opa),
            sig.size = isolate(input$sig.size),
            insig.size = isolate(input$insig.size),
            label.size = isolate(input$lab.size),
            webgl = isolate(input$webgl),
            webgl.ratio = isolate(input$webgl.ratio),
            show.counts = isolate(input$counts),
            show.hl.counts = isolate(input$hl.counts),
            counts.size = isolate(input$counts.size),
            highlight.featsets = isolate(input$hl.genesets),
            highlight.feats = highlight,
            featsets = robjects$genesets,
            highlight.feats.color = isolate(input$hl.genes.col),
            highlight.feats.size = isolate(input$hl.genes.size),
            highlight.feats.opac = isolate(input$hl.genes.opa),
            highlight.feats.linecolor = isolate(input$hl.genes.lcol),
            highlight.feats.linewidth = isolate(input$hl.genes.lw),
            highlight.feats.label = isolate(input$hl.genes.label),
            highlight.featsets.color = isolate(input$hl.genesets.col),
            highlight.featsets.size = isolate(input$hl.genesets.size),
            highlight.featsets.opac = isolate(input$hl.genesets.opa),
            highlight.featsets.linecolor = isolate(input$hl.genesets.lcol),
            highlight.featsets.linewidth = isolate(input$hl.genesets.lw),
            highlight.featsets.label = isolate(input$hl.genesets.label)
        )

        robjects$plot.gene1.vol <- fig

        if (robjects$click.obs$volc1$.suspended) {
            robjects$click.obs$volc1$resume()
        }

        if (robjects$dclick.obs$volc1$.suspended) {
            robjects$dclick.obs$volc1$resume()
        }

        fig
    })
    # nocov end

    # nocov start
    output$gene1.rank <- renderPlotly({
        req(robjects$set1.genes, input$gene.rankby)
        input$rank.update

        df <- robjects$set1.genes

        hov.info <- c("hit_type", "num", "goodsgrna")

        # Remove common essential genes if needed.
        if (isolate(input$rem.ess) & !is.null(df$essential)) {
            df <- df[!df$essential, ]
        }

        # Remove positive control genes if needed.
        if (isolate(input$rem.pos) & !is.null(df$Positive_Control)) {
            df <- df[!df$Positive_Control, ]
        }

        # Remove DepMap stuff if requested.
        if (!is.null(robjects$depmap.gene)) {
            if (isolate(input$dep.crispr.ess)) {
                df <- df[!df$DepMap_CRISPR_Essential, ]
            }

            if (isolate(input$dep.crispr.sel)) {
                df <- df[!df$DepMap_CRISPR_Selective, ]
            }

            if (isolate(input$dep.rnai.ess)) {
                df <- df[!df$DepMap_RNAi_Essential, ]
            }

            if (isolate(input$dep.rnai.sel)) {
                df <- df[!df$DepMap_RNAi_Selective, ]
            }
        }

        highlight <- NULL
        if (!is.null(isolate(input$hl.genes)) & isolate(input$hl.genes) != "") {
            highlight.feats <- strsplit(isolate(input$hl.genes), ",|\\s|,\\s")[[1]]
            highlight <- highlight.feats[highlight.feats != ""]
        }

        # Add common hits to highlight.
        if (isolate(input$highlight.common)) {
            highlight <- unique(c(robjects$common.hits, highlight))
        }

        fig <- plot_rank(
            res = df,
            ylim = list(isolate(input$rank.y.min), isolate(input$rank.y.max)),
            y.thresh = isolate(input$gene.lfc.th),
            y.lines = isolate(input$rank.fcline),
            sig.thresh = isolate(input$gene.fdr.th),
            h.id = robjects$h.id,
            h.id.suffix = "_rank1",
            sig.term = "FDR",
            rank.term = isolate(input$gene.rankby),
            rank.ascending = isolate(input$gene.rank.ascending),
            feat.term = "id",
            hover.info = c("hit_type", "goodsgrna"),
            fs = robjects$clicked.rank1,
            up.color = isolate(input$up.color),
            down.color = isolate(input$down.color),
            insig.color = isolate(input$insig.color),
            sig.opacity = isolate(input$sig.opa),
            insig.opacity = isolate(input$insig.opa),
            sig.size = isolate(input$sig.size),
            insig.size = isolate(input$insig.size),
            label.size = isolate(input$lab.size),
            webgl = isolate(input$webgl),
            webgl.ratio = isolate(input$webgl.ratio),
            show.counts = isolate(input$counts),
            show.hl.counts = isolate(input$hl.counts),
            counts.size = isolate(input$counts.size),
            highlight.featsets = isolate(input$hl.genesets),
            highlight.feats = highlight,
            featsets = robjects$genesets,
            highlight.feats.color = isolate(input$hl.genes.col),
            highlight.feats.size = isolate(input$hl.genes.size),
            highlight.feats.opac = isolate(input$hl.genes.opa),
            highlight.feats.linecolor = isolate(input$hl.genes.lcol),
            highlight.feats.linewidth = isolate(input$hl.genes.lw),
            highlight.feats.label = isolate(input$hl.genes.label),
            highlight.featsets.color = isolate(input$hl.genesets.col),
            highlight.featsets.size = isolate(input$hl.genesets.size),
            highlight.featsets.opac = isolate(input$hl.genesets.opa),
            highlight.featsets.linecolor = isolate(input$hl.genesets.lcol),
            highlight.featsets.linewidth = isolate(input$hl.genesets.lw),
            highlight.featsets.label = isolate(input$hl.genesets.label)
        )

        robjects$plot.gene1.rank <- fig

        if (robjects$click.obs$rank1$.suspended) {
            robjects$click.obs$rank1$resume()
        }

        if (robjects$dclick.obs$rank1$.suspended) {
            robjects$dclick.obs$rank1$resume()
        }

        fig
    })
    # nocov end

    # nocov start
    output$gene1.lawn <- renderPlotly({
        req(robjects$set1.genes)
        df <- robjects$set1.genes
        input$lawn.update

        hov.info <- c("hit_type", "num", "goodsgrna")

        # Remove common essential genes if needed.
        if (isolate(input$rem.ess) & !is.null(df$essential)) {
            df <- df[!df$essential, ]
        }

        # Remove positive control genes if needed.
        if (isolate(input$rem.pos) & !is.null(df$Positive_Control)) {
            df <- df[!df$Positive_Control, ]
        }

        # Remove DepMap stuff if requested.
        if (!is.null(robjects$depmap.gene)) {
            if (isolate(input$dep.crispr.ess)) {
                df <- df[!df$DepMap_CRISPR_Essential, ]
            }

            if (isolate(input$dep.crispr.sel)) {
                df <- df[!df$DepMap_CRISPR_Selective, ]
            }

            if (isolate(input$dep.rnai.ess)) {
                df <- df[!df$DepMap_RNAi_Essential, ]
            }

            if (isolate(input$dep.rnai.sel)) {
                df <- df[!df$DepMap_RNAi_Selective, ]
            }
        }

        highlight <- NULL
        if (!is.null(isolate(input$hl.genes)) & isolate(input$hl.genes) != "") {
            highlight.feats <- strsplit(isolate(input$hl.genes), ",|\\s|,\\s")[[1]]
            highlight <- highlight.feats[highlight.feats != ""]
        }

        # Add common hits to highlight.
        if (isolate(input$highlight.common)) {
            highlight <- unique(c(robjects$common.hits, highlight))
        }

        fig <- plot_lawn(
            res = df,
            ylim = isolate(input$lawn.y),
            fc.thresh = isolate(input$gene.lfc.th),
            sig.thresh = isolate(input$gene.fdr.th),
            sig.line = isolate(input$lawn.sigline),
            h.id = robjects$h.id,
            h.id.suffix = "_lawn1",
            sig.term = "FDR",
            lfc.term = "LFC",
            feat.term = "id",
            x.term = "RandomIndex",
            hover.info = hov.info,
            fs = robjects$clicked.lawn1,
            up.color = isolate(input$up.color),
            down.color = isolate(input$down.color),
            insig.color = isolate(input$insig.color),
            sig.opacity = isolate(input$sig.opa),
            insig.opacity = isolate(input$insig.opa),
            sig.size = isolate(input$sig.size),
            insig.size = isolate(input$insig.size),
            label.size = isolate(input$lab.size),
            webgl = isolate(input$webgl),
            webgl.ratio = isolate(input$webgl.ratio),
            show.counts = isolate(input$counts),
            show.hl.counts = isolate(input$hl.counts),
            counts.size = isolate(input$counts.size),
            highlight.featsets = isolate(input$hl.genesets),
            highlight.feats = highlight,
            featsets = robjects$genesets,
            highlight.feats.color = isolate(input$hl.genes.col),
            highlight.feats.size = isolate(input$hl.genes.size),
            highlight.feats.opac = isolate(input$hl.genes.opa),
            highlight.feats.linecolor = isolate(input$hl.genes.lcol),
            highlight.feats.linewidth = isolate(input$hl.genes.lw),
            highlight.feats.label = isolate(input$hl.genes.label),
            highlight.featsets.color = isolate(input$hl.genesets.col),
            highlight.featsets.size = isolate(input$hl.genesets.size),
            highlight.featsets.opac = isolate(input$hl.genesets.opa),
            highlight.featsets.linecolor = isolate(input$hl.genesets.lcol),
            highlight.featsets.linewidth = isolate(input$hl.genesets.lw),
            highlight.featsets.label = isolate(input$hl.genesets.label)
        )

        robjects$plot.gene1.lawn <- fig

        if (robjects$click.obs$lawn1$.suspended) {
            robjects$click.obs$lawn1$resume()
        }

        if (robjects$dclick.obs$lawn1$.suspended) {
            robjects$dclick.obs$lawn1$resume()
        }

        fig
    })
    # nocov end

    # nocov start
    output$gene2.summary <- renderDT(server = FALSE, {
        req(robjects$set2.genes)
        df <- robjects$set2.genes

        # Remove columns that are redundant or confusing.
        target <- which(names(robjects$set2.genes) %in% c(
            "neg|score", "neg|p-value", "neg|rank",
            "neg|lfc", "pos|score", "pos|p-value", "pos|rank",
            "pos|lfc", "RandomIndex", "Rank", "goodsgrna"
        )) - 1

        # Label overlapping hits between datasets if available.
        if (!is.null(robjects$common.hits)) {
            df$Overlap <- df$id %in% robjects$common.hits
        }

        DT::datatable(df,
            rownames = FALSE,
            filter = "top",
            extensions = c("Buttons"),
            caption = paste0(input$gene.sel2, " Gene Summary"),
            options = list(
                search = list(regex = TRUE),
                dom = "Blfrtip",
                buttons = c("copy", "csv", "excel", "pdf", "print"),
                pageLength = 10,
                columnDefs = list(list(visible = FALSE, targets = target))
            )
        ) %>% DT::formatStyle(0, target = "row", lineHeight = "50%")
    })
    # nocov end

    # nocov start
    output$gene2.vol <- renderPlotly({
        req(robjects$set2.genes)
        input$vol.update

        df <- robjects$set2.genes

        hov.info <- c("hit_type", "num", "goodsgrna")

        # Remove common essential genes if needed.
        if (isolate(input$rem.ess) & !is.null(df$essential)) {
            df <- df[!df$essential, ]
        }

        # Remove positive control genes if needed.
        if (isolate(input$rem.pos) & !is.null(df$Positive_Control)) {
            df <- df[!df$Positive_Control, ]
        }

        # Remove DepMap stuff if requested.
        if (!is.null(robjects$depmap.gene)) {
            if (isolate(input$dep.crispr.ess)) {
                df <- df[!df$DepMap_CRISPR_Essential, ]
            }

            if (isolate(input$dep.crispr.sel)) {
                df <- df[!df$DepMap_CRISPR_Selective, ]
            }

            if (isolate(input$dep.rnai.ess)) {
                df <- df[!df$DepMap_RNAi_Essential, ]
            }

            if (isolate(input$dep.rnai.sel)) {
                df <- df[!df$DepMap_RNAi_Selective, ]
            }
        }

        highlight <- NULL
        if (!is.null(isolate(input$hl.genes)) & isolate(input$hl.genes) != "") {
            highlight.feats <- strsplit(isolate(input$hl.genes), ",|\\s|,\\s")[[1]]
            highlight <- highlight.feats[highlight.feats != ""]
        }

        # Add common hits to highlight.
        if (isolate(input$highlight.common)) {
            highlight <- unique(c(robjects$common.hits, highlight))
        }

        fig <- plot_volcano(
            res = df,
            xlim = isolate(input$vol.x),
            ylim = isolate(input$vol.y),
            fc.thresh = isolate(input$gene.lfc.th),
            fc.lines = isolate(input$vol.fcline),
            sig.thresh = isolate(input$gene.fdr.th),
            sig.line = isolate(input$vol.sigline),
            h.id = robjects$h.id,
            h.id.suffix = "_volc2",
            sig.term = "FDR",
            lfc.term = "LFC",
            feat.term = "id",
            hover.info = hov.info,
            fs = robjects$clicked.volc2,
            up.color = isolate(input$up.color),
            down.color = isolate(input$down.color),
            insig.color = isolate(input$insig.color),
            sig.opacity = isolate(input$sig.opa),
            insig.opacity = isolate(input$insig.opa),
            sig.size = isolate(input$sig.size),
            insig.size = isolate(input$insig.size),
            label.size = isolate(input$lab.size),
            webgl = isolate(input$webgl),
            webgl.ratio = isolate(input$webgl.ratio),
            show.counts = isolate(input$counts),
            show.hl.counts = isolate(input$hl.counts),
            counts.size = isolate(input$counts.size),
            highlight.featsets = isolate(input$hl.genesets),
            highlight.feats = highlight,
            featsets = robjects$genesets,
            highlight.feats.color = isolate(input$hl.genes.col),
            highlight.feats.size = isolate(input$hl.genes.size),
            highlight.feats.opac = isolate(input$hl.genes.opa),
            highlight.feats.linecolor = isolate(input$hl.genes.lcol),
            highlight.feats.linewidth = isolate(input$hl.genes.lw),
            highlight.feats.label = isolate(input$hl.genes.label),
            highlight.featsets.color = isolate(input$hl.genesets.col),
            highlight.featsets.size = isolate(input$hl.genesets.size),
            highlight.featsets.opac = isolate(input$hl.genesets.opa),
            highlight.featsets.linecolor = isolate(input$hl.genesets.lcol),
            highlight.featsets.linewidth = isolate(input$hl.genesets.lw),
            highlight.featsets.label = isolate(input$hl.genesets.label)
        )

        robjects$plot.gene2.vol <- fig

        # Resume click observers if needed. Necessary to prevent warnings on load.
        if (robjects$click.obs$volc2$.suspended) {
            robjects$click.obs$volc2$resume()
        }

        if (robjects$dclick.obs$volc2$.suspended) {
            robjects$dclick.obs$volc2$resume()
        }

        fig
    })
    # nocov end

    # nocov start
    output$gene2.rank <- renderPlotly({
        req(robjects$set2.genes, input$gene.rankby)
        input$rank.update

        hov.info <- c("hit_type", "num", "goodsgrna")

        df <- robjects$set2.genes

        # Remove common essential genes if needed.
        if (isolate(input$rem.ess) & !is.null(df$essential)) {
            df <- df[!df$essential, ]
        }

        # Remove positive control genes if needed.
        if (isolate(input$rem.pos) & !is.null(df$Positive_Control)) {
            df <- df[!df$Positive_Control, ]
        }

        # Remove DepMap stuff if requested.
        if (!is.null(robjects$depmap.gene)) {
            if (isolate(input$dep.crispr.ess)) {
                df <- df[!df$DepMap_CRISPR_Essential, ]
            }

            if (isolate(input$dep.crispr.sel)) {
                df <- df[!df$DepMap_CRISPR_Selective, ]
            }

            if (isolate(input$dep.rnai.ess)) {
                df <- df[!df$DepMap_RNAi_Essential, ]
            }

            if (isolate(input$dep.rnai.sel)) {
                df <- df[!df$DepMap_RNAi_Selective, ]
            }
        }

        highlight <- NULL
        if (!is.null(isolate(input$hl.genes)) & isolate(input$hl.genes) != "") {
            highlight.feats <- strsplit(isolate(input$hl.genes), ",|\\s|,\\s")[[1]]
            highlight <- highlight.feats[highlight.feats != ""]
        }

        # Add common hits to highlight.
        if (isolate(input$highlight.common)) {
            highlight <- unique(c(robjects$common.hits, highlight))
        }

        fig <- plot_rank(
            res = df,
            ylim = list(isolate(input$rank.y.min), isolate(input$rank.y.max)),
            y.thresh = isolate(input$gene.lfc.th),
            y.lines = isolate(input$rank.fcline),
            sig.thresh = isolate(input$gene.fdr.th),
            h.id = robjects$h.id,
            h.id.suffix = "_rank2",
            sig.term = "FDR",
            rank.term = isolate(input$gene.rankby),
            rank.ascending = isolate(input$gene.rank.ascending),
            feat.term = "id",
            hover.info = hov.info,
            fs = robjects$clicked.rank2,
            up.color = isolate(input$up.color),
            down.color = isolate(input$down.color),
            insig.color = isolate(input$insig.color),
            sig.opacity = isolate(input$sig.opa),
            insig.opacity = isolate(input$insig.opa),
            sig.size = isolate(input$sig.size),
            insig.size = isolate(input$insig.size),
            label.size = isolate(input$lab.size),
            webgl = isolate(input$webgl),
            webgl.ratio = isolate(input$webgl.ratio),
            show.counts = isolate(input$counts),
            show.hl.counts = isolate(input$hl.counts),
            counts.size = isolate(input$counts.size),
            highlight.featsets = isolate(input$hl.genesets),
            highlight.feats = highlight,
            featsets = robjects$genesets,
            highlight.feats.color = isolate(input$hl.genes.col),
            highlight.feats.size = isolate(input$hl.genes.size),
            highlight.feats.opac = isolate(input$hl.genes.opa),
            highlight.feats.linecolor = isolate(input$hl.genes.lcol),
            highlight.feats.linewidth = isolate(input$hl.genes.lw),
            highlight.feats.label = isolate(input$hl.genes.label),
            highlight.featsets.color = isolate(input$hl.genesets.col),
            highlight.featsets.size = isolate(input$hl.genesets.size),
            highlight.featsets.opac = isolate(input$hl.genesets.opa),
            highlight.featsets.linecolor = isolate(input$hl.genesets.lcol),
            highlight.featsets.linewidth = isolate(input$hl.genesets.lw),
            highlight.featsets.label = isolate(input$hl.genesets.label)
        )

        robjects$plot.gene2.rank <- fig

        if (robjects$click.obs$rank2$.suspended) {
            robjects$click.obs$rank2$resume()
        }

        if (robjects$dclick.obs$rank2$.suspended) {
            robjects$dclick.obs$rank2$resume()
        }

        fig
    })
    # nocov end

    # nocov start
    output$gene2.lawn <- renderPlotly({
        req(robjects$set2.genes)
        input$lawn.update

        hov.info <- c("hit_type", "num", "goodsgrna")

        df <- robjects$set2.genes

        # Remove common essential genes if needed.
        if (isolate(input$rem.ess) & !is.null(df$essential)) {
            df <- df[!df$essential, ]
        }

        # Remove positive control genes if needed.
        if (isolate(input$rem.pos) & !is.null(df$Positive_Control)) {
            df <- df[!df$Positive_Control, ]
        }

        # Remove DepMap stuff if requested.
        if (!is.null(robjects$depmap.gene)) {
            if (isolate(input$dep.crispr.ess)) {
                df <- df[!df$DepMap_CRISPR_Essential, ]
            }

            if (isolate(input$dep.crispr.sel)) {
                df <- df[!df$DepMap_CRISPR_Selective, ]
            }

            if (isolate(input$dep.rnai.ess)) {
                df <- df[!df$DepMap_RNAi_Essential, ]
            }

            if (isolate(input$dep.rnai.sel)) {
                df <- df[!df$DepMap_RNAi_Selective, ]
            }
        }

        highlight <- NULL
        if (!is.null(isolate(input$hl.genes)) & isolate(input$hl.genes) != "") {
            highlight.feats <- strsplit(isolate(input$hl.genes), ",|\\s|,\\s")[[1]]
            highlight <- highlight.feats[highlight.feats != ""]
        }

        # Add common hits to highlight.
        if (isolate(input$highlight.common)) {
            highlight <- unique(c(robjects$common.hits, highlight))
        }

        fig <- plot_lawn(
            res = df,
            ylim = isolate(input$lawn.y),
            fc.thresh = isolate(input$gene.lfc.th),
            sig.thresh = isolate(input$gene.fdr.th),
            sig.line = isolate(input$lawn.sigline),
            h.id = robjects$h.id,
            h.id.suffix = "_lawn2",
            sig.term = "FDR",
            lfc.term = "LFC",
            feat.term = "id",
            x.term = "RandomIndex",
            hover.info = hov.info,
            fs = robjects$clicked.lawn2,
            up.color = isolate(input$up.color),
            down.color = isolate(input$down.color),
            insig.color = isolate(input$insig.color),
            sig.opacity = isolate(input$sig.opa),
            insig.opacity = isolate(input$insig.opa),
            sig.size = isolate(input$sig.size),
            insig.size = isolate(input$insig.size),
            label.size = isolate(input$lab.size),
            webgl = isolate(input$webgl),
            webgl.ratio = isolate(input$webgl.ratio),
            show.counts = isolate(input$counts),
            show.hl.counts = isolate(input$hl.counts),
            counts.size = isolate(input$counts.size),
            highlight.featsets = isolate(input$hl.genesets),
            highlight.feats = highlight,
            featsets = robjects$genesets,
            highlight.feats.color = isolate(input$hl.genes.col),
            highlight.feats.size = isolate(input$hl.genes.size),
            highlight.feats.opac = isolate(input$hl.genes.opa),
            highlight.feats.linecolor = isolate(input$hl.genes.lcol),
            highlight.feats.linewidth = isolate(input$hl.genes.lw),
            highlight.feats.label = isolate(input$hl.genes.label),
            highlight.featsets.color = isolate(input$hl.genesets.col),
            highlight.featsets.size = isolate(input$hl.genesets.size),
            highlight.featsets.opac = isolate(input$hl.genesets.opa),
            highlight.featsets.linecolor = isolate(input$hl.genesets.lcol),
            highlight.featsets.linewidth = isolate(input$hl.genesets.lw),
            highlight.featsets.label = isolate(input$hl.genesets.label)
        )

        robjects$plot.gene2.lawn <- fig

        if (robjects$click.obs$lawn2$.suspended) {
            robjects$click.obs$lawn2$resume()
        }

        if (robjects$dclick.obs$lawn2$.suspended) {
            robjects$dclick.obs$lawn2$resume()
        }

        fig
    })
    # nocov end
    invisible(NULL)
}
