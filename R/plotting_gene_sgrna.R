# For adding various lines to plots easily.
.vline <- function(x = 0, color = "red", width = 1, dash = "solid") {
    list(
        type = "line",
        y0 = 0,
        y1 = 1,
        yref = "paper",
        x0 = x,
        x1 = x,
        line = list(color = color, width = width, dash = dash)
    )
}

.hline <- function(y = 0, color = "blue", width = 1, dash = "solid") {
    list(
        type = "line",
        x0 = 0,
        x1 = 1,
        xref = "paper",
        y0 = y,
        y1 = y,
        line = list(color = color, width = width, dash = dash)
    )
}

.fitline <- function(df, color = "black", width = 0.75, dash = "solid") {
    list(
        type = "line",
        line = list(color = color, width = width, dash = dash),
        xref = "x",
        yref = "y",
        y0 = min(df$fv),
        y1 = max(df$fv),
        x0 = df$lfc.x[df$fv == min(df$fv)],
        x1 = df$lfc.x[df$fv == max(df$fv)]
    )
}


#' Create an interactive volcano plot
#'
#' Create an interactive volcano plot for data with fold change and significance terms.
#'
#' @param res Dataframe containing, at minimum, fold change and significance values.
#' @param xlim Positive numeric scalar indicating x-axis limits.
#'   The negative value will be used for the lower limit.
#' @param ylim Positive numeric scalar indicating y-axis limits.
#'   The negative value will be used for the lower limit.
#' @param fc.thresh Numeric scalar indicating the fold change threshold for coloring
#'   significant features.
#' @param fc.lines Logical indicating whether to add fold change threshold lines to the plot.
#' @param hover.info Character vector indicating which additional columns from \code{res}
#'   to include in the hover info.
#' @param sig.line Logical indicating whether to add a significance threshold line to the plot.
#' @param h.id Character scalar indicating the unique ID of the plotly object. Can usually be ignored,
#'   but should be used if multiple plots are being created in the same R session (e.g. Shiny app).
#' @param feat.term Character scalar indicating the column name of the feature IDs in \code{res}.
#' @param sig.term Character scalar indicating the column name of the significance values in \code{res}.
#' @param lfc.term Character scalar indicating the column name of the log fold change values in \code{res}.
#' @param down.color Character scalar indicating the color of down-regulated features.
#' @param up.color Character scalar indicating the color of up-regulated features.
#' @param insig.color Character scalar indicating the color of insignificant features.
#' @param sig.thresh Numeric scalar indicating the significance threshold for coloring significant features.
#' @param fs Dataframe containing coordinates and label information for points that should be labeled.
#'   Columns should be:
#'   \itemize{
#'    \item \code{x} - x coordinate of the point
#'    \item \code{y} - y coordinate of the point
#'    \item \code{customdata} - label to be displayed
#'   }
#' @param sig.size Numeric scalar indicating the size of significant feature points.
#' @param insig.size Numeric scalar indicating the size of insignificant feature points.
#' @param sig.opacity Numeric scalar indicating the opacity of significant feature points.
#' @param insig.opacity Numeric scalar indicating the opacity of insignificant feature points.
#' @param label.size Numeric scalar indicating the size of feature labels.
#' @param webgl Logical indicating whether to use WebGL for rendering the plot.
#' @param webgl.ratio Numeric scalar indicating the ratio of WebGL to HTML5 canvas rendering, increases resolution
#'   of saved plot when WebGL plotting is not used.
#' @param show.counts Logical indicating whether to show annotations for the number of features in the plot.
#' @param show.hl.counts Logical indicating whether to show annotations for the number of highlighted features in the plot.
#' @param counts.size Numeric scalar indicating the size of the feature counts labels.
#' @param highlight.featsets Character vector indicating which feature sets should be highlighted.
#' @param highlight.feats Character vector indicating which features should be highlighted.
#' @param featsets Named list of feature sets to be used for highlighting.
#' @param highlight.feats.color Character scalar indicating the color of highlighted features.
#' @param highlight.feats.size Numeric scalar indicating the size of highlighted features.
#' @param highlight.feats.opac Numeric scalar indicating the opacity of highlighted features.
#' @param highlight.feats.label Logical indicating whether to label highlighted features.
#' @param highlight.feats.linewidth Numeric scalar indicating the line width of highlighted features.
#' @param highlight.feats.linecolor Character scalar indicating the line color of highlighted features.
#' @param highlight.featsets.color Character scalar indicating the color of highlighted feature sets.
#' @param highlight.featsets.size Numeric scalar indicating the point size of highlighted feature sets.
#' @param highlight.featsets.opac Numeric scalar indicating the opacity of highlighted feature sets.
#' @param highlight.featsets.linewidth Numeric scalar indicating the line width of highlighted feature sets.
#' @param highlight.featsets.linecolor Character scalar indicating the line color of highlighted feature sets.
#' @param highlight.featsets.label Logical indicating whether to label highlighted feature sets.
#' @param h.id.suffix Character scalar indicating the suffix to be added to the plotly object ID.
#'
#' @return An interactive plotly volcano plot.
#'
#' @importFrom plotly plot_ly toWebGL layout config add_annotations
#'
#' @author Jared Andrews
#' @export
#' @examples
#' library(CRISPRball)
#' d1.genes <- read.delim(system.file("extdata", "esc1.gene_summary.txt",
#'     package = "CRISPRball"
#' ), check.names = FALSE)
#' plot.df <- gene_ingress(d1.genes,
#'     sig.thresh = 0.05, es.thresh = 0.5,
#'     es.col = "LFC", sig.col = "FDR"
#' )
#' plot_volcano(plot.df, feat.term = "id")
plot_volcano <- function(res,
                         xlim = 5,
                         ylim = 5,
                         fc.thresh = 0.5,
                         fc.lines = TRUE,
                         hover.info = NULL,
                         sig.line = TRUE,
                         h.id = "crispr",
                         feat.term = "rows",
                         sig.term = "FDR",
                         lfc.term = "LFC",
                         down.color = "#0026ff",
                         up.color = "#ff0000",
                         insig.color = "#A6A6A6",
                         sig.thresh = 0.05,
                         fs = NULL,
                         sig.size = 6,
                         insig.size = 5,
                         sig.opacity = 1,
                         insig.opacity = 0.5,
                         label.size = 10,
                         webgl = TRUE,
                         webgl.ratio = 7,
                         show.counts = TRUE,
                         show.hl.counts = TRUE,
                         counts.size = 8,
                         highlight.featsets = NULL,
                         highlight.feats = NULL,
                         featsets = NULL,
                         highlight.feats.color = "#E69F00",
                         highlight.feats.size = 7,
                         highlight.feats.opac = 1,
                         highlight.feats.linecolor = "#000000",
                         highlight.feats.linewidth = 1,
                         highlight.feats.label = TRUE,
                         highlight.featsets.color = "#009E73",
                         highlight.featsets.size = 7,
                         highlight.featsets.opac = 1,
                         highlight.featsets.linecolor = "#000000",
                         highlight.featsets.linewidth = 1,
                         highlight.featsets.label = FALSE,
                         h.id.suffix = "_volc") {
    # Check for required columns.
    if (!feat.term %in% colnames(res) & feat.term != "rows") {
        fig <- .empty_plot(paste0("Column '", feat.term, "' not found in dataframe."), plotly = TRUE)
    } else if (!sig.term %in% colnames(res)) {
        fig <- .empty_plot(paste0("Column '", sig.term, "' not found in dataframe."), plotly = TRUE)
    } else if (!lfc.term %in% colnames(res)) {
        fig <- .empty_plot(paste0("Column '", lfc.term, "' not found in dataframe."), plotly = TRUE)
    } else {
        # Styling.
        res$col <- rep(insig.color, nrow(res))
        res$cex <- rep(insig.size, nrow(res))
        res$order <- rep(0, nrow(res))
        res$lcol <- res$col
        res$lw <- 0
        res$opacity <- insig.opacity

        # Remove features with NA significance term (due to low expression, etc).
        res <- res[!is.na(res[[sig.term]]), ]

        # Get all gene IDs or symbols to be highlighted.
        highlight <- highlight.feats

        highlight.fs <- highlight.featsets
        if (!is.null(highlight.featsets)) {
            for (featset in highlight.featsets) {
                highlight.fs <- c(highlight.fs, featsets[[featset]])
            }
        }

        # Significance filter.
        up.feats <- res[[sig.term]] < sig.thresh & res[[lfc.term]] > fc.thresh
        res$col[up.feats] <- up.color
        res$cex[up.feats] <- sig.size
        res$order[up.feats] <- 1
        res$opacity[up.feats] <- sig.opacity

        dn.feats <- res[[sig.term]] < sig.thresh & res[[lfc.term]] < -fc.thresh
        res$col[dn.feats] <- down.color
        res$cex[dn.feats] <- sig.size
        res$order[dn.feats] <- 1
        res$opacity[dn.feats] <- sig.opacity

        res$x <- res[[lfc.term]]
        res$y <- -log10(res[[sig.term]])

        res$col[res$y < -log10(sig.thresh)] <- insig.color

        res$sh <- ifelse(res$y > ylim, "triangle-up-open",
            ifelse(res$x < -xlim, "triangle-left-open",
                ifelse(res$x > xlim, "triangle-right-open", 0)
            )
        )

        res$lw <- ifelse(res$sh != 0, 1, 0)

        res$y[res$y > ylim] <- ylim - 0.2
        res$x[res$x > xlim] <- xlim - 0.05
        res$x[res$x < -xlim] <- -xlim + 0.05
        if (feat.term == "rows") {
            res$feat <- rownames(res)
        } else {
            res$feat <- res[[feat.term]]
        }

        # Gene/geneset highlighting & labeling.
        n.fs.hl <- 0
        n.hl <- 0

        if (!is.null(highlight.fs)) {
            highlight.fs <- highlight.fs[highlight.fs %in% res$feat]
            n.fs.hl <- length(res$col[res$feat %in% highlight.fs])

            res$col[res$feat %in% highlight.fs] <- highlight.featsets.color
            res$cex[res$feat %in% highlight.fs] <- highlight.featsets.size
            res$opacity[res$feat %in% highlight.fs] <- highlight.featsets.opac
            res$lcol[res$feat %in% highlight.fs] <- highlight.featsets.linecolor
            res$lw[res$feat %in% highlight.fs] <- highlight.featsets.linewidth
            res$order[res$feat %in% highlight.fs] <- 2

            if (highlight.featsets.label) {
                add.fs <- data.frame(
                    x = res$x[res$feat %in% highlight.fs],
                    y = res$y[res$feat %in% highlight.fs],
                    customdata = res$feat[res$feat %in% highlight.fs]
                )

                if (!is.null(fs)) {
                    fs <- fs[, c("x", "y", "customdata")]
                    fs <- rbind(fs, add.fs)
                } else {
                    fs <- add.fs
                }
            }
        }

        # Want these to have precedence over the feature sets in case entries are in both.
        if (!is.null(highlight)) {
            highlight <- highlight[highlight %in% res$feat]
            if (length(highlight) > 0) {
                n.hl <- length(res$col[res$feat %in% highlight])

                res$col[res$feat %in% highlight] <- highlight.feats.color
                res$cex[res$feat %in% highlight] <- highlight.feats.size
                res$opacity[res$feat %in% highlight] <- highlight.feats.opac
                res$lcol[res$feat %in% highlight] <- highlight.feats.linecolor
                res$lw[res$feat %in% highlight] <- highlight.feats.linewidth
                res$order[res$feat %in% highlight] <- 3

                if (highlight.feats.label) {
                    add.fs <- data.frame(
                        x = res$x[res$feat %in% highlight],
                        y = res$y[res$feat %in% highlight],
                        customdata = res$feat[res$feat %in% highlight]
                    )

                    if (!is.null(fs)) {
                        fs <- fs[, c("x", "y", "customdata")]
                        fs <- rbind(fs, add.fs)
                    } else {
                        fs <- add.fs
                    }
                }
            }
        }

        res$hover.string <- paste0(
            "</br><b>", feat.term, ":</b> ", res$feat,
            "</br><b>", lfc.term, ":</b> ", format(round(res[[lfc.term]], 4), nsmall = 4),
            "</br><b>", sig.term, ":</b> ", format(round(res[[sig.term]], 6), nsmall = 6)
        )

        if (!is.null(hover.info)) {
            for (n in hover.info) {
                res$hover.string <- paste0(res$hover.string, "</br><b>", n, ":</b>", res[[n]])
            }
        }

        res <- as.data.frame(res)
        res <- res[order(res$order), ]

        # Get feature numbers.
        n.up.feats <- length(up.feats[up.feats == TRUE])
        n.dn.feats <- length(dn.feats[dn.feats == TRUE])
        n.feats <- nrow(res)

        # Add plot border.
        ay <- list(
            showline = TRUE,
            mirror = TRUE,
            linecolor = toRGB("black"),
            linewidth = 0.5,
            title = paste0("-log10(", sig.term, ")"),
            range = list(0, ylim),
            showgrid = FALSE,
            layer = "below traces",
            zeroline = FALSE,
            ticks = "outside",
            zerolinewidth = 0.5
        )

        ax <- list(
            showline = TRUE,
            mirror = TRUE,
            linecolor = toRGB("black"),
            linewidth = 0.5,
            title = lfc.term,
            range = list(-xlim, xlim),
            showgrid = FALSE,
            layer = "below traces",
            ticks = "outside",
            zerolinewidth = 0.5
        )

        # Create vertical and horizontal lines.
        fc.line1 <- NULL
        fc.line2 <- NULL

        sig.hline <- NULL
        if (sig.line) {
            sig.hline <- .hline(y = -log10(sig.thresh), color = "#000000", width = 1, dash = "longdash")
        }

        if (fc.thresh != 0 & fc.lines) {
            fc.line1 <- .vline(x = fc.thresh, color = "#000000", width = 1, dash = "longdash")
            fc.line2 <- .vline(x = -fc.thresh, color = "#000000", width = 1, dash = "longdash")
        }

        # Figure generation.
        fig <- plot_ly(res,
            x = ~x,
            y = ~y,
            customdata = ~feat,
            type = "scatter",
            mode = "markers",
            marker = list(
                color = ~col,
                size = ~cex,
                symbol = ~sh,
                line = list(color = ~lcol, width = ~lw),
                opacity = ~opacity
            ),
            text = ~hover.string,
            hoverinfo = "text",
            source = paste0(h.id, h.id.suffix)
        ) %>%
            config(
                edits = list(
                    annotationPosition = TRUE,
                    annotationTail = TRUE
                ),
                toImageButtonOptions = list(format = "svg"),
                displaylogo = FALSE,
                plotGlPixelRatio = webgl.ratio
            )

        if (!is.null(fs)) {
            fig <- fig %>%
                layout(
                    xaxis = ax,
                    yaxis = ay,
                    showlegend = FALSE,
                    shapes = list(sig.hline, fc.line1, fc.line2),
                    hoverlabel = list(font = list(size = 10))
                ) %>%
                add_annotations(
                    x = fs$x, y = fs$y, text = fs$customdata,
                    font = list(size = label.size, family = "Arial"), arrowside = "none"
                )
        } else {
            fig <- fig %>% layout(
                xaxis = ax,
                yaxis = ay,
                showlegend = FALSE,
                shapes = list(sig.hline, fc.line1, fc.line2),
                hoverlabel = list(font = list(size = 10))
            )
        }

        # Feature count annotations.
        if (show.counts) {
            fig <- fig %>%
                add_annotations(
                    x = 1,
                    y = 1,
                    xref = "paper",
                    yref = "paper",
                    text = paste0(
                        "Up features: ", n.up.feats,
                        "\nDown features: ", n.dn.feats,
                        "\nTotal features: ", n.feats
                    ),
                    showarrow = FALSE,
                    font = list(size = counts.size)
                )
        }

        if (show.hl.counts) {
            fig <- fig %>%
                add_annotations(
                    x = 0,
                    y = 1,
                    xref = "paper",
                    yref = "paper",
                    text = paste0(
                        "Set features: ", n.fs.hl,
                        "\nHighlighted features: ", n.hl
                    ),
                    showarrow = FALSE,
                    font = list(size = counts.size)
                )
        }

        if (webgl) {
            fig <- fig %>% toWebGL()
        }
    }

    fig
}


#' Create an interactive rank plot
#'
#' Create an interactive rank plot for data with fold change, significance terms, and rank.
#'
#' @inheritParams plot_volcano
#' @param res Dataframe containing, at minimum, significance values and something to rank by (LFC, RRA score, betas, etc).
#' @param ylim Numeric vector of length two for the y-axis limits.
#' @param y.thresh Numeric scalar used as the y-axis threshold for point coloring.
#'   The negative of this value is also used as the threshold.
#' @param y.lines Logical as for whether or not to show horizontal lines at \code{y.thresh}.
#' @param rank.term Character scalar for the term to rank by from \code{res}. This will be used as the y-axis.
#' @param rank.ascending Boolean indicating whether or not the rank should be ascending.
#'
#' @return An interactive plotly rank plot.
#'
#' @importFrom plotly plot_ly toWebGL layout config add_annotations
#'
#' @author Jared Andrews
#' @export
#' @examples
#' library(CRISPRball)
#' d1.genes <- read.delim(system.file("extdata", "esc1.gene_summary.txt",
#'     package = "CRISPRball"
#' ), check.names = FALSE)
#' plot.df <- gene_ingress(d1.genes,
#'     sig.thresh = 0.05, es.thresh = 0.5,
#'     es.col = "LFC", sig.col = "FDR"
#' )
#' plot_rank(plot.df, feat.term = "id")
plot_rank <- function(res,
                      ylim = c(-10, 10),
                      y.thresh = 0.5,
                      y.lines = TRUE,
                      hover.info = NULL,
                      h.id = "crispr",
                      feat.term = "rows",
                      sig.term = "FDR",
                      rank.term = "LFC",
                      rank.ascending = TRUE,
                      down.color = "#0026ff",
                      up.color = "#ff0000",
                      insig.color = "#A6A6A6",
                      sig.thresh = 0.05,
                      fs = NULL,
                      sig.size = 6,
                      insig.size = 5,
                      sig.opacity = 1,
                      insig.opacity = 0.5,
                      label.size = 10,
                      webgl = TRUE,
                      webgl.ratio = 7,
                      show.counts = TRUE,
                      show.hl.counts = TRUE,
                      counts.size = TRUE,
                      highlight.featsets = NULL,
                      highlight.feats = NULL,
                      featsets = NULL,
                      highlight.feats.color = "#E69F00",
                      highlight.feats.size = 7,
                      highlight.feats.opac = 1,
                      highlight.feats.linecolor = "#000000",
                      highlight.feats.linewidth = 1,
                      highlight.feats.label = TRUE,
                      highlight.featsets.color = "#009E73",
                      highlight.featsets.size = 7,
                      highlight.featsets.opac = 1,
                      highlight.featsets.linecolor = "#000000",
                      highlight.featsets.linewidth = 1,
                      highlight.featsets.label = FALSE,
                      h.id.suffix = "_volc") {
    # Check for required columns.
    if (!rank.term %in% colnames(res)) {
        fig <- .empty_plot(paste0("Column '", rank.term, "' not found in dataframe."), plotly = TRUE)
    } else if (!is.null(sig.term) & !sig.term %in% colnames(res)) {
        fig <- .empty_plot(paste0("Column '", sig.term, "' not found in dataframe."), plotly = TRUE)
    } else if (!feat.term %in% colnames(res) & feat.term != "rows") {
        fig <- .empty_plot(paste0("Column '", feat.term, "' not found in dataframe."), plotly = TRUE)
    } else {
        # Styling.
        res$col <- rep(insig.color, nrow(res))
        res$cex <- rep(insig.size, nrow(res))
        res$order <- rep(0, nrow(res))
        res$lcol <- res$col
        res$lw <- 0
        res$opacity <- insig.opacity

        # Remove features with NA significance term (due to low expression, etc).
        if (!is.null(sig.term)) {
            res <- res[!is.na(res[[sig.term]]), ]
        }

        # Get all feature IDs to be highlighted.
        highlight <- highlight.feats

        highlight.fs <- NULL
        if (!is.null(highlight.featsets)) {
            for (featset in highlight.featsets) {
                highlight.fs <- c(highlight.fs, featsets[[featset]])
            }
        }

        # Significance filter, if provided.
        if (!is.null(sig.term) & !is.null(sig.thresh)) {
            up.feats <- res[[sig.term]] < sig.thresh & res[[rank.term]] > y.thresh
            dn.feats <- res[[sig.term]] < sig.thresh & res[[rank.term]] < -y.thresh
        } else {
            up.feats <- res[[rank.term]] > y.thresh
            dn.feats <- res[[rank.term]] < -y.thresh
        }

        res$col[up.feats] <- up.color
        res$cex[up.feats] <- sig.size
        res$order[up.feats] <- 1
        res$opacity[up.feats] <- sig.opacity

        res$col[dn.feats] <- down.color
        res$cex[dn.feats] <- sig.size
        res$order[dn.feats] <- 1
        res$opacity[dn.feats] <- sig.opacity

        # Add rank column by rank.term
        if (rank.ascending) {
            res$x <- rank(res[[rank.term]], ties.method = "first")
        } else {
            res$x <- rank(-res[[rank.term]], ties.method = "first")
        }

        res$y <- res[[rank.term]]

        res$col[res$y < y.thresh & res$y > -y.thresh] <- insig.color

        res$sh <- ifelse(res$y > ylim[[2]], "triangle-up-open",
            ifelse(res$y < ylim[[1]], "triangle-down-open", 0)
        )

        res$lw <- ifelse(res$sh != 0, 1, 0)

        # Get feature identifier for hover/labeling.
        if (feat.term == "rows") {
            res$feat <- rownames(res)
        } else {
            res$feat <- res[[feat.term]]
        }

        # Gene/geneset highlighting.
        n.fs.hl <- 0
        n.hl <- 0

        if (!is.null(highlight.fs)) {
            highlight.fs <- highlight.fs[highlight.fs %in% res$feat]
            n.fs.hl <- length(res$col[res$feat %in% highlight.fs])

            res$col[res$feat %in% highlight.fs] <- highlight.featsets.color
            res$cex[res$feat %in% highlight.fs] <- highlight.featsets.size
            res$opacity[res$feat %in% highlight.fs] <- highlight.featsets.opac
            res$lcol[res$feat %in% highlight.fs] <- highlight.featsets.linecolor
            res$lw[res$feat %in% highlight.fs] <- highlight.featsets.linewidth
            res$order[res$feat %in% highlight.fs] <- 2

            if (highlight.featsets.label) {
                add.fs <- data.frame(
                    x = res$x[res$feat %in% highlight.fs],
                    y = res$y[res$feat %in% highlight.fs],
                    customdata = res$feat[res$feat %in% highlight.fs]
                )

                if (!is.null(fs)) {
                    fs <- fs[, c("x", "y", "customdata")]
                    fs <- rbind(fs, add.fs)
                } else {
                    fs <- add.fs
                }
            }
        }

        # Want these to have precedence over the feature sets in case entries are in both.
        if (!is.null(highlight)) {
            highlight <- highlight[highlight %in% res$feat]
            if (length(highlight) > 0) {
                n.hl <- length(res$col[res$feat %in% highlight])

                res$col[res$feat %in% highlight] <- highlight.feats.color
                res$cex[res$feat %in% highlight] <- highlight.feats.size
                res$opacity[res$feat %in% highlight] <- highlight.feats.opac
                res$lcol[res$feat %in% highlight] <- highlight.feats.linecolor
                res$lw[res$feat %in% highlight] <- highlight.feats.linewidth
                res$order[res$feat %in% highlight] <- 3

                if (highlight.feats.label) {
                    add.fs <- data.frame(
                        x = res$x[res$feat %in% highlight],
                        y = res$y[res$feat %in% highlight],
                        customdata = res$feat[res$feat %in% highlight]
                    )

                    if (!is.null(fs)) {
                        fs <- fs[, c("x", "y", "customdata")]
                        fs <- rbind(fs, add.fs)
                    } else {
                        fs <- add.fs
                    }
                }
            }
        }

        res$hover.string <- paste(
            "</br><b>", feat.term, ":</b> ", res$feat,
            "</br><b>Rank:</b> ", res$x,
            "</br><b>", rank.term, ":</b> ", format(round(res[[rank.term]], 4), nsmall = 4),
            "</br><b>", sig.term, ":</b> ", format(round(res[[sig.term]], 6), nsmall = 6)
        )

        if (!is.null(hover.info)) {
            for (n in hover.info) {
                res$hover.string <- paste0(res$hover.string, "</br><b>", n, ":</b> ", res[[n]])
            }
        }

        res <- as.data.frame(res)
        res <- res[order(res$order), ]

        # Get feature numbers.
        n.up.feats <- length(up.feats[up.feats == TRUE])
        n.dn.feats <- length(dn.feats[dn.feats == TRUE])
        n.feats <- nrow(res)

        # Add plot border.
        ay <- list(
            showline = TRUE,
            mirror = TRUE,
            linecolor = toRGB("black"),
            linewidth = 0.5,
            title = rank.term,
            range = ylim,
            showgrid = FALSE,
            layer = "below traces",
            zeroline = TRUE,
            ticks = "outside",
            zerolinewidth = 0.5
        )

        ax <- list(
            showline = TRUE,
            mirror = TRUE,
            linecolor = toRGB("black"),
            linewidth = 0.5,
            title = "Rank",
            range = list(-(0.03 * nrow(res)), nrow(res) + (0.03 * nrow(res))),
            showgrid = FALSE,
            layer = "below traces",
            ticks = "outside",
            zerolinewidth = 0.5,
            zeroline = FALSE
        )

        # Create horizontal lines.
        y.line1 <- NULL
        y.line2 <- NULL

        if (y.thresh != 0 & y.lines) {
            y.line1 <- .hline(y = y.thresh, color = "#000000", width = 1, dash = "longdash")
            y.line2 <- .hline(y = -y.thresh, color = "#000000", width = 1, dash = "longdash")
        }

        # Figure generation.
        fig <- plot_ly(res,
            x = ~x,
            y = ~y,
            customdata = ~feat,
            type = "scatter",
            mode = "markers",
            marker = list(
                color = ~col,
                size = ~cex,
                symbol = ~sh,
                line = list(color = ~lcol, width = ~lw),
                opacity = ~opacity
            ),
            text = ~hover.string,
            hoverinfo = "text",
            source = paste0(h.id, h.id.suffix)
        ) %>%
            config(
                edits = list(
                    annotationPosition = TRUE,
                    annotationTail = TRUE
                ),
                toImageButtonOptions = list(format = "svg"),
                displaylogo = FALSE,
                plotGlPixelRatio = webgl.ratio
            )

        if (!is.null(fs)) {
            fig <- fig %>%
                layout(
                    xaxis = ax,
                    yaxis = ay,
                    showlegend = FALSE,
                    shapes = list(y.line1, y.line2)
                ) %>%
                add_annotations(
                    x = fs$x,
                    y = fs$y,
                    text = fs$customdata,
                    font = list(size = label.size, family = "Arial"),
                    arrowside = "none"
                )
        } else {
            fig <- fig %>% layout(
                xaxis = ax,
                yaxis = ay,
                showlegend = FALSE,
                shapes = list(y.line1, y.line2)
            )
        }

        # Feature count annotations.
        if (show.counts) {
            fig <- fig %>%
                add_annotations(
                    x = 1,
                    y = 0,
                    xref = "paper",
                    yref = "paper",
                    text = paste0(
                        "Up features: ", n.up.feats,
                        "\nDown features: ", n.dn.feats,
                        "\nTotal features: ", n.feats
                    ),
                    showarrow = FALSE,
                    font = list(size = counts.size)
                )
        }

        if (show.hl.counts) {
            fig <- fig %>%
                add_annotations(
                    x = 0,
                    y = 1,
                    xref = "paper",
                    yref = "paper",
                    text = paste0(
                        "Set features: ", n.fs.hl,
                        "\nHighlighted features: ", n.hl
                    ),
                    showarrow = FALSE,
                    font = list(size = counts.size)
                )
        }

        if (webgl) {
            fig <- fig %>% toWebGL()
        }
    }

    fig
}


#' Create an interactive lawn plot
#'
#' Create an interactive lawn plot for data with significance values.
#' Typically, this plot is randomly ordered along the x-axis,
#' but the user is free to order it by any term in \code{res} that they'd like.
#'
#' @inheritParams plot_volcano
#' @param res Dataframe containing, at minimum, significance values and a term by which the x-axis can be ordered.
#' @param x.term Character scalar for the x-axis term from \code{res} to be plotted.
#'
#' @return An interactive plotly rank plot.
#'
#' @importFrom plotly plot_ly toWebGL layout config add_annotations ggplotly
#'
#' @author Jared Andrews
#' @export
#' @examples
#' library(CRISPRball)
#' d1.genes <- read.delim(system.file("extdata", "esc1.gene_summary.txt",
#'     package = "CRISPRball"
#' ), check.names = FALSE)
#' plot.df <- gene_ingress(d1.genes,
#'     sig.thresh = 0.05, es.thresh = 0.5,
#'     es.col = "LFC", sig.col = "FDR"
#' )
#' plot_lawn(plot.df, feat.term = "id")
plot_lawn <- function(res,
                      ylim = 5,
                      fc.thresh = 0.5,
                      hover.info = NULL,
                      sig.line = TRUE,
                      h.id = "crispr",
                      feat.term = "rows",
                      x.term = "RandomIndex",
                      sig.term = "FDR",
                      lfc.term = "LFC",
                      down.color = "#0026ff",
                      up.color = "#ff0000",
                      insig.color = "#A6A6A6",
                      sig.thresh = 0.05,
                      fs = NULL,
                      sig.size = 6,
                      insig.size = 5,
                      sig.opacity = 1,
                      insig.opacity = 0.5,
                      label.size = 10,
                      webgl = TRUE,
                      webgl.ratio = 7,
                      show.counts = TRUE,
                      show.hl.counts = TRUE,
                      counts.size = TRUE,
                      highlight.featsets = NULL,
                      highlight.feats = NULL,
                      featsets = NULL,
                      highlight.feats.color = "#E69F00",
                      highlight.feats.size = 7,
                      highlight.feats.opac = 1,
                      highlight.feats.linecolor = "#000000",
                      highlight.feats.linewidth = 1,
                      highlight.feats.label = TRUE,
                      highlight.featsets.color = "#009E73",
                      highlight.featsets.size = 7,
                      highlight.featsets.opac = 1,
                      highlight.featsets.linecolor = "#000000",
                      highlight.featsets.linewidth = 1,
                      highlight.featsets.label = FALSE,
                      h.id.suffix = "_lawn") {
    # Check for required columns.
    if (!x.term %in% colnames(res)) {
        fig <- .empty_plot(paste0("Column '", x.term, "' not found in dataframe."), plotly = TRUE)
    } else if (!lfc.term %in% colnames(res)) {
        fig <- .empty_plot(paste0("Column '", lfc.term, "' not found in dataframe."), plotly = TRUE)
    } else if (!is.null(sig.term) & !sig.term %in% colnames(res)) {
        fig <- .empty_plot(paste0("Column '", sig.term, "' not found in dataframe."), plotly = TRUE)
    } else if (!feat.term %in% colnames(res) & feat.term != "rows") {
        fig <- .empty_plot(paste0("Column '", feat.term, "' not found in dataframe."), plotly = TRUE)
    } else {
        # Styling.
        res$col <- rep(insig.color, nrow(res))
        res$cex <- rep(insig.size, nrow(res))
        res$order <- rep(0, nrow(res))
        res$lcol <- res$col
        res$lw <- 0
        res$opacity <- insig.opacity

        # Remove features with NA significance term (due to low expression, etc).
        res <- res[!is.na(res[[sig.term]]), ]

        # Get all gene IDs or symbols to be highlighted.
        highlight <- highlight.feats

        highlight.fs <- NULL
        if (!is.null(highlight.featsets)) {
            for (featset in highlight.featsets) {
                highlight.fs <- c(highlight.fs, featsets[[featset]])
            }
        }

        # Significance filter.
        up.feats <- res[[sig.term]] < sig.thresh & res[[lfc.term]] > fc.thresh
        res$col[up.feats] <- up.color
        res$cex[up.feats] <- sig.size
        res$order[up.feats] <- 1
        res$opacity[up.feats] <- sig.opacity

        dn.feats <- res[[sig.term]] < sig.thresh & res[[lfc.term]] < -fc.thresh
        res$col[dn.feats] <- down.color
        res$cex[dn.feats] <- sig.size
        res$order[dn.feats] <- 1
        res$opacity[dn.feats] <- sig.opacity

        res$x <- res[[x.term]]
        res$y <- -log10(res[[sig.term]])

        res$col[res$y < -log10(sig.thresh)] <- insig.color

        res$sh <- ifelse(res$y > ylim, "triangle-up-open", 0)

        res$lw <- ifelse(res$sh != 0, 1, 0)

        res$y[res$y > ylim] <- ylim - 0.2

        if (feat.term == "rows") {
            res$feat <- rownames(res)
        } else {
            res$feat <- res[[feat.term]]
        }

        # Gene/geneset highlighting.
        n.fs.hl <- 0
        n.hl <- 0

        if (!is.null(highlight.fs)) {
            highlight.fs <- highlight.fs[highlight.fs %in% res$feat]
            n.fs.hl <- length(res$col[res$feat %in% highlight.fs])

            res$col[res$feat %in% highlight.fs] <- highlight.featsets.color
            res$cex[res$feat %in% highlight.fs] <- highlight.featsets.size
            res$opacity[res$feat %in% highlight.fs] <- highlight.featsets.opac
            res$lcol[res$feat %in% highlight.fs] <- highlight.featsets.linecolor
            res$lw[res$feat %in% highlight.fs] <- highlight.featsets.linewidth
            res$order[res$feat %in% highlight.fs] <- 2

            if (highlight.featsets.label) {
                add.fs <- data.frame(
                    x = res$x[res$feat %in% highlight.fs],
                    y = res$y[res$feat %in% highlight.fs],
                    customdata = res$feat[res$feat %in% highlight.fs]
                )

                if (!is.null(fs)) {
                    fs <- fs[, c("x", "y", "customdata")]
                    fs <- rbind(fs, add.fs)
                } else {
                    fs <- add.fs
                }
            }
        }

        # Want these to have precedence over the feature sets in case entries are in both.
        if (!is.null(highlight)) {
            highlight <- highlight[highlight %in% res$feat]
            if (length(highlight) > 0) {
                n.hl <- length(res$col[res$feat %in% highlight])

                res$col[res$feat %in% highlight] <- highlight.feats.color
                res$cex[res$feat %in% highlight] <- highlight.feats.size
                res$opacity[res$feat %in% highlight] <- highlight.feats.opac
                res$lcol[res$feat %in% highlight] <- highlight.feats.linecolor
                res$lw[res$feat %in% highlight] <- highlight.feats.linewidth
                res$order[res$feat %in% highlight] <- 3

                if (highlight.feats.label) {
                    add.fs <- data.frame(
                        x = res$x[res$feat %in% highlight],
                        y = res$y[res$feat %in% highlight],
                        customdata = res$feat[res$feat %in% highlight]
                    )

                    if (!is.null(fs)) {
                        fs <- fs[, c("x", "y", "customdata")]
                        fs <- rbind(fs, add.fs)
                    } else {
                        fs <- add.fs
                    }
                }
            }
        }

        res$hover.string <- paste0(
            "</br><b>", feat.term, ":</b> ", res$feat,
            "</br><b>", lfc.term, ":</b> ", format(round(res[[lfc.term]], 4), nsmall = 4),
            "</br><b>", sig.term, ":</b> ", format(round(res[[sig.term]], 6), nsmall = 6)
        )

        if (!is.null(hover.info)) {
            for (n in hover.info) {
                res$hover.string <- paste0(res$hover.string, "</br><b>", n, ":</b> ", res[[n]])
            }
        }

        res <- as.data.frame(res)
        res <- res[order(res$order), ]

        # Get feature numbers.
        n.up.feats <- length(up.feats[up.feats == TRUE])
        n.dn.feats <- length(dn.feats[dn.feats == TRUE])
        n.feats <- nrow(res)

        # Add plot border.
        ay <- list(
            showline = TRUE,
            mirror = TRUE,
            linecolor = toRGB("black"),
            linewidth = 0.5,
            title = paste0("-log10(", sig.term, ")"),
            range = list(0, ylim),
            showgrid = FALSE,
            layer = "below traces",
            zeroline = FALSE,
            ticks = "outside",
            zerolinewidth = 0.5
        )

        ax <- list(
            showline = TRUE,
            mirror = TRUE,
            linecolor = toRGB("black"),
            title = "Features",
            linewidth = 0.5,
            title = lfc.term,
            showgrid = FALSE,
            zeroline = FALSE,
            layer = "below traces",
            ticks = "outside",
            zerolinewidth = 0.5
        )

        # Create vertical and horizontal lines.
        sig.hline <- NULL
        if (sig.line) {
            sig.hline <- .hline(y = -log10(sig.thresh), color = "#000000", width = 1, dash = "longdash")
        }

        # Figure generation.
        fig <- plot_ly(res,
            x = ~x,
            y = ~y,
            customdata = ~feat,
            type = "scatter",
            mode = "markers",
            marker = list(
                color = ~col,
                size = ~cex,
                symbol = ~sh,
                line = list(color = ~lcol, width = ~lw),
                opacity = ~opacity
            ),
            text = ~hover.string,
            hoverinfo = "text",
            source = paste0(h.id, h.id.suffix)
        ) %>%
            config(
                edits = list(
                    annotationPosition = TRUE,
                    annotationTail = TRUE
                ),
                toImageButtonOptions = list(format = "svg"),
                displaylogo = FALSE,
                plotGlPixelRatio = webgl.ratio
            )

        if (!is.null(fs)) {
            fig <- fig %>%
                layout(
                    xaxis = ax,
                    yaxis = ay,
                    showlegend = FALSE,
                    shapes = list(sig.hline)
                ) %>%
                add_annotations(
                    x = fs$x, y = fs$y, text = fs$customdata,
                    font = list(size = label.size, family = "Arial"), arrowside = "none"
                )
        } else {
            fig <- fig %>% layout(
                xaxis = ax,
                yaxis = ay,
                showlegend = FALSE,
                shapes = list(sig.hline)
            )
        }

        # Feature count annotations.
        if (show.counts) {
            fig <- fig %>%
                add_annotations(
                    x = 1,
                    y = 1,
                    xref = "paper",
                    yref = "paper",
                    text = paste0(
                        "Up features: ", n.up.feats,
                        "\nDown features: ", n.dn.feats,
                        "\nTotal features: ", n.feats
                    ),
                    showarrow = FALSE,
                    font = list(size = counts.size)
                )
        }

        if (show.hl.counts) {
            fig <- fig %>%
                add_annotations(
                    x = 0,
                    y = 1,
                    xref = "paper",
                    yref = "paper",
                    text = paste0(
                        "Set features: ", n.fs.hl,
                        "\nHighlighted features: ", n.hl
                    ),
                    showarrow = FALSE,
                    font = list(size = counts.size)
                )
        }

        if (webgl) {
            fig <- fig %>% toWebGL()
        }
    }

    fig
}


#' Create a sgRNA pair plot
#'
#' This function creates a pair plot of sgRNA counts for control and treatment conditions.
#'
#' @param df A data frame containing sgRNA counts for control and treatment.
#' @return A plotly object.
#'
#' @author Jared Andrews
#' @rdname INTERNAL_make_sgrna_pairplot
.make_sgrna_pairplot <- function(df) {
    gene <- df$Gene[1]

    df <- data.frame(
        group = c(rep("control", nrow(df)), rep("treatment", nrow(df))),
        counts = c(df$control_mean, df$treat_mean), id = rep(df$sgrna, 2)
    )

    df$hover.string <- paste0(
        "</br><b>Control counts:</b> ", df$counts[df$group == "control"],
        "</br><b>Treatment counts:</b> ", df$counts[df$group == "treatment"],
        "</br><b>sgRNA:</b> ", df$id
    )

    plot_ly(df,
        x = ~group,
        y = ~ counts + 1,
        split = ~id,
        type = "scatter",
        mode = "lines+markers",
        text = ~hover.string,
        hoverinfo = "text"
    ) %>%
        layout(
            showlegend = FALSE, title = paste0(gene, " sgRNAs"),
            yaxis = list(
                range = c(log10(0.8), log10(max(df$counts + 100))),
                type = "log",
                rangemode = "tozero",
                zerolinecolor = "black",
                ticks = "outside",
                showline = TRUE,
                mirror = TRUE,
                zerolinewidth = 2,
                gridcolor = "#ffff",
                title = "Normalized Counts + 1"
            ),
            xaxis = list(
                ticks = "outside",
                showline = TRUE,
                mirror = TRUE,
                title = "",
                showgrid = FALSE
            )
        ) %>%
        config(
            toImageButtonOptions = list(format = "svg"),
            displaylogo = FALSE,
            plotGlPixelRatio = 7
        )
}
