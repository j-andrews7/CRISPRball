#' Create DepMap dataframe for plotting
#'
#' Based on the requested data type, this function will create a dataframe
#' from the DepMap database that can be used for plotting.
#'
#' @param gene Character scalar of gene name.
#' @param data.type Character scalar of data type to retrieve. One of "dependency",
#'   "crispr", "rnai", "cn", or "ccle_tpm".
#' @param depmap.meta data.frame of DepMap cell line metadata, as stored in the 'meta' table
#'   of the SQLite database built by \code{\link{build_depmap_db}}.
#' @param depmap.pool pool connection to DepMap SQLite database built with \code{\link{build_depmap_db}}.
#' @return data.frame containing appropriate DepMap data for plotting.
#'
#' @export
#' @seealso \code{\link{plot_depmap_lineages}}
#'
#' @author Jared Andrews
#' @examples
#' \donttest{
#' library(CRISPRball)
#' build_depmap_db()
#' pool <- pool::dbPool(RSQLite::SQLite(), dbname = "depmap_db.sqlite")
#' depmap.meta <- pool::dbGetQuery(pool, "SELECT * FROM 'meta'")
#' depmap.gene <- pool::dbGetQuery(pool, "SELECT * FROM 'gene.summary'")
#'
#' df <- get_depmap_plot_data(
#'     gene = "CDK2", data.type = "crispr",
#'     depmap.meta = depmap.meta, depmap.pool = pool
#' )
#' }
get_depmap_plot_data <- function(gene, data.type, depmap.meta, depmap.pool) {
    # Get appropriate plot stuff based on datatype.
    switch(data.type,
        dependency = {
            h.text <- "Dependency"
            colname <- "dependency"

            df.c <- pool::dbGetQuery(depmap.pool,
                'SELECT * FROM "crispr" WHERE "gene_name" == (:x)',
                params = list(x = gene)
            )

            df.r <- pool::dbGetQuery(depmap.pool,
                'SELECT * FROM "rnai" WHERE "gene_name" == (:x)',
                params = list(x = gene)
            )

            df <- data.frame()

            if (nrow(df.c) > 0) {
                df.c$dataset <- "CRISPR"
                df <- df.c
            }

            if (nrow(df.r) > 0) {
                df.r$dataset <- "RNAi"

                if (nrow(df) > 0) {
                    df <- rbind(df, df.r)
                } else {
                    df <- df.r
                }
            }
        },
        crispr = {
            h.text <- "Dependency"
            colname <- "dependency"
        },
        rnai = {
            h.text <- "Dependency"
            colname <- "dependency"
        },
        cn = {
            h.text <- "log2(Copy Number)"
            colname <- "log_copy_number"
        },
        ccle_tpm = {
            h.text <- "log2(TPM+1)"
            colname <- "rna_expression"
        }
    )

    # Get correct table.
    if (data.type != "dependency") {
        query <- sprintf('SELECT * FROM "%s" WHERE "gene_name" == (:x)', data.type)
        df <- pool::dbGetQuery(depmap.pool, query, params = list(x = gene))
    }

    if (!is.null(df) && nrow(df) > 0) {
        df$cell_line_name <- depmap.meta$cell_line_name[match(df$depmap_id, depmap.meta$depmap_id)]
        df$primary_disease <- depmap.meta$primary_disease[match(df$depmap_id, depmap.meta$depmap_id)]
        df$lineage <- depmap.meta$lineage[match(df$depmap_id, depmap.meta$depmap_id)]
        df$lineage_subtype <- depmap.meta$lineage_subtype[match(df$depmap_id, depmap.meta$depmap_id)]

        df$hover.string <- paste0(
            "</br><b>Cell Line:</b> ", df$cell_line_name,
            "</br><b>", h.text, ":</b> ", format(round(df[[colname]], 3), nsmall = 3),
            "</br><b>Lineage:</b> ", df$lineage,
            "</br><b>Disease:</b> ", df$primary_disease
        )
    } else {
        df <- NULL
    }

    return(df)
}


#' Plot gene dependency information from DepMap CRISPR and RNAi tables
#'
#' @param df data.frame containing information for a single gene
#'   as returned by \code{\link{get_depmap_plot_data}}.
#' @param crispr.color Character scalar for CRISPR trace color as hexcode.
#' @param rnai.color Character scalar for RNAi trace color as hexcode.
#' @param depline Boolean indicating whether to show the dependency threshold line.
#' @param plot.grid Boolean indicating whether to plot gridlines.
#' @return plotly object
#'
#' @importFrom ggplot2 ggplot theme_bw theme scale_color_manual scale_fill_manual geom_density geom_rug xlab ylab aes_string element_blank geom_vline
#' @importFrom plotly ggplotly layout config %>%
#'
#' @seealso \code{\link{get_depmap_plot_data}}
#'
#' @export
#' @author Jared Andrews
#' @examples
#' library(CRISPRball)
#' data(depmap_22q1_crispr_rnai)
#' plot_depmap_dependency(depmap_22q1_crispr_rnai)
plot_depmap_dependency <- function(df,
                                   crispr.color = "#3584B5",
                                   rnai.color = "#52288E",
                                   depline = TRUE,
                                   plot.grid = FALSE) {
    # Plot construction.
    if (!is.null(df) && nrow(df) > 0) {
        gg <- ggplot() +
            geom_density(data = df, aes(
                x = .data[["dependency"]], color = .data[["dataset"]],
                fill = .data[["dataset"]]
            ), alpha = 0.6) +
            geom_rug(data = df[df$dataset == "CRISPR", ], aes(
                x = .data[["dependency"]],
                color = .data[["dataset"]], text = .data[["hover.string"]]
            ), outside = FALSE) +
            geom_rug(data = df[df$dataset == "RNAi", ], aes(
                x = .data[["dependency"]],
                color = .data[["dataset"]], text = .data[["hover.string"]]
            ), sides = "t") +
            ylab("") +
            xlab("") +
            theme_bw() +
            scale_color_manual(
                values = c(crispr.color, rnai.color),
                breaks = c("CRISPR", "RNAi")
            ) +
            scale_fill_manual(
                values = c(crispr.color, rnai.color),
                breaks = c("CRISPR", "RNAi")
            ) +
            geom_vline(xintercept = 0) +
            theme(legend.position = "none")

        if (!plot.grid) {
            gg <- gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }

        if (depline) {
            gg <- gg + geom_vline(xintercept = -1, color = "red", linetype = "dashed")
        }

        gg <- ggplotly(gg, tooltip = "text") %>%
            layout(
                xaxis = list(
                    title = "Gene Effect"
                ),
                yaxis = list(
                    title = "Density"
                )
            )

        gg %>%
            config(
                edits = list(
                    annotationPosition = TRUE,
                    annotationTail = TRUE
                ),
                toImageButtonOptions = list(format = "svg"),
                displaylogo = FALSE,
                plotGlPixelRatio = 7
            )
    } else {
        .empty_plot("Gene not found in DepMap.", plotly = TRUE)
    }
}


#' Plot gene expression information from DepMap, mostly from CCLE
#'
#' @inheritParams plot_depmap_dependency
#' @param color Character scalar for trace color.
#' @return plotly object
#'
#' @importFrom ggplot2 ggplot theme_bw theme scale_color_manual scale_fill_manual geom_density geom_rug xlab ylab aes_string element_blank
#' @importFrom plotly ggplotly layout config %>%
#'
#' @seealso \code{\link{get_depmap_plot_data}}
#'
#' @export
#' @author Jared Andrews
#' @examples
#' library(CRISPRball)
#' data(depmap_22q1_TPM)
#' plot_depmap_expression(depmap_22q1_TPM)
plot_depmap_expression <- function(df,
                                   color = "#7B8CB2",
                                   plot.grid = FALSE) {
    if (!is.null(df) && nrow(df) > 0) {
        df$color <- color

        gg <- ggplot(show.legend = FALSE) +
            geom_density(data = df, aes(
                x = .data[["rna_expression"]],
                color = .data[["color"]], fill = .data[["color"]]
            )) +
            geom_rug(data = df, aes(
                x = .data[["rna_expression"]],
                color = .data[["color"]], text = .data[["hover.string"]],
                fill = .data[["color"]]
            ), outside = FALSE) +
            ylab("Density") +
            xlab("log2(TPM+1)") +
            theme_bw() +
            scale_color_manual(values = c(df$color), breaks = c(df$color)) +
            scale_fill_manual(values = c(df$color), breaks = c(df$color)) +
            theme(legend.position = "none")

        if (!plot.grid) {
            gg <- gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }

        gg <- ggplotly(gg, tooltip = "text")

        gg %>%
            config(
                edits = list(
                    annotationPosition = TRUE,
                    annotationTail = TRUE
                ),
                toImageButtonOptions = list(format = "svg"),
                displaylogo = FALSE,
                plotGlPixelRatio = 7
            )
    } else {
        .empty_plot("Gene not found in DepMap.", plotly = TRUE)
    }
}


#' Plot gene copy number information from DepMap
#'
#' @inheritParams plot_depmap_dependency
#' @param color Character scalar for trace color.
#' @return plotly object
#'
#' @importFrom ggplot2 ggplot theme_bw theme scale_color_manual scale_fill_manual geom_density geom_rug xlab ylab aes_string element_blank
#' @importFrom plotly ggplotly layout config %>%
#'
#' @seealso \code{\link{get_depmap_plot_data}}
#'
#' @export
#' @author Jared Andrews
#' @examples
#' library(CRISPRball)
#' data(depmap_22q1_cn)
#' plot_depmap_cn(depmap_22q1_cn)
plot_depmap_cn <- function(df,
                           color = "#CEA3CB",
                           plot.grid = FALSE) {
    if (!is.null(df) && nrow(df) > 0) {
        df$color <- color

        gg <- ggplot(show.legend = FALSE) +
            geom_density(data = df, aes(
                x = .data[["log_copy_number"]],
                color = .data[["color"]], fill = .data[["color"]]
            )) +
            geom_rug(data = df, aes(
                x = .data[["log_copy_number"]],
                color = .data[["color"]], text = .data[["hover.string"]],
                fill = .data[["color"]]
            ), outside = FALSE) +
            ylab("Density") +
            xlab("log2(Copy Number)") +
            theme_bw() +
            scale_color_manual(values = c(df$color), breaks = c(df$color)) +
            scale_fill_manual(values = c(df$color), breaks = c(df$color)) +
            theme(legend.position = "none")

        if (!plot.grid) {
            gg <- gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }

        gg <- ggplotly(gg, tooltip = "text")

        gg %>%
            config(
                edits = list(
                    annotationPosition = TRUE,
                    annotationTail = TRUE
                ),
                toImageButtonOptions = list(format = "svg"),
                displaylogo = FALSE,
                plotGlPixelRatio = 7
            )
    } else {
        .empty_plot("Gene not found in DepMap.", plotly = TRUE)
    }
}


#' Plot selected information across lineages from DepMap.
#'
#' @inheritParams plot_depmap_dependency
#' @param plot.val Character scalar of column name to plot values from.
#' @param group.by Character scalar of column name to group by.
#' @param lineage Character scalar of lineage for which
#'   to plot sub-lineage data.
#' @param label.size Numeric scaler for axis label size.
#' @param pt.size Numeric scalar for point size.
#' @param pt.color Character scalar for point color.
#' @param boxplot.fill Character scalar for boxplot fill color.
#' @param boxplot.line.color Character scalar for boxplot line color.
#' @return plotly object
#'
#' @importFrom ggplot2 ggplot theme_bw theme scale_color_manual scale_fill_manual geom_density geom_rug xlab ylab aes_string element_blank
#' @importFrom plotly ggplotly layout config %>%
#'
#' @seealso \code{\link{get_depmap_plot_data}}
#'
#' @export
#' @author Jared Andrews
#' @examples
#' library(CRISPRball)
#' data("depmap_22q1_rnai")
#' plot_depmap_lineages(df = depmap_22q1_rnai, plot.val = "dependency", group.by = "lineage")
plot_depmap_lineages <- function(df,
                                 plot.val,
                                 group.by,
                                 lineage = NULL,
                                 depline = TRUE,
                                 label.size = 12,
                                 pt.size = 5,
                                 pt.color = "#56B4E9",
                                 boxplot.fill = "#E2E2E2",
                                 boxplot.line.color = "#000000") {
    if (!is.null(df) && nrow(df) > 0) {
        # Get correct lineage.
        if (!is.null(lineage)) {
            df <- df[df$lineage == lineage, ]
        }

        # Get counts in each group and add to labels.
        ylabs <- paste0(names(table(df[[group.by]])), " (", table(df[[group.by]]), ")")

        # Add plot border.
        ay <- list(
            showline = TRUE,
            mirror = TRUE,
            linecolor = toRGB("black"),
            linewidth = 0.5,
            showgrid = FALSE,
            layer = "below traces",
            zeroline = FALSE,
            ticks = "outside",
            zerolinewidth = 0.5,
            tickvals = names(table(df[[group.by]])),
            ticktext = ylabs,
            tickfont = list(size = label.size),

            # Fix extra whitespace at top and bottom of plot
            range = ~ c(-1, length(unique(df[[group.by]])) - 0.5)
        )

        ax <- list(
            showline = TRUE,
            mirror = TRUE,
            linecolor = toRGB("black"),
            linewidth = 0.5,
            zeroline = TRUE,
            showgrid = FALSE,
            layer = "below traces",
            ticks = "outside",
            zerolinewidth = 0.5
        )

        fig <- plot_ly(df,
            x = as.formula(paste0("~", plot.val)),
            y = as.formula(paste0("~", group.by)),
            fillcolor = boxplot.fill,
            color = I(boxplot.line.color),
            type = "box",
            boxpoints = FALSE,
            alpha = 1
        )

        fig <- fig %>% add_trace(
            type = "scatter",
            x = as.formula(paste0("~", plot.val)),
            y = as.formula(paste0("~", group.by)),
            mode = "markers",
            text = ~hover.string,
            hoverinfo = "text",
            marker = list(
                color = pt.color,
                size = pt.size
            )
        )

        if (depline & plot.val == "dependency") {
            dline <- .vline(x = -1, dash = "longdash", width = 1)
        } else {
            dline <- NULL
        }

        fig <- fig %>% layout(
            showlegend = FALSE,
            shapes = list(dline),
            xaxis = ax,
            yaxis = ay
        )

        fig %>%
            config(
                edits = list(
                    annotationPosition = TRUE,
                    annotationTail = TRUE
                ),
                toImageButtonOptions = list(format = "svg"),
                displaylogo = FALSE,
                plotGlPixelRatio = 7
            )
    } else {
        .empty_plot("Gene not found in DepMap.", plotly = TRUE)
    }
}
