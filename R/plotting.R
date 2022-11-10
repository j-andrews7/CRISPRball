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

# Generic volcano plot function.
# TODO: Add defaults, document, and export.
.make_volcano <- function(res, xlim, ylim, fc.thresh, fc.lines, hover.info = NULL,
                          sig.line, h.id, feat.term, sig.term, lfc.term, down.color, up.color,
                          insig.color, sig.thresh = 0.05, fs = NULL, sig.size, insig.size,
                          sig.opacity, insig.opacity, label.size, webgl, webgl.ratio, show.counts,
                          show.hl.counts, counts.size, highlight.featsets, highlight.feats, featsets,
                          highlight.feats.color, highlight.feats.size, highlight.feats.opac,
                          highlight.feats.linecolor, highlight.feats.linewidth,
                          highlight.featsets.color, highlight.featsets.size, highlight.featsets.opac,
                          highlight.featsets.linecolor, highlight.featsets.linewidth, h.id.suffix = "_volc") {

  # Styling.
  res$col <- rep(insig.color, nrow(res))
  res$cex <- rep(insig.size, nrow(res))
  res$order <- rep(0, nrow(res))
  res$lcol <- res$col
  res$lw <- 0
  res$opacity <- insig.opacity

  # Remove features with NA significance term (due to low expression, etc).
  res <- res[!is.na(res[[sig.term]]),]

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
                          ifelse(res$x > xlim, "triangle-right-open", 0)))

  res$lw <- ifelse(res$sh != 0, 1, 0)

  res$y[res$y > ylim] <- ylim - 0.2
  res$x[res$x > xlim] <- xlim - 0.05
  res$x[res$x < -xlim] <- -xlim + 0.05
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
  }

  # Want these to have precedence over the feature sets in case entries are in both.
  if (!is.null(highlight)) {
    highlight <- highlight[highlight %in% res$feat]
    n.hl <- length(res$col[res$feat %in% highlight])

    res$col[res$feat %in% highlight] <- highlight.feats.color
    res$cex[res$feat %in% highlight] <- highlight.feats.size
    res$opacity[res$feat %in% highlight] <- highlight.feats.opac
    res$lcol[res$feat %in% highlight] <- highlight.feats.linecolor
    res$lw[res$feat %in% highlight] <- highlight.feats.linewidth
    res$order[res$feat %in% highlight] <- 3
  }

  res$hover.string <- paste0("</br><b>", feat.term, ":</b> ", res$feat,
                            "</br><b>", lfc.term, ":</b> ", format(round(res[[lfc.term]], 4), nsmall = 4),
                            "</br><b>", sig.term, ":</b> ", format(round(res[[sig.term]], 6), nsmall = 6))

  if (!is.null(hover.info)) {
    for (n in hover.info) {
      res$hover.string <- paste0(res$hover.string, "</br><b>", n, ":</b>", res[[n]])
    }
  }

  res <- as.data.frame(res)
  res <- res[order(res$order),]

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
  if(sig.line) {
    sig.hline <- .hline(y = -log10(sig.thresh), color = "#000000", width = 1, dash = "longdash")
  }

  if (fc.thresh != 0 & fc.lines) {
    fc.line1 <- .vline(x = fc.thresh, color = "#000000", width = 1, dash = "longdash")
    fc.line2 <- .vline(x = -fc.thresh, color = "#000000", width = 1, dash = "longdash")
  }

  # Figure generation.
  fig <- plot_ly(res, x = ~x,
                 y = ~y,
                 customdata = ~feat,
                 type = "scatter",
                 mode = "markers",
                 marker = list(color = ~col,
                               size = ~cex,
                               symbol = ~sh,
                               line = list(color = ~lcol, width = ~lw),
                               opacity = ~opacity),
                 text = ~hover.string,
                 hoverinfo = "text",
                 source = paste0(h.id, h.id.suffix)) %>%
    config(edits = list(annotationPosition = TRUE,
                        annotationTail = TRUE),
           toImageButtonOptions = list(format = "svg"),
           displaylogo = FALSE,
           plotGlPixelRatio = webgl.ratio)

  if (!is.null(fs)) {
    fig <- fig %>%
      layout(xaxis = ax,
             yaxis = ay,
             showlegend = FALSE,
             shapes = list(sig.hline, fc.line1, fc.line2),
             hoverlabel = list(font=list(size=10))) %>%
      add_annotations(x = fs$x, y = fs$y, text = fs$customdata,
                      font = list(size = label.size, family = "Arial"), arrowside = "none")
  } else {
    fig <- fig %>% layout(xaxis = ax,
                          yaxis = ay,
                          showlegend = FALSE,
                          shapes = list(sig.hline, fc.line1, fc.line2),
                          hoverlabel = list(font=list(size=10)))
  }

  # Feature count annotations.
  if (show.counts) {
    fig <- fig %>%
      add_annotations(
        x= 1,
        y= 1,
        xref = "paper",
        yref = "paper",
        text = paste0("Up features: ", n.up.feats,
                      "\nDown features: ", n.dn.feats,
                      "\nTotal features: ", n.feats),
        showarrow = FALSE,
        font = list(size = counts.size)
      )
  }

  if (show.hl.counts) {
    fig <- fig %>%
      add_annotations(
        x= 0,
        y= 1,
        xref = "paper",
        yref = "paper",
        text = paste0("Set features: ", n.fs.hl,
                      "\nHighlighted features: ", n.hl),
        showarrow = FALSE,
        font = list(size = counts.size)
      )
  }

  if (webgl) {
    fig <- fig %>% toWebGL()
  }

  fig
}


# Generic rank plot function.
# TODO: Add defaults, document, and export.
.make_rank <- function(df, ylim, y.thresh, y.lines, hover.info = NULL,
                          h.id, feat.term, sig.term, y.term, x.term, down.color, up.color,
                          insig.color, sig.thresh = 0.05, fs = NULL, sig.size, insig.size,
                          sig.opacity, insig.opacity, label.size, webgl, webgl.ratio, show.counts,
                          show.hl.counts, counts.size, highlight.featsets, highlight.feats, featsets,
                          highlight.feats.color, highlight.feats.size, highlight.feats.opac,
                          highlight.feats.linecolor, highlight.feats.linewidth,
                          highlight.featsets.color, highlight.featsets.size, highlight.featsets.opac,
                          highlight.featsets.linecolor, highlight.featsets.linewidth, h.id.suffix = "_rank") {
  
  # Styling.
  df$col <- rep(insig.color, nrow(df))
  df$cex <- rep(insig.size, nrow(df))
  df$order <- rep(0, nrow(df))
  df$lcol <- df$col
  df$lw <- 0
  df$opacity <- insig.opacity

  # Remove features with NA significance term (due to low expression, etc).
  if (!is.null(sig.term)){
    df <- df[!is.na(df[[sig.term]]),]
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
    up.feats <- df[[sig.term]] < sig.thresh & df[[y.term]] > y.thresh
  } else {
    up.feats <- df[[y.term]] > y.thresh
  }

  df$col[up.feats] <- up.color
  df$cex[up.feats] <- sig.size
  df$order[up.feats] <- 1
  df$opacity[up.feats] <- sig.opacity

  if (!is.null(sig.term) & !is.null(sig.thresh)) {
    dn.feats <- df[[sig.term]] < sig.thresh & df[[y.term]] < -y.thresh
  } else {
    dn.feats <- df[[y.term]] < -y.thresh
  }

  df$col[dn.feats] <- down.color
  df$cex[dn.feats] <- sig.size
  df$order[dn.feats] <- 1
  df$opacity[dn.feats] <- sig.opacity

  df$x <- df[[x.term]]
  df$y <- df[[y.term]]

  df$col[df$y < y.thresh & df$y > -y.thresh] <- insig.color

  df$sh <- ifelse(df$y > ylim[[2]], "triangle-up-open", ifelse(df$y < ylim[[1]], "triangle-down-open",0))

  df$lw <- ifelse(df$sh != 0, 1, 0)

  # Get feature identifier for hover/labeling.
  if (feat.term == "rows") {
    df$feat <- rownames(df)
  } else {
    df$feat <- df[[feat.term]]
  }

  # Gene/geneset highlighting.
  n.fs.hl <- 0
  n.hl <- 0

  if (!is.null(highlight.fs)) {
    highlight.fs <- highlight.fs[highlight.fs %in% df$feat]
    n.fs.hl <- length(df$col[df$feat %in% highlight.fs])

    df$col[df$feat %in% highlight.fs] <- highlight.featsets.color
    df$cex[df$feat %in% highlight.fs] <- highlight.featsets.size
    df$opacity[df$feat %in% highlight.fs] <- highlight.featsets.opac
    df$lcol[df$feat %in% highlight.fs] <- highlight.featsets.linecolor
    df$lw[df$feat %in% highlight.fs] <- highlight.featsets.linewidth
    df$order[df$feat %in% highlight.fs] <- 2
  }

  # Want these to have precedence over the feature sets in case entries are in both.
  if (!is.null(highlight)) {
    highlight <- highlight[highlight %in% df$feat]
    n.hl <- length(df$col[df$feat %in% highlight])

    df$col[df$feat %in% highlight] <- highlight.feats.color
    df$cex[df$feat %in% highlight] <- highlight.feats.size
    df$opacity[df$feat %in% highlight] <- highlight.feats.opac
    df$lcol[df$feat %in% highlight] <- highlight.feats.linecolor
    df$lw[df$feat %in% highlight] <- highlight.feats.linewidth
    df$order[df$feat %in% highlight] <- 3
  }

  df$hover.string <- paste("</br><b>", feat.term, ":</b> ", df$feat,
                           "</br><b>", x.term, ":</b> ", df[[x.term]],
                            "</br><b>", y.term, ":</b> ", format(round(df[[y.term]], 4), nsmall = 4),
                            "</br><b>", sig.term, ":</b> ", format(round(df[[sig.term]], 6), nsmall = 6))

  if (!is.null(hover.info)) {
    for (n in hover.info) {
      df$hover.string <- paste0(df$hover.string, "</br><b>", n, ":</b> ", df[[n]])
    }
  }

  df <- as.data.frame(df)
  df <- df[order(df$order),]

  # Get feature numbers.
  n.up.feats <- length(up.feats[up.feats == TRUE])
  n.dn.feats <- length(dn.feats[dn.feats == TRUE])
  n.feats <- nrow(df)

  # Add plot border.
  ay <- list(
    showline = TRUE,
    mirror = TRUE,
    linecolor = toRGB("black"),
    linewidth = 0.5,
    title = y.term,
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
    title = x.term,
    range = list(-(0.03 * nrow(df)), nrow(df) + (0.03 * nrow(df))),
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
  fig <- plot_ly(df, x = ~x,
                 y = ~y,
                 customdata = ~feat,
                 type = "scatter",
                 mode = "markers",
                 marker = list(color = ~col,
                               size = ~cex,
                               symbol = ~sh,
                               line = list(color = ~lcol, width = ~lw),
                               opacity = ~opacity),
                 text = ~hover.string,
                 hoverinfo = "text",
                 source = paste0(h.id, h.id.suffix)) %>%
    config(edits = list(annotationPosition = TRUE,
                        annotationTail = TRUE),
           toImageButtonOptions = list(format = "svg"),
           displaylogo = FALSE,
           plotGlPixelRatio = webgl.ratio)

  if (!is.null(fs)) {
    fig <- fig %>%
      layout(xaxis = ax,
             yaxis = ay,
             showlegend = FALSE,
             shapes = list(y.line1, y.line2)) %>%
      add_annotations(x = fs$x, y = fs$y, text = fs$customdata,
                      font = list(size = label.size, family = "Arial"), arrowside = "none")
  } else {
    fig <- fig %>% layout(xaxis = ax,
                          yaxis = ay,
                          showlegend = FALSE,
                          shapes = list(y.line1, y.line2))
  }

  # Feature count annotations.
  if (show.counts) {
    fig <- fig %>%
      add_annotations(
        x= 1,
        y= 0,
        xref = "paper",
        yref = "paper",
        text = paste0("Up features: ", n.up.feats,
                      "\nDown features: ", n.dn.feats,
                      "\nTotal features: ", n.feats),
        showarrow = FALSE,
        font = list(size = counts.size)
      )
  }

  if (show.hl.counts) {
    fig <- fig %>%
      add_annotations(
        x= 0,
        y= 1,
        xref = "paper",
        yref = "paper",
        text = paste0("Set features: ", n.fs.hl,
                      "\nHighlighted features: ", n.hl),
        showarrow = FALSE,
        font = list(size = counts.size)
      )
  }

  if (webgl) {
    fig <- fig %>% toWebGL()
  }

  fig
}

# Generic lawn plot function.
# TODO: Add defaults, document, and export.
.make_lawn <- function(res, ylim, fc.thresh, hover.info = NULL,
                          sig.line, h.id, feat.term, x.term, sig.term, lfc.term, down.color, up.color,
                          insig.color, sig.thresh = 0.05, fs = NULL, sig.size, insig.size,
                          sig.opacity, insig.opacity, label.size, webgl, webgl.ratio, show.counts,
                          show.hl.counts, counts.size, highlight.featsets, highlight.feats, featsets,
                          highlight.feats.color, highlight.feats.size, highlight.feats.opac,
                          highlight.feats.linecolor, highlight.feats.linewidth,
                          highlight.featsets.color, highlight.featsets.size, highlight.featsets.opac,
                          highlight.featsets.linecolor, highlight.featsets.linewidth, h.id.suffix = "_lawn") {

  # Styling.
  res$col <- rep(insig.color, nrow(res))
  res$cex <- rep(insig.size, nrow(res))
  res$order <- rep(0, nrow(res))
  res$lcol <- res$col
  res$lw <- 0
  res$opacity <- insig.opacity

  # Remove features with NA significance term (due to low expression, etc).
  res <- res[!is.na(res[[sig.term]]),]

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
  }

  # Want these to have precedence over the feature sets in case entries are in both.
  if (!is.null(highlight)) {
    highlight <- highlight[highlight %in% res$feat]
    n.hl <- length(res$col[res$feat %in% highlight])

    res$col[res$feat %in% highlight] <- highlight.feats.color
    res$cex[res$feat %in% highlight] <- highlight.feats.size
    res$opacity[res$feat %in% highlight] <- highlight.feats.opac
    res$lcol[res$feat %in% highlight] <- highlight.feats.linecolor
    res$lw[res$feat %in% highlight] <- highlight.feats.linewidth
    res$order[res$feat %in% highlight] <- 3
  }

  res$hover.string <- paste0("</br><b>", feat.term, ":</b> ", res$feat,
                             "</br><b>", lfc.term, ":</b> ", format(round(res[[lfc.term]], 4), nsmall = 4),
                             "</br><b>", sig.term, ":</b> ", format(round(res[[sig.term]], 6), nsmall = 6))

  if (!is.null(hover.info)) {
    for (n in hover.info) {
      res$hover.string <- paste0(res$hover.string, "</br><b>", n, ":</b> ", res[[n]])
    }
  }

  res <- as.data.frame(res)
  res <- res[order(res$order),]

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
  if(sig.line) {
    sig.hline <- .hline(y = -log10(sig.thresh), color = "#000000", width = 1, dash = "longdash")
  }

  # Figure generation.
  fig <- plot_ly(res, x = ~x,
                 y = ~y,
                 customdata = ~feat,
                 type = "scatter",
                 mode = "markers",
                 marker = list(color = ~col,
                               size = ~cex,
                               symbol = ~sh,
                               line = list(color = ~lcol, width = ~lw),
                               opacity = ~opacity),
                 text = ~hover.string,
                 hoverinfo = "text",
                 source = paste0(h.id, h.id.suffix)) %>%
    config(edits = list(annotationPosition = TRUE,
                        annotationTail = TRUE),
           toImageButtonOptions = list(format = "svg"),
           displaylogo = FALSE,
           plotGlPixelRatio = webgl.ratio)

  if (!is.null(fs)) {
    fig <- fig %>%
      layout(xaxis = ax,
             yaxis = ay,
             showlegend = FALSE,
             shapes = list(sig.hline)) %>%
      add_annotations(x = fs$x, y = fs$y, text = fs$customdata,
                      font = list(size = label.size, family = "Arial"), arrowside = "none")
  } else {
    fig <- fig %>% layout(xaxis = ax,
                          yaxis = ay,
                          showlegend = FALSE,
                          shapes = list(sig.hline))
  }

  # Feature count annotations.
  if (show.counts) {
    fig <- fig %>%
      add_annotations(
        x= 1,
        y= 1,
        xref = "paper",
        yref = "paper",
        text = paste0("Up features: ", n.up.feats,
                      "\nDown features: ", n.dn.feats,
                      "\nTotal features: ", n.feats),
        showarrow = FALSE,
        font = list(size = counts.size)
      )
  }

  if (show.hl.counts) {
    fig <- fig %>%
      add_annotations(
        x= 0,
        y= 1,
        xref = "paper",
        yref = "paper",
        text = paste0("Set features: ", n.fs.hl,
                      "\nHighlighted features: ", n.hl),
        showarrow = FALSE,
        font = list(size = counts.size)
      )
  }

  if (webgl) {
    fig <- fig %>% toWebGL()
  }

  fig
}

# sgRNA pair plot.
.make_sgrna_pairplot <- function(df) {
  gene <- df$Gene[1]

  df <- data.frame(group = c(rep("control", nrow(df)), rep("treatment", nrow(df))), 
                   counts = c(df$control_mean, df$treat_mean), id = rep(df$sgrna, 2))
  
  df$hover.string <- paste0("</br><b>Control counts:</b> ", df$counts[df$group == "control"],
                            "</br><b>Treatment counts:</b> ", df$counts[df$group == "treatment"],
                            "</br><b>sgRNA:</b> ", df$id)
  
  plot_ly(df,
          x = ~group,
          y = ~counts + 1,
          split = ~id,
          type = "scatter",
          mode = "lines+markers",
          text = ~hover.string,
          hoverinfo = "text") %>%
    layout(showlegend = FALSE, title = paste0(gene, " sgRNAs"),
           yaxis = list(range = c(log10(0.8), log10(max(df$counts+100))),
                        type = "log",
                        rangemode = "tozero",
                        zerolinecolor = "black",
                        ticks = "outside",
                        showline = TRUE,
                        mirror = TRUE,
                        zerolinewidth = 2,
                        gridcolor = "#ffff",
                        title = "Normalized Counts + 1"),
           xaxis = list(ticks = "outside",
                        showline = TRUE,
                        mirror = TRUE,
                        title = "",
                        showgrid = FALSE)) %>%
    config(toImageButtonOptions = list(format = "svg"),
           displaylogo = FALSE,
           plotGlPixelRatio = 7)
}


#' Plot gene dependency information from DepMap CRISPR and RNAi tables
#' 
#' @param gene Character scalar for gene symbol.
#' @param depmap.meta data.frame of DepMap cell line metadata, as stored in the 'meta' table 
#'   of the SQLite database built by \code{\link{build_depmap_db}}.
#' @param crispr.color Character scalar for CRISPR trace color.
#' @param rnai.color Character scalar for RNAi trace color.
#' @param depline Boolean indicating whether to show the dependency line.
#' @param depmap.pool pool connection to DepMap SQLite database built with \code{\link{build_depmap_db}}.
#' @param plot.grid Boolean indicating whether to plot gridlines.
#' @return plotly object
#'   
#' @export
#' @author Jared Andrews  
plot_depmap_dependency <- function(gene, depmap.meta, crispr.color, 
                                   rnai.color, depline, depmap.pool, plot.grid) {
  
  # Data pull and aggregation.
  df.c <- pool::dbGetQuery(depmap.pool, 'SELECT * FROM "crispr" WHERE "gene_name" == (:x)', params = list(x = gene))
  df.r <- pool::dbGetQuery(depmap.pool, 'SELECT * FROM "rnai" WHERE "gene_name" == (:x)', params = list(x = gene))

  df <- data.frame()

  if (nrow(df.c) > 0) {
    df.c$dataset <- "CRISPR"
    df <- df.c
  }
  
  if (nrow(df.r) > 0) {
    df.r$dataset <- "RNAi"

    if(nrow(df) > 0) {
      df <- rbind(df, df.r)
    } else {
      df <- df.r
    }
  }

  # Plot construction.
  if (nrow(df) > 0) {
    df$cell_line_name <- depmap.meta$cell_line_name[match(df$depmap_id, depmap.meta$depmap_id)]
    df$primary_disease <- depmap.meta$primary_disease[match(df$depmap_id, depmap.meta$depmap_id)]
    df$lineage <- depmap.meta$lineage[match(df$depmap_id, depmap.meta$depmap_id)]
    df$lineage_subtype <- depmap.meta$lineage_subtype[match(df$depmap_id, depmap.meta$depmap_id)]
    
    df$hover.string <- paste0("</br><b>Cell Line:</b> ", df$cell_line_name,
                              "</br><b>Gene Effect:</b> ", format(round(df$dependency, 3), nsmall = 3),
                              "</br><b>Lineage:</b> ", df$lineage,
                              "</br><b>Disease:</b> ", df$primary_disease)

    gg <- ggplot() +
      geom_density(data = df, aes_string(x="dependency", color="dataset", fill="dataset"), alpha = 0.6) +
      geom_rug(data = df[df$dataset == "CRISPR",], aes_string(x="dependency", color="dataset", text="hover.string"), outside = FALSE) +
      geom_rug(data = df[df$dataset == "RNAi",], aes_string(x="dependency", color="dataset", text="hover.string"), sides = "t") +
      ylab("") +
      xlab("") +
      theme_bw() +
      scale_color_manual(values=c(crispr.color, rnai.color), 
                         breaks = c("CRISPR", "RNAi")) +
      scale_fill_manual(values=c(crispr.color, rnai.color), 
                        breaks = c("CRISPR", "RNAi")) +
      geom_vline(xintercept = 0) + theme(legend.position="none")

    if (!plot.grid) {
      gg <- gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }

    if (depline) {
      gg <- gg + geom_vline(xintercept = -1, color = "red", linetype = "dashed")
    }

    gg <- ggplotly(gg, tooltip = "text") %>%
      layout(
        xaxis = list(
          title="Gene Effect"
        ),
        yaxis = list(
          title="Density"))
    
    gg %>%
      config(edits = list(annotationPosition = TRUE,
                          annotationTail = TRUE),
             toImageButtonOptions = list(format = "svg"),
             displaylogo = FALSE,
             plotGlPixelRatio = 7)
  } else {
    .plot_gene_not_found(gene)
  }
}


#' Plot gene expression information from DepMap, mostly from CCLE
#' 
#' @inheritParams plot_depmap_dependency
#' @param color Character scalar for trace color.
#' 
#' @return plotly object
#' 
#' @export
#' @author Jared Andrews  
plot_depmap_expression <- function(gene, depmap.meta, depmap.pool, color, plot.grid) {
  
  df <- pool::dbGetQuery(depmap.pool, 'SELECT * FROM "ccle_tpm" WHERE "gene_name" == (:x)', params = list(x = gene))
  
  if (nrow(df) > 0) {
    df$cell_line_name <- depmap.meta$cell_line_name[match(df$depmap_id, depmap.meta$depmap_id)]
    df$primary_disease <- depmap.meta$primary_disease[match(df$depmap_id, depmap.meta$depmap_id)]
    df$lineage <- depmap.meta$lineage[match(df$depmap_id, depmap.meta$depmap_id)]
    df$lineage_subtype <- depmap.meta$lineage_subtype[match(df$depmap_id, depmap.meta$depmap_id)]
    
    df$hover.string <- paste0("</br><b>Cell Line:</b> ", df$cell_line_name,
                              "</br><b>log2(TPM+1):</b> ", format(round(df$rna_expression, 3), nsmall = 3),
                              "</br><b>Lineage:</b> ", df$lineage,
                              "</br><b>Disease:</b> ", df$primary_disease)
    df$color <- color
    
    gg <- ggplot(show.legend = FALSE) +
      geom_density(data = df, aes_string(x="rna_expression", color="color", fill="color")) +
      geom_rug(data = df, aes_string(x="rna_expression", color="color", text="hover.string", fill="color"), outside = FALSE) +
      ylab("Density") +
      xlab("log2(TPM+1)") +
      theme_bw() +
      scale_color_manual(values=c(df$color), breaks = c(df$color)) +
      scale_fill_manual(values=c(df$color), breaks = c(df$color)) + theme(legend.position="none")
    
    if (!plot.grid) {
      gg <- gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }

    gg <- ggplotly(gg, tooltip = "text")
    
    gg %>%
      config(edits = list(annotationPosition = TRUE,
                          annotationTail = TRUE),
             toImageButtonOptions = list(format = "svg"),
             displaylogo = FALSE,
             plotGlPixelRatio = 7)
  } else {
    .plot_gene_not_found(gene)
  }
}


#' Plot gene copy number information from DepMap, mostly from CCLE
#' 
#' @inheritParams plot_depmap_dependency
#' @param color Character scalar for trace color.
#' @return plotly object
#'   
#' @export
#' @author Jared Andrews  
plot_depmap_cn <- function(gene, depmap.meta, depmap.pool, color, plot.grid) {
  
  df <- pool::dbGetQuery(depmap.pool, 'SELECT * FROM "cn" WHERE "gene_name" == (:x)', params = list(x = gene))
  
  if (nrow(df) > 0) {
    df$cell_line_name <- depmap.meta$cell_line_name[match(df$depmap_id, depmap.meta$depmap_id)]
    df$primary_disease <- depmap.meta$primary_disease[match(df$depmap_id, depmap.meta$depmap_id)]
    df$lineage <- depmap.meta$lineage[match(df$depmap_id, depmap.meta$depmap_id)]
    df$lineage_subtype <- depmap.meta$lineage_subtype[match(df$depmap_id, depmap.meta$depmap_id)]
    
    df$hover.string <- paste0("</br><b>Cell Line:</b> ", df$cell_line_name,
                              "</br><b>Copy Number (log2):</b> ", format(round(df$log_copy_number, 3), nsmall = 3),
                              "</br><b>Lineage:</b> ", df$lineage,
                              "</br><b>Disease:</b> ", df$primary_disease)
    df$color <- color
    
    gg <- ggplot(show.legend = FALSE) +
      geom_density(data = df, aes_string(x="log_copy_number", color="color", fill="color")) +
      geom_rug(data = df, aes_string(x="log_copy_number", color="color", text="hover.string", fill="color"), outside = FALSE) +
      ylab("Density") +
      xlab("log2(Copy Number)") +
      theme_bw() +
      scale_color_manual(values=c(df$color), breaks = c(df$color)) +
      scale_fill_manual(values=c(df$color), breaks = c(df$color)) + theme(legend.position="none")
    
    if (!plot.grid) {
      gg <- gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }

    gg <- ggplotly(gg, tooltip = "text")
    
    gg %>%
      config(edits = list(annotationPosition = TRUE,
                          annotationTail = TRUE),
             toImageButtonOptions = list(format = "svg"),
             displaylogo = FALSE,
             plotGlPixelRatio = 7)
  } else {
    .plot_gene_not_found(gene)
  }
}

#' Plot selected information across lineages from DepMap.
#'
#' @inheritParams plot_depmap_dependency
#' @param data.type One of "crispr", "rnai", "cn", or "ccle_tpm".
#' @param group.by Character scalar of column name to group by.
#' @param lineage Character scalar of lineage for which to plot sub-lineage data.
#' @param label.size Numeric scaler for axis label size.
#' @param pt.size Numeric scalar for point size.
#' @param pt.color Character scalar for point color.
#' @param boxplot.fill Character scalar for boxplot fill color.
#' @param boxplot.line.color Character scalar for boxplot line color.
#' @return plotly object
#'
#' @export
#' @author Jared Andrews
plot_depmap_lineages <- function(gene, data.type, group.by, depmap.meta, depmap.pool, 
                                 lineage = NULL,
                                 depline = TRUE, label.size = 12, pt.size = 5, pt.color = "#56B4E9", 
                                 boxplot.fill = "#E2E2E2", 
                                 boxplot.line.color = "#000000") {
  
  # Get correct table.
  query <- sprintf('SELECT * FROM "%s" WHERE "gene_name" == (:x)', data.type)
  df <- pool::dbGetQuery(depmap.pool, query, params = list(x = gene))

  # Get appropriate plot stuff based on datatype.
  switch(data.type,
         crispr={
           h.text <- "Dependency"
           colname <- "dependency"
         },
         rnai={
           h.text <- "Dependency"
           colname <- "dependency"
         },
         cn={
           h.text <- "log2(Copy Number)"
           colname <- "log_copy_number"
         },
         ccle_tpm={
           h.text <- "log2(TPM+1)"
           colname <- "rna_expression"
         })
  
  if (nrow(df) > 0) {
    df$cell_line_name <- depmap.meta$cell_line_name[match(df$depmap_id, depmap.meta$depmap_id)]
    df$primary_disease <- depmap.meta$primary_disease[match(df$depmap_id, depmap.meta$depmap_id)]
    df$lineage <- depmap.meta$lineage[match(df$depmap_id, depmap.meta$depmap_id)]
    df$lineage_subtype <- depmap.meta$lineage_subtype[match(df$depmap_id, depmap.meta$depmap_id)]

    # Get correct lineage.
    if (!is.null(lineage)) {
      df <- df[df$lineage == lineage,]
    }

    df$hover.string <- paste0("</br><b>Cell Line:</b> ", df$cell_line_name,
                              "</br><b>", h.text, ":</b> ", format(round(df[[colname]], 3), nsmall = 3),
                              "</br><b>Lineage:</b> ", df$lineage,
                              "</br><b>Disease:</b> ", df$primary_disease)

    # Get counts in each group and add to labels.
    ylabs <- paste0(names(table(df[[group.by]]))," (",table(df[[group.by]]),")")

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
                   x = as.formula(paste0("~", colname)), 
                   y = as.formula(paste0("~", group.by)), 
                   fillcolor = boxplot.fill,
                   color = I(boxplot.line.color),
                   type = "box", 
                   boxpoints = FALSE, 
                   alpha = 1)

    fig <- fig %>% add_trace(type = "scatter",
                            x = as.formula(paste0("~", colname)), 
                            y = as.formula(paste0("~", group.by)), 
                            mode = "markers", 
                            text = ~hover.string,
                            hoverinfo = "text",
                            marker = list(color = pt.color,
                                          size = pt.size))

    if (depline & data.type %in% c("crispr", "rnai")) {
      dline <- .vline(x = -1, dash = "longdash", width = 1)
    } else {
      dline <- NULL
    }

    fig <- fig %>% layout(showlegend = FALSE,
                          shapes = list(dline),
                          xaxis = ax,
                          yaxis = ay)

    fig %>%
      config(edits = list(annotationPosition = TRUE,
                          annotationTail = TRUE),
             toImageButtonOptions = list(format = "svg"),
             displaylogo = FALSE,
             plotGlPixelRatio = 7)
  } else {
    .plot_gene_not_found(gene)
  }
}

#' Plot text indicating that gene was not found
#' 
#' @param gene Character scalar of gene identifier.
#' @author Jared Andrews
#' @rdname INTERNAL_plot_gene_not_found
.plot_gene_not_found <- function(gene) {
  # Just plots text for when a gene isn't found in depmap.
  fig <- plot_ly()
  fig <- fig %>%
    add_trace(
      mode = "text",
      text = paste0(gene, " not found in DepMap."),
      type = "scattergl",
      textfont = list(
        size = 20
      ),
      x = 2,
      y = 2
    )
  
  fig <- fig %>%
    layout(
      xaxis = list(
        range = c(0, 4),
        showline = FALSE,
        zeroline = FALSE,
        showgrid = FALSE,
        showticklabels = FALSE
      ),
      yaxis = list(
        range = c(0, 4),
        showline = FALSE,
        zeroline = FALSE,
        showgrid = FALSE,
        showticklabels = FALSE
      )
    )
  
  fig %>% style(hoverinfo = 'none')
}
