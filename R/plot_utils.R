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

  # LFC filter.
  if(fc.thresh > 0) {
    fc.threshed <- abs(res[[lfc.term]]) < fc.thresh
    res$col[fc.threshed] <- insig.color
    res$cex[fc.threshed] <- insig.size
    res$order[fc.threshed] <- 0
  }

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

  # LFC filter.
  if(y.thresh > 0) {
    y.threshed <- abs(df[[y.term]]) < y.thresh
    df$col[y.threshed] <- insig.color
    df$cex[y.threshed] <- insig.size
    df$order[y.threshed] <- 0
  }

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

  # LFC filter.
  if(fc.thresh > 0) {
    fc.threshed <- abs(res[[lfc.term]]) < fc.thresh
    res$col[fc.threshed] <- insig.color
    res$cex[fc.threshed] <- insig.size
    res$order[fc.threshed] <- 0
  }

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
