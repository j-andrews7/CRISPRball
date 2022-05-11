# Generate easier columns for plotting for various data summaries.
.gene_ingress <- function(df, sig.thresh, lfc.thresh, positive.ctrl.genes = NULL, essential.genes = NULL, depmap.genes = NULL) {

  if (!is.null(essential.genes)) {
    df$essential <- df$id %in% essential.genes
  }

  if (!is.null(positive.ctrl.genes)) {
    df$Positive_Control <- df$id %in% positive.ctrl.genes
  }

  if (!is.null(depmap.genes)) {
    df$DepMap_CRISPR_Essential <- df$id %in% depmap.genes$gene_name[depmap.genes$dataset == "Chronos_Combined" &
                                                                 depmap.genes$common_essential == "True"]
    df$DepMap_CRISPR_Selective <- df$id %in% depmap.genes$gene_name[depmap.genes$dataset == "Chronos_Combined" &
                                                                 depmap.genes$strongly_selective == "True"]

    df$DepMap_RNAi_Essential <- df$id %in% depmap.genes$gene_name[depmap.genes$dataset == "RNAi_merged" &
                                                               depmap.genes$common_essential == "True"]
    df$DepMap_RNAi_Selective <- df$id %in% depmap.genes$gene_name[depmap.genes$dataset == "RNAi_merged" &
                                                               depmap.genes$strongly_selective == "True"]
  }

  df$LFC <- as.numeric(df$`neg|lfc`)

  df$RRAscore <- apply(df, 1, function(x) {
    x <- split(unname(x),names(x))
    as.numeric(ifelse(as.numeric(x$`neg|score`) < as.numeric(x$`pos|score`), x$`neg|score`, x$`pos|score`))
  })

  df$FDR <- apply(df, 1, function(x) {
    x <- split(unname(x),names(x))
    as.numeric(ifelse(as.numeric(x$`neg|fdr`) < as.numeric(x$`pos|fdr`), x$`neg|fdr`, x$`pos|fdr`))
  })

  df$hit_type <- apply(df, 1, function(x) {
    x <- split(unname(x),names(x))
    ifelse(as.numeric(x$`neg|fdr`) < sig.thresh & as.numeric(x$LFC) < -as.numeric(lfc.thresh), "neg",
           ifelse(as.numeric(x$`pos|fdr`) < sig.thresh & as.numeric(x$LFC) > as.numeric(lfc.thresh), "pos", NA))
  })

  df$pval <- apply(df, 1, function(x) {
    x <- split(unname(x),names(x))
    as.numeric(ifelse(as.numeric(x$`neg|p-value`) < as.numeric(x$`pos|p-value`), x$`neg|p-value`, x$`pos|p-value`))
  })

  df$goodsgrna <- apply(df, 1, function(x) {
    x <- split(unname(x),names(x))
    as.integer(ifelse(as.numeric(x$`neg|score`) < as.numeric(x$`pos|score`), x$`neg|goodsgrna`, x$`pos|goodsgrna`))
  })

  df$Rank <- rank(df$LFC)

  df$RandomIndex <- sample(1:nrow(df), nrow(df))

  df
}

# sgRNA pair plot.
.make_sgrna_pairplot <- function(df) {
  gene <- df$Gene[1]
  df <- data.frame(group = c(rep("control", nrow(df)), rep("treatment", nrow(df))), counts = c(df$control_count, df$treatment_count), id = rep(df$sgrna, 2))
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
