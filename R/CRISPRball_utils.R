.utils.js <- "
shinyjs.disableTab = function(name) {
  var tab = $('.nav.navbar-nav li a[data-value=\"' + name + '\"]');
  tab.bind('click.tab', function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass('disabled');
}

shinyjs.enableTab = function(name) {
  var tab = $('.nav.navbar-nav li a[data-value=\"' + name + '\"]');
  tab.unbind('click.tab');
  tab.removeClass('disabled');
}
"

# .count_norm_ingress <- function(fileList) {
#  out <- read.delim(fileList$datapath)
#  # names(out) <- sapply(fileList$name, gsub, pattern='.count_normalized.txt', replacement='', fixed=TRUE)
#  # browser()
#  return(out)
# }
# 
# .count_summ_ingress <- function(fileList) {
#  out <- unlist(lapply(fileList$datapath, read.table, sep = '\t', header = TRUE))
#  # names(out) <- sapply(fileList$name, gsub, pattern='.countsummary.txt', replacement='', fixed=TRUE)
#  return(out)
# }

# TODO: These still need testing/tweaks.
.gene_summ_ingress <- function(fileList) {
  out <- lapply(fileList$datapath, read.delim, check.names = FALSE)
  names(out) <- sapply(fileList$name, gsub, pattern='.gene_summary.txt', replacement='', fixed=TRUE)
  return(out)
}

.sgrna_summ_ingress <- function(fileList) {
  out <- lapply(fileList$datapath, read.delim, check.names = FALSE)
  names(out) <- sapply(fileList$name, gsub, pattern='.sgrna_summary.txt', replacement='', fixed=TRUE)
  return(out)
}

# Generate easier columns for plotting for various data summaries.
.gene_ingress <- function(df, sig.thresh, lfc.thresh, positive.ctrl.genes = NULL, essential.genes = NULL, depmap.genes = NULL) {

  if (!is.null(essential.genes)) {
    df$essential <- df$id %in% essential.genes
  }

  if (!is.null(positive.ctrl.genes)) {
    df$Positive_Control <- df$id %in% positive.ctrl.genes
  }

  if (!is.null(depmap.genes)) {
    df$DepMap_CRISPR_Essential <- df$id %in%
      depmap.genes$gene_name[depmap.genes$dataset %in%
                               c("Chronos_Combined", "Chronos_Score", "Chronos_Achilles") &
                               depmap.genes$common_essential == TRUE]

    df$DepMap_CRISPR_Selective <- df$id %in%
      depmap.genes$gene_name[depmap.genes$dataset %in%
                               c("Chronos_Combined", "Chronos_Score", "Chronos_Achilles") &
                               depmap.genes$strongly_selective == TRUE]

    df$DepMap_RNAi_Essential <- df$id %in% depmap.genes$gene_name[depmap.genes$dataset == "RNAi_merged" &
                                                                    depmap.genes$common_essential == TRUE]
    df$DepMap_RNAi_Selective <- df$id %in% depmap.genes$gene_name[depmap.genes$dataset == "RNAi_merged" &
                                                                    depmap.genes$strongly_selective == TRUE]
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
