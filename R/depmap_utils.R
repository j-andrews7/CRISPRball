#' Get essential/selective gene information from DepMap summary table.
#' 
#' @param gene Character scalar for gene symbol.
#' @param depmap.summary data.frame containing DepMap summary information.
#' @return Named list containing RNAi and CRISPR named lists containing dataset information for the provided gene, 
#'   if available. If the gene is not found in the summary data.frame, the \code{avail} element for the RNAi and CRISPR lists
#'   will be set to \code{FALSE}.
#'   
#' @export
#' @author Jared Andrews  
get_depmap_essentiality <- function(gene, depmap.summary) {

  crispr = list(avail = FALSE)
  rnai = list(avail = FALSE)
    
  if (gene %in% depmap.summary$gene_name) {
    df <- depmap.summary[depmap.summary$gene_name == gene,]

    if ("RNAi_merged" %in% df$dataset) {
      rnai$avail <- TRUE
      rnai$dataset <- "RNAi_merged"
      rnai$dep_lines <- df$dependent_cell_lines[df$dataset == rnai$dataset]
      rnai$total_lines <- df$cell_lines_with_data[df$dataset == rnai$dataset]
      
      if (df$common_essential[df$gene_name == gene & df$dataset == rnai$dataset] == 1) {
        rnai$label <- "COMMON ESSENTIAL"
      } else if (df$strongly_selective[df$gene_name == gene & df$dataset == rnai$dataset] == 1) {
        rnai$label <- "STRONGLY SELECTIVE"
      } 
    }
    
    # Check for various potential CRISPR datasets.
    if ("Chronos_Combined" %in% df$dataset) {
      crispr$dataset <- "Chronos_Combined"
    } else if ("Chronos_Achilles" %in% df$dataset) {
      crispr$dataset <- "Chronos_Achilles"
    } else if ("Chronos_Score" %in% df$dataset) {
      crispr$dataset <- "Chronos_Score"
    }
    
    if (!is.null(crispr$dataset)) {
      crispr$avail <- TRUE
      crispr$dep_lines <- df$dependent_cell_lines[df$dataset == crispr$dataset]
      crispr$total_lines <- df$cell_lines_with_data[df$dataset == crispr$dataset]
      
      if (df$common_essential[df$gene_name == gene & df$dataset == crispr$dataset] == 1) {
        crispr$label <- "COMMON ESSENTIAL"
      } else if (df$strongly_selective[df$gene_name == gene & df$dataset == crispr$dataset] == 1) {
        crispr$label <- "STRONGLY SELECTIVE"
      } 
    }
  }
  
  outs <- list(crispr = crispr, rnai = rnai)
  return(outs)
}


#' Plot gene dependency information from DepMap CRISPR and RNAi tables.
#' 
#' @param gene Character scalar for gene symbol.
#' @param depmap.meta data.frame of DepMap cell line metadata, as stored in the 'meta' table 
#'   of the SQLite database built by \code{\link{build_depmap_db}}.
#' @param depmap.pool pool connection to DepMap SQLite database built with \code{\link{build_depmap_db}}.
#' @return plotly object
#'   
#' @export
#' @author Jared Andrews  
plot_depmap_dependency <- function(gene, depmap.meta, depmap.pool) {
  
  df.c <- pool::dbGetQuery(depmap.pool, 'SELECT * FROM "crispr" WHERE "gene_name" == (:x)', params = list(x = gene))
  df.r <- pool::dbGetQuery(depmap.pool, 'SELECT * FROM "rnai" WHERE "gene_name" == (:x)', params = list(x = gene))
  
  df <- NULL
  
  if (nrow(df.c) > 0) {
    df.c$dataset <- "CRISPR"
    df <- df.c
  }
  
  if (nrow(df.r) > 0) {
    df.r$dataset <- "RNAi"
    
    if(!is.null(df)) {
      df <- rbind(df, df.r)
    } else {
      df <- df.r
    }
  }
  
  df$cell_line_name <- depmap.meta$cell_line_name[df$depmap_id == depmap.meta$depmap_id]
  df$primary_disease <- depmap.meta$primary_disease[df$depmap_id == depmap.meta$depmap_id]
  df$lineage <- depmap.meta$lineage[df$depmap_id == depmap.meta$depmap_id]
  df$lineage_subtype <- depmap.meta$lineage_subtype[df$depmap_id == depmap.meta$depmap_id]
  
  gg <- ggplot() +  
    geom_density(data = df, aes(x=dependency, color=dataset, fill=dataset), alpha = 0.6) + 
    geom_rug(data = df[df$dataset == "CRISPR",], aes(x=dependency, color=dataset), outside = FALSE) + 
    geom_rug(data = df[df$dataset == "RNAi",], aes(x=dependency, color=dataset), sides = "t") + 
    ylab("") + 
    xlab("") + 
    theme_bw() + 
    scale_color_manual(values=c("#3584B5", "#52288E"), breaks = c("CRISPR", "RNAi")) + 
    scale_fill_manual(values=c("#3584B5", "#52288E"), breaks = c("CRISPR", "RNAi")) +
    geom_vline(xintercept = 0) + 
    geom_vline(xintercept = -1, color = "red", linetype = "dashed")
    
  gg <- ggplotly(gg) %>% 
    layout(xaxis = list(
             title="Gene Effect"),   
           yaxis = list(   
             title="Density"))
  
  gg %>%
    config(edits = list(annotationPosition = TRUE,
                        annotationTail = TRUE),
           toImageButtonOptions = list(format = "svg"),
           displaylogo = FALSE,
           plotGlPixelRatio = webgl.ratio)
}
