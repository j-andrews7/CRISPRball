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

#' Generate dependency summary info tagList
#' @param dep.info Named list containing summary CRISPR and RNAi info.
#' @param dep.release Character scalar for DepMap release as returned by \code{\link[depmap]{depmap_release}}.
#' @param crispr.color Character scalar for color to use for CRISPR title.
#' @param rnai.color Character scalar for color to use for RNAi title.
#' 
#' @return TagList containing dependency summary information.
.make_dependency_tag <- function(dep.info, dep.release, crispr.color, rnai.color) {  
  cinfo <- "N/A"
  if (dep.info$crispr$avail) {
    cinfo <- paste0(dep.info$crispr$dep_lines, "/", dep.info$crispr$total_lines)
  }
  
  rinfo <- "N/A"
  if (dep.info$rnai$avail) {
    rinfo <- paste0(dep.info$rnai$dep_lines, "/", dep.info$rnai$total_lines)
  }
  
  c.lab <- NULL
  r.lab <- NULL
  
  cess.mess <- c("A gene which, in a large, pan-cancer screen, ranks in the top X most ",
                 "depleting genes in at least 90% of cell lines. X is chosen empirically ",
                 "using the minimum of the distribution of gene ranks in their 90th percentile ",
                 "least depleting lines.")
  
  ssel.mess <- c("A gene whose dependency is at least 100 times more likely to have been sampled",
                 "from a skewed distribution than a normal distribution.")
  
  if (!is.null(dep.info$crispr$label)) {
    outpop <- if (dep.info$crispr$label == "COMMON ESSENTIAL") cess.mess else ssel.mess
    c.lab <- span(strong(dep.info$crispr$label), 
                  popify(icon("info-circle", style="font-size: 12px"), dep.info$crispr$label,
                        outpop, placement = "right", trigger = c("hover", "click"), options = list(container = "body")), br(), 
                  style = paste0("background: ", crispr.color,"; color: #ffffff; border-radius: 5px; padding: 3px;"))
  }
  
  if (!is.null(dep.info$rnai$label)) {
    outpop <- if (dep.info$rnai$label == "COMMON ESSENTIAL") cess.mess else ssel.mess
    r.lab <- span(strong(dep.info$rnai$label), 
                  popify(icon("info-circle", style="font-size: 12px"), dep.info$rnai$label,
                         outpop, placement = "right", trigger = c("hover", "click"), options = list(container = "body")), br(), 
                  style = paste0("background: ", rnai.color, "; color: #ffffff; border-radius: 5px; padding: 3px;"))
  }
  
  out <- tagList(div(span(strong(paste0("CRISPR (DepMap ", dep.release, ", ", dep.info$crispr$dataset, "): ", cinfo)), 
                          style = paste0("color: ", crispr.color, ";")), style = "margin-bottom: 7px;"),
                 c.lab,
                 div(span(strong(paste0("RNAi (DepMap ", dep.release, ", ", dep.info$rnai$dataset, "): ", rinfo)), 
                          style = paste0("color: ", rnai.color, ";")), style = "margin-bottom: 7px; margin-top: 8px"),
                 r.lab)
  
  return(out)
}

.error_if_no_mygene <- function() {
  if (!requireNamespace("mygene", quietly = TRUE)) {
    stop("'mygene' installation required to display gene information.")
  }
}

#' Generate gene tagList via mygene API
#' @param gene Character scalar for gene symbol.
#' 
#' @return TagList containing dependency summary information.
.make_gene_tag <- function(gene) { 
  .error_if_no_mygene()
  info <- mygene::query(gene, fields = "all", size = 1)

  if (length(info) > 0) {
    info <- info$hits
    out <- tagList(splitLayout(span(strong("Gene: "), info$symbol), span(strong("Aliases: "), paste0(unlist(info$alias), collapse = ", "))),
                   splitLayout(span(strong("Position: "), paste0(info$genomic_pos$chr, ":", info$genomic_pos$start, "-", info$genomic_pos$end)),
                               span(strong("Gene type: "), info$type_of_gene)),
                   splitLayout(span(strong("Entrez: "), a(info$entrezgene, href = paste0("https://www.ncbi.nlm.nih.gov/gene/", info$entrezgene))),
                               span(strong("Ensembl: "), a(info$ensembl$gene, href = paste0("http://www.ensembl.org/id/", info$ensembl$gene)))),
                   div(br(), span(strong("Summary: "), info$summary))
                   )
  } else {
    out <- tagList(div(span("Unable to find gene information.")))
  }

  return(out)
}
