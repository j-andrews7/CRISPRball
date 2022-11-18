#' Create a tabPanel for the Gene tab
#' 
#' Create a \code{\link{tabPanel}} with UI elements for the Gene tab.
#' 
#' @param dataset.choices A character vector containing the dataset names used for selecting a dataset for plotting.
#' @param genesets A named list containing character vectors of genes, used in selecting pre-chosen sets of genes to be highlighted.
#' 
#' @return
#' A \code{\link{tabPanel}} with UI elements for the Gene tab.
#' 
#' @author Jared Andrews
#' 
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom shinyBS tipify
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyjqui jqui_resizable
#' @importFrom InteractiveComplexHeatmap InteractiveComplexHeatmapOutput
#' @importFrom plotly plotlyOutput
#' @importFrom colourpicker colourInput
#' @importFrom shinyBS tipify popify bsCollapse bsCollapsePanel
#' 
#' @rdname INTERNAL_create_tab_gene

.create_tab_gene <- function(dataset.choices, genesets) {
    tabPanel(
      title = "Gene (Overview)",
      id = "gene-overview",
      sidebarLayout(
       sidebarPanel(
         width = 2,
         h4("Plot Controls"),
         hr(),
         div(
           fluidRow(
             column(6,
                    tipify(selectizeInput("gene.sel1", "Dataset 1:", choices = dataset.choices),
                           "Dataset shown in top row.", "right", options = list(container = "body")),
                    numericInput("gene.fdr.th", "FDR threshold:",
                                 min = 0, max = 1, step = 0.01, value = 0.05)
             ),
             column(6,
                    tipify(selectizeInput("gene.sel2", "Dataset 2:", choices = dataset.choices),
                           "Dataset shown in bottom row.", "right", options = list(container = "body")),
                    numericInput("gene.lfc.th", "log2FC threshold:",
                                 min = 0, max = Inf, step = 0.05, value = 0.5)
             )
           ),
           fluidRow(
             column(12,
                    tipify(prettyCheckbox("rem.ess", label = "Remove essential genes", value = FALSE,
                                          animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                           "Remove essential genes if any provided to function.", "right", options = list(container = "body")),
                    tipify(prettyCheckbox("dep.crispr.ess", label = "Remove DepMap CRISPR essential genes", value = FALSE,
                                          animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                           "Remove DepMap Chronos Combined, Score, and Achilles common essential genes from latest release.",
                           "right", options = list(container = "body")),
                    tipify(prettyCheckbox("dep.rnai.ess", label = "Remove DepMap RNAi essential genes", value = FALSE,
                                          animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                           "Remove RNAi common essential genes from latest DepMap release.", "right", options = list(container = "body")),
                    tipify(prettyCheckbox("dep.crispr.sel", label = "Remove DepMap CRISPR selective genes", value = FALSE,
                                          animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                           "Remove DepMap Chronos Combined, Score, and Achilles strongly selective genes from latest release.", "right",
                           options = list(container = "body")),
                    tipify(prettyCheckbox("dep.rnai.sel", label = "Remove DepMap RNAi selective genes", value = FALSE,
                                          animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                           "Remove DepMap RNAi strongly selective genes from latest release.", "right", options = list(container = "body")),
                    tipify(prettyCheckbox("rem.pos", label = "Remove positive control genes", value = FALSE,
                                          animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                           "Remove positive control genes if provided to function.", "right", options = list(container = "body")),
                    tipify(prettyCheckbox("highlight.common", label = "Highlight common hits", value = FALSE,
                                          animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                           "Highlight common hits between datasets.", "right", options = list(container = "body"))
             )
           ),
           style = "background-color: #FFFFFF; padding: 3px; margin-bottom: 3px; border: 1px solid #bce8f1; "),
         bsCollapse(open = "com.settings",
                    bsCollapsePanel(
                      title = span(icon("plus"), "Common Plot Settings"), value = "com.settings", style = "info",
                      fluidRow(
                        column(width = 6,
                               tipify(colourInput("down.color", "Down colour", value = "#0026ff"),
                                      "Color of negatively selected genes.", "right", options = list(container = "body")),
                               tipify(colourInput("up.color", "Up colour", value = "red"),
                                      "Color of positively selected genes.", "right", options = list(container = "body")),
                               tipify(colourInput("insig.color", "Insig colour", value = "#A6A6A6"),
                                      "Color of insignificant genes.", "right", options = list(container = "body")),
                               tipify(numericInput("sig.opa", label = "Sig opacity:", value = 1, step = 0.05, min = 0),
                                      "Opacity of significantly selected genes.", "right", options = list(container = "body"))
                        ),
                        column(width = 6,
                               tipify(numericInput("sig.size", label = "Sig pt size:", value = 6, step = 0.1, min = 0),
                                      "Point size of significantly selected genes.", "right", options = list(container = "body")),
                               tipify(numericInput("lab.size", label = "Label size:", value = 10, step = 0.5, min = 1),
                                      "Label size of significantly selected genes.", "right", options = list(container = "body")),
                               tipify(numericInput("insig.opa", label = "Insig opacity:", value = 0.5, step = 0.05, min = 0),
                                      "Opacity of insignificant genes.", "right", options = list(container = "body")),
                               tipify(numericInput("insig.size", label = "Insig pt size:", value = 5, step = 0.1, min = 0),
                                      "Point size of insignificant genes.", "right", options = list(container = "body"))
                        )
                      ),
                      splitLayout(
                        tipify(prettyCheckbox("counts", label = "Show counts", TRUE, bigger = TRUE,
                                              animation = "smooth", status = "success",
                                              icon = icon("check"), width = "100%"),
                               "Show hit counts on plot.", "right", options = list(container = "body")),
                        tipify(prettyCheckbox("webgl", label = "Use webGL", TRUE, bigger = TRUE,
                                              animation = "smooth", status = "success",
                                              icon = icon("check"), width = "100%"),
                               "Use webGL for plot generation (faster to update, sometimes has visual artifacts).", "right", options = list(container = "body"))
                      ),
                      splitLayout(
                        tipify(prettyCheckbox("hl.counts", label = "Show highlight counts", FALSE, bigger = TRUE,
                                              animation = "smooth", status = "success",
                                              icon = icon("check"), width = "100%"),
                               "Show highlighted gene counts on plot.", "right", options = list(container = "body"))
                      ),
                      splitLayout(
                        tipify(numericInput("counts.size", label = "Counts size:", value = 8, step = 0.1, min = 0),
                               "Size of counts text.", "right", options = list(container = "body")),
                        tipify(numericInput("webgl.ratio", label = "webGL pixel ratio:", value = 7, step = 0.1, min = 1),
                               "webGL rasterization ratio. Recommend leaving this alone for high-res rasterization.", "right", options = list(container = "body"))
                      )
                    ),
                    bsCollapsePanel(
                      title = span(icon("plus"), "Volcano Plot Settings"), value = "vol.settings", style = "info",
                      fluidRow(
                        column(width = 6,
                               numericInput("vol.x", label = "x-axis limits:", value = 5, step = 0.1, min = 0.1)
                        ),
                        column(width = 6,
                               numericInput("vol.y", label = "y-axis limits:", value = 5, step = 0.5, min = 1)
                        )
                      ),
                      splitLayout(
                        prettyCheckbox("vol.fcline", label = "Show FC threshold", value = TRUE,
                                       animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                        prettyCheckbox("vol.sigline", label = "Show Sig. threshold", value = TRUE,
                                       animation = "smooth", status = "success", bigger = TRUE, icon = icon("check"))
                      ),
                      div(actionButton("vol.update", "Update Volcano Plots"), align = "center")
                    ),
                    bsCollapsePanel(
                      title = span(icon("plus"), "Rank Plot Settings"), value = "rank.settings", style = "info",
                      fluidRow(
                        column(width = 6,
                               numericInput("rank.y.max", label = "y-axis max:", value = 10, step = 0.5),
                               prettyCheckbox("rank.fcline", label = "Show FC threshold", value = TRUE,
                                              animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                        ),
                        column(width = 6,
                               numericInput("rank.y.min", label = "y-axis min:", value = -10, step = 0.5, min = 1),
                        )
                      ),
                      div(actionButton("rank.update", "Update Rank Plots"), align = "center")
                    ),
                    bsCollapsePanel(
                      title = span(icon("plus"), "Lawn Plot Settings"), value = "lawn.settings", style = "info",
                      splitLayout(
                        prettyCheckbox("lawn.sigline", label = "Show Sig. threshold", value = TRUE,
                                       animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                        numericInput("lawn.y", label = "y-axis limits:", value = 5, step = 0.5, min = 1)
                      ),
                      div(actionButton("lawn.update", "Update Lawn Plots"), align = "center")
                    ),
                    bsCollapsePanel(title = span(icon("plus"), "Highlight Gene(sets)"), value = "highlight.settings", style = "info",
                                    tipify(textAreaInput("hl.genes", "Highlight Genes:", value = "", rows = 4,
                                                         placeholder = "Enter space, comma, or newline delimited genes"),
                                           "Genes to highlight on plots.", "right", options = list(container = "body")),
                                    tipify(pickerInput("hl.genesets", "Highlight Genesets:", choices = c("", names(genesets)),
                                                       multiple = TRUE, options = list(`live-search` = TRUE, `actions-box` = TRUE)),
                                           "If provided, genesets available to highlight on plots.", "right", options = list(container = "body")),
                                    fluidRow(
                                      column(6,
                                             tipify(numericInput("hl.genes.opa", label = "Genes opacity:", value = 1, step = 0.05, min = 0),
                                                    "Opacity of highlighted genes.", "right", options = list(container = "body")),
                                             tipify(numericInput("hl.genes.size", label = "Genes pt size:", value = 7, step = 0.1, min = 0),
                                                    "Point size of highlighted genes.", "right", options = list(container = "body")),
                                             tipify(numericInput("hl.genes.lw", label = "Genes border width:", value = 1, step = 0.05, min = 0),
                                                    "Border width of highlighted genes.", "right", options = list(container = "body")),
                                             tipify(colourInput("hl.genes.col", "Genes color:", value = "#E69F00"),
                                                    "Fill color of highlighted genes.", "right", options = list(container = "body")),
                                             tipify(colourInput("hl.genes.lcol", "Genes border:", value = "#000000"),
                                                    "Border color of highlighted genes.", "right", options = list(container = "body")))
                                      ,
                                      column(6,
                                             tipify(numericInput("hl.genesets.opa", label = "Sets opacity:", value = 1, step = 0.05, min = 0),
                                                    "Opacity of genes in highlighted geneset(s).", "right", options = list(container = "body")),
                                             tipify(numericInput("hl.genesets.size", label = "Sets pt size:", value = 7, step = 0.1, min = 0),
                                                    "Point size of genes in highlighted geneset(s).", "right", options = list(container = "body")),
                                             tipify(numericInput("hl.genesets.lw", label = "Sets border width:", value = 1, step = 0.05, min = 0),
                                                    "Border width genes in of highlighted geneset(s).", "right", options = list(container = "body")),
                                             tipify(colourInput("hl.genesets.col", "Sets color:", value = "#009E73"),
                                                    "Fill color of genes in highlighted geneset(s).", "right", options = list(container = "body")),
                                             tipify(colourInput("hl.genesets.lcol", "Sets border:", value = "#000000"),
                                                    "Border color of genes in highlighted geneset(s).", "right", options = list(container = "body")))
                                    )
                    )
         ),
         div(actionButton("gene.update", "Update Plots"), align = "center")
       ),
       mainPanel(
         width = 10,
         fluidRow(
           column(width = 4,
                  span(popify(icon("circle-info", style="font-size: 20px"), title = "Volcano Plot",
                              c("This volcano plot shows the log2 fold change on the x-axis and the -log10(FDR) value on the y-axis. ",
                                "Thresholds are adjustable. Gene labels can be added (or removed) by clicking on a point ",
                                "and can be moved by clicking and dragging the label. The plot is fully customizable with the settings on the left. ",
                                "Click and drag to zoom in. Hover over a point for additional info.",
                                "Genes with full sgRNA depletion tend to all have the same significance value, forming a shelf-like max y-axis value."),
                              placement = "bottom", trigger = "hover", options = list(container = "body")),
                       withSpinner(jqui_resizable(plotlyOutput("gene1.vol"))))
           ),
           column(width = 4,
                  span(popify(icon("circle-info", style="font-size: 20px"), title = "Rank Plot",
                              c("This rank plot shows the log2 fold change on the y-axis and the gene rank on the x-axis. ",
                                "Thresholds are adjustable. Gene labels can be added (or removed) by clicking on a point ",
                                "and can be moved by clicking and dragging the label. The plot is fully customizable with the settings on the left. ",
                                "Click and drag to zoom in. Hover over a point for additional info."),
                              placement = "bottom", trigger = "hover", options = list(container = "body")),
                       withSpinner(jqui_resizable(plotlyOutput("gene1.rank"))))
           ),
           column(width = 4,
                  span(popify(icon("circle-info", style="font-size: 20px"), title = "Lawn Plot",
                              c("This plot shows the log2 fold change on the y-axis and the genes randomly ordered on the x-axis. ",
                                "Thresholds are adjustable. Gene labels can be added (or removed) by clicking on a point ",
                                "and can be moved by clicking and dragging the label. The plot is fully customizable with the settings on the left. ",
                                "Click and drag to zoom in. Hover over a point for additional info."),
                              placement = "bottom", trigger = "hover", options = list(container = "body")),
                       withSpinner(jqui_resizable(plotlyOutput("gene1.lawn"))))
           )
         ),
         hr(),
         fluidRow(
           column(width = 4,
                  withSpinner(jqui_resizable(plotlyOutput("gene2.vol")))
           ),
           column(width = 4,
                  withSpinner(jqui_resizable(plotlyOutput("gene2.rank")))
           ),
           column(width = 4,
                  withSpinner(jqui_resizable(plotlyOutput("gene2.lawn")))
           )
         )
       )
      )
    )
}

tab_gene_summary <- tabPanel(
  title = "Gene Summary Tables",
  id = "gene-summ",
  br(),
  div(DT::dataTableOutput("gene1.summary"), style = "font-size:80%;"),
  br(),
  div(DT::dataTableOutput("gene2.summary"), style = "font-size:80%;")
)