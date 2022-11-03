.create_tab_depmap <- function(sgrna.data, depmap.meta) {
  tabPanel(
    title = "DepMap",
    id = "depmap",
    sidebarLayout(
      sidebarPanel(
        width = 2,
        h4("Plot Controls"),
        hr(),
        div(
          fluidRow(
            column(12,
                  pickerInput("depmap.gene", "Choose gene:", choices = unique(c(sgrna.data[[1]]$Gene)),
                              multiple = FALSE, options = list(`live-search` = TRUE, `actions-box` = TRUE))
            )
          ),
          style = "background-color: #FFFFFF; padding: 3px; margin-bottom: 3px; border: 1px solid #bce8f1; "),
        bsCollapse(
          open = "dm.dep.settings",
          bsCollapsePanel(
            title = span(icon("plus"), "Dependency Plot Settings"), value = "dm.dep.settings", style = "info",
            fluidRow(
              column(width = 6,
                      tipify(colourInput("dep.crispr.color", "CRISPR color:", value = "#3584B5"),
                            "Fill color of CRISPR rug and density plots.", "right", options = list(container = "body")),
                      tipify(prettyCheckbox("dep.plot.grid", label = "Show gridlines", value = TRUE,
                            animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                            "Plot gridlines.", "right", options = list(container = "body"))
              ),
              column(width = 6,
                      tipify(colourInput("dep.rnai.color", "RNAi color:", value = "#52288E"),
                            "Fill color of RNAi rug and density plots.", "right", options = list(container = "body")),
                      tipify(prettyCheckbox("dep.depline", label = "Show dep threshold", value = TRUE,
                                    animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                                    "Plot vertical dependency line threshold.", "right", options = list(container = "body"))
              )
            ),
            div(actionButton("dm.dep.update", "Update Dependency Plot"), align = "center")
          ),
          bsCollapsePanel(
            title = span(icon("plus"), "Expression Plot Settings"), value = "dm.exp.settings", style = "info",
            fluidRow(
              column(width = 12,
                      tipify(prettyCheckbox("exp.plot.grid", label = "Show gridlines", value = TRUE,
                            animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                            "Plot gridlines.", "right", options = list(container = "body")),
                      tipify(colourInput("exp.color", "Expression color:", value = "#7B8CB2"),
                            "Fill color of expression rug and density plot.", "right", options = list(container = "body"))
              )
            ),
            div(actionButton("dm.exp.update", "Update Expression Plot"), align = "center")
          ),
          bsCollapsePanel(
            title = span(icon("plus"), "Copy Number Plot Settings"), value = "dm.cn.settings", style = "info",
            fluidRow(
              column(width = 12,
                      tipify(prettyCheckbox("cn.plot.grid", label = "Show gridlines", value = TRUE,
                              animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                              "Plot gridlines.", "right", options = list(container = "body")),
                      tipify(colourInput("cn.color", "Copy number color:", value = "#CEA3CB"),
                            "Fill color of copy number rug and density plot.", "right", options = list(container = "body"))
              )
            ),
            div(actionButton("dm.cn.update", "Update Copy Number Plot"), align = "center")
          ),
          bsCollapsePanel(
            title = span(icon("plus"), "Characterization Plot Settings"), value = "dm.char.settings", style = "info",
            fluidRow(
              column(width = 6,
                    pickerInput("lin.group", "Group by:", 
                                choices = c("lineage", "primary_disease", "lineage_subtype"),
                                multiple = FALSE),
                    tipify(colourInput("lin.box.color", "Boxplot line color:", value = "#000000"),
                            "Color of boxplot lines.", "right", options = list(container = "body")),
                    numericInput("lin.pt.size", "Point size:",
                                  min = 0, step = 0.1, value = 5),
                    tipify(numericInput("lin.label.size", "Label font size:",
                                        min = 0, step = 0.1, value = 12),
                            "Font size of labels. Useful for fitting more labels into subtype plots", 
                              "right", options = list(container = "body")),
                    tipify(prettyCheckbox("lin.depline", label = "Show dep threshold", value = TRUE,
                            animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                            "Plot vertical dependency line threshold (if 'crispr' or 'rnai' dataset selected).", 
                            "right", options = list(container = "body"))
              ),
              column(6,
                    pickerInput("lin.data", "Choose dataset:", 
                                choices = c("crispr", "rnai", "cn", "ccle_tpm"),
                                multiple = FALSE),
                    tipify(colourInput("lin.box.fill", "Boxplot fill color:", value = "#E2E2E2"),
                            "Color of boxplot fill.", "right", options = list(container = "body")),
                    tipify(colourInput("lin.pt.color", "Point color:", value = "#56B4E9"),
                            "Color of points.", "right", options = list(container = "body")),
              )
            ),
            div(actionButton("dm.lineage.update", "Update Lineage Plot"), align = "center")
          ),
          bsCollapsePanel(
            title = span(icon("plus"), "Characterization Subplot Settings"), value = "dm.subchar.settings", style = "info",
            fluidRow(
              column(width = 6,
                    tipify(selectInput("sub.lineage", "Lineage:", choices = unique(depmap.meta$lineage)),
                            "Choose lineage for which to plot sublineages.", "right", options = list(container = "body")),
                    numericInput("sub.pt.size", "Point size:",
                                  min = 0, step = 0.1, value = 5),
                    tipify(numericInput("sub.label.size", "Label font size:",
                                        min = 0, step = 0.1, value = 12),
                            "Font size of labels. Useful for fitting more labels into subtype plots", 
                              "right", options = list(container = "body")),
                    tipify(prettyCheckbox("sub.depline", label = "Show dep threshold", value = TRUE,
                            animation = "smooth", status = "success", bigger = TRUE, icon = icon("check")),
                            "Plot vertical dependency line threshold (if 'crispr' or 'rnai' dataset selected).", 
                            "right", options = list(container = "body"))
              ),
              column(6,
                    tipify(colourInput("sub.box.color", "Boxplot line color:", value = "#000000"),
                            "Color of boxplot lines.", "right", options = list(container = "body")),
                    tipify(colourInput("sub.box.fill", "Boxplot fill color:", value = "#E2E2E2"),
                            "Color of boxplot fill.", "right", options = list(container = "body")),
                    tipify(colourInput("sub.pt.color", "Point color:", value = "#56B4E9"),
                            "Color of points.", "right", options = list(container = "body")),
              )
            ),
            div(actionButton("dm.sublineage.update", "Update Sublineage Plot"), align = "center")
          )
        )
      ),
      mainPanel(
        width = 10,
        fluidRow(
          column(width = 4,
                span(h3("Dependent Cell Lines", popify(icon("circle-info", style="font-size: 20px"), "Dependent Cell Lines",
                                                        c("This plot shows DepMap dependency scores for the selected gene. ",
                                                          "A cell line is considered dependent if it has a probability of dependency ",
                                                          "greater than 50%. <br><br>",
                                                          "Probabilities of dependency are calculated for each gene score in a cell ",
                                                          "line as the probability that score arises from the distribution of essential ",
                                                          "gene scores rather than nonessential gene scores.",
                                                          "See the <a href=https://www.biorxiv.org/content/10.1101/720243v1>DepMap preprint</a> ",
                                                          "for more info. <br><br>",
                                                          "<b>Gene Effect</b><br>",
                                                          "Outcome from <a href=https://www.nature.com/articles/s41467-018-06916-5>DEMETER2</a>",
                                                          " or <a href=https://www.biorxiv.org/content/10.1101/2021.02.25.432728v1>Chronos</a>. ",
                                                          "A lower score means that a gene is more likely to be dependent in a given cell line. ",
                                                          "A score of 0 is equivalent to a gene that is not essential whereas a score of -1 corresponds ",
                                                          "to the median of all common essential genes."),
                                                        placement = "bottom", trigger = c("hover", "click"), options = list(container = "body")), .noWS="outside"),
                      uiOutput("depmap.deplines"),
                      withSpinner(jqui_resizable(plotlyOutput("depmap.essplot", height = 250)))),
        
                span(h3("Expression", popify(icon("circle-info", style="font-size: 20px"), "Gene Expression",
                                              c("RNASeq files are aligned with STAR and quantified with RSEM, then TPM-normalized. ",
                                                "Reported values are log2(TPM+1)."),
                                              placement = "bottom", trigger = "hover", options = list(container = "body")), .noWS="outside"),
                      withSpinner(jqui_resizable(plotlyOutput("depmap.expplot", height = 200)))),
                span(h3("Copy Number", popify(icon("circle-info", style="font-size: 20px"), "Copy Number",
                                              c("The <a href=https://forum.depmap.org/t/what-is-relative-copy-number-copy-number-ratio/104/2 target=_blank>relative ",
                                                "copy number</a> pipeline used varies by cell line. For around 1000 lines, Sanger WES data ",
                                                "was used, while for around 700 lines, Broad WES data was used. The remaining lines use SNP ",
                                                "array data as explained in <a href=https://doi.org/10.1038/s41586-019-1186-3 target=_blank rel=noopener>",
                                                "10.1038/s41586-019-1186-3</a>. See <a href=https://doi.org/10.1101/720243 target=_blank ",
                                                "rel=noopener>10.1101/720243</a> for details on how CN source is chosen per line. Lines with ",
                                                "WES data were processed through GATK using PONs from TCGA without matched normals and transformed by log2(x+1)."),
                                              placement = "bottom", trigger = c("hover", "click"), options = list(container = "body")), .noWS="outside"),
                      withSpinner(jqui_resizable(plotlyOutput("depmap.cnplot", height = 200))))
          ),
          column(width = 4,
                span(h3("Characterization", popify(icon("circle-info", style="font-size: 20px"), "Characterization",
                                              c("This plot shows the distribution of values for the selected gene from the chosen dataset. ",
                                                "<br><br><b>crispr</b> will display CRISPR dependency scores from Public, Chronos. <br><br><b>rnai</b> will display ",
                                                "RNAi perturbation data from Achilles+DRIVE+Marcotte, DEMETER2. ",
                                                "<br><br><b>cn</b> will display copy number data from Sanger WES or Broad WES.",
                                                "<br><br><b>ccle_tpm</b> will display gene expression data from CCLE in log2(TPM+1).<br><br>",
                                                "The grouping can be changed between lineage, disease, and lineage subtypes."),
                                              placement = "bottom", trigger = c("hover", "click"), options = list(container = "body")), .noWS="outside"),
                      withSpinner(jqui_resizable(plotlyOutput("depmap.lineages", height = 800))))
          ),
          column(width = 4,
                span(h3("Sublineages", popify(icon("circle-info", style="font-size: 20px"), "Sublineage",
                                              c("This plot shows the distribution of values for the selected gene for the sublineages of the chosen lineage. ",
                                                "<br><br><b>crispr</b> will display CRISPR dependency scores from Public, Chronos. <br><br><b>rnai</b> will display ",
                                                "RNAi perturbation data from Achilles+DRIVE+Marcotte, DEMETER2. ",
                                                "<br><br><b>cn</b> will display copy number data from Sanger WES or Broad WES.",
                                                "<br><br><b>ccle_tpm</b> will display gene expression data from CCLE in log2(TPM+1).<br><br>",
                                                "The grouping can be changed between lineage subtypes."),
                                              placement = "bottom", trigger = c("hover", "click"), options = list(container = "body")), .noWS="outside"),
                      withSpinner(jqui_resizable(plotlyOutput("depmap.sublineage", height = 300)))),
                      br(),
                wellPanel(span(h3("Gene Info",popify(icon("circle-info", style="font-size: 20px"), "Gene Info",
                            c("Gene info and accessions."),
                            placement = "bottom", trigger = "hover", options = list(container = "body")), .noWS="outside"),
                      withSpinner(uiOutput("depmap.geneinfo"))))
          )
        )
      )
    )
  )
}