#' Create a tabPanel for the data upload tab
#'
#' Create a \code{\link{tabPanel}} with UI elements for the data upload tab.
#'
#' @return
#' A \code{\link{tabPanel}} with UI elements for the data upload tab.
#'
#' @author Jared Andrews
#'
#' @importFrom shiny wellPanel fluidRow column hr div numericInput selectizeInput h3 h2 a strong tabPanel span actionButton fileInput br hr code icon
#'
#' @rdname INTERNAL_create_tab_data_upload
.create_tab_data_upload <- function() {
    # nocov start
    tabPanel(
        title = "Data Upload & Usage",
        id = "upload",
        fluidRow(
            column(
                6,
                wellPanel(
                    div(
                        class = "white",
                        h3("Count Files"),
                        fileInput("countNormFile", "Choose Normalized Count File", multiple = FALSE, accept = ".txt"),
                        fileInput("countSummary", "Choose Count Summary File", multiple = FALSE, accept = ".txt"),
                    )
                )
            ),
            column(
                6,
                wellPanel(
                    div(
                        class = "white",
                        h3("Gene and sgRNA Summary Files"),
                        fileInput("geneSummaryFiles", "Choose Gene Summary File(s)", multiple = TRUE, accept = ".txt"),
                        fileInput("sgrnaSummaryFiles", "Choose sgRNA Summary File(s)", multiple = TRUE, accept = ".txt"),
                    )
                )
            )
        ),
        fluidRow(
            column(
                12,
                wellPanel(
                    div(
                        class = "white",
                        fluidRow(
                            column(
                                12,
                                h2("Getting Started with CRISPRball"),
                                strong("CRISPRball"), " is a Shiny/web application for visualizing MAGeCK RRA or MLE analysis results.",
                                "For more information on running the analysis, see the ",
                                a("MAGeCK Analysis User Guide.", href = "https://sourceforge.net/p/mageck/wiki/Home/"),
                                "This tool complements the ",
                                a("MAGeCKFlute package", href = "https://www.bioconductor.org/packages/release/bioc/html/MAGeCKFlute.html"),
                                " by focusing on an interactive experience and the simple production of publication-quality figures.",
                                br(),
                                br(),
                                hr()
                            ),
                            column(
                                6,
                                h3("Data Upload"),
                                span("The ", strong("Data Upload"), "tab allows you to upload a normalized counts file, a count summary file, and gene and sgRNA summary files."),
                                span(
                                    "The normalized counts and count summary files should be tab-delimited text files returned by ", code("mageck count"),
                                    " with ", code("*.countsummary.txt"), "and ", code("*.count_norm.txt"), " extensions, respectively. ",
                                    "The gene and sgRNA summary files should be tab-delimited text files with ", code("*.gene_summary.txt"), " and ", code("*.sgrna_summary.txt"), "  extensions, ",
                                    "multiple sets of which can be uploaded. The prefixes of the gene summary and sgrna summary files should match, ",
                                    "e.g. P7.gene_summary.txt and P7.sgrna_summary.txt. MAGeCK MLE outputs will be appropriately parsed upon upload into ",
                                    "individual tables for each comparison."
                                ),
                                br(),
                                br(),
                                br(),
                                h3("Plot Content, Gene(set) Highlighting, and Labeling"),
                                span(
                                    "All plots have explanations of their content that can be accessed by hovering and/or clicking the ",
                                    icon("circle-info", style = "font-size: 14px"), " icon above each plot."
                                ),
                                br(),
                                br(),
                                span(
                                    "In the", strong("Gene (Overview)"), " tab, labels can be added on click to each point in the plot. These labels can also be ",
                                    "dragged around the plot to change their position. Double clicking the plot will remove all labels. ",
                                    "Individual labels can be removed by clicking the same point. ",
                                    "The other text annotations can also be re-positioned by dragging. ",
                                    "Genes or genesets (if provided) can be highlighted on the plots using the ",
                                    strong("Highlight Gene(sets)"), " inputs in the sidebar. ",
                                    "Highlighted gene(set)s can also be automatically labeled."
                                ),
                                br(),
                                br(),
                                br(),
                                h3("Adjusting Plot Settings/Appearance"),
                                "Nearly all plots have settings that can be adjusted by changing the input in the sidebar for each tab. ",
                                "Hovering over an input field will display a tooltip with a description of the input and how it affects the plot.",
                                span("After tweaking the settings for a given, click the ", strong("Update Plot"), " button to update the plot.")
                            ),
                            column(
                                6,
                                h3("Downloading Data & Plots"),
                                span(
                                    "All data in tables can be easily downloaded via the buttons above the table. Filtering will be applied. ",
                                    "Plots can be downloaded as SVG files by clicking the ", strong("Download Plot"),
                                    " button above the plot visible on hover. Standalone, interactive versions of plots can be downloaded by clicking ",
                                    "the associated ", strong("Download Plot"), " button above/below the plot. Aesthetics and labels will be retained."
                                ),
                                br(),
                                br(),
                                br(),
                                h3("Comparing Datasets"),
                                span(
                                    "In the", strong("Comparisons"), " tab, hits from multiple datasets or lines can be compared to each other. ",
                                    "The Upset plots displays the number of hits that are shared between datasets. ",
                                    "The datasets to compare can be selected in the sidebar. ",
                                    "Clicking on the connected points in the Upset plot will display the corresponding genes in the right-side table. "
                                ),
                                br(),
                                br(),
                                br(),
                                h3("DepMap Integration"),
                                span(
                                    "If the DepMap database is provided, the ", strong("DepMap"), " tab will contain plots for the DepMap data. ",
                                    "This can be useful to explore potential ",
                                    "hits for selectivity, expression, and copy number in a given lineage or disease. ",
                                    "Common essential genes from DepMap's CRISPR and/or RNAi screens can be removed from the gene plots by ",
                                    "checking the appropriate box in the sidebar. ",
                                    br(),
                                    br(),
                                    "See the ", a("DepMap site", href = "http://depmap.org/"), " for more information, additional tools, and ",
                                    "the appropriate ", a("publication to cite", href = "http://depmap.org/publications/"), " if you utilize this data."
                                )
                            )
                        )
                    )
                )
            )
        )
    )
    # nocov end
}
