.tab_sgrna <- function(sgrna.choices, sgrna.gene) {
    tabPanel(
      title = "sgRNA",
      id = "sgrna",
      sidebarLayout(
       sidebarPanel(
         width = 2,
         h4("Plot Controls"),
         hr(),
         div(
           fluidRow(
             column(6,
                    tipify(selectizeInput("sgrna.sel1", "Dataset 1:", choices = sgrna.choices),
                           "Dataset shown in top row.", "right", options = list(container = "body"))
             ),
             column(6,
                    tipify(selectizeInput("sgrna.sel2", "Dataset 2:", choices = sgrna.choices),
                           "Dataset shown in bottom row.", "right", options = list(container = "body"))
             )
           ),
           fluidRow(
             column(12,
                    pickerInput("sgrna.gene", "Choose gene:", choices = sgrna.gene,
                                multiple = FALSE, options = list(`live-search` = TRUE, `actions-box` = TRUE))
             )
           ),
           style = "background-color: #FFFFFF; padding: 3px; margin-bottom: 3px; border: 1px solid #bce8f1; "),
       ),
       mainPanel(
         width = 10,
         fluidRow(
           column(width = 2,
                  span(popify(icon("circle-info", style="font-size: 20px"), title = "Counts Plot",
                              c("This rank plot shows the normalized counts for each sgRNA for the ",
                                "specified gene across the samples that make up the dataset comparision. ",
                                "Hover over a point for additional info."),
                              placement = "bottom", trigger = "hover", options = list(container = "body")),
                       withSpinner(jqui_resizable(plotlyOutput("sgrna1.counts"))))
           ),
           column(width = 4,
                  span(popify(icon("circle-info", style="font-size: 20px"), title = "Rank Plot",
                              c("This rank plot shows the log2 fold change on the y-axis and the sgRNA rank on the x-axis. ",
                                "sgRNAs for the selected gene will be highlighted. ",
                                "Click and drag to zoom in. Hover over a point for additional info."),
                              placement = "bottom", trigger = "hover", options = list(container = "body")),
                       withSpinner(jqui_resizable(plotlyOutput("sgrna1.rank"))))
           ),
           column(width = 6,
                  jqui_resizable(div(DT::dataTableOutput("sgrna1.detail"), style = "font-size:80%;"))
           )
         ),
         hr(),
         fluidRow(
           column(width = 2,
                  withSpinner(jqui_resizable(plotlyOutput("sgrna2.counts")))
           ),
           column(width = 4,
                  withSpinner(jqui_resizable(plotlyOutput("sgrna2.rank")))
           ),
           column(width = 6,
                  jqui_resizable(div(DT::dataTableOutput("sgrna2.detail"), style = "font-size:80%;"))
           )
         )
       )
      )
    )
}