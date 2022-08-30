tab_about <- tabPanel(
      title = "About",
      id = "about",
      fluidRow(
        column(width = 12,
          h2("About CRISPRball"),
          hr(),
          HTML("<p>CRISPRball was developed by <a href='https://github.com/j-andrews7' target=_blank>Jared Andrews</a> ",
            "in the Department of Developmental Neurobiology and <a href='https://github.com/jake-steele' target=_blank>Jake Steele</a> ",
            "in the Center for Advanced Genome Engineering (CAGE) at St. Jude Children's Research Hospital.</p>"),
          br(),
          HTML("<p>CRISPRball is released under the <a href='https://github.com/j-andrews7' target=_blank>MIT license</a> and should ",
            "be used only for research purposes. The CRISPRball package is developed and available on ",
            "<a href='https://github.com/j-andrews7/CRISPRball' target=_blank>Github</a>.</p>"),
        )  
      )
    )