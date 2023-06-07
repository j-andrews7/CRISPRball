#' Create a tabPanel for the about tab
#'
#' Create a \code{\link{tabPanel}} with UI elements for the about tab.
#'
#' @return
#' A \code{\link{tabPanel}} with UI elements for the about tab.
#'
#' @author Jared Andrews
#'
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#'
#' @rdname INTERNAL_create_tab_about
.create_tab_about <- function() {
    # nocov start
    tabPanel(
        title = "About",
        id = "about",
        fluidRow(
            column(
                4,
                a(img(src = "logo/CRISPRball_Hex.png", height = "700"), href = "https://bioconductor.org/packages/CRISPRball"),
            ),
            column(
                8,
                h2("About CRISPRball"),
                hr(),
                HTML(
                    "<p>CRISPRball was developed by <a href='https://github.com/j-andrews7' target=_blank>Jared Andrews</a> ",
                    "in the Department of Developmental Neurobiology and <a href='https://github.com/jake-steele' target=_blank>Jake Steele</a> ",
                    "in the Center for Advanced Genome Engineering (CAGE) at St. Jude Children's Research Hospital.</p>",
                    "<p>CRISPRball is released under the <a href='https://github.com/j-andrews7/CRISPRball/blob/main/LICENSE' target=_blank>MIT license</a> and should ",
                    "be used only for research purposes. The CRISPRball package is developed and available on ",
                    "<a href='https://github.com/j-andrews7/CRISPRball' target=_blank>Github</a>.</p>"
                )
            )
        )
    )
    # nocov end
}
