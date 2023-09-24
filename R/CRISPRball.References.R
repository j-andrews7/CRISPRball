#' CRISPRball
#'
#' @docType package
#' @name CRISPRball
#' @description This package was designed to make the visualization, interpretation, and exploration of CRISPR screen analyses 
#'   accessible for both analysts and bench scientists. It makes high-quality, publication-ready figures simple to create and customize.
#'   It also provides a simple interface to DepMap data for comparison and filtering.
#' @details Most users will be interested in the main function to create a Shiny application (\code{\link{CRISPRball}}),
#'   though there are certain plotting functions that may be of interest to developers (\code{\link{plot_volcano}}, \code{\link{plot_rank}}, \code{\link{plot_lawn}}, etc).
#'   Most included plotting functions produce a ggplot object by default and have few required arguments.
#'   Many additional arguments are available for customization to generate complex, publication-ready figures.
#' 
#' This package supplements the \href{https://www.bioconductor.org/packages/release/bioc/html/MAGeCKFlute.html}{MAGeCKFlute} package, adding
#'   additional functionality, visualizations, and a Shiny interface.
#'
#' To report bugs, suggest new features, or ask for help, the best method is to create an issue on the github, \href{https://github.com/j-andrews7/CRISPRball}{here}, 
#' or the bioconductor support site (be sure to tag 'CRISPRball' so that I get a notification!), \href{https://support.bioconductor.org/}{here}
NULL