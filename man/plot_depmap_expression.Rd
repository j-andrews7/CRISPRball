% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_depmap.R
\name{plot_depmap_expression}
\alias{plot_depmap_expression}
\title{Plot gene expression information from DepMap, mostly from CCLE}
\usage{
plot_depmap_expression(df, color = "#7B8CB2", plot.grid = FALSE)
}
\arguments{
\item{df}{data.frame containing information for a single gene
as returned by \code{\link{get_depmap_plot_data}}.}

\item{color}{Character scalar for trace color.}

\item{plot.grid}{Boolean indicating whether to plot gridlines.}
}
\value{
plotly object
}
\description{
Plot gene expression information from DepMap, mostly from CCLE
}
\examples{
library(CRISPRball)
data(depmap_22q1_TPM)
plot_depmap_expression(depmap_22q1_TPM)
}
\seealso{
\code{\link{get_depmap_plot_data}}
}
\author{
Jared Andrews
}
