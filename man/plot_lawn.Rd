% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_gene_sgrna.R
\name{plot_lawn}
\alias{plot_lawn}
\title{Create an interactive lawn plot}
\usage{
plot_lawn(
  res,
  ylim = 5,
  fc.thresh = 0.5,
  hover.info = NULL,
  sig.line = TRUE,
  h.id = "crispr",
  feat.term = "rows",
  x.term = "RandomIndex",
  sig.term = "fdr",
  lfc.term = "LFC",
  down.color = "#0026ff",
  up.color = "#ff0000",
  insig.color = "#A6A6A6",
  sig.thresh = 0.05,
  fs = NULL,
  sig.size = 6,
  insig.size = 5,
  sig.opacity = 1,
  insig.opacity = 0.5,
  label.size = 10,
  webgl = TRUE,
  webgl.ratio = 7,
  show.counts = TRUE,
  show.hl.counts = TRUE,
  counts.size = TRUE,
  highlight.featsets = NULL,
  highlight.feats = NULL,
  featsets = NULL,
  highlight.feats.color = "#E69F00",
  highlight.feats.size = 7,
  highlight.feats.opac = 1,
  highlight.feats.linecolor = "#000000",
  highlight.feats.linewidth = 1,
  highlight.feats.label = TRUE,
  highlight.featsets.color = "#009E73",
  highlight.featsets.size = 7,
  highlight.featsets.opac = 1,
  highlight.featsets.linecolor = "#000000",
  highlight.featsets.linewidth = 1,
  highlight.featsets.label = FALSE,
  h.id.suffix = "_lawn"
)
}
\arguments{
\item{res}{Dataframe containing, at minimum, significance values and a term by which the x-axis can be ordered.}

\item{ylim}{Positive numeric scalar indicating y-axis limits.
The negative value will be used for the lower limit.}

\item{fc.thresh}{Numeric scalar indicating the fold change threshold for coloring
significant features.}

\item{hover.info}{Character vector indicating which additional columns from \code{res}
to include in the hover info.}

\item{sig.line}{Logical indicating whether to add a significance threshold line to the plot.}

\item{h.id}{Character scalar indicating the unique ID of the plotly object. Can usually be ignored,
but should be used if multiple plots are being created in the same R session (e.g. Shiny app).}

\item{feat.term}{Character scalar indicating the column name of the feature IDs in \code{res}.}

\item{x.term}{Character scalar for the x-axis term from \code{res} to be plotted.}

\item{sig.term}{Character scalar indicating the column name of the significance values in \code{res}.}

\item{lfc.term}{Character scalar indicating the column name of the log fold change values in \code{res}.}

\item{down.color}{Character scalar indicating the color of down-regulated features.}

\item{up.color}{Character scalar indicating the color of up-regulated features.}

\item{insig.color}{Character scalar indicating the color of insignificant features.}

\item{sig.thresh}{Numeric scalar indicating the significance threshold for coloring significant features.}

\item{fs}{Dataframe containing coordinates and label information for points that should be labeled.
Columns should be:
\itemize{
\item \code{x} - x coordinate of the point
\item \code{y} - y coordinate of the point
\item \code{customdata} - label to be displayed
}}

\item{sig.size}{Numeric scalar indicating the size of significant feature points.}

\item{insig.size}{Numeric scalar indicating the size of insignificant feature points.}

\item{sig.opacity}{Numeric scalar indicating the opacity of significant feature points.}

\item{insig.opacity}{Numeric scalar indicating the opacity of insignificant feature points.}

\item{label.size}{Numeric scalar indicating the size of feature labels.}

\item{webgl}{Logical indicating whether to use WebGL for rendering the plot.}

\item{webgl.ratio}{Numeric scalar indicating the ratio of WebGL to HTML5 canvas rendering, increases resolution
of saved plot when WebGL plotting is not used.}

\item{show.counts}{Logical indicating whether to show annotations for the number of features in the plot.}

\item{show.hl.counts}{Logical indicating whether to show annotations for the number of highlighted features in the plot.}

\item{counts.size}{Numeric scalar indicating the size of the feature counts labels.}

\item{highlight.featsets}{Character vector indicating which feature sets should be highlighted.}

\item{highlight.feats}{Character vector indicating which features should be highlighted.}

\item{featsets}{Named list of feature sets to be used for highlighting.}

\item{highlight.feats.color}{Character scalar indicating the color of highlighted features.}

\item{highlight.feats.size}{Numeric scalar indicating the size of highlighted features.}

\item{highlight.feats.opac}{Numeric scalar indicating the opacity of highlighted features.}

\item{highlight.feats.linecolor}{Character scalar indicating the line color of highlighted features.}

\item{highlight.feats.linewidth}{Numeric scalar indicating the line width of highlighted features.}

\item{highlight.feats.label}{Logical indicating whether to label highlighted features.}

\item{highlight.featsets.color}{Character scalar indicating the color of highlighted feature sets.}

\item{highlight.featsets.size}{Numeric scalar indicating the point size of highlighted feature sets.}

\item{highlight.featsets.opac}{Numeric scalar indicating the opacity of highlighted feature sets.}

\item{highlight.featsets.linecolor}{Character scalar indicating the line color of highlighted feature sets.}

\item{highlight.featsets.linewidth}{Numeric scalar indicating the line width of highlighted feature sets.}

\item{highlight.featsets.label}{Logical indicating whether to label highlighted feature sets.}

\item{h.id.suffix}{Character scalar indicating the suffix to be added to the plotly object ID.}
}
\value{
An interactive plotly rank plot.
}
\description{
Create an interactive lawn plot for data with significance values.
Typically, this plot is randomly ordered along the x-axis,
but the user is free to order it by any term in \code{res} that they'd like.
}
\examples{
library(CRISPRball)
d1.genes <- read.delim(system.file("extdata", "esc1.gene_summary.txt",
    package = "CRISPRball"
), check.names = FALSE)
plot.df <- gene_ingress(d1.genes,
    sig.thresh = 0.05, es.thresh = 0.5,
    es.col = "LFC", sig.col = "fdr"
)
plot_lawn(plot.df, feat.term = "id")
}
\author{
Jared Andrews
}
