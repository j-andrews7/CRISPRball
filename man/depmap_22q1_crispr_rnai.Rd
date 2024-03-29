% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{depmap_22q1_crispr_rnai}
\alias{depmap_22q1_crispr_rnai}
\title{DepMap CRISPR & RNAi screen data}
\format{
\subsection{\code{depmap_22q1_crispr_rnai}}{

A dataframe with 1782 rows and 9 columns:
\describe{
\item{depmap_id}{Depmap cell line ID}
\item{dependency}{Depmap dependency score}
\item{entrez_id}{Numeric entrez gene ID}
\item{gene_name}{Gene symbol}
\item{dataset}{Screen type, either CRISPR or RNAi}
\item{cell_line_name}{Cell line name}
\item{primary_disease}{Cell line disease origin}
\item{lineage}{Cell line lineage}
\item{lineage_subtype}{Cell line lineage subtype, usually a malignancy}
\item{hover.string}{String for hover text when plotting}
}
}
}
\source{
\url{https://depmap.org/portal/download/all/}

\url{https://bioconductor.org/packages/release/data/experiment/html/depmap.html}
}
\usage{
data(depmap_22q1_crispr_rnai)
}
\value{
A dataframe containing DepMap 22Q1 CRISPR and RNAi data.
}
\description{
The DepMap CRISPR & RNAi screen data for CDK2.
}
\keyword{datasets}
