% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{depmap_22q1_TPM}
\alias{depmap_22q1_TPM}
\title{DepMap expression data}
\format{
\subsection{\code{depmap_22q1_TPM}}{

A dataframe with 1393 rows and 9 columns:
\describe{
\item{depmap_id}{Depmap cell line ID}
\item{rna_expression}{log2(TPM+1) expression}
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
data(depmap_22q1_TPM)
}
\value{
A dataframe containing DepMap 22Q1 expression data.
}
\description{
The DepMap expression data for CDK2.
}
\keyword{datasets}
