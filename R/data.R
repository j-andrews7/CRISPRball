#' DepMap RNAi screen data
#'
#' The DepMap RNAi screen data for CDK2.
#'
#' @format ## `depmap_22q1_rnai`
#' A data frame with 712 rows and 9 columns:
#' \describe{
#'   \item{dependency}{Depmap dependency score}
#'   \item{entrez_id}{Numeric entrez gene ID}
#'   \item{gene_name}{Gene symbol}
#'   \item{depmap_id}{Depmap cell line ID}
#'   \item{cell_line_name}{Cell line name}
#'   \item{primary_disease}{Cell line disease origin}
#'   \item{lineage}{Cell line lineage}
#'   \item{lineage_subtype}{Cell line lineage subtype, usually a malignancy}
#'   \item{hover.string}{String for hover text when plotting}
#' }
#' @source <https://depmap.org/portal/download/all/>
#' @source <https://bioconductor.org/packages/release/data/experiment/html/depmap.html>
#' @usage
#' data(depmap_22q1_rnai)
"depmap_22q1_rnai"


#' DepMap CRISPR screen data
#'
#' The DepMap CRISPR screen data for CDK2.
#'
#' @format ## `depmap_22q1_crispr`
#' A data frame with 1070 rows and 9 columns:
#' \describe{
#'   \item{depmap_id}{Depmap cell line ID}
#'   \item{dependency}{Depmap dependency score}
#'   \item{entrez_id}{Numeric entrez gene ID}
#'   \item{gene_name}{Gene symbol}
#'   \item{cell_line_name}{Cell line name}
#'   \item{primary_disease}{Cell line disease origin}
#'   \item{lineage}{Cell line lineage}
#'   \item{lineage_subtype}{Cell line lineage subtype, usually a malignancy}
#'   \item{hover.string}{String for hover text when plotting}
#' }
#' @source <https://depmap.org/portal/download/all/>
#' @source <https://bioconductor.org/packages/release/data/experiment/html/depmap.html>
#' @usage
#' data(depmap_22q1_crispr)
"depmap_22q1_crispr"


#' DepMap CRISPR & RNAi screen data
#'
#' The DepMap CRISPR & RNAi screen data for CDK2.
#'
#' @format ## `depmap_22q1_crispr_rnai`
#' A data frame with 1782 rows and 9 columns:
#' \describe{
#'   \item{depmap_id}{Depmap cell line ID}
#'   \item{dependency}{Depmap dependency score}
#'   \item{entrez_id}{Numeric entrez gene ID}
#'   \item{gene_name}{Gene symbol}
#'   \item{dataset}{Screen type, either CRISPR or RNAi}
#'   \item{cell_line_name}{Cell line name}
#'   \item{primary_disease}{Cell line disease origin}
#'   \item{lineage}{Cell line lineage}
#'   \item{lineage_subtype}{Cell line lineage subtype, usually a malignancy}
#'   \item{hover.string}{String for hover text when plotting}
#' }
#' @source <https://depmap.org/portal/download/all/>
#' @source <https://bioconductor.org/packages/release/data/experiment/html/depmap.html>
#' @usage
#' data(depmap_22q1_crispr_rnai)
"depmap_22q1_crispr_rnai"


#' DepMap expression data
#'
#' The DepMap expression data for CDK2.
#'
#' @format ## `depmap_22q1_TPM`
#' A data frame with 1393 rows and 9 columns:
#' \describe{
#'   \item{depmap_id}{Depmap cell line ID}
#'   \item{rna_expression}{log2(TPM+1) expression}
#'   \item{entrez_id}{Numeric entrez gene ID}
#'   \item{gene_name}{Gene symbol}
#'   \item{dataset}{Screen type, either CRISPR or RNAi}
#'   \item{cell_line_name}{Cell line name}
#'   \item{primary_disease}{Cell line disease origin}
#'   \item{lineage}{Cell line lineage}
#'   \item{lineage_subtype}{Cell line lineage subtype, usually a malignancy}
#'   \item{hover.string}{String for hover text when plotting}
#' }
#' @source <https://depmap.org/portal/download/all/>
#' @source <https://bioconductor.org/packages/release/data/experiment/html/depmap.html>
#' @usage
#' data(depmap_22q1_TPM)
"depmap_22q1_TPM"


#' DepMap copy number data
#'
#' The DepMap copy number data for CDK2.
#'
#' @format ## `depmap_22q1_TPM`
#' A data frame with 1754 rows and 9 columns:
#' \describe{
#'   \item{depmap_id}{Depmap cell line ID}
#'   \item{log_copy_number}{log2 copy number}
#'   \item{entrez_id}{Numeric entrez gene ID}
#'   \item{gene_name}{Gene symbol}
#'   \item{dataset}{Screen type, either CRISPR or RNAi}
#'   \item{cell_line_name}{Cell line name}
#'   \item{primary_disease}{Cell line disease origin}
#'   \item{lineage}{Cell line lineage}
#'   \item{lineage_subtype}{Cell line lineage subtype, usually a malignancy}
#'   \item{hover.string}{String for hover text when plotting}
#' }
#' @source <https://depmap.org/portal/download/all/>
#' @source <https://bioconductor.org/packages/release/data/experiment/html/depmap.html>
#' @usage
#' data(depmap_22q1_cn)
"depmap_22q1_cn"
