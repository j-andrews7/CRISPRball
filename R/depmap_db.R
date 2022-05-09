.error_if_no_depmap <- function() {
  if (!requireNamespace("depmap", quietly = TRUE)) {
    stop("'depmap' installation required for using depmap data.")
  }
}

.error_if_no_rsqlite <- function() {
  if (!requireNamespace("RSQLite", quietly = TRUE)) {
    stop("'RSQLite' installation required for building and using depmap database.")
  }
}

.error_if_no_pool <- function() {
  if (!requireNamespace("pool", quietly = TRUE)) {
    stop("'pool' installation required for building and using depmap database.")
  }
}


#' Build SQLite database of DepMap data
#' 
#' 
#' @return Name of SQLite database containing DepMap data.
#' 
#' @export
build_depmap_db <- function() {
  .error_if_no_depmap()
  .error_if_no_pool()
  .error_if_no_rsqlite()
  
  pool <- pool::dbPool(RSQLite::SQLite(), dbname = "depmap_db.sqlite")
  
  # Get depmap data and make table in database.
  rnai <- depmap::depmap_rnai()
  rnai$gene <- NULL
  rnai$cell_line <- NULL
  pool::dbWriteTable(pool, "rnai", rnai, overwrite = TRUE, append = FALSE)
  rm(rnai)
  
  crispr <- depmap::depmap_crispr()
  crispr$gene <- NULL
  crispr$cell_line <- NULL
  pool::dbWriteTable(pool, "crispr", crispr, overwrite = TRUE, append = FALSE)
  rm(crispr)
  
  cn <- depmap::depmap_copyNumber()
  cn$gene <- NULL
  cn$cell_line <- NULL
  pool::dbWriteTable(pool, "cn", cn, overwrite = TRUE, append = FALSE)
  rm(cn)
  
  ccle_tpm <- depmap::depmap_TPM()
  ccle_tpm$gene <- NULL
  ccle_tpm$cell_line <- NULL
  pool::dbWriteTable(pool, "ccle_tpm", ccle_tpm, overwrite = TRUE, append = FALSE)
  rm(ccle_tpm)
  
  meta <- depmap::depmap_metadata()
  pool::dbWriteTable(pool, "meta", as.data.frame(meta), overwrite = TRUE, append = FALSE)
  
  drug <- depmap::depmap_drug_sensitivity()
  drug$gene <- NULL
  drug$cell_line <- NULL
  drug$smiles <- NULL
  pool::dbWriteTable(pool, "drug", drug, overwrite = TRUE, append = FALSE)
  
  pool::poolClose(pool)
  
  return("depmap_db.sqlite")
}
