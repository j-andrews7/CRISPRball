# Example data for tests
stopifnot(
    require(pool),
    require(RSQLite)
)

# plot_bar
count.summary <- data.frame(
    Label = c("Sample1", "Sample2", "Sample3"),
    GiniIndex = c(0.2, 0.5, 0.3),
    Zerocounts = c(10, 15, 12)
)

count.summary.empty <- data.frame()

# plot_hist
mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), nrow = 5, ncol = 2)
colnames(mat) <- c("col1", "col2")

mat.one.col <- matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1)
colnames(mat.one.col) <- c("col1")

mat.empty <- matrix(nrow = 0, ncol = 0)

# plot_correlation_heatmap
mat.heatmap <- matrix(c(1, 0.8, 0.6, 0.8, 1, 0.7, 0.6, 0.7, 1), nrow = 3, ncol = 3)
mat.one.col.heatmap <- matrix(c(1), nrow = 1, ncol = 1)

# plot_pca_biplot
# Create a sample matrix and pca object with metadata
mat.pca <- matrix(
    rexp(10 * 20, rate = 0.1),
    ncol = 10
)
rownames(mat.pca) <- paste0("gene", seq_len(nrow(mat.pca)))
colnames(mat.pca) <- paste0("sample", seq_len(ncol(mat.pca)))

metadata <- data.frame(row.names = colnames(mat.pca))
metadata$Group <- rep(NA, ncol(mat.pca))
metadata$Group[seq(1, 10, 2)] <- "A"
metadata$Group[seq(2, 10, 2)] <- "B"

pca.res <- pca(mat.pca, metadata = metadata, removeVar = 0.1)

# get_depmap_plot_data
# Build a SQLite database for testing
con <- dbConnect(RSQLite::SQLite(), ":memory:")
dbWriteTable(con, "meta", data.frame(
    depmap_id = 1:2,
    cell_line_name = c("Cell1", "Cell2"),
    primary_disease = c("Disease1", "Disease2"),
    lineage = c("Lineage1", "Lineage2"),
    lineage_subtype = c("Subtype1", "Subtype2")
))

dbWriteTable(con, "crispr", data.frame(
    depmap_id = 1:2,
    gene_name = c("Gene1", "Gene2"),
    dependency = c(0.1, 0.2)
))

dbWriteTable(con, "rnai", data.frame(
    depmap_id = 1:2,
    gene_name = c("Gene1", "Gene2"),
    dependency = c(0.1, 0.2)
))

dbWriteTable(con, "cn", data.frame(
    depmap_id = 1:2,
    gene_name = c("Gene1", "Gene2"),
    log_copy_number = c(1, 0)
))

dbWriteTable(con, "ccle_tpm", data.frame(
    depmap_id = 1:2,
    gene_name = c("Gene1", "Gene2"),
    rna_expression = c(10, 20)
))

depmap.meta <- dbGetQuery(con, "SELECT * FROM 'meta'")

# plot_depmap_dependency
df <- data.frame(
    dependency = c(0.5, 0.2, 0.3, 0.4, 0.6),
    dataset = c("CRISPR", "RNAi", "CRISPR", "RNAi", "CRISPR"),
    hover.string = c("string1", "string2", "string3", "string4", "string5")
)

df.empty <- data.frame()

# plot_depmap_expression
df.expression <- data.frame(
    rna_expression = c(0.5, 0.2, 0.3, 0.4, 0.6),
    hover.string = c("string1", "string2", "string3", "string4", "string5")
)

# plot_depmap_cn
df.cn <- data.frame(
    log_copy_number = c(1, 1.2, 1.3, 1.4, 0.6),
    hover.string = c("string1", "string2", "string3", "string4", "string5")
)

# plot_depmap_lineages
df.lineages <- data.frame(
    dependency = c(0.5, 0.2, 0.3, 0.4, 0.6),
    lineage = c("A", "B", "A", "B", "A"),
    hover.string = c("string1", "string2", "string3", "string4", "string5")
)

# plot_lawn
df.lawn <- data.frame(
    RandomIndex = runif(50),
    LFC = rnorm(50),
    FDR = runif(50),
    id = paste0("Gene", 1:50)
)
