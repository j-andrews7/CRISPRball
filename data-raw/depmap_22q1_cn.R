library(depmap)

gene <- "CDK2"
df <- depmap::depmap_copyNumber()
df$gene <- NULL
df$cell_line <- NULL

depmap.meta <- depmap::depmap_metadata()

df <- df[df$gene_name == gene, ]

df$cell_line_name <- depmap.meta$cell_line_name[match(df$depmap_id, depmap.meta$depmap_id)]
df$primary_disease <- depmap.meta$primary_disease[match(df$depmap_id, depmap.meta$depmap_id)]
df$lineage <- depmap.meta$lineage[match(df$depmap_id, depmap.meta$depmap_id)]
df$lineage_subtype <- depmap.meta$lineage_subtype[match(df$depmap_id, depmap.meta$depmap_id)]

df$hover.string <- paste0(
    "</br><b>Cell Line:</b> ", df$cell_line_name,
    "</br><b>log2(Copy Number):</b> ", format(round(df[["log_copy_number"]], 3), nsmall = 3),
    "</br><b>Lineage:</b> ", df$lineage,
    "</br><b>Disease:</b> ", df$primary_disease
)

depmap_22q1_cn <- df

write_csv(depmap_22q1_cn, "data-raw/depmap_22q1_cn.csv")

usethis::use_data(depmap_22q1_cn, overwrite = TRUE)
