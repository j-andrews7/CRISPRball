library(depmap)

gene <- "CDK2"

df.c <- depmap::depmap_crispr()
df.c$gene <- NULL
df.c$cell_line <- NULL
df.c <- df.c[df.c$gene_name == gene, ]

df.r <- depmap::depmap_rnai()
df.r$gene <- NULL
df.r$cell_line <- NULL
df.r <- df.r[df.r$gene_name == gene, ]

df <- data.frame()

if (nrow(df.c) > 0) {
    df.c$dataset <- "CRISPR"
    df <- df.c
}

if (nrow(df.r) > 0) {
    df.r$dataset <- "RNAi"

    if (nrow(df) > 0) {
        df <- rbind(df, df.r)
    } else {
        df <- df.r
    }
}

depmap.meta <- depmap::depmap_metadata()

df$cell_line_name <- depmap.meta$cell_line_name[match(df$depmap_id, depmap.meta$depmap_id)]
df$primary_disease <- depmap.meta$primary_disease[match(df$depmap_id, depmap.meta$depmap_id)]
df$lineage <- depmap.meta$lineage[match(df$depmap_id, depmap.meta$depmap_id)]
df$lineage_subtype <- depmap.meta$lineage_subtype[match(df$depmap_id, depmap.meta$depmap_id)]

df$hover.string <- paste0(
    "</br><b>Cell Line:</b> ", df$cell_line_name,
    "</br><b>Dependency:</b> ", format(round(df[["dependency"]], 3), nsmall = 3),
    "</br><b>Lineage:</b> ", df$lineage,
    "</br><b>Disease:</b> ", df$primary_disease
)

depmap_22q1_crispr_rnai <- df

write_csv(depmap_22q1_crispr_rnai, "data-raw/depmap_22q1_crispr_rnai.csv")

usethis::use_data(depmap_22q1_crispr_rnai, overwrite = TRUE)
