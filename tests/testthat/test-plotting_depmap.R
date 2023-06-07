library(pool)
library(RSQLite)
library(plotly)

## get_depmap_plot_data
# Test basic functionality
test_that("get_depmap_plot_data returns a dataframe", {
    df <- get_depmap_plot_data("Gene1", "dependency", depmap.meta, con)
    expect_s3_class(df, "data.frame")

    df <- get_depmap_plot_data("Gene1", "crispr", depmap.meta, con)
    expect_s3_class(df, "data.frame")

    df <- get_depmap_plot_data("Gene1", "rnai", depmap.meta, con)
    expect_s3_class(df, "data.frame")

    df <- get_depmap_plot_data("Gene1", "cn", depmap.meta, con)
    expect_s3_class(df, "data.frame")

    df <- get_depmap_plot_data("Gene1", "ccle_tpm", depmap.meta, con)
    expect_s3_class(df, "data.frame")
})

# Test for correct input gene_name
test_that("get_depmap_plot_data returns correct gene data", {
    df <- get_depmap_plot_data("Gene1", "dependency", depmap.meta, con)
    expect_equal(df$gene_name[1], "Gene1")

    df <- get_depmap_plot_data("Gene1", "crispr", depmap.meta, con)
    expect_equal(df$gene_name[1], "Gene1")

    df <- get_depmap_plot_data("Gene1", "rnai", depmap.meta, con)
    expect_equal(df$gene_name[1], "Gene1")

    df <- get_depmap_plot_data("Gene1", "cn", depmap.meta, con)
    expect_equal(df$gene_name[1], "Gene1")

    df <- get_depmap_plot_data("Gene1", "ccle_tpm", depmap.meta, con)
    expect_equal(df$gene_name[1], "Gene1")
})

# Test for correct input data.type
test_that("get_depmap_plot_data returns expected columns", {
    df <- get_depmap_plot_data("Gene1", "dependency", depmap.meta, con)
    expect_named(df, c(
        "gene_name", "depmap_id", "dependency", "dataset", "hover.string",
        "cell_line_name", "primary_disease", "lineage", "lineage_subtype"
    ),
    ignore.order = TRUE
    )

    df <- get_depmap_plot_data("Gene1", "crispr", depmap.meta, con)
    expect_named(df, c(
        "gene_name", "depmap_id", "dependency", "hover.string",
        "cell_line_name", "primary_disease", "lineage", "lineage_subtype"
    ),
    ignore.order = TRUE
    )

    df <- get_depmap_plot_data("Gene1", "rnai", depmap.meta, con)
    expect_named(df, c(
        "gene_name", "depmap_id", "dependency", "hover.string",
        "cell_line_name", "primary_disease", "lineage", "lineage_subtype"
    ),
    ignore.order = TRUE
    )

    df <- get_depmap_plot_data("Gene1", "cn", depmap.meta, con)
    expect_named(df, c(
        "gene_name", "depmap_id", "log_copy_number", "hover.string",
        "cell_line_name", "primary_disease", "lineage", "lineage_subtype"
    ),
    ignore.order = TRUE
    )

    df <- get_depmap_plot_data("Gene1", "ccle_tpm", depmap.meta, con)
    expect_named(df, c(
        "gene_name", "depmap_id", "rna_expression", "hover.string",
        "cell_line_name", "primary_disease", "lineage", "lineage_subtype"
    ),
    ignore.order = TRUE
    )
})

# Test for correct data values
test_that("get_depmap_plot_data returns expected value", {
    df <- get_depmap_plot_data("Gene1", "dependency", depmap.meta, con)
    expect_equal(df$dependency[1], 0.1)

    df <- get_depmap_plot_data("Gene1", "crispr", depmap.meta, con)
    expect_equal(df$dependency[1], 0.1)

    df <- get_depmap_plot_data("Gene1", "rnai", depmap.meta, con)
    expect_equal(df$dependency[1], 0.1)

    df <- get_depmap_plot_data("Gene1", "cn", depmap.meta, con)
    expect_equal(df$log_copy_number[1], 1)

    df <- get_depmap_plot_data("Gene1", "ccle_tpm", depmap.meta, con)
    expect_equal(df$rna_expression[1], 10)
})

# Test for correct cell_line_name
test_that("get_depmap_plot_data returns correct cell line name", {
    df <- get_depmap_plot_data("Gene1", "dependency", depmap.meta, con)
    expect_equal(df$cell_line_name[1], "Cell1")

    df <- get_depmap_plot_data("Gene1", "crispr", depmap.meta, con)
    expect_equal(df$cell_line_name[1], "Cell1")

    df <- get_depmap_plot_data("Gene1", "rnai", depmap.meta, con)
    expect_equal(df$cell_line_name[1], "Cell1")

    df <- get_depmap_plot_data("Gene1", "cn", depmap.meta, con)
    expect_equal(df$cell_line_name[1], "Cell1")

    df <- get_depmap_plot_data("Gene1", "ccle_tpm", depmap.meta, con)
    expect_equal(df$cell_line_name[1], "Cell1")
})

# Test for NULL output if data is not available
test_that("get_depmap_plot_data returns NULL if no data is available", {
    df <- get_depmap_plot_data("Gene3", "dependency", depmap.meta, con)
    expect_null(df)

    df <- get_depmap_plot_data("Gene3", "crispr", depmap.meta, con)
    expect_null(df)

    df <- get_depmap_plot_data("Gene3", "rnai", depmap.meta, con)
    expect_null(df)

    df <- get_depmap_plot_data("Gene3", "cn", depmap.meta, con)
    expect_null(df)

    df <- get_depmap_plot_data("Gene3", "ccle_tpm", depmap.meta, con)
    expect_null(df)
})


## plot_depmap_dependency
# Test that function returns a plotly object
test_that("plot_depmap_dependency returns a plotly object", {
    suppressWarnings(plot <- plot_depmap_dependency(df))
    expect_s3_class(plot, "plotly")
})

# Test that function can handle empty dataframes
test_that("plot_depmap_dependency handles empty dataframes", {
    plot <- plot_depmap_dependency(df.empty)
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Gene not found in DepMap.")
})

# Test that plot contains the correct traces (CRISPR and RNAi)
test_that("plot_depmap_dependency contains correct traces", {
    suppressWarnings(plot <- plot_depmap_dependency(df))
    traces <- plot$x$data
    expect_true(any(sapply(traces, function(x) x$name == "CRISPR")))
    expect_true(any(sapply(traces, function(x) x$name == "RNAi")))
})

# Test the correct color assignment
test_that("plot_depmap_dependency assigns correct colors", {
    suppressWarnings(plot <- plot_depmap_dependency(df, crispr.color = "#FF0000", rnai.color = "#2f00ff"))
    traces <- plot$x$data
    crispr_trace <- traces[[which(sapply(traces, function(x) x$name == "CRISPR"))]]
    rnai_trace <- traces[[which(sapply(traces, function(x) x$name == "RNAi"))]]
    expect_equal(crispr_trace$line$color, "rgba(255,0,0,0.6)")
    expect_equal(rnai_trace$line$color, "rgba(47,0,255,0.6)")
})

# Test the dependency threshold line presence
test_that("plot_depmap_dependency shows dependency threshold line when depline = TRUE", {
    suppressWarnings(plot <- plot_depmap_dependency(df, depline = TRUE))
    shapes <- plot$x$layout$shapes
    expect_true(any(sapply(shapes, function(x) x$xref == "paper" && x$line$color == "rgba(51,51,51,1)")))
})

# Test the dependency threshold line absence
test_that("plot_depmap_dependency does not show dependency threshold line when depline = FALSE", {
    suppressWarnings(plot <- plot_depmap_dependency(df, depline = FALSE))
    shapes <- plot$x$layout$shapes
    expect_false(any(sapply(shapes, function(x) x$xref == "x" && x$line$color == "rgba(51,51,51,1)")))
})

# Test that the function can handle NULL input
test_that("plot_depmap_dependency handles NULL input", {
    suppressWarnings(plot <- plot_depmap_dependency(NULL))
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Gene not found in DepMap.")
})


## plot_depmap_expression
# Test that function returns a plotly object
test_that("plot_depmap_expression returns a plotly object", {
    suppressWarnings(plot <- plot_depmap_expression(df.expression))
    expect_s3_class(plot, "plotly")
})

# Test that function can handle empty dataframes
test_that("plot_depmap_expression handles empty dataframes", {
    suppressWarnings(plot <- plot_depmap_expression(df.empty))
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Gene not found in DepMap.")
})

# Test the correct color assignment
test_that("plot_depmap_expression assigns correct colors", {
    suppressWarnings(plot <- plot_depmap_expression(df.expression, color = "#FF0000"))
    traces <- plot$x$data
    expect_equal(traces[[1]]$line$color, "rgba(255,0,0,1)")
})

# Test that the function can handle NULL input
test_that("plot_depmap_expression handles NULL input", {
    suppressWarnings(plot <- plot_depmap_expression(NULL))
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Gene not found in DepMap.")
})


## plot_depmap_cn
# Test that function returns a plotly object
test_that("plot_depmap_cn returns a plotly object", {
    suppressWarnings(plot <- plot_depmap_cn(df.cn))
    expect_s3_class(plot, "plotly")
})

# Test that function can handle empty dataframes
test_that("plot_depmap_cn handles empty dataframes", {
    df.empty <- data.frame()
    suppressWarnings(plot <- plot_depmap_cn(df.empty))
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Gene not found in DepMap.")
})

# Test the correct color assignment
test_that("plot_depmap_cn assigns correct colors", {
    suppressWarnings(plot <- plot_depmap_cn(df.cn, color = "#FF0000"))
    traces <- plot$x$data
    expect_equal(traces[[1]]$line$color, "rgba(255,0,0,1)")
})

# Test that the function can handle NULL input
test_that("plot_depmap_cn handles NULL input", {
    suppressWarnings(plot <- plot_depmap_cn(NULL))
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Gene not found in DepMap.")
})


## plot_depmap_lineages
# Test that function returns a plotly object
test_that("plot_depmap_lineages returns a plotly object", {
    suppressWarnings(plot <- plot_depmap_lineages(df.lineages, plot.val = "dependency", group.by = "lineage"))
    expect_s3_class(plot, "plotly")
})

# Test that function can handle empty dataframes
test_that("plot_depmap_lineages handles empty dataframes", {
    df.empty <- data.frame()
    suppressWarnings(plot <- plot_depmap_lineages(df.empty, plot.val = "dependency", group.by = "lineage"))
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Gene not found in DepMap.")
})

# Test that function handles NULL lineage
test_that("plot_depmap_lineages handles NULL lineage", {
    suppressWarnings(plot <- plot_depmap_lineages(df.lineages, plot.val = "dependency", group.by = "lineage", lineage = NULL))
    suppressWarnings(pp <- plotly_build(plot))
    traces <- pp$x$data
    expect_equal(length(traces), 2) # Should be box and scatter plots
})

# Test that function filters data by specified lineage
test_that("plot_depmap_lineages filters data by specified lineage", {
    suppressWarnings(plot <- plot_depmap_lineages(df.lineages, plot.val = "dependency", group.by = "lineage", lineage = "A"))
    suppressWarnings(pp <- plotly_build(plot))
    traces <- pp$x$data
    expect_equal(length(traces), 2) # Should be box and scatter plots
    # Verify data for each trace
    expect_equal(length(traces[[1]]$x), 3) # 3 points in lineage "A"
    expect_equal(length(traces[[2]]$x), 3) # 3 points in lineage "A"
})

# Test that the function can handle NULL input
test_that("plot_depmap_lineages handles NULL input", {
    suppressWarnings(plot <- plot_depmap_lineages(NULL, plot.val = "dependency", group.by = "lineage"))
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Gene not found in DepMap.")
})
