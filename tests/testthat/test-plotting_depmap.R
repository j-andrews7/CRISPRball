library(pool)
library(RSQLite)
library(plotly)

test_that("get_depmap_plot_data function", {
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
})


test_that("plot_depmap_dependency function", {

    # Create sample dataframe
    df <- data.frame(
        dependency = c(0.5, 0.2, 0.3, 0.4, 0.6),
        dataset = c("CRISPR", "RNAi", "CRISPR", "RNAi", "CRISPR"),
        hover.string = c("string1", "string2", "string3", "string4", "string5")
    )

    # Test that function returns a plotly object
    test_that("plot_depmap_dependency returns a plotly object", {
        plot <- plot_depmap_dependency(df)
        expect_is(plot, "plotly")
    })

    # Test that function can handle empty dataframes
    test_that("plot_depmap_dependency handles empty dataframes", {
        df.empty <- data.frame()
        plot <- plot_depmap_dependency(df.empty)
        expect_equal(plot$x$data[[1]]$text[1], "Gene not found in DepMap.")
    })

    # Test that plot contains the correct traces (CRISPR and RNAi)
    test_that("plot_depmap_dependency contains correct traces", {
        plot <- plot_depmap_dependency(df)
        traces <- plot$x$data
        expect_true(any(sapply(traces, function(x) x$name == "CRISPR")))
        expect_true(any(sapply(traces, function(x) x$name == "RNAi")))
    })

    # Test the correct color assignment
    test_that("plot_depmap_dependency assigns correct colors", {
        plot <- plot_depmap_dependency(df, crispr.color = "#FF0000", rnai.color = "#2f00ff")
        traces <- plot$x$data
        crispr_trace <- traces[[which(sapply(traces, function(x) x$name == "CRISPR"))]]
        rnai_trace <- traces[[which(sapply(traces, function(x) x$name == "RNAi"))]]
        expect_equal(crispr_trace$line$color, "#FF0000")
        expect_equal(rnai_trace$line$color, "#2f00ff")
    })

    # Test the dependency threshold line presence
    test_that("plot_depmap_dependency shows dependency threshold line when depline = TRUE", {
        plot <- plot_depmap_dependency(df, depline = TRUE)
        shapes <- plot$x$layout$shapes
        expect_true(any(sapply(shapes, function(x) x$xref == "x" && x$line$color == "red")))
    })

    # Test the dependency threshold line absence
    test_that("plot_depmap_dependency does not show dependency threshold line when depline = FALSE", {
        plot <- plot_depmap_dependency(df, depline = FALSE)
        shapes <- plot$x$layout$shapes
        expect_false(any(sapply(shapes, function(x) x$xref == "x" && x$line$color == "red")))
    })

    # Test that the function can handle NULL input
    test_that("plot_depmap_dependency handles NULL input", {
        plot <- plot_depmap_dependency(NULL)
        expect_equal(plot$x$data[[1]]$text[1], "Gene not found in DepMap.")
    })
})
