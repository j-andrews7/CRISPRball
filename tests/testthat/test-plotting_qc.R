test_that("plot_bar function", {
    # Create sample dataframe
    count.summary <- data.frame(
        Label = c("Sample1", "Sample2", "Sample3"),
        GiniIndex = c(0.2, 0.5, 0.3),
        Zerocounts = c(10, 15, 12)
    )

    # Test that function returns a plotly object
    test_that("plot_bar returns a plotly object", {
        plot <- plot_bar(count.summary, x = "Label", y = "GiniIndex")
        expect_s3_class(plot, "plotly")
    })

    # Test that function can handle different y parameters
    test_that("plot_bar handles different y parameters", {
        plot <- plot_bar(count.summary, x = "Label", y = "Zerocounts")
        expect_s3_class(plot, "plotly")
    })

    # Test that function can handle empty dataframes
    test_that("plot_bar handles empty dataframes", {
        count.summary.empty <- data.frame()
        plot <- plot_bar(count.summary.empty, x = "Label", y = "GiniIndex")
        expect_s3_class(plot, "plotly")
        expect_equal(plot$x$layoutAttrs[[2]]$title$text, "count.summary is empty.")
    })

    # Test that the function can handle NULL input
    test_that("plot_bar handles NULL input", {
        plot <- plot_bar(NULL, x = "Label", y = "GiniIndex")
        expect_s3_class(plot, "plotly")
        expect_equal(plot$x$layoutAttrs[[2]]$title$text, "count.summary is NULL.")
    })
})


test_that("plot_hist function", {
    # Create a sample matrix
    mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), nrow = 5, ncol = 2)
    colnames(mat) <- c("col1", "col2")

    # Test that function returns a plotly object
    test_that("plot_hist returns a plotly object", {
        plot <- plot_hist(mat)
        expect_s3_class(plot, "plotly")
    })

    # Test that function can handle matrix with one column
    test_that("plot_hist handles matrix with one column", {
        mat_one_col <- matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1)
        colnames(mat_one_col) <- c("col1")
        plot <- plot_hist(mat_one_col)
        expect_s3_class(plot, "plotly")
    })

    # Test that function handles empty matrix
    test_that("plot_hist handles empty matrix", {
        mat_empty <- matrix(nrow = 0, ncol = 0)
        plot <- plot_hist(mat_empty)
        expect_s3_class(plot, "plotly")
        expect_equal(plot$x$layoutAttrs[[2]]$title$text, "mat is empty.")
    })

    # Test that the function can handle NULL input
    test_that("plot_hist handles NULL input", {
        plot <- plot_hist(NULL)
        expect_s3_class(plot, "plotly")
        expect_equal(plot$x$layoutAttrs[[2]]$title$text, "mat is NULL.")
    })
})


test_that("plot_correlation_heatmap function", {
    # Create a sample correlation matrix
    mat <- matrix(c(1, 0.8, 0.6, 0.8, 1, 0.7, 0.6, 0.7, 1), nrow = 3, ncol = 3)

    # Test that function returns a Heatmap object
    test_that("plot_correlation_heatmap returns a Heatmap object", {
        heatmap <- plot_correlation_heatmap(mat)
        expect_s4_class(heatmap, "Heatmap")
    })

    # Test that function can handle matrix with one column
    test_that("plot_correlation_heatmap handles matrix with one column", {
        mat_one_col <- matrix(c(1), nrow = 1, ncol = 1)
        heatmap <- plot_correlation_heatmap(mat_one_col)
        expect_s4_class(heatmap, "Heatmap")
    })

    # Test that function handles empty matrix
    test_that("plot_correlation_heatmap handles empty matrix", {
        mat_empty <- matrix(nrow = 0, ncol = 0)
        heatmap <- plot_correlation_heatmap(mat_empty)
        expect_s4_class(heatmap, "Heatmap")
    })

    # Test that the function can handle NULL input
    test_that("plot_correlation_heatmap handles NULL input", {
        # Expect error as the function cannot handle NULL input
        expect_error(plot_correlation_heatmap(NULL), "mat is NULL.")
    })
})


test_that("plot_pca_biplot function", {
    # Create a sample matrix and pca object with metadata
    mat <- matrix(
        rexp(10 * 20, rate = 0.1),
        ncol = 10
    )
    rownames(mat) <- paste0("gene", seq_len(nrow(mat)))
    colnames(mat) <- paste0("sample", seq_len(ncol(mat)))

    metadata <- data.frame(row.names = colnames(mat))
    metadata$Group <- rep(NA, ncol(mat))
    metadata$Group[seq(1, 10, 2)] <- "A"
    metadata$Group[seq(2, 10, 2)] <- "B"

    pca.res <- pca(mat, metadata = metadata, removeVar = 0.1)

    # Test that function returns a plotly object
    test_that("plot_pca_biplot returns a plotly object", {
        plot <- plot_pca_biplot(pca.res)
        expect_s3_class(plot, "plotly")
    })

    # Test that function can handle NULL pca object
    test_that("plot_pca_biplot handles NULL pca object", {
        expect_error(plot_pca_biplot(NULL))
    })

    # Test that function can handle pca object with no metadata
    test_that("plot_pca_biplot handles pca object with no metadata", {
        pca.res$metadata <- NULL
        plot <- suppressWarnings(plot_pca_biplot(pca.res, color.by = "Group"))
        expect_s3_class(plot, "plotly")
        # Expect a warning as the function can handle this case, but should warn the user
        expect_warning(plot_pca_biplot(pca.res, color.by = "Group"))
    })

    # Test that the function correctly handles valid `color.by` parameter
    test_that("plot_pca_biplot correctly handles valid `color.by` parameter", {
        plot <- plot_pca_biplot(pca.res, color.by = "Group")
        expect_s3_class(plot, "plotly")
    })

    # Test that the function correctly handles invalid `color.by` parameter
    test_that("plot_pca_biplot correctly handles invalid `color.by` parameter", {
        # Expect error as the column does not exist in the metadata
        expect_error(plot_pca_biplot(pca.res, color.by = "InvalidColumnName"))
    })

    # Test that the function correctly handles `show.loadings` parameter
    test_that("plot_pca_biplot correctly handles `show.loadings` parameter", {
        plot <- plot_pca_biplot(pca.res, show.loadings = TRUE)
        expect_s3_class(plot, "plotly")
    })

    # Test that the function correctly handles `dim.z` parameter
    test_that("plot_pca_biplot correctly handles `dim.z` parameter", {
        plot <- plot_pca_biplot(pca.res, dim.z = "PC3")
        expect_s3_class(plot, "plotly")
    })
})
