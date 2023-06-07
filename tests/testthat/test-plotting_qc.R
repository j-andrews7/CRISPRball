## plot_bar
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


## plot_hist
# Test that function returns a plotly object
test_that("plot_hist returns a plotly object", {
    plot <- plot_hist(mat)
    expect_s3_class(plot, "plotly")
})

# Test that function can handle matrix with one column
test_that("plot_hist handles matrix with one column", {
    plot <- plot_hist(mat.one.col)
    expect_s3_class(plot, "plotly")
})

# Test that function handles empty matrix
test_that("plot_hist handles empty matrix", {
    plot <- plot_hist(mat.empty)
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "mat is empty.")
})

# Test that the function can handle NULL input
test_that("plot_hist handles NULL input", {
    plot <- plot_hist(NULL)
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "mat is NULL.")
})


# plot_correlation_heatmap
# Test that function returns a Heatmap object
test_that("plot_correlation_heatmap returns a Heatmap object", {
    heatmap <- plot_correlation_heatmap(mat.heatmap)
    expect_s4_class(heatmap, "Heatmap")
})

# Test that function can handle matrix with one column
test_that("plot_correlation_heatmap handles matrix with one column", {
    heatmap <- plot_correlation_heatmap(mat.one.col.heatmap)
    expect_s4_class(heatmap, "Heatmap")
})

# Test that function handles empty matrix
test_that("plot_correlation_heatmap handles empty matrix", {
    heatmap <- plot_correlation_heatmap(mat.empty)
    expect_s4_class(heatmap, "Heatmap")
})

# Test that the function can handle NULL input
test_that("plot_correlation_heatmap handles NULL input", {
    # Expect error as the function cannot handle NULL input
    expect_error(plot_correlation_heatmap(NULL), "mat is NULL.")
})


## plot_pca_biplot
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
