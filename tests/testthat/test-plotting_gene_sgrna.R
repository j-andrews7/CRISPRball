## plot_volcano
test_that("plot_volcano returns expected output with default parameters", {
    # generate a random dataframe for testing
    res <- data.frame(
        gene = sample(LETTERS, 100, replace = TRUE),
        LFC = rnorm(100, mean = 0, sd = 1),
        FDR = runif(100, min = 0, max = 1)
    )

    plot <- plot_volcano(res)

    # check that the output is a plotly object
    expect_s3_class(plot, "plotly")
})

test_that("plot_volcano handles different column names correctly", {
    # generate a random dataframe for testing with different column names
    res <- data.frame(
        gene = sample(LETTERS, 100, replace = TRUE),
        logFC = rnorm(100, mean = 0, sd = 1),
        p_value = runif(100, min = 0, max = 1)
    )

    plot <- plot_volcano(res, lfc.term = "logFC", sig.term = "p_value", feat.term = "gene")

    # check that the output is a plotly object
    expect_s3_class(plot, "plotly")
})

test_that("plot_volcano throws an error with missing columns", {
    # generate a random dataframe for testing with missing necessary columns
    res <- data.frame(
        gene = sample(LETTERS, 100, replace = TRUE),
        logFC = rnorm(100, mean = 0, sd = 1),
        p_value = runif(100, min = 0, max = 1)
    )

    plot <- plot_volcano(res, sig.term = "p_value")

    # check that the output is a plotly object
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Column 'LFC' not found in dataframe.")

    plot <- plot_volcano(res, sig.term = "p_value", feat.term = "genes")

    # check that the output is a plotly object
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Column 'genes' not found in dataframe.")
})


## plot_rank
test_that("plot_rank runs without error with minimum parameters", {
    test_df <- data.frame(
        Rank = c(1, 2, 3, 4, 5),
        LFC = c(-0.5, 0.5, -0.6, 0.6, 0),
        FDR = c(0.04, 0.06, 0.03, 0.07, 0.01)
    )
    expect_s3_class(plot_rank(test_df), "plotly")
})

test_that("plot_rank returns error when necessary columns are not present", {
    incorrect_df <- data.frame(
        test1 = c(1, 2, 3, 4, 5),
        test2 = c(-0.5, 0.5, -0.6, 0.6, 0),
        test3 = c(0.04, 0.06, 0.03, 0.07, 0.01)
    )

    plot <- plot_rank(incorrect_df)
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Column 'Rank' not found in dataframe.")
})

test_that("plot_rank handles NULL significance term gracefully", {
    test_df <- data.frame(
        Rank = c(1, 2, 3, 4, 5),
        LFC = c(-0.5, 0.5, -0.6, 0.6, 0),
        FDR = c(NA, NA, 0.03, 0.07, NA)
    )
    expect_s3_class(plot_rank(test_df), "plotly")
})


## plot_lawn
test_that("plot_lawn runs without errors with required inputs", {
    plot <- plot_lawn(
        res = dummy.data,
        x.term = "RandomIndex",
        lfc.term = "LFC",
        sig.term = "FDR",
        feat.term = "id"
    )
    expect_s3_class(plot, "plotly")
})

test_that("plot_lawn fails gracefully with non-existent x.term", {
    plot <- plot_lawn(
        res = dummy.data,
        x.term = "NonExistentColumn",
        lfc.term = "LFC",
        sig.term = "FDR",
        feat.term = "id"
    )
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Column 'NonExistentColumn' not found in dataframe.")
})

test_that("plot_lawn fails gracefully with non-existent lfc.term", {
    plot <- plot_lawn(
        res = dummy.data,
        x.term = "RandomIndex",
        lfc.term = "NonExistentColumn",
        sig.term = "FDR",
        feat.term = "id"
    )
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Column 'NonExistentColumn' not found in dataframe.")
})

test_that("plot_lawn fails gracefully with non-existent sig.term", {
    plot <- plot_lawn(
        res = dummy.data,
        x.term = "RandomIndex",
        lfc.term = "LFC",
        sig.term = "NonExistentColumn",
        feat.term = "id"
    )
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Column 'NonExistentColumn' not found in dataframe.")
})

test_that("plot_lawn fails gracefully with non-existent feat.term", {
    plot <- plot_lawn(
        res = dummy.data,
        x.term = "RandomIndex",
        lfc.term = "LFC",
        sig.term = "FDR",
        feat.term = "NonExistentColumn"
    )
    expect_s3_class(plot, "plotly")
    expect_equal(plot$x$layoutAttrs[[2]]$title$text, "Column 'NonExistentColumn' not found in dataframe.")
})
