test_that("DepMap database is retrieved properly", {
    db <- tempfile()
    build_depmap_db(file = db)

    expect_true(file.exists(db))

    pool <- pool::dbPool(RSQLite::SQLite(), dbname = db)
    depmap.gene <- pool::dbGetQuery(pool, "SELECT * FROM 'gene.summary'")

    essentials <- get_depmap_essentiality(gene = "CDK2", depmap.summary = depmap.gene)
    expect_equal(length(essentials), 2)
    expect_equal(names(essentials), c("crispr", "rnai"))

    depmap.meta <- pool::dbGetQuery(pool, "SELECT * FROM 'meta'")
    expect_true(is(depmap.meta, "data.frame"))

    df <- get_depmap_plot_data(
        gene = "CDK2", data.type = "dependency",
        depmap.meta = depmap.meta, depmap.pool = pool
    )
    expect_equal(ncol(df), 10)
    expect_true(nrow(df) > 0)

    df <- get_depmap_plot_data(
        gene = "CDK2", data.type = "crispr",
        depmap.meta = depmap.meta, depmap.pool = pool
    )
    expect_equal(ncol(df), 9)
    expect_true(nrow(df) > 0)

    df <- get_depmap_plot_data(
        gene = "CDK2", data.type = "rnai",
        depmap.meta = depmap.meta, depmap.pool = pool
    )
    expect_equal(ncol(df), 9)
    expect_true(nrow(df) > 0)

    df <- get_depmap_plot_data(
        gene = "CDK2", data.type = "cn",
        depmap.meta = depmap.meta, depmap.pool = pool
    )
    expect_equal(ncol(df), 9)
    expect_true(nrow(df) > 0)

    df <- get_depmap_plot_data(
        gene = "CDK2", data.type = "ccle_tpm",
        depmap.meta = depmap.meta, depmap.pool = pool
    )
    expect_equal(ncol(df), 9)
    expect_true(nrow(df) > 0)

    pool::poolClose(pool)
})
