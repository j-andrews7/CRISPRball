test_that("get_depmap_essentiality returns expected output", {
    # Create a mock depmap.summary dataframe
    depmap.summary <- data.frame(
        gene_name = c("A", "B", "C", "A", "B", "C"),
        dataset = c(
            "RNAi_merged", "RNAi_merged", "RNAi_merged",
            "Chronos_Combined", "Chronos_Combined", "Chronos_Combined"
        ),
        dependent_cell_lines = c(10, 20, 30, 40, 50, 60),
        cell_lines_with_data = c(100, 200, 300, 400, 500, 600),
        common_essential = c(1, 0, 1, 0, 1, 0),
        strongly_selective = c(0, 1, 0, 1, 0, 1)
    )

    # Test gene "A"
    res_A <- get_depmap_essentiality("A", depmap.summary)

    expect_equal(res_A$crispr$avail, TRUE)
    expect_equal(res_A$crispr$dataset, "Chronos_Combined")
    expect_equal(res_A$crispr$dep_lines, 40)
    expect_equal(res_A$crispr$total_lines, 400)
    expect_equal(res_A$crispr$label, "STRONGLY SELECTIVE")

    expect_equal(res_A$rnai$avail, TRUE)
    expect_equal(res_A$rnai$dataset, "RNAi_merged")
    expect_equal(res_A$rnai$dep_lines, 10)
    expect_equal(res_A$rnai$total_lines, 100)
    expect_equal(res_A$rnai$label, "COMMON ESSENTIAL")

    # Test gene "B"
    res_B <- get_depmap_essentiality("B", depmap.summary)

    expect_equal(res_B$crispr$avail, TRUE)
    expect_equal(res_B$crispr$dataset, "Chronos_Combined")
    expect_equal(res_B$crispr$dep_lines, 50)
    expect_equal(res_B$crispr$total_lines, 500)
    expect_equal(res_B$crispr$label, "COMMON ESSENTIAL")

    expect_equal(res_B$rnai$avail, TRUE)
    expect_equal(res_B$rnai$dataset, "RNAi_merged")
    expect_equal(res_B$rnai$dep_lines, 20)
    expect_equal(res_B$rnai$total_lines, 200)
    expect_equal(res_B$rnai$label, "STRONGLY SELECTIVE")

    # Test gene "Z" that does not exist in depmap.summary
    res_Z <- get_depmap_essentiality("Z", depmap.summary)

    expect_equal(res_Z$crispr$avail, FALSE)
    expect_equal(res_Z$rnai$avail, FALSE)
})


test_that(".make_dependency_tag generates correct output", {
    dep.info <- list(
        crispr = list(
            avail = TRUE, dep_lines = 10,
            total_lines = 20, dataset = "dataset1",
            label = "COMMON ESSENTIAL"
        ),
        rnai = list(
            avail = TRUE, dep_lines = 30,
            total_lines = 40, dataset = "dataset2",
            label = "STRONGLY SELECTIVE"
        )
    )
    dep.release <- "Release1"
    crispr.color <- "#ff0000"
    rnai.color <- "#00ff00"

    res <- .make_dependency_tag(dep.info, dep.release, crispr.color, rnai.color)

    # Check that the result is a tagList
    expect_s3_class(res, "shiny.tag.list")

    # Check the CRISPR section
    crispr_tag <- res[[1]]
    expect_equal(crispr_tag$name, "div")
    expect_equal(crispr_tag$attribs$style, "margin-bottom: 7px;")

    crispr_span <- crispr_tag$children[[1]]
    expect_equal(crispr_span$name, "span")
    expect_equal(crispr_span$attribs$style, paste0("color: ", crispr.color, ";"))

    crispr_strong <- crispr_span$children[[1]]
    expect_equal(crispr_strong$name, "strong")
    expect_equal(
        crispr_strong$children[[1]],
        paste0(
            "CRISPR (DepMap ", dep.release, ", ",
            dep.info$crispr$dataset, "): ", dep.info$crispr$dep_lines,
            "/", dep.info$crispr$total_lines
        )
    )

    # Check the RNAi section
    rnai_tag <- res[[3]]
    expect_equal(rnai_tag$name, "div")
    expect_equal(rnai_tag$attribs$style, "margin-bottom: 7px; margin-top: 8px")

    rnai_span <- rnai_tag$children[[1]]
    expect_equal(rnai_span$name, "span")
    expect_equal(rnai_span$attribs$style, paste0("color: ", rnai.color, ";"))

    rnai_strong <- rnai_span$children[[1]]
    expect_equal(rnai_strong$name, "strong")
    expect_equal(
        rnai_strong$children[[1]],
        paste0(
            "RNAi (DepMap ", dep.release, ", ",
            dep.info$rnai$dataset, "): ", dep.info$rnai$dep_lines,
            "/", dep.info$rnai$total_lines
        )
    )
})


test_that(".make_gene_tag generates correct output", {
    gene <- "BRCA1"

    res <- .make_gene_tag(gene)

    # Check that the result is a tagList
    expect_s3_class(res, "shiny.tag.list")

    # Check that the result has the expected length
    expect_equal(length(res), 4)

    # Check that the first element contains the gene name
    expect_true(grepl(gene, as.character(res[[1]])))
})


test_that(".make_gene_tag handles non-existent genes", {
    gene <- "nonexistentgene"

    res <- .make_gene_tag(gene)

    # Check that the result is a tagList
    expect_s3_class(res, "shiny.tag.list")

    # Check that the result has the expected length
    expect_equal(length(res), 1)

    # Check that the message is as expected
    expect_true(grepl(
        "Unable to find gene information.",
        as.character(res[[1]])
    ))
})
