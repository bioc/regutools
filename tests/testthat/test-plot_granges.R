context("plot_granges")
test_that("output with plot_map = FALSE is a GRanges", {
    fileToDb <- file.path(tempdir(), "regulondb_sqlite3.db")
    if (!file.exists(fileToDb))
        download_database(tempdir())
    regdb <-
        regulondb(
            fileToDb,
            organism = "chr",
            genome_version = "prueba",
            database_version = "prueba"
        )
    reg_result <-  plot_GRanges(regdb, from = 5000, to = 10000, plot_map = FALSE)
    expect_equivalent(class(reg_result), "GRanges")
})


test_that("function with plot_map = NULL causes error", {
    expect_error(plot_GRanges(e_coli_regulondb, from = 5000, to = 10000, plot_map = NULL))
})

test_that("non-valid genomic elements causes error ", {
    expect_error(plot_GRanges(regdb, from = 5000, to - 10000, plot_map = FALSE, elements = c("gene", "promotr")) )
})
