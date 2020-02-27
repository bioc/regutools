context("plot_dna_objects")
test_that("plot_dna_objects output is a GRanges object", {
    ## Connect to the RegulonDB database if necessary
    if (!exists('regulondb_conn'))
        regulondb_conn <- connect_database()

    ## Build a regulondb object
    regdb <-
        regulondb(
            database_conn = regulondb_conn,
            organism = "chr",
            genome_version = "prueba",
            database_version = "prueba"
        )
    reg_result <-
        plot_dna_objects(regdb,
            from = 5000,
            to = 10000,)
    expect_equivalent(class(reg_result), "GRanges")
})


test_that("non-valid genomic elements causes error ", {
    expect_error(plot_dna_objects(
        regdb,
        from = 5000,
        to - 10000,
        elements = c("gene", "promotr")
    ))
})
