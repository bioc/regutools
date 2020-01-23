context("get_gene_regulators")
test_that("function get_gene_regulators works as expected", {
    fileToDb <- file.path(tempdir(), "regulondb_sqlite3.db")
    if (!file.exists(fileToDb))
        download_database(tempdir())
    regdb <-
        regulondb(
            fileToDb,
            organism = "prueba",
            genome_version = "prueba",
            database_version = "prueba"
        )

    expect_s4_class(get_gene_regulators(regdb, genes = "araC"),
        "regulondb_result")
    expect_s4_class(
        get_gene_regulators(regdb, genes = "araC", output.type = "GENE"),
        "regulondb_result"
    )
    expect_s4_class(get_gene_regulators(regdb, genes = "araC", format = "table"),
        "regulondb_result")
    expect_s4_class(get_gene_regulators(regdb, genes = "araC", format = "onerow"),
        "regulondb_result")

    expect_error(get_gene_regulators(regdb, genes = "araC", format = ""), "")
})
