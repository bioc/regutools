test_that("get_gene_synonyms works ", {
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

    expect_s4_class(get_gene_synonyms(regdb, "araC", from = "name"), "regulondb_result")
    expect_error(get_gene_synonyms(regdb, 2, from = "name"))
    expect_error(get_gene_synonyms(regdb, "araC", from = "namez"))
    expect_error(get_gene_synonyms(regdb, "araC", from = "name", to = "namez"))


})
