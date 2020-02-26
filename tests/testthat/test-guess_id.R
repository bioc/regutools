test_that("guess_id works ", {
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

    expect_error(guess_id(1, regdb))
    expect_error(guess_id(c("araC","araG"), regdb))
    expect_length(guess_id("araC", regdb), 1)
    expect_length(guess_id("b0064", regdb), 1)
    expect_length(guess_id("ECK120000050", regdb), 1)

})
