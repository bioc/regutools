context("build_condition")

test_that("existing_partial_match returns an expected value", {
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

    # Type character
    expect_type(
        existing_partial_match(
            filters = list(name = c("rec")),
            partialmatch = "name",
            operator = NULL
        ),
        "character"
    )
    # Length 1
    expect_length(existing_partial_match(
        filters = list(name = c("ara")),
        partialmatch = "name",
        operator = NULL
    ),
        1)


})
