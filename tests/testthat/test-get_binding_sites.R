context("Binding sites")
test_that("Function to get transcription factor binding sites works as expected", {
    fileToDb <- file.path(tempdir(), "regulondb_sqlite3.db")
    if(!file.exists(fileToDb))
        download_database(tempdir())
    regdb <-
        regulondb(
            fileToDb,
            organism = "prueba",
            genome_version = "prueba",
            database_version = "prueba" )
    rs <- get_binding_sites( regdb, transcription_factor = "AraC" )
    expect_s4_class( rs, "GRanges" )
    rs <- get_binding_sites( regdb, transcription_factor = "AraC",
                             output_format="Biostrings" )
    expect_s4_class( rs, "DNAStringSet" )
    expect_error( get_binding_sites( regdb, transcription_factor = "AraC",
                             output_format="other" ) )
    expect_null( suppressWarnings(
        get_binding_sites( regdb, transcription_factor = "thisisnotagene" )) )
} )
