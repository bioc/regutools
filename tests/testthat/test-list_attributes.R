context("list_attributes")
test_that("list_attributes works as expected ", {
  download_database( tempdir() )
  fileToDb <- file.path(tempdir(), "regulondb_sqlite3.db")
  regdb <- regulondb(fileToDb, organism="prueba", genome_version="prueba", database_version="prueba")

  # Returns a character vector
  expect_type(list_attributes(regdb,dataset = "DNA_OBJECTS"), "character")

  # Test that function fails if not enough arguments are given
  expect_error(list_attributes(),"")
  expect_error(list_attributes(regdb),"")

})