context("list_datasets")
test_that("function list_dataset works as expected", {
  download_database( tempdir() )
  fileToDb <- file.path(tempdir(), "regulondb_sqlite3.db")
  regdb <- regulondb(fileToDb, organism="prueba", genome_version="prueba", database_version="prueba")

  # Returns a character vector
  expect_type(list_datasets(regdb), "character")

  # Test that function fails if not database is provided
  expect_error(list_datasets(), "")
})
