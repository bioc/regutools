context("get_attributes")
test_that( "Function get attributes works as expected", {
  download_database( tempdir() )
  fileToDb <- file.path(tempdir(), "regulondb_sqlite3.db")
  regdb <- regulondb(fileToDb, organism="prueba", genome_version="prueba", database_version="prueba")
  query <- get_dataset( regdb,
                        dataset="GENE",
                        attributes=c("posleft", "posright","name", "strand"),
                        filters=list(name=c("araC","crp","lacI") ) )
  expect_s4_class( query, "regulondb_result" )
} )
