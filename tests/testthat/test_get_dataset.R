context("get_dataset")
test_that( "Function get dataset works as expected", {
  download_database( tempdir() )
  fileToDb <- file.path(tempdir(), "regulondb_sqlite3.db")
  regdb <- regulondb(fileToDb, organism="prueba", genome_version="prueba", database_version="prueba")
  query <- get_dataset( regdb,
                        dataset="GENE",
                        attributes=c("posleft", "posright","name", "strand"),
                        filters=list(name=c("araC","crp","lacI") ) )
  expect_s4_class( query, "regulondb_result" )
  expect_error(
    get_dataset( regdb, dataset="GENE", attributes=c("posleft", "posright","name", "strand"),
                 filters=list(name=c("araC","crp","lacI") ), output_format="ss" ) )
  expect_s4_class(
    get_dataset( regdb, dataset="GENE", attributes=c("posleft", "posright","name", "strand"),
               filters=list(name=c("araC","crp","lacI") ), output_format="GRanges" ), "GRanges" )
  expect_s4_class(
    get_dataset( regdb, dataset="GENE",
               filters=list(name=c("araC","crp","lacI") ),
               attributes=c("posleft", "posright","name", "strand", "dna_sequence"),
               output_format="DNAStringSet"), "DNAStringSet" )
  expect_s4_class(
    get_dataset( regdb, dataset="GENE",
                 filters=list(name=c("araC","crp","lacI") ),
                 attributes=c("posleft", "posright","name", "strand", "product_sequence"),
                 output_format="BStringSet"), "BStringSet" )
} )

