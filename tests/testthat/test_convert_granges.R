context("GRanges_convert")
test_that( "Results from regulondb queries can be converted to GRanges", {
  download_database( tempdir() )
  fileToDb <- file.path(tempdir(), "regulondb_sqlite3.db")
  regdb <- regulondb(fileToDb, organism="prueba", genome_version="prueba", database_version="prueba")
  query <- get_dataset( regdb,
               dataset="GENE",
               attributes=c("posleft", "posright","name", "strand"),
               filters=list(name=c("araC","crp","lacI") ) )
  expect_s4_class( convert_to_granges( query ), "GRanges" )
  expect_error( convert_to_granges( query[,c("posleft", "name")] ) )
  expect_s4_class( convert_to_granges( query[,c("posleft", "posright")] ), "GRanges" )
  query <- get_dataset( regulondb=regdb,
               dataset="DNA_OBJECTS" )
  expect_s4_class( convert_to_granges( query ), "GRanges" )
  expect_error( convert_to_granges( as.data.frame( query ) ) )
  query <- get_dataset( regulondb=regdb, dataset="OPERON" )
  expect_s4_class( convert_to_granges( query ), "GRanges" )
  query <- get_dataset( regulondb=regdb, dataset="RI" )
  expect_s4_class( convert_to_granges( query ), "GRanges" )
  query <- get_dataset( regulondb=regdb, dataset="TF",
                        attributes="total_regulated_tf",
                        filters=list(name=c("Fis", "Rob", "TyrR") ) )
  expect_error( convert_to_granges( query ) )
} )
