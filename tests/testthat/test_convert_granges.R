context("GRanges_convert")
test_that( "Results from regulondb queries can be converted to GRanges", {
  download_database( tempdir() )
  fileToDb <- file.path(tempdir(), "regulondb_sqlite3.db")
  regdb <- regulondb(fileToDb, organism="prueba", genome_version="prueba", database_version="prueba")
  query <- get_dataset( regdb,
               dataset="GENE",
               attributes=c("posleft", "posright","name", "strand"),
               filters=list(name=c("araC","crp","lacI") ) )
  query

  expect_s4_class( convert_to_granges( query ), "GRanges")
  expect_error( convert_to_granges( query[,c("posleft", "name")] ) )

  regulondb_result <- query
  keep <- !(is.na(regulondb_result$posleft) | is.na(regulondb_result$posright))
  grdata <- with( regulondb_result[which(keep),],
                  GRanges( attributes( regulondb_result )$organism,
                           IRanges(start=posleft, end=posright)) )

} )
