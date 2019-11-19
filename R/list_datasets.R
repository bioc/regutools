#' @title List available datasets in RegulonDB database
#' @description  This function returns a vector of all available tables from a regulondb class.
#' @param regulondb A regulondb class.
#' @keywords data retrieval, datasets, database,
#' @examples
#' download_database(getwd())
#' ecoli_regulondb <-build_regulondb( database_path=file.path(getwd(), "regulondb_sqlite3.db"), organism="E.coli", database_version="1", genome_version="1" )
#' list_datasets( ecoli_regulondb )
#' @export
list_datasets <- function( regulondb ){
  stopifnot(validObject(regulondb))
  result <- dbListTables( regulondb )
  return( result )
}
