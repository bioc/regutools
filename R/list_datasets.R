#' @title List available datasets in RegulonDB database
#' @description  This function returns a vector of all available tables from a
#' regulondb class.
#' @param regulondb A regulondb class.
#' @keywords data retrieval, datasets, database,
#' @examples
#'
#' ## Download the database if necessary
#' if(!file.exists(file.path(tempdir(), 'regulondb-sqlite3.db'))) {
#'     download_database(tempdir())
#' }
#'
#' ## Build the regulon db object
#' e_coli_regulondb <-
#'     regulondb(
#'         database_path = file.path(tempdir(), "regulondb_sqlite3.db"),
#'         organism = "E.coli",
#'         database_version = "1",
#'         genome_version = "1"
#'     )
#'
#' ## List the available datasets
#' list_datasets(e_coli_regulondb)
#'
#' @export
list_datasets <- function(regulondb) {
    stopifnot(validObject(regulondb))
    dbListTables(regulondb)
}
