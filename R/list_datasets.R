#' @title List available datasets in RegulonDB database
#' @description  This function returns a vector of all available tables from a
#' regulondb class.
#' @param regulondb A regulondb class.
#' @keywords data retrieval, datasets, database,
#' @examples
#'
#' ## Connect to the RegulonDB database if necessary
#' if(!exists('regulondb_conn')) regulondb_conn <- connect_database()
#'
#' ## Build the regulon db object
#' e_coli_regulondb <-
#'     regulondb(
#'         database_conn = regulondb_conn,
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
