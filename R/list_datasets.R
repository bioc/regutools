#' @title List available datasets in RegulonDB database
#' @description  This function returns a vector of all available datasets. No parameters are provided.
#' @keywords data retrieval, datasets, database,
#' @author Carmina Barberena Jonás, Jesús Emiliano Sotelo Fonseca, Josá Alquicira Hernández
#' @examples
#' list_datasets()
#' @export

list_datasets <- function( regulonData ){
  # Connect to database

  # List tables
  result <- dbListTables(regulon)
  dbDisconnect(regulon)
  return(result)
}
