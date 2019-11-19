#' @title List attributes/fields from a dataset/table
#' @description List all attributes and their description of a dataset from RegulonDB. The result of this function may
#' be used as parameter 'values' in \code{\link{list_attributes}} function.
#' @author Carmina Barberena Jonás, Jesús Emiliano Sotelo Fonseca, José Alquicira Hernández, Joselyn Chavez
#' @keywords data retrieval, attributes,
#' @param dataset Dataset of interest. The name should correspond to a table of the database.
#' @return A character vector with the field names.
#' @examples
#' download_database(getwd())
#' ecoli_regulondb <-build_regulondb( database_path=file.path(getwd(), "regulondb_sqlite3.db"), organism="E.coli", database_version="1", genome_version="1" )
#' list_attributes(ecoli_regulondb, "TF")
#' list_attributes(ecoli_regulondb, "OPERON")
#' @export

list_attributes <- function( regulondb, dataset#, description=FALSE
                             ){
  if(missing(dataset))
    stop("Parameter 'dataset' is missing, please specify\n")
  stopifnot( validObject( regulondb ) )
#  if( !description ){
    dbListFields( regulondb, dataset )
  # }else{
  #   query <- paste0("SELECT attribute, description FROM REGULONDB_OBJECTS WHERE table_name = '", dataset, "'")
  #   dbGetQuery( regulondb, query )
  # }
}
