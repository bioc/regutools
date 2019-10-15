#' @title List attributes from a dataset
#' @description List all attributes and their description of a dataset from RegulonDB. The result of this function may
#' be used as parameter 'values' in \code{\link{GetAttr}} function.
#' @author Carmina Barberena Jonás, Jesús Emiliano Sotelo Fonseca, José Alquicira Hernández, Joselyn Chavez
#' @keywords data retrieval, attributes,
#' @param dataset Dataset of interest
#' @param coments Default FALSE. If TRUE, displays the description of each attribute
#' @return A data frame with two columns:
#' \itemize{
#' \item \code{attribute}. Name of the attribute
#' \item \code{description}. Description of attribute
#' }
#' @examples
#' list_attributes("TF")
#' list_attributes("OPERON", comments = TRUE)
#' @export

list_attributes <- function(dataset, comments=FALSE){

  # Validate mart
  if(!all(dataset %in% list_datasets())){
    cat("Invalid Dataset. Available datasets:\n")
    cat(paste0(list_datasets(), "\n"))
    stop("Please check list_datasets() function.")
  }

  # Query REGULONDB_OBJECTS table
  if (comments){
  query <- paste0("SELECT attribute, description FROM REGULONDB_OBJECTS WHERE table_name = '", dataset, "'")
  }else{
    query <- paste0("SELECT attribute FROM REGULONDB_OBJECTS WHERE table_name = '", dataset, "'")
  }
  # Connect to database
  regulon <- dbConnect(SQLite(),
                       system.file("extdata", "regulondb_sqlite3.db", package = "regutools"))

  # Retrieve data
  result <- dbGetQuery(regulon, query)
  dbDisconnect(regulon)

  return(result)
}
