#' @title Extract data from RegulonDB
#' @description This function retrieves data from RegulonDB. Attributes from datasets can be selected and filtered.
#' @param regulondb A regulondb object.
#' @param dataset Dataset of interest. Use the function list_datasets for an overview of valid datasets.
#' @param attributes Vector of attributes to be retrieved.
#' @param filters List of filters to be used. The names should correspond to the attribute and the values correspond to the condition for selection.
#' @param and Logical argument. If FALSE, filters will be considered under the "OR" operator
#' @param partialmatch name of the condition(s) with a string pattern for full or partial match in the query
#' @param interval the filters whose values will be considered as interval
#' @keywords data retrieval, attributes, filters,
#' @author Carmina Barberena Jonas, Jesús Emiliano Sotelo Fonseca, José Alquicira Hernández, Joselyn Chávez
#' @examples
#' # Obtain all the information from the "GENE" dataset
#' get_dataset(dataset="GENE")
#'
#' # Get the attributes posright and name from the "GENE" dataset
#' get_dataset(dataset="GENE", attributes=c("posright","name") )
#'
#' # From "GENE" dataset, get the gene name, strand, posright, product name and id of all genes regulated
#'   with name like "ara", strand as "forward" with a position rigth between 2000 and 40000
#'
#' get_dataset(dataset = "GENE",
#'        attributes = c("name", "strand", "posright", "product_name", "id"),
#'        filters = list(name=c("ara"),
#'                       strand=c("forward"),
#'                       posright=c("2000","40000")),
#'        and=TRUE,
#'        partialmatch = "name",
#'        interval= "posright")
#' @export
get_dataset <- function(regulondb, dataset = NULL, attributes = NULL, filters = NULL, and = TRUE, interval=NULL, partialmatch=NULL){

  # Validate if attributes is a vector
  if(!is.null(attributes) & (!is.vector(attributes))){
      stop("Parameter 'attributes' must be a vector", call.= FALSE)
  }

  # Validate dataset
  if(is.null(dataset)) {stop("Non dataset provided", call.= FALSE)}
  if(!all(dataset %in% list_datasets(regulondb))){
    stop("Invalid dataset. See valid datasets in list_datasets()", call.= FALSE)
  }

  # Validate attributes
  if(!all(attributes %in% list_attributes(regulondb, dataset))){
    non.existing.attrs.index <- attributes %in% list_attributes(regulondb, dataset)
    non.existing.attrs <- attributes[!non.existing.attrs.index]
    stop("Attribute(s) ", paste0('"',paste(non.existing.attrs, collapse = ", "), '"'),
         " do not exist. See valid attributes in list_attributes()", call.= FALSE)
  }

  # Validate partialmatch

  if(!all(partialmatch %in% list_attributes(regulondb, dataset))) {
    non.existing.attrs.index <- partialmatch %in% list_attributes(regulondb, dataset)
    non.existing.attrs <- partialmatch[!non.existing.attrs.index]
    stop("Partialmatch ", paste0('"',paste(non.existing.attrs, collapse = ", "), '"'),
         " do not exist.", call.= FALSE)
  }

  if(!all(partialmatch %in% names(filters))) {
    non.existing.attrs.index <- partialmatch %in% names(filters)
    non.existing.attrs <- partialmatch[!non.existing.attrs.index]
    stop("Partialmatch ", paste0('"',paste(non.existing.attrs, collapse = ", "), '"'),
         " not defined in 'filters' ", call.= FALSE)
  }

  # Sets logical operator
  if(and){ operator <- "AND"
  }else{operator <- "OR"}

  if(is.null(filters) & is.null(attributes)){
    query <- paste0("SELECT * FROM ", dataset, ";")
  }else if (is.null(attributes) & !is.null(filters) ) {
    cond <- buil_condition(dataset, filters, operator, interval, partialmatch)
    query <- paste0("SELECT * FROM ", dataset, " WHERE ", cond, ";")
  }else if (!is.null(attributes) & is.null(filters)){
    query <- paste0("SELECT ", paste(attributes, collapse=" , ")," FROM ", dataset, ";")
  } else {
    cond <- buil_condition(dataset, filters, operator, interval, partialmatch)
    query <- paste("SELECT ", paste(attributes, collapse = " , "), "FROM ", dataset, " WHERE ", cond , ";") #Construct query
  }
  # Connect to database
  # regulon <- dbConnect(SQLite(), system.file("extdata", "regulondb_sqlite3.db", package = "regutools"))
  # Retrieve data
  result <- dbGetQuery( regulondb, query )
  #dbDisconnect(regulon)

  #Check if results exist
  if(!nrow(result)){
    stop("Your query returned no results.", call.= FALSE)
  }
  return(result)
}
