#' @title Extract data from RegulonDB
#' @description This function retrieves data from RegulonDB. Attributes from datasets can be selected and filtered.
#' @param regulondb A regulondb object.
#' @param dataset Dataset of interest. Use the function list_datasets for an overview of valid datasets.
#' @param attributes Vector of attributes to be retrieved.
#' @param filters List of filters to be used. The names should correspond to the attribute and the values correspond to the condition for selection.
#' @param and Logical argument. If FALSE, filters will be considered under the "OR" operator
#' @param partialmatch name of the condition(s) with a string pattern for full or partial match in the query
#' @param interval the filters whose values will be considered as interval
#' @author Carmina Barberena Jonas, Jesús Emiliano Sotelo Fonseca, José Alquicira Hernández, Joselyn Chávez
#' @return A data frame.
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
    cond <- build_condition(regulondb, dataset, filters, operator, interval, partialmatch)
    query <- paste0("SELECT * FROM ", dataset, " WHERE ", cond, ";")
  }else if (!is.null(attributes) & is.null(filters)){
    query <- paste0("SELECT ", paste(attributes, collapse=" , ")," FROM ", dataset, ";")
  } else {
    cond <- build_condition(regulondb, dataset, filters, operator, interval, partialmatch)
    query <- paste("SELECT ", paste(attributes, collapse = " , "), "FROM ", dataset, " WHERE ", cond , ";") #Construct query
  }
  result <- dbGetQuery( regulondb, query )

  #Check if results exist
  if(!nrow(result)){
    warning("Your query returned no results.")
  }

  result <- new( "regulondb_result",
                 DataFrame(result), organism=regulondb@organism,
       genome_version=regulondb@genome_version,
       database_version=regulondb@database_version,
       dataset=dataset )
  return(result)
}

#' @title Function to convert output of regulondb queries to Bioconductor S4 objects
#' @description This function
#' @param regulondb_result A data frame that resulted from a call to the function get_dataset.
#' @keywords data retrieval, attributes, filters,
#' @author Alejandro Reyes
#' @importFrom GenomicRanges GRanges
#' @export
convert_to_granges <- function( regulondb_result ){
  stopifnot( is( regulondb_result, "data.frame" ) )
  if( !all(c("source", "dataset") %in% names( attributes( regulondb_result ) ) ) ){
    stop( "Attributes 'source' and/or 'dataset' were not found in the input.
         Is the input the result of get_dataset?",
         call.=FALSE )
  }
  dataset <- attributes( regulondb_result )$dataset
  if( dataset %in% c("GENE", "DNA_OBJECTS") ){
    posLeft <- "posleft"
    posRight <- "posright"
  }else if( dataset == "" ){
  }
  if( all( c( posLeft, posRight ) %in% names( regulondb_result ) ) ){
    keep <- !(is.na(regulondb_result$posleft) | is.na(regulondb_result$posright))
    grdata <- with( regulondb_result[which(keep),],
          GRanges( attributes( regulondb_result )$organism,
                   IRanges(start=posleft, end=posright)) )
    if( "strand" %in% names(regulondb_result) ){
      strand(grdata) <- ifelse( regulondb_result$strand[which(keep)] == "forward", "+", "-" )
    }
    mcols(grdata) <- DataFrame(
      regulondb_result[keep,!colnames( regulondb_result ) %in% c(posLeft, posRight, "strand"),drop=FALSE])
    if(sum(!keep)> 0)
      warning( sprintf("Dropped %s entries since NAs were found in the genomic coordinates", sum(keep)) )
    grdata
  }else{
    stop( sprintf( "Not enough information to convert into a GRanges object. Please make sure that the input the following columns: \n\t%s",
      paste(c(posLeft, posRight), collapse="\n\t") ), call.=FALSE )
  }
  grdata
}
