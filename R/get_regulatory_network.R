#' @title Return complete regulatory network.
#' @description This function retrieves all the regulation networks in regulonDB between TF-TF, GENE-GENE or TF-GENE depending on the parameter 'type'
#' @param regulondb A regulondb object
#' @param regulator Name of TF or gene that acts as regulator. If NULL, the function retrieves all existent networks in the regulonDB.
#' @param type "TF-GENE", "TF-TF", "GENE-GENE"
#' @param cytograph If TRUE, displays network in Cytoscape. This option requires previous instalation and launch of Cytoscape.
#'
#' @keywords regulation retrieval, TF, networks,
#' @author Carmina Barberena Jonas, Jesús Emiliano Sotelo Fonseca, José Alquicira Hernández, Joselyn Chávez
#' @examples
#' # Retrieve regulation of 'araC'
#'
#' get_regulatory_network(regulator = "AraC", type="TF-GENE")
#'
#' # Retrieve all GENE-GENE networks
#'
#' get_regulatory_network(type = "GENE-GENE")
#'
#' # Retrieve TF-GENE network of AraC and display in Cytoscape
#'
#' get_regulatory_network(regulator = "AraC", type = "TF-GENE", cytograph = TRUE)
#' @export

get_regulatory_network <- function(regulondb, regulator = NULL, type="TF-GENE", cytograph = FALSE){
  #Check type parameter
  if(! type %in% c("GENE-GENE","TF-GENE","TF-TF")){
    stop("Parameter 'type' must be TF-GENE, TF-TF, or GENE-GENE.",call.=FALSE)
  }

  #Get Network
  if (!is.null(regulator)) {
    network <- as.data.frame(get_dataset(regulondb, attributes=c("regulator_name","regulated_name","effect"),
                                             filters=list(
                                               "network_type"= type,
                                               "regulator_name" = regulator),
                                             dataset="NETWORK") )
  } else {
    network <- as.data.frame(get_dataset(regulondb, attributes=c("regulator_name","regulated_name","effect"),
                     filters=list("network_type"=c(type)),
                     dataset="NETWORK") ) }

  if (cytograph) { #CytoScape conection
    tryCatch(cytoscapePing(), error = function(e) {stop("Launch Cytoscape before running GetNetwork()", call. = FALSE)})
    colnames(network) <- c("source", "target", "interaction")
    my_new_network <- createNetworkFromDataFrames(edges =  network)
    setVisualStyle('Sample1')
    setEdgeColorMapping('interaction', c('activator', 'repressor','dual'), c('#339900', '#CC3300', '#0033CC'))
  } else {
    colnames(network)<-c("regulator","regulated","effect")

    #Change effect to +, - and +/-
    network$effect<-sub(pattern="activator",replacement="+",x=network$effect)
    network$effect<-sub(pattern="repressor",replacement="-",x=network$effect)
    network$effect<-sub(pattern="dual",replacement="+/-",x=network$effect)

    return(network)
    }

}


