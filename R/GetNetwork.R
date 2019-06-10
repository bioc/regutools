#' @title Return complete regulatory network.
#' @description This function retrieves all the regulation networks in regulonDB between TF-TF, GENE-GENE or TF-GENE depending on the parameter 'type'
#' @param type "TF-GENE", "TF-TF", "GENE-GENE"
#' @keywords regulation retrieval, TF, networks,
#' @author Carmina Barberena Jonas, Jesús Emiliano Sotelo Fonseca, José Alquicira Hernández.
#' @examples
#' # Retrieve regulation of 'araC'
#'
#' GetNetwork(type="TF-GENE")
#' @export

GetNetwork<-function(type="TF-GENE"){
  #Check type parameter
  if(! type %in% c("GENE-GENE","TF-GENE","TF-TF")){
    stop("Parameter 'type' must be TF-GENE, TF-TF, or GENE-GENE.",call.=FALSE)
  }

  #Get Network
  network<-GetAttr(attributes=c("regulator_name","regulated_name","effect"),
                   filters=list("network_type"=c(type)),
                   dataset="NETWORK")
  colnames(network)<-c("regulator","gene","effect")

  #Change effect to +, - and +/-
  network$effect<-sub(pattern="activator",replacement="+",x=network$effect)
  network$effect<-sub(pattern="repressor",replacement="-",x=network$effect)
  network$effect<-sub(pattern="dual",replacement="+/-",x=network$effect)

  return(network)
}


