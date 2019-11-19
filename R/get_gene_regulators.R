#' @title Get TFs or genes that regulate the genes of interest
#' @description Given a list of genes (name, bnumber or GI),
#' get all transcription factors or genes that regulate them.
#' The effect of regulators over the gene of interest can be positive (+), negative (-) or dual (+/-)
#' @param genes Vector of genes (name, bnumber or GI).
#' @param format Output format: multirow, onerow, table
#' @param output.type How regulators will be represented: "TF"/"GENE"
#' @keywords regulation retrieval, TFs, networks,
#' @author Carmina Barberena Jonas, Jesús Emiliano Sotelo Fonseca, José Alquicira Hernández, Joselyn Chávez
#' @examples
#' # Get Transcription factors that regulate araC in one row
#' get_gene_regulators(genes = c("araC"),output.type = "TF",format = "onerow")
#'
#' # Get genes that regulate araC in table format
#' get_gene_regulators(genes = c("araC"),output.type = "GENE",format = "table")
#' @export

get_gene_regulators <- function(genes,format="multirow",output.type="TF"){
  #Check genes parameter class
  if(! class(genes) %in% c("vector","list","character")){
    stop("Parameter 'genes' must be a character vector or list.",call.=FALSE)
  }
  #Check format parameter
  if(! format %in% c("multirow","onerow", "table")){
    stop("Parameter 'format' must be multirow, onerow, or table.",call.=FALSE)
  }
  #Check output.type
  if(! output.type %in% c("TF","GENE")){
    stop("Parameter 'output.type' must be either 'TF' or 'GENE'",call.=FALSE)
  }

  #Convert GIs to gene names
  equivalence_table<- get_dataset(dataset="GENE", attributes=c("id","name","bnumber","gi"))
  genes<-lapply(as.list(genes),function(gene){
    if(grepl("ECK12",gene)){
      return(equivalence_table[equivalence_table[,"id"]==gene,"name"])
    }
    if(grepl("^b[0-9]", gene)){
      return(equivalence_table[equivalence_table[,"bnumber"] %in% as.character(gene) , "name"])
      }
    ifelse(grepl("[a-z]",gene), return(gene), return(equivalence_table[equivalence_table[,"gi"] %in% as.character(gene) , "name"]))
  })

  #Checks for type of network
  if (output.type == "TF"){
    network.type <- "TF-GENE"
  } else if (output.type == "GENE"){
    network.type <- "GENE-GENE"
  }

  #Retrieve data from NETWORK table
  regulation <- get_dataset(attributes=c("regulated_name","regulator_name","effect"),
                        filters=list("regulated_name"=genes, "network_type"=network.type),
                        dataset="NETWORK")
  colnames(regulation)<-c("genes","regulators","effect")

  #Change effect to +, - and +/-
  regulation$effect<-sub(pattern="activator",replacement="+",x=regulation$effect)
  regulation$effect<-sub(pattern="repressor",replacement="-",x=regulation$effect)
  regulation$effect<-sub(pattern="dual",replacement="+/-",x=regulation$effect)

  #Format output
  #Multirow
  if(format=="multirow"){
    #Add internal attribute "format" to use in GetSummary function.
    attributes(regulation)$format <- format
    return(regulation)

  #Onerow
  } else if (format=="onerow"){
    regulation<-lapply(as.list(genes), function(x){
      genereg<-regulation[regulation[,"genes"]==x,]
      genereg<- paste(paste(genereg$regulators, genereg$effect,sep="(", collapse="), "),")",sep="")
    })
    regulation<-unlist(regulation)
    genes<-unlist(genes)
    regulation<-data.frame(genes,regulation)
    colnames(regulation)<-c("genes","regulators")

    #Add internal attribute "format" to use in GetSummary function.
    attributes(regulation)$format <- format
    return(regulation)

  #Table
  } else if (format=="table"){
    #Empty dataframe
    rtable<-data.frame(matrix(nrow=length(genes),ncol=length(c(unique(regulation$regulators)))))
    colnames(rtable)<-unique(regulation$regulators)
    rownames(rtable)<-genes

    #Fill dataframe with regulation
    for(i in 1:dim(regulation)[1]){
      rtable[regulation[i,1],regulation[i,2]]<-regulation[i,3]
    }
    regulation<-rtable

    #Add internal attribute "format" to use in GetSummary function.
    attributes(regulation)$format <- format
    return(regulation)

  }
}

