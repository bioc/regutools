#' @title Return summary of gene regulation.
#' @description This function takes the output of \code{\link{GetGeneRegulation}} with format multirow,
#' onerow or table, or a vector with genes and retrieves information about the TFs and their regulated genes
#' @param geneRegulation Result from \code{\link{GetGeneRegulation}} or vector of genes
#' @return A data frame with the following columns:
#' \itemize{
#' \item The name or gene of TF
#' \item Regulated Genes per TF
#' \item Percent of regulated genes per TF
#' \item positive, negative or dual regulation
#' \item Name(s) of regulated genes
#' }
#' @keywords regulation retrieval, summary, networks,
#' @author Carmina Barberena Jonas, Jesús Emiliano Sotelo Fonseca, José Alquicira Hernández, Joselyn Chávez
#' @examples
#' # Retrieve summary of  araC's regulation
#'
#' araC_regulation<- GetGeneRegulation(genes = c("araC"),
#'                                     format = "multirow",
#'                                     output.type = "TF")
#' GetSummary(araC_regulation)
#'
#' # Retrieve summary of genes 'araC' and 'modB'
#'
#' GetSummary(geneRegulation = c("araC", "modB"))
#' GetSummary(geneRegulation = c("ECK120000050", "modB"))
#' @export

GetSummary<-function(geneRegulation){
  regulation <- geneRegulation

  #Regulation of a gene list
  if(class(regulation)%in%c("character","vector")){
    regulation <- GetGeneRegulation(genes=regulation)
  }

  #Checks if user changed the attribute format in the output data.frame from GetGeneRegulation.
  if( ! attributes(regulation)$format %in% c("multirow","onerow","table")){
    stop("The input data.frame attribute 'format' must be multirow, onerow, or table.",call.=FALSE)
  }

  # Conver onerow to multirow
  if(attributes(regulation)$format == "onerow"){
    regulation <- GetGeneRegulation(genes = as.character(regulation$genes))
  }

  # Convert table to multirow
  if(attributes(regulation)$format == "table") {
    regulation$genes <- rownames(regulation)
    regulation <- melt(regulation, id = "genes")
    colnames(regulation) <- c("genes", "regulators", "effect")
    regulation <- regulation[complete.cases(regulation), ]
  }

  TF_counts<-data.frame(table(regulation$regulators), stringsAsFactors=FALSE)

  summary<- apply(TF_counts, 1, function(x){

    #Get rows for a specific TF
    TF_data <-regulation[regulation[,"regulators"]==x[1],]

    #Count regulated effects
    effect <- table(TF_data$effect)
    #Adds regulation +, - or +/- in case their frequency is 0.
    if(! "+"%in%names(effect)){
      effect <-c(effect, "+"=0)
    }
    if(! "-"%in%names(effect)){
      effect <-c(effect, "-"=0)
    }
    if(! "+/-"%in%names(effect)){
      effect <-c(effect, "+/-"=0)
    }

    #Concatenates row to include the summary information for each TF
    summary_row<-c(TF=x[1], #TF name
                   regulated_number=x[2],#Number of genes regulated per TF
                   regulated_percentage= (as.numeric(x[2])/length(regulation$genes))*100,#Percent of genes in query regulated per TF
                   activator= effect["+"],#Frequency of activation interactions
                   repressor= effect["-"],#Frequency of repression interactions
                   dual=effect["+/-"] ,#Frequency of dual interactions
                   regulated= paste(TF_data$genes, collapse=", ")#List of genes regulated per TF
    )

    return(summary_row)

  })

  summary <- data.frame(t(summary)) #Format summary as a data.frame
  colnames(summary) <- c("TF", "Regulated_genes_per_TF", "Percent","+","-","+/-","Regulated_genes")

  return(summary)
}



