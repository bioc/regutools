#' @title Get the binding sites for a Transcription Factor (TF)
#' @description Retrieve the binding sites and genome location for a given transcription factor.
#' @author José Alquicira Hernández, Jacques van Helden, Joselyn Chávez
#' @param TF name of the transcription factor of interest
#' @param seq.format Default: \code{table}. Supported: \code{fasta}, \code{wconsensus}.
#' @return A dataframe with the following columns corresponding to the TFBSs associated with the TF.
#' If seq.format = "fasta", returns a string with the TFBS sequences with '\n' breakline character.
#' If the transcription factor does not exist, returns NA.
#' \itemize{
#' \item ID (identifier)
#' \item left position
#' \item right position
#' \item strand
#' \item sequence
#' }
#' @examples
#' # Get the binding sites for AraC
#'
#' GetBindingSites(TF = "AraC")
#'
#' # Get the fasta binding sites sequences for AraC
#'
#' GetBindingSites(TF = "AraC", seq.format = "fasta")
#' @export

GetBindingSites <- function(TF, seq.format = "table", evidence = NULL) {

  if (!is.null(evidence)) {
    tfbs.raw <- tryCatch({
      GetAttr(dataset = "TF",
              attributes = c("name", "tfbs_unique", "evidence_reference"),
              filters = list(name = TF))
    },error = function(cond) return(NULL))
  } else {
    tfbs.raw <- tryCatch({
      GetAttr(dataset = "TF",
              attributes = c("name", "tfbs_unique"),
              filters = list(name = TF))
    },error = function(cond) return(NULL))
  }



  if(is.null(tfbs.raw)){
    return(NA)
  }else{
    # convert raw info into a table with 1 row per TFBS and 1 col per attribute
    tfbs.table <-  as.data.frame(strsplit(x = tfbs.raw$tfbs_unique, split = " "), col.names = "V1")
    tfbs.table <- t(as.data.frame(strsplit(x = as.character(tfbs.table$V1), split =  "\t")))
    if (!is.null(evidence)) {
      tfbs.table <- rbind(tfbs.table,as.data.frame(tfbs.raw$evidence_reference))
      colnames(tfbs.table) <- c("ID", "left", "right", "strand", "sequence", "evidence_reference")
    } else {
      colnames(tfbs.table) <- c("ID", "left", "right", "strand", "sequence")
    }
    rownames(tfbs.table) <- NULL
    tfbs.table <- as.data.frame(tfbs.table)
    if (seq.format == "table") {return(tfbs.table)}
    if (seq.format == "fasta") {
      header <- paste0(">", tfbs.table$ID, "_", tfbs.table$left, ":", tfbs.table$right, ":", tfbs.table$strand)
      seq.string <- paste0(header,"\n",as.vector(tfbs$sequence))
      return(seq.string)
    } else { stop("seq.format must be 'table' or 'fasta' ", call. = FALSE)}

  }
}

## TO DO:
## - add chromosome (for future multi-organisms RegulonDB)
## - add evidence
