#' @title Extract sequences from a set of Transcription Factor Binding Sites (TFBSs)
#' @description Returns the sequences for a set of Transcription Factor Binding Sites.
#' @author José Alquicira Hernández, Jacques van Helden, Joselyn Chávez.
#' @param tfbs A collection of sites previously retrieved from RegulonDB with the function \code{\link{GetSitesForTF}}
#' @param seq.format Default: \code{fasta} output format. Supported: \code{fasta}, \code{wconsensus}.
#' @return A string with the TFBS sequences with "\n" breakline character.
#' @examples
#' ## Get the binding sites for AraC and extract their sequence in fasta format
#' tfbs <- GetSitesForTF(TF = "AraC")
#' tfbs.seq <- GetSiteSequences(tfbs)
#' cat(tfbs.seq)
#' @export


GetSiteSequences <- function(tfbs, seq.format = "fasta") {

  if(!(is.data.frame(tfbs) & all(names(tfbs) %in% c("ID", "left", "right", "strand", "sequence")))){
    stop("Must provide a valid collection of sites retrieved with 'GetSitesForTF' function")
  }

  if (seq.format == "fasta") {

    header <- paste0(">", tfbs$ID, "_", tfbs$left, ":", tfbs$right, ":", tfbs$strand)
    seq.string <- paste0(header,"\n",as.vector(tfbs$sequence))
  } else {

    stop(paste(seq.format, " is not a valid format for GetSiteSequences()."))

  }

  return(seq.string)

}
